use libc::{c_double, c_int};
use libffi::high::{Closure10, Closure5, Closure6, Closure7};
use ndarray::ArrayViewMut2;
use once_cell::sync::Lazy;
use std::slice;
use std::sync::Mutex;

use super::low;

static FLAG: Lazy<Mutex<()>> = Lazy::new(|| Mutex::new(()));

pub trait Callback {
    fn fcn<'b>(&self, t: f64, y: &'b [f64], f: &'b mut [f64]);
    fn jac<'b>(&self, t: f64, y: &'b [f64], pd: ArrayViewMut2<'b, f64>);
    fn mas<'b>(&self, mas: ArrayViewMut2<'b, f64>);
    fn solout<'b, CONTRA>(
        &self,
        nr: usize,
        t_old: f64,
        t: f64,
        y: &'b [f64],
        contra: &'b CONTRA,
    ) -> i32
    where
        CONTRA: 'b + Fn(f64) -> Vec<f64>;
    fn jacobian_type(&self) -> CallbackType;
    fn mass_type(&self) -> CallbackType;
    fn use_solout(&self) -> bool;
}

#[derive(Debug, Clone)]
pub struct Control {
    pub atol: f64,
    pub rtol: f64,
    pub max_steps: usize,
}

impl Default for Control {
    fn default() -> Self {
        Self {
            atol: 1e-8,
            rtol: 1e-8,
            max_steps: 500,
        }
    }
}

pub enum CallbackType {
    FullMatrix {
        supplied: bool,
    },
    BandedMatrix {
        supplied: bool,
        mu: usize,
        ml: usize,
    },
}

impl CallbackType {
    fn settings(&self, n: i32) -> (i32, i32, i32) {
        match self {
            Self::FullMatrix { supplied } => {
                if *supplied {
                    (1, n, n)
                } else {
                    (0, n, n)
                }
            }
            Self::BandedMatrix { supplied, mu, ml } => {
                if *supplied {
                    (1, *mu as i32, *ml as i32)
                } else {
                    (0, *mu as i32, *ml as i32)
                }
            }
        }
    }

    fn l(&self, n: i32) -> i32 {
        match self {
            Self::FullMatrix { supplied: _ } => n as i32,
            Self::BandedMatrix {
                supplied: _,
                mu,
                ml,
            } => (*mu + *ml + 1) as i32,
        }
    }

    fn le(&self, n: i32) -> i32 {
        match self {
            Self::FullMatrix { supplied: _ } => n as i32,
            Self::BandedMatrix {
                supplied: _,
                mu,
                ml,
            } => 2 * (*mu + *ml + 1) as i32,
        }
    }
}

pub fn real_work_space(
    n: usize,
    jac_type: &CallbackType,
    mas_type: &CallbackType,
    _option: &Control,
) -> Vec<f64> {
    let nsmax = 7;
    let ljac = jac_type.l(n as i32);
    let lmas = mas_type.l(n as i32);
    let le = jac_type.le(n as i32);
    let rwork = vec![0_f64; n * (ljac + lmas + nsmax * le + 3 * nsmax + 3) as usize + 20];

    rwork
}

pub fn integer_work_space(
    n: usize,
    _jac_type: &CallbackType,
    _mas_type: &CallbackType,
    option: &Control,
) -> Vec<i32> {
    let nsmax = 7;
    let mut iwork = vec![0_i32; ((2 + (nsmax - 1) / 2) * n + 20) as usize];

    iwork[1] = option.max_steps as i32;

    iwork
}

pub fn radau<'a, 'b, T: Callback>(
    callback: &'a T,
    y: &mut [f64],
    mut t0: f64,
    t1: f64,
    h: f64,
    rtol: &[f64],
    atol: &[f64],
    jac_type: &CallbackType,
    mas_type: &CallbackType,
    use_solout: bool,
    rwork: &mut [f64],
    iwork: &mut [i32],
) -> i32
where
{
    let n = y.len() as i32;
    let (ijac, mujac, mljac) = jac_type.settings(n);
    let (imas, mumas, mlmas) = mas_type.settings(n);
    let iout = if use_solout { 1_i32 } else { 0_i32 };

    let fcn = |n: *const c_int,
               x_ptr: *const c_double,
               y_ptr: *const c_double,
               f_ptr: *mut c_double,
               _: *mut c_double,
               _: *mut c_int| {
        let (f, y, t) = unsafe {
            (
                slice::from_raw_parts_mut(f_ptr, *n as usize),
                slice::from_raw_parts(y_ptr, *n as usize),
                *x_ptr,
            )
        };
        callback.fcn(t, y, f);
    };
    let closure = Closure6::new(&fcn);
    let fcn = closure.code_ptr();

    let jac = |n: *const c_int,
               t_ptr: *const c_double,
               y_ptr: *const c_double,
               pd_ptr: *mut c_double,
               nrow_pd: *const c_int,
               _: *mut c_double,
               _: *mut c_double| {
        let (nrow, pd, y, t) = unsafe {
            debug_assert!(*nrow_pd == *n || *nrow_pd == 1 + mujac + mljac);
            let y = slice::from_raw_parts(y_ptr, *n as usize);
            let t = *t_ptr;
            let nrow = *nrow_pd as usize;
            let pd = slice::from_raw_parts_mut(pd_ptr, y.len() * nrow);
            (nrow, pd, y, t)
        };
        let mut pd = ArrayViewMut2::from_shape((y.len(), nrow), pd).unwrap();
        pd.swap_axes(0, 1); // make pd fortran-ordered
        callback.jac(t, y, pd);
    };
    let closure2 = Closure7::new(&jac);
    let jac = closure2.code_ptr();

    let mas = |n: *const c_int,
               m_ptr: *mut c_double,
               nrow_m: *const c_int,
               _: *mut c_double,
               _: *mut c_int| {
        let (m, n, nrow) = unsafe {
            debug_assert!(*nrow_m == *n || *nrow_m == 1 + mumas + mlmas);
            let nrow = *nrow_m as usize;
            let n = *n as usize;
            let m = slice::from_raw_parts_mut(m_ptr, n * nrow);
            (m, n, nrow)
        };
        let mut m = ArrayViewMut2::from_shape((n, nrow), m).unwrap();
        m.swap_axes(0, 1); // make pd fortran-ordered
        callback.mas(m);
    };
    let closure3 = Closure5::new(&mas);
    let mas = closure3.code_ptr();

    let solout = |nr: *const c_int,
                  xold: *const c_double,
                  x: *const c_double,
                  y_ptr: *const c_double,
                  cont_ptr: *const c_double,
                  lrc: *const c_int,
                  n_ptr: *const c_int,
                  _: *mut c_double,
                  _: *mut c_int,
                  irtrn: *mut c_int| {
        let (nr, xold, x, y) = unsafe {
            let y = slice::from_raw_parts(y_ptr, *n_ptr as usize);
            (*nr as usize, *xold, *x, y)
        };
        let contra = |s: f64| {
            let mut interp = Vec::with_capacity(y.len());
            for i in 0_i32..n {
                let val = unsafe { low::contra_(&i, &s, cont_ptr, lrc) };
                interp.push(val);
            }
            interp
        };

        let ret = callback.solout(nr, xold, x, y, &contra);
        unsafe {
            *irtrn = ret;
        }
    };
    let closure4 = Closure10::new(&solout);
    let solout = closure4.code_ptr();

    let itol = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let mut rpar = vec![0.; 0];
    let mut ipar = vec![0; 0];
    let mut idid = 0_i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::radau_(
            &n,
            *fcn,
            &mut t0,
            y.as_mut_ptr(),
            &t1,
            &h,
            rtol.as_ptr(),
            atol.as_ptr(),
            &itol,
            *jac,
            &ijac,
            &mljac,
            &mujac,
            *mas,
            &imas,
            &mlmas,
            &mumas,
            *solout,
            &iout,
            rwork.as_mut_ptr(),
            &lrw,
            iwork.as_mut_ptr(),
            &liw,
            rpar.as_mut_ptr(),
            ipar.as_mut_ptr(),
            &mut idid,
        );
    }
    idid
}
