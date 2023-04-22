use libc::{c_double, c_int};
use ndarray::ArrayViewMut2;
use std::slice;

use super::low;

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

extern "C" fn fcn<T: Callback>(
    n: *const c_int,
    x_ptr: *const c_double,
    y_ptr: *const c_double,
    f_ptr: *mut c_double,
    _: *mut c_double,
    callback_ptr: *mut c_int,
) {
    let (f, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(f_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *x_ptr,
        )
    };
    let callback = unsafe { &*(callback_ptr as *const T) };
    callback.fcn(t, y, f);
}

extern "C" fn jac<T: Callback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    pd_ptr: *mut c_double,
    nrow_pd: *const c_int,
    _: *mut c_double,
    callback_ptr: *mut c_double,
) {
    let (nrow, pd, y, t) = unsafe {
        let y = slice::from_raw_parts(y_ptr, *n as usize);
        let t = *t_ptr;
        let nrow = *nrow_pd as usize;
        let pd = slice::from_raw_parts_mut(pd_ptr, y.len() * nrow);
        (nrow, pd, y, t)
    };
    let callback = unsafe { &*(callback_ptr as *const T) };
    let mut pd = ArrayViewMut2::from_shape((y.len(), nrow), pd).unwrap();
    pd.swap_axes(0, 1); // make pd fortran-ordered
    callback.jac(t, y, pd);
}

extern "C" fn mas<T: Callback>(
    n: *const c_int,
    m_ptr: *mut c_double,
    nrow_m: *const c_int,
    _: *mut c_double,
    callback_ptr: *mut c_int,
) {
    let (m, n, nrow) = unsafe {
        let nrow = *nrow_m as usize;
        let n = *n as usize;
        let m = slice::from_raw_parts_mut(m_ptr, n * nrow);
        (m, n, nrow)
    };
    let callback = unsafe { &*(callback_ptr as *const T) };
    let mut m = ArrayViewMut2::from_shape((n, nrow), m).unwrap();
    m.swap_axes(0, 1); // make pd fortran-ordered
    callback.mas(m);
}

extern "C" fn solout<T: Callback>(
    nr: *const c_int,
    xold: *const c_double,
    x: *const c_double,
    y_ptr: *const c_double,
    cont_ptr: *const c_double,
    lrc: *const c_int,
    n_ptr: *const c_int,
    _: *mut c_double,
    callback_ptr: *mut c_int,
    irtrn: *mut c_int,
) {
    let (nr, xold, x, y) = unsafe {
        let y = slice::from_raw_parts(y_ptr, *n_ptr as usize);
        (*nr as usize, *xold, *x, y)
    };
    let callback = unsafe { &*(callback_ptr as *const T) };
    let contra = |s: f64| {
        let mut interp = Vec::with_capacity(y.len());
        for i in 0_i32..y.len() as i32 {
            let val = unsafe { low::contra_(&i, &s, cont_ptr, lrc) };
            interp.push(val);
        }
        interp
    };

    let ret = callback.solout(nr, xold, x, y, &contra);
    unsafe {
        *irtrn = ret;
    }
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

    let itol = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let rpar = std::ptr::null_mut();
    let mut idid = 0_i32;

    unsafe {
        let ipar = std::mem::transmute::<&T, *mut c_int>(callback);
        low::radau_(
            &n,
            fcn::<T>,
            &mut t0,
            y.as_mut_ptr(),
            &t1,
            &h,
            rtol.as_ptr(),
            atol.as_ptr(),
            &itol,
            jac::<T>,
            &ijac,
            &mljac,
            &mujac,
            mas::<T>,
            &imas,
            &mlmas,
            &mumas,
            solout::<T>,
            &iout,
            rwork.as_mut_ptr(),
            &lrw,
            iwork.as_mut_ptr(),
            &liw,
            rpar,
            ipar,
            &mut idid,
        );
    }
    idid
}
