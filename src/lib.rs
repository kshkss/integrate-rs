use libc::{c_double, c_int};
use libffi::high::Closure4;
use libffi::high::Closure7;
use std::slice;

#[link(name = "gfortran")]
extern "C" {
    /// Call `DLSODE` subroutine from ODEPACK
    ///
    /// For info on passed arguments look inside ODEPACK.
    pub fn dlsode_(
        f: extern "C" fn(*const c_int, *const c_double, *mut c_double, *mut c_double),
        neq: &c_int,
        y: *mut c_double,
        t: &mut c_double,
        tout: &c_double,
        itol: &c_int,
        rtol: &c_double,
        atol: &c_double,
        itask: &c_int,
        istate: &mut c_int,
        iopt: &c_int,
        rwork: *mut c_double,
        lrw: &c_int,
        iwork: *mut c_int,
        liw: &c_int,
        jac: extern "C" fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *const c_int,
            *const c_int,
            *mut c_double,
            *const c_int,
        ),
        mf: &c_int,
    );
}

/// A dummy function to pass to `dlsode_` in case the user does not want to specify a Jacobian.
pub extern "C" fn fake_jacobian(
    _neq: *const c_int,
    _t: *const c_double,
    _y: *const c_double,
    _ml: *const c_int,
    _mu: *const c_int,
    _pd: *mut c_double,
    _nr: *const c_int,
) {
}

use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

use ndarray::prelude::*;

#[derive(Clone, Copy)]
enum MethodFlag {
    EstimateFullJacobian,
    FullJacobian,
    EstimateBandedJacobian(usize, usize),
    BandedJacobian(usize, usize),
}

pub struct Jacobian<'a> {
    mf: MethodFlag,
    udf: Box<
        dyn 'a
            + Fn(
                *const c_int,
                *const c_double,
                *const c_double,
                *const c_int,
                *const c_int,
                *mut c_double,
                *const c_int,
            ),
    >,
}

impl<'a> Jacobian<'a> {
    fn new<G>(mf: MethodFlag, g: G) -> Self
    where
        G: 'a
            + Fn(
                *const c_int,
                *const c_double,
                *const c_double,
                *const c_int,
                *const c_int,
                *mut c_double,
                *const c_int,
            ),
    {
        let udf = Box::new(g);
        Self { mf, udf }
    }

    fn method_flag(&self) -> c_int {
        use MethodFlag::*;
        match self.mf {
            EstimateFullJacobian => 22 as c_int,
            EstimateBandedJacobian(_, _) => 25 as c_int,
            FullJacobian => 21,
            BandedJacobian(_, _) => 24 as c_int,
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<c_double> {
        use MethodFlag::*;
        match self.mf {
            EstimateFullJacobian | FullJacobian => {
                vec![0. as c_double; 22 + 9 * n_eq + n_eq * n_eq]
            }
            EstimateBandedJacobian(ml, mu) | BandedJacobian(ml, mu) => {
                vec![0. as c_double; 22 + 10 * n_eq + (2 * ml + mu) * n_eq]
            }
        }
    }

    fn integer_work_space(&self, n_eq: usize) -> Vec<c_int> {
        use MethodFlag::*;
        match self.mf {
            EstimateFullJacobian | FullJacobian => vec![0 as c_int; 22 + n_eq],
            EstimateBandedJacobian(ml, mu) | BandedJacobian(ml, mu) => {
                let mut iwork = vec![0 as c_int; 22 + n_eq];
                iwork[0] = ml as c_int;
                iwork[1] = mu as c_int;
                iwork
            }
        }
    }
}

pub struct Lsode<'a> {
    dydt: Box<dyn 'a + Fn(&[f64], f64) -> Vec<f64>>,
    jacobian: Jacobian<'a>,
}

impl<'a> Lsode<'a> {
    pub fn new<F>(dydt: F) -> Self
    where
        F: 'a + Fn(&[f64], f64) -> Vec<f64>,
    {
        let g = |_neq: *const c_int,
                 _t: *const c_double,
                 _y: *const c_double,
                 _ml: *const c_int,
                 _mu: *const c_int,
                 _pd: *mut c_double,
                 _nr: *const c_int| {};
        Self {
            dydt: Box::new(dydt),
            jacobian: Jacobian {
                mf: MethodFlag::EstimateFullJacobian,
                udf: Box::new(g),
            },
        }
    }

    pub fn estimate_banded_jacobian(self, ml: usize, mu: usize) -> Self {
        let g = |_neq: *const c_int,
                 _t: *const c_double,
                 _y: *const c_double,
                 _ml: *const c_int,
                 _mu: *const c_int,
                 _pd: *mut c_double,
                 _nr: *const c_int| {};
        Lsode {
            jacobian: Jacobian {
                mf: MethodFlag::EstimateBandedJacobian(ml, mu),
                udf: Box::new(g),
            },
            ..self
        }
    }

    pub fn with_full_jacobian<G>(self, udf: G) -> Self
    where
        G: 'a + Fn(&[f64], f64) -> Array2<f64>,
    {
        let g = move |n: *const c_int,
                      t_ptr: *const c_double,
                      y_ptr: *const c_double,
                      _ml: *const c_int,
                      _mu: *const c_int,
                      dy_ptr: *mut c_double,
                      _nrow: *const c_int| {
            let n = unsafe { *n as usize };
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, n * n),
                    slice::from_raw_parts(y_ptr, n),
                    *t_ptr,
                )
            };
            let mut dy = ArrayViewMut2::<f64>::from_shape((n, n), dy).expect("somthing wrong");
            dy.swap_axes(0, 1); // make dy fortran-ordered
            let dy_new = udf(y, t);
            dy.assign(&dy_new);
        };
        Lsode {
            jacobian: Jacobian {
                mf: MethodFlag::FullJacobian,
                udf: Box::new(g),
            },
            ..self
        }
    }

    pub fn with_banded_jacobian<G>(self, ml: usize, mu: usize, udf: G) -> Self
    where
        G: 'a + Fn(&[f64], f64) -> Vec<Array1<f64>>,
    {
        let g = move |n: *const c_int,
                      t_ptr: *const c_double,
                      y_ptr: *const c_double,
                      _ml: *const c_int,
                      _mu: *const c_int,
                      dy_ptr: *mut c_double,
                      _nrow: *const c_int| {
            let n = unsafe { *n as usize };
            let ms: usize = todo!();
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, ms),
                    slice::from_raw_parts(y_ptr, n),
                    *t_ptr,
                )
            };
            let mut dy = ArrayViewMut1::<f64>::from(dy);
            let dy_new = udf(y, t);
            /*
            dy.assign(&dy_new);
            */
            todo!()
        };
        todo!();
        Lsode {
            jacobian: Jacobian {
                mf: MethodFlag::BandedJacobian(ml, mu),
                udf: Box::new(g),
            },
            ..self
        }
    }

    /// Solves system of ODEs for times in `t_dense`.
    /// First time in `t_dense` has to be the initial time.
    ///
    /// Each equation in the system of ODEs has the form:
    ///
    /// > *dy/dt = f(y, t)*
    ///
    /// The function expects the function *f* as the first argument `rhs`.
    /// Initial state is given in `y0`.
    ///
    /// # Example
    ///
    /// ```
    /// let y0 = [1.0];
    /// let ts = vec![0.0, 1.0];
    /// let f = |y: &[f64], t: f64| {
    ///     let mut dy = vec![0.0];
    ///     dy[0] = t * y[0];
    ///     dy
    ///     };
    /// let sol = lsode::Lsode::new(f).solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// assert!((sol[1][0] - y0[0]*0.5_f64.exp()).abs() < 1e-3, "error too large");
    /// ```
    pub fn solve(&self, y0: &[f64], t: &[f64], atol: f64, rtol: f64) -> Vec<Vec<f64>> {
        let f = |n: *const c_int,
                 t_ptr: *const c_double,
                 y_ptr: *mut c_double,
                 dy_ptr: *mut c_double| {
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, *n as usize),
                    slice::from_raw_parts(y_ptr, *n as usize),
                    *t_ptr,
                )
            };
            let dy_new = (self.dydt)(y, t);
            for (dest, &deriv) in dy.iter_mut().zip(dy_new.iter()) {
                *dest = deriv;
            }
        };
        let closure = Closure4::new(&f);
        let call = closure.code_ptr();

        let jacobian_closure = Closure7::new(&self.jacobian.udf);
        let call_jacobian = jacobian_closure.code_ptr();

        let mut y: Vec<f64> = y0.to_vec();
        let n = y0.len();
        let mut t0 = t[0];

        let itol = 1;
        let itask = 1;
        let iopt = 0;
        let mut istate = 1;
        let mf = self.jacobian.method_flag();

        let mut rwork = self.jacobian.real_work_space(n);
        let mut iwork = self.jacobian.integer_work_space(n);
        let lrw = rwork.len();
        let liw = iwork.len();

        let mut result = Vec::new();

        let _lock = FLAG.lock().unwrap();
        for &tout in t.iter() {
            unsafe {
                dlsode_(
                    *call,
                    &(n as i32),
                    y.as_mut_ptr(),
                    &mut t0,
                    &tout,
                    &itol,
                    &rtol,
                    &atol,
                    &itask,
                    &mut istate,
                    &iopt,
                    rwork.as_mut_ptr(),
                    &(lrw as i32),
                    iwork.as_mut_ptr(),
                    &(liw as i32),
                    *call_jacobian,
                    &mf,
                );
            }

            result.push(y.clone());
        }
        result
    }
}

/// Generates vector from `start` to `stop` with constant  spacing.
///
/// # Example
///
/// ```
/// assert_eq!(lsode::linspace(0.0, 2.0, 3), vec![0.0, 1.0, 2.0]);
/// ```
pub fn linspace(start: f64, stop: f64, n_points: usize) -> Vec<f64> {
    (0..n_points)
        .map(|i| start + i as f64 * (stop - start) / (n_points - 1) as f64)
        .collect()
}

#[cfg(test)]
mod tests {
    use ndarray::prelude::*;

    #[test]
    fn ndarray_swap_axis() {
        let a = array![[0, 1], [2, 3]];
        assert_eq!(a[[0, 1]], 1);
        assert_eq!(a[[1, 0]], 2);
        for (ai, c) in a.as_slice().unwrap().iter().zip(&[0, 1, 2, 3]) {
            assert_eq!(ai, c);
        }

        let mut b = array![[0, 0], [0, 0]];
        b.assign(&a);
        assert_eq!(b[[0, 1]], 1);
        assert_eq!(b[[1, 0]], 2);
        for (bi, c) in b.as_slice().unwrap().iter().zip(&[0, 1, 2, 3]) {
            assert_eq!(bi, c);
        }

        b.swap_axes(0, 1);
        b.assign(&a);
        assert_eq!(b[[0, 1]], 1);
        assert_eq!(b[[1, 0]], 2);
        b.swap_axes(0, 1);
        for (bi, c) in b.as_slice().unwrap().iter().zip(&[0, 2, 1, 3]) {
            assert_eq!(bi, c);
        }
    }
}
