use libc::{c_double, c_int};
use libffi::high::Closure4;
use libffi::high::Closure7;
use std::slice;

pub use ndarray;

pub mod odepack;

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
pub mod dlsode;
use dlsode::Lsode;
mod dlsodes;
use dlsodes::Lsodes;

use once_cell::sync::Lazy;
use std::sync::Mutex;

static FLAG: Lazy<Mutex<()>> = Lazy::new(|| {
    Mutex::<()>::new(())
});

use ndarray::prelude::*;

pub struct Adams<'a> {
    odepack: Lsode<'a>,
}

impl<'a> Adams<'a> {
    /// Solves system of ODEs for times in `t`.
    /// First time in `t` has to be the initial time.
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
    /// let sol = integrate::Adams::new(f).solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// assert!((sol[1][0] - y0[0]*0.5_f64.exp()).abs() < 1e-3, "error too large");
    /// ```
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
        let odepack = Lsode {
            mf: dlsode::Generator::None,
            udf: Box::new(g),
            dydt: Box::new(dydt),
            max_steps: 500,
        };
        Self { odepack }
    }

    pub fn max_steps(&mut self, n: usize) {
        self.odepack.max_steps = n;
    }

    pub fn solve(&self, y0: &[f64], t: &[f64], atol: f64, rtol: f64) -> Vec<Vec<f64>> {
        self.odepack.solve(y0, t, atol, rtol)
    }
}

pub enum BDF<'a> {
    Dense(Lsode<'a>),
    Sparse(Lsodes<'a>),
}

impl<'a> BDF<'a> {
    /// Solves system of ODEs for times in `t`.
    /// First time in `t` has to be the initial time.
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
    /// let sol = integrate::BDF::new(f).solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// assert!((sol[1][0] - y0[0]*0.5_f64.exp()).abs() < 1e-3, "error too large");
    /// ```
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
        Self::Dense(Lsode {
            dydt: Box::new(dydt),
            mf: dlsode::Generator::InternalFull,
            udf: Box::new(g),
            max_steps: 500,
        })
    }

    pub fn max_steps(self, max_steps: usize) -> Self {
        use BDF::*;
        match self {
            Dense(odepack) => Dense(Lsode {
                max_steps,
                ..odepack
            }),
            Sparse(odepack) => Sparse(Lsodes {
                max_steps,
                ..odepack
            }),
        }
    }

    /// Default
    ///
    pub fn gen_full_jacobian(self) -> Self {
        let g = |_neq: *const c_int,
                 _t: *const c_double,
                 _y: *const c_double,
                 _ml: *const c_int,
                 _mu: *const c_int,
                 _pd: *mut c_double,
                 _nr: *const c_int| {};
        use BDF::*;
        match self {
            Dense(odepack) => Dense(Lsode {
                mf: dlsode::Generator::InternalFull,
                udf: Box::new(g),
                ..odepack
            }),
            Sparse(odepack) => Dense(Lsode {
                mf: dlsode::Generator::InternalFull,
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
        }
    }

    /// # Example
    ///
    /// ```
    /// extern crate approx;
    /// use ndarray::prelude::*;
    ///
    /// let y0 = [1., 0.];
    /// let ts = Array1::linspace(0., 1., 10).to_vec();
    /// let f = |y: &[f64], t: f64| vec![
    ///     998. * y[0] + 1998. * y[1],
    ///     -999. * y[0] - 1999. * y[1],
    ///     ];
    /// let g = |_y: &[f64], _t: f64| array![
    ///     [998., 1998.],
    ///     [-999., -1999.],
    ///     ];
    /// let sol = integrate::BDF::new(f)
    ///     .gen_full_jacobian_by(g)
    ///     .solve(&y0, &ts, 1e-6, 1e-8);
    ///
    /// for (y, t) in sol.iter().zip(ts) {
    ///     approx::assert_abs_diff_eq!(y[0], 2. * (-t).exp() - (-1000. * t).exp(), epsilon = 1e-3);
    ///     approx::assert_abs_diff_eq!(y[1], -(-t).exp() + (-1000. * t).exp(), epsilon = 1e-3);
    /// }
    /// ```
    pub fn gen_full_jacobian_by<G>(self, udf: G) -> Self
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
        use BDF::*;
        match self {
            Dense(odepack) => Dense(Lsode {
                mf: dlsode::Generator::UserSuppliedFull,
                udf: Box::new(g),
                ..odepack
            }),
            Sparse(odepack) => Dense(Lsode {
                mf: dlsode::Generator::UserSuppliedFull,
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
        }
    }

    /// # Example
    ///
    /// ```
    /// extern crate approx;
    ///
    /// let y0 = [1., 1.];
    /// let ts = vec![0., 1.];
    /// let f = |y: &[f64], t: f64| vec![t * y[0], (t + 2.) * y[1]];
    /// let sol = integrate::BDF::new(f)
    ///     .gen_banded_jacobian(0, 0)
    ///     .solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// approx::assert_abs_diff_eq!(sol[1][0], y0[0]*0.5_f64.exp(), epsilon = 1e-3);
    /// ```
    pub fn gen_banded_jacobian(self, ml: usize, mu: usize) -> Self {
        let g = |_neq: *const c_int,
                 _t: *const c_double,
                 _y: *const c_double,
                 _ml: *const c_int,
                 _mu: *const c_int,
                 _pd: *mut c_double,
                 _nr: *const c_int| {};
        use BDF::*;
        match self {
            Dense(odepack) => Dense(Lsode {
                mf: dlsode::Generator::InternalBanded(ml, mu),
                udf: Box::new(g),
                ..odepack
            }),
            Sparse(odepack) => Dense(Lsode {
                mf: dlsode::Generator::InternalBanded(ml, mu),
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
        }
    }

    ///
    ///
    ///     The matrix a is stored in `ab` using the matrix diagonal ordered form::
    ///         ab[mu + i - j][j] == a[i,j]
    ///     Example of `ab` (shape of a is (6,6), `u` =1, `l` =2)::
    ///         *    a01  a12  a23  a34  a45
    ///         a00  a11  a22  a33  a44  a55
    ///         a10  a21  a32  a43  a54   *
    ///         a20  a31  a42  a53   *    *
    pub fn gen_banded_jacobian_by<G>(self, ml: usize, mu: usize, udf: G) -> Self
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
            let ms: usize = (ml + mu + 1) * n;
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, ms),
                    slice::from_raw_parts(y_ptr, n),
                    *t_ptr,
                )
            };
            let mut dy = ArrayViewMut2::<f64>::from_shape((ml + mu + 1, n), dy).unwrap();
            dy.swap_axes(0, 1); // make dy fortran-ordered
            let dy_new = udf(y, t);
            for (mut dy, dy_new) in dy.axis_iter_mut(Axis(0)).zip(dy_new.iter()) {
                dy.assign(dy_new);
            }
        };
        use BDF::*;
        match self {
            Dense(odepack) => Dense(Lsode {
                mf: dlsode::Generator::UserSuppliedBanded(ml, mu),
                udf: Box::new(g),
                ..odepack
            }),
            Sparse(odepack) => Dense(Lsode {
                mf: dlsode::Generator::UserSuppliedBanded(ml, mu),
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
        }
    }

    /// # Example
    ///
    /// ```
    /// extern crate approx;
    ///
    /// const RK1: f64 = 0.1;
    /// const RK2: f64 = 10.;
    /// const RK3: f64 = 50.;
    /// const RK4: f64 = 2.5;
    /// const RK5: f64 = 0.1;
    /// const RK6: f64 = 10.;
    /// const RK7: f64 = 50.;
    /// const RK8: f64 = 2.5;
    /// const RK9:f64 = 50.;
    /// const RK10:f64 = 5.;
    /// const RK11:f64 = 50.;
    /// const RK12:f64 = 50.;
    /// const RK13:f64 = 50.;
    /// const RK14:f64 = 30.;
    /// const RK15:f64 = 100.;
    /// const RK16:f64 = 2.5;
    /// const RK17:f64 = 100.;
    /// const RK18:f64 = 2.5;
    /// const RK19:f64 = 50.;
    /// const RK20:f64 = 50.;
    /// let f = |y: &[f64], t: f64| {
    /// vec![ -RK1 * y[0],
    ///  RK1 * y[0] + RK11 * RK14 * y[3] + RK19 * RK14 * y[4] - RK3 * y[1] * y[2] - RK15 * y[1] *
    ///  y[11] - RK2 * y[1],
    ///  RK2 * y[1] - RK5 * y[2] - RK3 * y[1] * y[2] - RK7 * y[9] * y[2] + RK11 * RK14 * y[3] + RK12
    ///  * RK14 * y[5],
    ///  RK3 * y[1] * y[2] - RK11 * RK14 * y[3] - RK4 * y[3],
    ///  RK15 * y[1] * y[11] - RK19 * RK14 * y[4] - RK16 * y[4],
    ///  RK7 * y[9] * y[2] - RK12 * RK14 * y[5] - RK8 * y[5],
    ///  RK17 * y[9] * y[11] - RK20 * RK14 * y[6] - RK18 * y[6],
    ///  RK9 * y[9] - RK13 * RK14 * y[7] - RK10 * y[7],
    ///  RK4 * y[3] + RK16 * y[4] + RK8 * y[5] + RK18 * y[6],
    ///  RK5 * y[2] + RK12 * RK14 * y[5] + RK20 * RK14 * y[6] + RK13 * RK14 * y[7] - RK7 * y[9] *
    ///  y[2] - RK17 * y[9] * y[11] - RK6 * y[9] - RK9 * y[9],
    ///  RK10 * y[7],
    ///  RK6 * y[9] + RK19 * RK14 * y[4] + RK20 * RK14 * y[6] - RK15 * y[1] * y[11] - RK17 * y[9]
    ///  * y[11],
    ///  ]
    /// };
    /// let y0 = [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,];
    /// let ts = vec![0., 0.1];
    /// let sol = integrate::BDF::new(f)
    ///     .gen_sparse_jacobian(58)
    ///     .solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// for (&result, &gt) in sol[1].iter().zip(
    ///        &[ 9.90050e-01,  6.28228e-03,  3.65313e-03,  7.51934e-07,
    ///           1.12167e-09,  1.18458e-09,  1.77291e-12,  3.26476e-07,
    ///           5.46720e-08,  9.99500e-06,  4.48483e-08,  2.76398e-06,
    ///         ])
    /// {
    ///     approx::assert_abs_diff_eq!(result, gt, epsilon = 1e-6);
    /// }
    /// ```
    pub fn gen_sparse_jacobian(self, max_nnz: usize) -> Self {
        let g = |_neq: *const c_int,
                 _t: *const c_double,
                 _y: *const c_double,
                 _j: *const c_int,
                 _ian: *const c_int,
                 _jan: *const c_int,
                 _pd: *mut c_double| {};
        use BDF::*;
        match self {
            Dense(odepack) => Sparse(Lsodes {
                mf: dlsodes::Generator::InternalSparse(max_nnz),
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
            Sparse(odepack) => Sparse(Lsodes {
                mf: dlsodes::Generator::InternalSparse(max_nnz),
                udf: Box::new(g),
                ..odepack
            }),
        }
    }

    /// # Example
    ///
    /// ```
    /// extern crate approx;
    ///
    /// const RK1: f64 = 0.1;
    /// const RK2: f64 = 10.;
    /// const RK3: f64 = 50.;
    /// const RK4: f64 = 2.5;
    /// const RK5: f64 = 0.1;
    /// const RK6: f64 = 10.;
    /// const RK7: f64 = 50.;
    /// const RK8: f64 = 2.5;
    /// const RK9:f64 = 50.;
    /// const RK10:f64 = 5.;
    /// const RK11:f64 = 50.;
    /// const RK12:f64 = 50.;
    /// const RK13:f64 = 50.;
    /// const RK14:f64 = 30.;
    /// const RK15:f64 = 100.;
    /// const RK16:f64 = 2.5;
    /// const RK17:f64 = 100.;
    /// const RK18:f64 = 2.5;
    /// const RK19:f64 = 50.;
    /// const RK20:f64 = 50.;
    /// let f = |y: &[f64], t: f64| {
    /// vec![ -RK1 * y[0],
    ///  RK1 * y[0] + RK11 * RK14 * y[3] + RK19 * RK14 * y[4] - RK3 * y[1] * y[2] - RK15 * y[1] *
    ///  y[11] - RK2 * y[1],
    ///  RK2 * y[1] - RK5 * y[2] - RK3 * y[1] * y[2] - RK7 * y[9] * y[2] + RK11 * RK14 * y[3] + RK12
    ///  * RK14 * y[5],
    ///  RK3 * y[1] * y[2] - RK11 * RK14 * y[3] - RK4 * y[3],
    ///  RK15 * y[1] * y[11] - RK19 * RK14 * y[4] - RK16 * y[4],
    ///  RK7 * y[9] * y[2] - RK12 * RK14 * y[5] - RK8 * y[5],
    ///  RK17 * y[9] * y[11] - RK20 * RK14 * y[6] - RK18 * y[6],
    ///  RK9 * y[9] - RK13 * RK14 * y[7] - RK10 * y[7],
    ///  RK4 * y[3] + RK16 * y[4] + RK8 * y[5] + RK18 * y[6],
    ///  RK5 * y[2] + RK12 * RK14 * y[5] + RK20 * RK14 * y[6] + RK13 * RK14 * y[7] - RK7 * y[9] *
    ///  y[2] - RK17 * y[9] * y[11] - RK6 * y[9] - RK9 * y[9],
    ///  RK10 * y[7],
    ///  RK6 * y[9] + RK19 * RK14 * y[4] + RK20 * RK14 * y[6] - RK15 * y[1] * y[11] - RK17 * y[9]
    ///  * y[11],
    ///  ]
    /// };
    ///
    /// let jac = |y: &[f64], t: f64, j: usize| {
    ///     let mut jac = vec![0.; y.len()];
    ///     match j {
    ///         0 => {
    ///             jac[0] = -RK1;
    ///             jac[1] = RK1;
    ///         },
    ///         1 => {
    ///             jac[1] = - RK3 * y[2] -RK15 * y[11] - RK2;
    ///             jac[2] = RK2 - RK3 * y[2];
    ///             jac[3] = RK3 * y[2];
    ///             jac[4] = RK15 * y[11];
    ///             jac[11] = -RK15 * y[11];
    ///         },
    ///         2 => {
    ///             jac[1] = - RK3 * y[1];
    ///             jac[2] = - RK5 - RK3 * y[1] - RK7 * y[9];
    ///             jac[3] = RK3 * y[1];
    ///             jac[5] = RK7 * y[9];
    ///             jac[9] = RK5 - RK7 * y[9];
    ///         },
    ///         3 => {
    ///             jac[1] = RK11 * RK14;
    ///             jac[2] = RK11 * RK14;
    ///             jac[3] = - RK11 * RK14 - RK4;
    ///             jac[8] = RK4;
    ///         },
    ///         4 => {
    ///             jac[1] = RK19 * RK14;
    ///             jac[4] = -RK19 * RK14 - RK16;
    ///             jac[8] = RK16;
    ///             jac[11] = RK19 * RK14;
    ///         },
    ///         5 => {
    ///             jac[2] = RK12 * RK14;
    ///             jac[5] = -RK12 * RK14 -RK8;
    ///             jac[8] = RK8;
    ///             jac[9] = RK12 * RK14;
    ///         },
    ///         6 => {
    ///             jac[6] = -RK20 * RK14 - RK18;
    ///             jac[8] = RK18;
    ///             jac[9] = RK20 * RK14;
    ///             jac[11] = RK20 * RK14;
    ///         },
    ///         7 => {
    ///             jac[7] = - RK13 * RK14 - RK10;
    ///             jac[9] = RK13 * RK14;
    ///             jac[10] = RK10;
    ///         },
    ///         8 => {
    ///         },
    ///         9 => {
    ///             jac[2] = -RK7 * y[2];
    ///             jac[5] = RK7 * y[2];
    ///             jac[6] = RK17 * y[11];
    ///             jac[7] = RK9;
    ///             jac[9] = -RK7 * y[2] - RK17 * y[11] - RK6 - RK9;
    ///             jac[11] = RK6 - RK17 * y[11];
    ///         },
    ///         10 => {
    ///         },
    ///         11 => {
    ///             jac[1] = - RK15 * y[1];
    ///             jac[4] = RK15 * y[1];
    ///             jac[6] = RK17 * y[9];
    ///             jac[9] = - RK17 * y[9];
    ///             jac[11] = -RK15 * y[1] - RK17 * y[9];
    ///         },
    ///         _ => {
    ///             panic!();
    ///         }
    ///     }
    ///     jac
    /// };
    /// let y0 = [1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,];
    /// let ts = vec![0., 0.1];
    /// let sol = integrate::BDF::new(f)
    ///     .gen_sparse_jacobian_by(58, jac)
    ///     .solve(&y0, &ts, 1e-6, 1e-6);
    ///
    /// for (&result, &gt) in sol[1].iter().zip(
    ///        &[ 9.90050e-01,  6.28228e-03,  3.65313e-03,  7.51934e-07,
    ///           1.12167e-09,  1.18458e-09,  1.77291e-12,  3.26476e-07,
    ///           5.46720e-08,  9.99500e-06,  4.48483e-08,  2.76398e-06,
    ///         ])
    /// {
    ///     approx::assert_abs_diff_eq!(result, gt, epsilon = 1e-6);
    /// }
    /// ```
    pub fn gen_sparse_jacobian_by(
        self,
        max_nnz: usize,
        udf: impl 'a + Fn(&[f64], f64, usize) -> Vec<f64>,
    ) -> Self {
        let g = move |neq: *const c_int,
                      t_ptr: *const c_double,
                      y_ptr: *const c_double,
                      j_ptr: *const c_int,
                      _ian: *const c_int,
                      _jan: *const c_int,
                      pd_ptr: *mut c_double| {
            let n = unsafe { *neq as usize };
            let j = unsafe { *j_ptr as usize };
            let (pd, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(pd_ptr, n),
                    slice::from_raw_parts(y_ptr, n),
                    *t_ptr,
                )
            };
            let mut pd = ArrayViewMut1::<f64>::from(pd);
            let pd_new = udf(y, t, j - 1);
            pd.assign(&ArrayView1::from(&pd_new));
        };
        use BDF::*;
        match self {
            Dense(odepack) => Sparse(Lsodes {
                mf: dlsodes::Generator::UserSuppliedSparse(max_nnz),
                udf: Box::new(g),
                dydt: odepack.dydt,
                max_steps: odepack.max_steps,
            }),
            Sparse(odepack) => Sparse(Lsodes {
                mf: dlsodes::Generator::UserSuppliedSparse(max_nnz),
                udf: Box::new(g),
                ..odepack
            }),
        }
    }

    pub fn solve(&self, y0: &[f64], t: &[f64], atol: f64, rtol: f64) -> Vec<Vec<f64>> {
        match self {
            Self::Dense(odepack) => odepack.solve(y0, t, atol, rtol),
            Self::Sparse(odepack) => odepack.solve(y0, t, atol, rtol),
        }
    }
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
