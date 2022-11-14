use super::mid;
use ndarray::prelude::*;
use std::rc::Rc;

enum Jac<'a> {
    NoJac,
    UserSuppliedFull {
        jac: Box<dyn 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)>,
    },
    InternalFull,
    UserSuppliedBanded {
        ml: usize,
        mu: usize,
        jac: Box<dyn 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)>,
    },
    InternalBanded {
        ml: usize,
        mu: usize,
    },
    UserSuppliedSparse {
        max_nnz: usize,
        jac: Box<dyn 'a + Fn(f64, &[f64], usize, &mut [f64])>,
    },
    InternalSparse {
        max_nnz: usize,
    },
}

impl<'a> Jac<'a> {
    fn method_flag(&self) -> i32 {
        use Jac::*;
        match self {
            NoJac => 10,
            InternalFull => 22,
            InternalBanded { .. } => 25,
            InternalSparse { .. } => 222,
            UserSuppliedFull { .. } => 21,
            UserSuppliedBanded { .. } => 24,
            UserSuppliedSparse { .. } => 121,
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<f64> {
        use std::mem::size_of;
        use Jac::*;
        match self {
            NoJac => vec![0_f64; 20 + 16 * n_eq],
            InternalFull | UserSuppliedFull { .. } => {
                vec![0_f64; 22 + 9 * n_eq + n_eq * n_eq]
            }
            InternalBanded { ml, mu } | UserSuppliedBanded { ml, mu, .. } => {
                vec![0_f64; 22 + 10 * n_eq + (2 * ml + mu) * n_eq]
            }
            InternalSparse { max_nnz } => {
                const LENRAT: usize = size_of::<f64>() / size_of::<i32>();
                let lwm = 2 * max_nnz + 2 * n_eq + (max_nnz + 10 * n_eq) / LENRAT;
                vec![0_f64; 20 + 9 * n_eq + lwm]
            }
            UserSuppliedSparse { max_nnz, .. } => {
                const LENRAT: usize = size_of::<f64>() / size_of::<i32>();
                let lwm = 2 * max_nnz + 2 * n_eq + (max_nnz + 9 * n_eq) / LENRAT;
                vec![0_f64; 20 + 9 * n_eq + lwm]
            }
        }
    }

    fn integer_work_space(&self, n_eq: usize, option: &Control) -> Vec<i32> {
        use Jac::*;
        match self {
            NoJac => {
                let mut iwork = vec![0_i32; 20];
                iwork[5] = option.max_steps as i32;
                iwork
            }
            InternalFull | UserSuppliedFull { .. } => {
                let mut iwork = vec![0_i32; 22 + n_eq];
                iwork[5] = option.max_steps as i32;
                iwork
            }
            InternalBanded { ml, mu } | UserSuppliedBanded { ml, mu, .. } => {
                let mut iwork = vec![0_i32; 22 + n_eq];
                iwork[0] = *ml as i32;
                iwork[1] = *mu as i32;
                iwork[5] = option.max_steps as i32;
                iwork
            }
            InternalSparse { max_nnz } | UserSuppliedSparse { max_nnz, .. } => {
                let mut iwork = vec![0_i32; 30];
                iwork[5] = option.max_steps as i32;
                iwork[17] = *max_nnz as i32;
                iwork
            }
        }
    }
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

pub struct Adam<'a> {
    f: Rc<dyn 'a + Fn(f64, &[f64], &mut [f64])>,
    option: Control,
}

impl<'a> Adam<'a> {
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
    pub fn new(f: impl 'a + Fn(f64, &[f64], &mut [f64]), option: Control) -> Self {
        Self {
            f: Rc::new(f),
            option,
        }
    }

    pub fn solve(&self, y0: &[f64], t: (f64, f64)) -> Vec<f64> {
        let lsode = BDF {
            f: Rc::clone(&self.f),
            jac: Jac::NoJac,
            option: self.option.clone(),
        };
        lsode.solve(y0, t)
    }
}

pub struct BDF<'a> {
    f: Rc<dyn 'a + Fn(f64, &[f64], &mut [f64])>,
    jac: Jac<'a>,
    option: Control,
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
    pub fn new(f: impl 'a + Fn(f64, &[f64], &mut [f64]), option: Control) -> Self {
        Self {
            f: Rc::new(f),
            jac: Jac::InternalFull,
            option,
        }
    }

    pub fn gen_full_jacobian(self) -> Self {
        Self {
            jac: Jac::InternalFull,
            ..self
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
    pub fn gen_full_jacobian_by(self, udf: impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)) -> Self {
        Self {
            jac: Jac::UserSuppliedFull { jac: Box::new(udf) },
            ..self
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
        Self {
            jac: Jac::InternalBanded { ml, mu },
            ..self
        }
    }

    /// The matrix a is stored in `ab` using the matrix diagonal ordered form::
    ///   ab[mu + i - j][j] == a[i,j]
    /// Example of `ab` (shape of a is (6,6), `u` =1, `l` =2)::
    ///  *   a01  a12  a23  a34  a45
    /// a00  a11  a22  a33  a44  a55
    /// a10  a21  a32  a43  a54   *
    /// a20  a31  a42  a53   *    *
    pub fn gen_banded_jacobian_by(
        self,
        ml: usize,
        mu: usize,
        udf: impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>),
    ) -> Self {
        Self {
            jac: Jac::UserSuppliedBanded {
                ml,
                mu,
                jac: Box::new(udf),
            },
            ..self
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
        Self {
            jac: Jac::InternalSparse { max_nnz },
            ..self
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
        udf: impl 'a + Fn(f64, &[f64], usize, &mut [f64]),
    ) -> Self {
        Self {
            jac: Jac::UserSuppliedSparse {
                max_nnz,
                jac: Box::new(udf),
            },
            ..self
        }
    }

    pub fn solve(&self, y0: &[f64], t: (f64, f64)) -> Vec<f64> {
        let mut y = y0.to_vec();
        let mf = self.jac.method_flag();
        let mut rwork = self.jac.real_work_space(y0.len());
        let mut iwork = self.jac.integer_work_space(y0.len(), &self.option);

        match self.jac {
            Jac::NoJac | Jac::InternalFull | Jac::InternalBanded { .. } => {
                let jac = |_t, _y, _pd| {
                    unreachable!();
                };
                mid::dlsode(
                    self.f.as_ref(),
                    &jac,
                    &mut y,
                    t.0,
                    t.1,
                    self.option.rtol,
                    self.option.atol,
                    &mut rwork,
                    &mut iwork,
                    mf,
                )
            }

            Jac::InternalSparse { .. } => {
                let jac = |_t, _y, _j, _pd| {
                    unreachable!();
                };
                mid::dlsodes(
                    self.f.as_ref(),
                    &jac,
                    &mut y,
                    t.0,
                    t.1,
                    self.option.rtol,
                    self.option.atol,
                    &mut rwork,
                    &mut iwork,
                    mf,
                )
            }

            Jac::UserSuppliedFull { ref jac } | Jac::UserSuppliedBanded { ref jac, .. } => {
                let jac = |t: f64, y: &[f64], pd: &mut [f64]| {
                    assert_eq!(0, pd.len() % y.len());
                    let mut pd = ArrayViewMut2::from_shape((y.len(), pd.len() / y.len()), pd)
                        .expect("size mismatch");
                    pd.swap_axes(0, 1); // make pd fortran-ordered
                    jac(t, y, pd);
                };
                mid::dlsode(
                    self.f.as_ref(),
                    &jac,
                    &mut y,
                    t.0,
                    t.1,
                    self.option.rtol,
                    self.option.atol,
                    &mut rwork,
                    &mut iwork,
                    mf,
                )
            }

            Jac::UserSuppliedSparse { ref jac, .. } => mid::dlsodes(
                self.f.as_ref(),
                &jac,
                &mut y,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            ),
        }

        y
    }
}
