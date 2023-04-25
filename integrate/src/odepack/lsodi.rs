//! Solvers that treat ODE and DAE systems in the linearly implicit form *A(t,y) dy/dt = f(t,y)*
//!
//! A is a square matrix ans allowed to be
//! singular, in which case the system is a differential-algebraic
//! equation (DAE) system.  In that case, the user must be very careful
//! to supply a well-posed problem with consistent initial conditions.

use super::mid;
use super::Control;
use ndarray::prelude::*;

impl<F, ADDA, JAC> mid::LsodiCallback for (&F, &ADDA, &JAC)
where
    F: Fn(f64, &[f64], &[f64], &mut [f64]),
    ADDA: Fn(f64, &[f64], ArrayViewMut2<f64>),
    JAC: Fn(f64, &[f64], &[f64], ArrayViewMut2<f64>),
{
    fn residual(&self, t: f64, y: &[f64], s: &[f64], res: &mut [f64]) {
        (self.0)(t, y, s, res);
    }

    fn adda(&self, t: f64, y: &[f64], pd: &mut [f64]) {
        let mut pd =
            ArrayViewMut2::from_shape((y.len(), pd.len() / y.len()), pd).expect("size mismatch");
        pd.swap_axes(0, 1); // make pd fortran-ordered
        (self.1)(t, y, pd);
    }

    fn jac(&self, t: f64, y: &[f64], s: &[f64], pd: &mut [f64]) {
        let mut pd =
            ArrayViewMut2::from_shape((y.len(), pd.len() / y.len()), pd).expect("size mismatch");
        pd.swap_axes(0, 1); // make pd fortran-ordered
        (self.2)(t, y, s, pd);
    }
}

impl<F, ADDA, JAC> mid::LsodisCallback for (&F, &ADDA, &JAC)
where
    F: Fn(f64, &[f64], &[f64], &mut [f64]),
    ADDA: Fn(f64, &[f64], usize, &mut [f64]),
    JAC: Fn(f64, &[f64], &[f64], usize, &mut [f64]),
{
    fn residual(&self, t: f64, y: &[f64], s: &[f64], res: &mut [f64]) {
        (self.0)(t, y, s, res);
    }

    fn adda(&self, t: f64, y: &[f64], j: usize, pd: &mut [f64]) {
        (self.1)(t, y, j, pd);
    }

    fn jac(&self, t: f64, y: &[f64], s: &[f64], j: usize, pd: &mut [f64]) {
        (self.2)(t, y, s, j, pd);
    }
}

/// Solver for systems whose jacobian matrix and A are full matrix.
///
/// # Example
///
/// ```
/// use integrate::odepack::*;
/// use integrate::odepack::lsodi::*;
/// use ndarray::prelude::*;
///
/// let res = |t: f64, y: &[f64], s: &[f64], r: &mut [f64]| {
///     r[0] = -0.4 * y[0] + 1e4 * y[1] * y[2] - s[0];
///     r[1] = 0.4 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1] - s[1];
///     r[2] = y[0] + y[1] + y[2] - 1.
/// };
/// let adda = |t: f64, y: &[f64], mut pd: ArrayViewMut2<f64>| {
///     pd[[0,0]] = pd[[0,0]] + 1.;
///     pd[[1,1]] = pd[[1,1]] + 1.;
/// };
///
/// let y = vec![1., 0., 0.,];
/// let dy = vec![-0.04, 0.04, 0.];
/// let t = vec![0., 4e-1, 4e0, 4e1, 4e2, 4e3, 4e4, 4e5, 4e6, 4e7, 4e8, 4e9, 4e10];
///
/// let results = Lsodi::new(&res, &adda, Control::default()).run(&t, &y, &dy);
/// ```
pub struct LsodiFullJacobian<'a> {
    option: Control,
    residual: &'a (dyn 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
    adda: &'a (dyn 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
    jac: Option<&'a (dyn 'a + Fn(f64, &[f64], &[f64], ArrayViewMut2<f64>))>,
}

impl<'a> LsodiFullJacobian<'a> {
    /// Constructor.
    ///
    /// # Arguments
    ///
    /// * `residual` - Closure represents residual function, f(t,y) - A(t,y) s
    /// * `adda` - Closure represents A(t,y)
    /// * `jac` - Closure represents jacobian matrix of **residual function**
    pub fn new(
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
        jac: &'a (impl 'a + Fn(f64, &[f64], &[f64], ArrayViewMut2<f64>)),
        option: Control,
    ) -> Self {
        Self {
            residual,
            adda,
            jac: Some(jac),
            option,
        }
    }
    fn is_gen_internally(&self) -> bool {
        self.jac.is_none()
    }

    fn method_flag(&self) -> i32 {
        if self.is_gen_internally() {
            22
        } else {
            21
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<f64> {
        let rwork = vec![0.; 22 + 9 * n_eq + n_eq * n_eq];
        rwork
    }

    fn integer_work_space(&self, n_eq: usize, control: &Control) -> Vec<i32> {
        let mut iwork = vec![0; 20 + n_eq];
        iwork[5] = control.max_steps as i32;
        iwork
    }

    /// Solve the system.
    ///
    /// # Arguments
    ///
    /// * `t` - Values of independent variable. Two values are needed at least.
    /// The first one is initial value and the last one is the value at the end.
    /// * `y0` - Initial values of y
    /// * `dy0` - Initial values of dy/dt
    pub fn run(&self, t: &[f64], y0: &[f64], dy0: &[f64]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut ys = Vec::with_capacity(t.len());
        let mut dys = Vec::with_capacity(t.len());
        let mut y = y0.to_vec();
        let mut dy = dy0.to_vec();
        let mut t0 = t[0];
        for &t1 in &t[1..] {
            ys.push(y.to_owned());
            dys.push(dy.to_owned());
            self.step((t0, t1), &mut y, &mut dy);
            t0 = t1;
        }
        ys.push(y);
        dys.push(dy);
        (ys, dys)
    }

    fn step(&self, t: (f64, f64), y: &mut [f64], dy: &mut [f64]) {
        let mf = self.method_flag();
        let mut rwork = self.real_work_space(y.len());
        let mut iwork = self.integer_work_space(y.len(), &self.option);

        if let Some(ref jac) = self.jac {
            mid::dlsodi(
                &(&self.residual, &self.adda, &jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            )
        } else {
            let jac = |_t: f64, _y: &[f64], _s: &[f64], _pd: ArrayViewMut2<f64>| {
                unreachable!();
            };
            mid::dlsodi(
                &(&self.residual, &self.adda, &jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            )
        }
    }
}

/// Solver for systems whose jacobian matrix and A are banded matrix.
pub struct LsodiBandedJacobian<'a> {
    option: Control,
    ml: usize,
    mu: usize,
    residual: &'a (dyn 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
    adda: &'a (dyn 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
    jac: Option<&'a (dyn 'a + Fn(f64, &[f64], &[f64], ArrayViewMut2<f64>))>,
}

impl<'a> LsodiBandedJacobian<'a> {
    /// Constructor.
    ///
    /// # Arguments
    ///
    /// * `ml` - Width of the lower parts of jacobian matrix and A, excluding the main diagonal.
    /// * `mu` - Width of the upper parts of jacobian matrix and A, excluding the main diagonal.
    /// * `residual` - Closure represents residual function, f(t,y) - A(t,y) s
    /// * `adda` - Closure represents A(t,y)
    /// * `jac` - Closure represents jacobian matrix of **residual function**
    pub fn new(
        ml: usize,
        mu: usize,
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
        jac: &'a (impl 'a + Fn(f64, &[f64], &[f64], ArrayViewMut2<f64>)),
        option: Control,
    ) -> Self {
        Self {
            ml,
            mu,
            residual,
            adda,
            jac: Some(jac),
            option,
        }
    }

    fn is_gen_internally(&self) -> bool {
        self.jac.is_none()
    }

    fn method_flag(&self) -> i32 {
        if self.is_gen_internally() {
            25
        } else {
            24
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<f64> {
        vec![0.; 22 + 10 * n_eq + (2 * self.ml + self.mu) * n_eq]
    }

    fn integer_work_space(&self, n_eq: usize, control: &Control) -> Vec<i32> {
        let mut iwork = vec![0; 20 + n_eq];
        iwork[0] = self.ml as i32;
        iwork[1] = self.mu as i32;
        iwork[5] = control.max_steps as i32;
        iwork
    }

    /// Solve the system.
    ///
    /// # Arguments
    ///
    /// * `t` - Values of independent variable. Two values are needed at least.
    /// The first one is initial value and the last one is the value at the end.
    /// * `y0` - Initial values of y
    /// * `dy0` - Initial values of dy/dt
    pub fn run(&self, t: &[f64], y0: &[f64], dy0: &[f64]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut ys = Vec::with_capacity(t.len());
        let mut dys = Vec::with_capacity(t.len());
        let mut y = y0.to_vec();
        let mut dy = dy0.to_vec();
        let mut t0 = t[0];
        for &t1 in &t[1..] {
            ys.push(y.to_owned());
            dys.push(dy.to_owned());
            self.step((t0, t1), &mut y, &mut dy);
            t0 = t1;
        }
        ys.push(y);
        dys.push(dy);
        (ys, dys)
    }

    fn step(&self, t: (f64, f64), y: &mut [f64], dy: &mut [f64]) {
        let mf = self.method_flag();
        let mut rwork = self.real_work_space(y.len());
        let mut iwork = self.integer_work_space(y.len(), &self.option);

        if let Some(ref jac) = self.jac {
            mid::dlsodi(
                &(&self.residual, &self.adda, &jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            )
        } else {
            let jac = |_t: f64, _y: &[f64], _s: &[f64], _pd: ArrayViewMut2<f64>| {
                unreachable!();
            };
            mid::dlsodi(
                &(&self.residual, &self.adda, &jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            )
        }
    }
}

/// Solver for systems whose jacobian matrix and A are sparse matrix.
pub struct LsodiSparseJacobian<'a> {
    option: Control,
    max_nnz: usize,
    residual: &'a (dyn 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
    adda: &'a (dyn 'a + Fn(f64, &[f64], usize, &mut [f64])),
    jac: Option<&'a (dyn 'a + Fn(f64, &[f64], &[f64], usize, &mut [f64]))>,
}

impl<'a> LsodiSparseJacobian<'a> {
    /// Constructor.
    ///
    /// # Arguments
    ///
    /// * `max_nnz` - Maximum number of non-zero elements in jacobian matrix and A.
    /// * `residual` - Closure represents residual function, f(t,y) - A(t,y) s
    /// * `adda` - Closure represents A(t,y)
    /// * `jac` - Closure represents jacobian matrix of **residual function**
    pub fn new(
        max_nnz: usize,
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], usize, &mut [f64])),
        jac: &'a (impl 'a + Fn(f64, &[f64], &[f64], usize, &mut [f64])),
        option: Control,
    ) -> Self {
        Self {
            max_nnz,
            residual,
            adda,
            jac: Some(jac),
            option,
        }
    }

    fn is_gen_internally(&self) -> bool {
        self.jac.is_none()
    }

    fn method_flag(&self) -> i32 {
        if self.is_gen_internally() {
            222
        } else {
            121
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<f64> {
        use std::mem::size_of;
        let lenrat = size_of::<f64>() / size_of::<i32>();

        if self.is_gen_internally() {
            vec![0.; 2 * self.max_nnz + 2 * n_eq + (self.max_nnz + 10 * n_eq) / lenrat]
        } else {
            vec![0.; 2 * self.max_nnz + 2 * n_eq + (self.max_nnz + 9 * n_eq) / lenrat]
        }
    }

    fn integer_work_space(&self, _n_eq: usize, control: &Control) -> Vec<i32> {
        let mut iwork = vec![0; 30];
        iwork[5] = control.max_steps as i32;
        iwork[18] = self.max_nnz as i32;
        iwork
    }

    /// Solve the system.
    ///
    /// # Arguments
    ///
    /// * `t` - Values of independent variable. Two values are needed at least.
    /// The first one is initial value and the last one is the value at the end.
    /// * `y0` - Initial values of y
    /// * `dy0` - Initial values of dy/dt
    pub fn run(&self, t: &[f64], y0: &[f64], dy0: &[f64]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        let mut ys = Vec::with_capacity(t.len());
        let mut dys = Vec::with_capacity(t.len());
        let mut y = y0.to_vec();
        let mut dy = dy0.to_vec();
        let mut t0 = t[0];
        for &t1 in &t[1..] {
            ys.push(y.to_owned());
            dys.push(dy.to_owned());
            self.step((t0, t1), &mut y, &mut dy);
            t0 = t1;
        }
        ys.push(y);
        dys.push(dy);
        (ys, dys)
    }

    fn step(&self, t: (f64, f64), y: &mut [f64], dy: &mut [f64]) {
        let mf = self.method_flag();
        let mut rwork = self.real_work_space(y.len());
        let mut iwork = self.integer_work_space(y.len(), &self.option);

        if let Some(ref jac) = self.jac {
            mid::dlsodis(
                &(&self.residual, &self.adda, jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            );
        } else {
            let jac = |_t: f64, _y: &[f64], _s: &[f64], _j: usize, _pd: &mut [f64]| {
                unreachable!();
            };
            mid::dlsodis(
                &(&self.residual, &self.adda, &jac),
                y,
                dy,
                t.0,
                t.1,
                self.option.rtol,
                self.option.atol,
                &mut rwork,
                &mut iwork,
                mf,
            );
        }
    }
}

pub struct Lsodi;

impl Lsodi {
    /// ```
    /// use integrate::odepack::*;
    /// use integrate::odepack::lsodi::*;
    /// use ndarray::prelude::*;
    ///
    /// let res = |t: f64, y: &[f64], s: &[f64], r: &mut [f64]| {
    ///     r[0] = -0.4 * y[0] + 1e4 * y[1] * y[2] - s[0];
    ///     r[1] = 0.4 * y[0] - 1e4 * y[1] * y[2] - 3e7 * y[1] * y[1] - s[1];
    ///     r[2] = y[0] + y[1] + y[2] - 1.
    /// };
    /// let adda = |t: f64, y: &[f64], mut pd: ArrayViewMut2<f64>| {
    ///     pd[[0,0]] = pd[[0,0]] + 1.;
    ///     pd[[1,1]] = pd[[1,1]] + 1.;
    /// };
    ///
    /// let y = vec![1., 0., 0.,];
    /// let dy = vec![-0.04, 0.04, 0.];
    /// let t = vec![0., 4e-1, 4e0, 4e1, 4e2, 4e3, 4e4, 4e5, 4e6, 4e7, 4e8, 4e9, 4e10];
    ///
    /// let results = Lsodi::new(&res, &adda, Control::default()).run(&t, &y, &dy);
    /// ```
    pub fn new<'a>(
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
        control: Control,
    ) -> LsodiFullJacobian<'a> {
        LsodiFullJacobian::<'a> {
            residual,
            option: control,
            adda,
            jac: None,
        }
    }

    pub fn new_with_banded_jacobian<'a>(
        ml: usize,
        mu: usize,
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], ArrayViewMut2<f64>)),
        control: Control,
    ) -> LsodiBandedJacobian<'a> {
        LsodiBandedJacobian::<'a> {
            residual,
            option: control,
            ml,
            mu,
            adda,
            jac: None,
        }
    }

    pub fn new_with_sparse_jacobian<'a>(
        max_nnz: usize,
        residual: &'a (impl 'a + Fn(f64, &[f64], &[f64], &mut [f64])),
        adda: &'a (impl 'a + Fn(f64, &[f64], usize, &mut [f64])),
        control: Control,
    ) -> LsodiSparseJacobian<'a> {
        LsodiSparseJacobian::<'a> {
            residual,
            option: control,
            max_nnz,
            adda,
            jac: None,
        }
    }
}
