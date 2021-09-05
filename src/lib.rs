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
    static ref flag: Mutex<()> = Mutex::<()>::new(());
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
/// let f = |y: &[f64], t: &f64| {
///     let mut dy = vec![0.0];
///     dy[0] = *t * y[0];
///     dy
///     };
/// let sol = lsode::solve_ode(f, &y0, ts, 1e-6, 1e-6);
///
/// assert!((sol[1][0] - y0[0]*0.5_f64.exp()).abs() < 1e-3, "error too large");
/// ```
pub fn solve_ode<F>(rhs: F, y0: &[f64], t_dense: Vec<f64>, atol: f64, rtol: f64) -> Vec<Vec<f64>>
where
    F: Fn(&[f64], &f64) -> Vec<f64>,
{
    let f =
        |n: *const c_int, t_ptr: *const c_double, y_ptr: *mut c_double, dy_ptr: *mut c_double| {
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, *n as usize),
                    slice::from_raw_parts(y_ptr, *n as usize),
                    *t_ptr,
                )
            };
            let dy_new = rhs(y, &t);
            for (i, deriv) in dy_new.iter().enumerate() {
                dy[i] = *deriv
            }
        };
    let closure = Closure4::new(&f);
    let call = closure.code_ptr();

    let mut y: Vec<f64> = y0.to_vec();
    let n = y0.len();
    let mut t = t_dense[0];

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;
    let mf = 22;

    let lrw = 22 + 9 * n + n * n;
    let liw = 20 + n;
    let mut rwork: Vec<f64> = (0..lrw).map(|_i| 0.0 as c_double).collect();
    let mut iwork: Vec<i32> = (0..liw).map(|_i| 0 as c_int).collect();

    let mut result = Vec::new();

    let _lock = flag.lock().unwrap();
    for tout in t_dense {
        unsafe {
            dlsode_(
                *call,
                &(n as i32),
                y.as_mut_ptr(),
                &mut t,
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
                fake_jacobian,
                &mf,
            );
        }

        result.push(y.clone());
    }
    result
}

use ndarray::prelude::*;

pub fn solve_ode_with_jacobian<F, G>(
    rhs: F,
    jacobian: G,
    y0: &[f64],
    t_dense: Vec<f64>,
    atol: f64,
    rtol: f64,
) -> Vec<Vec<f64>>
where
    F: Fn(&[f64], &f64) -> Vec<f64>,
    G: Fn(&[f64], &f64) -> Array2<f64>,
{
    let f =
        |n: *const c_int, t_ptr: *const c_double, y_ptr: *mut c_double, dy_ptr: *mut c_double| {
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, *n as usize),
                    slice::from_raw_parts(y_ptr, *n as usize),
                    *t_ptr,
                )
            };
            let dy_new = rhs(y, &t);
            for (i, deriv) in dy_new.iter().enumerate() {
                dy[i] = *deriv
            }
        };
    let closure = Closure4::new(&f);
    let call = closure.code_ptr();

    let g = |n: *const c_int,
             t_ptr: *const c_double,
             y_ptr: *const c_double,
             _ml: *const c_int,
             _mu: *const c_int,
             dy_ptr: *mut c_double,
             _nrow: *const c_int| {
        let (dy, y, t) = unsafe {
            (
                slice::from_raw_parts_mut(dy_ptr, *n as usize * *n as usize),
                slice::from_raw_parts(y_ptr, *n as usize),
                *t_ptr,
            )
        };
        let mut dy = ArrayViewMut1::<f64>::from(dy);
        let dy_new = jacobian(y, &t);
        for (dest, jac) in dy.iter_mut().zip(dy_new.t().iter()) {
            *dest = *jac;
        }
    };

    let jacobian_closure = Closure7::new(&g);
    let call_jacobian = jacobian_closure.code_ptr();

    let mut y: Vec<f64> = y0.to_vec();
    let n = y0.len();
    let mut t = t_dense[0];

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;
    let mf = 21; // Stiff (BDF) method, user-supplied full Jacobian.

    let lrw = 22 + 9 * n + n * n;
    let liw = 20 + n;
    let mut rwork: Vec<f64> = (0..lrw).map(|_i| 0.0 as c_double).collect();
    let mut iwork: Vec<i32> = (0..liw).map(|_i| 0 as c_int).collect();

    let mut result = Vec::new();

    for tout in t_dense {
        unsafe {
            dlsode_(
                *call,
                &(n as i32),
                y.as_mut_ptr(),
                &mut t,
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
