use libc::{c_double, c_int};
use libffi::high::Closure4;
use libffi::high::Closure7;
use std::slice;

#[link(name = "gfortran")]
extern "C" {
    /// Call `DLSODES` subroutine from ODEPACK
    ///
    /// For info on passed arguments look inside ODEPACK.
    pub fn dlsodes_(
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
            *const c_int,    // neq
            *const c_double, // t
            *const c_double, // y
            *const c_int,    // j
            *const c_int,    // ian
            *const c_int,    // jan
            *mut c_double,   //pdj
        ),
        mf: &c_int,
    );
}

use super::FLAG;

#[derive(Debug, Clone)]
pub enum Generator {
    InternalSparse(usize),
    UserSuppliedSparse(usize),
}

type JacobianGenerator<'a> = Box<
    dyn 'a
        + Fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *const c_int,
            *const c_int,
            *const c_int,
            *mut c_double,
        ),
>;

pub struct Lsodes<'a> {
    pub dydt: Box<dyn 'a + Fn(&[f64], f64) -> Vec<f64>>,
    pub mf: Generator,
    pub udf: JacobianGenerator<'a>,
}

impl<'a> Lsodes<'a> {
    fn method_flag(&self) -> c_int {
        use Generator::*;
        match self.mf {
            InternalSparse(_) => 222,
            UserSuppliedSparse(_) => 121,
        }
    }

    fn real_work_space(&self, n_eq: usize) -> Vec<c_double> {
        use Generator::*;
        match self.mf {
            InternalSparse(nnz) => {
                const LENRAT: usize = 2;
                let lwm = 2 * nnz + 2 * n_eq + (nnz + 10 * n_eq) / LENRAT;
                vec![0_f64; 20 + 9 * n_eq + lwm]
            }
            UserSuppliedSparse(nnz) => {
                const LENRAT: usize = 2;
                let lwm = 2 * nnz + 2 * n_eq + (nnz + 9 * n_eq) / LENRAT;
                vec![0_f64; 20 + 9 * n_eq + lwm]
            }
        }
    }

    fn integer_work_space(&self, n_eq: usize) -> Vec<c_int> {
        use Generator::*;
        match self.mf {
            InternalSparse(_) | UserSuppliedSparse(_) => vec![0_i32; 30],
        }
    }

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

        let mut y: Vec<f64> = y0.to_vec();
        let n = y0.len();
        let mut t0 = t[0];

        let itol = 1;
        let itask = 1;
        let iopt = 0;
        let mut istate = 1;
        let mf = self.method_flag();

        let mut rwork = self.real_work_space(n);
        let mut iwork = self.integer_work_space(n);
        let lrw = rwork.len();
        let liw = iwork.len();

        let jacobian_closure = Closure7::new(&self.udf);
        let call_jacobian = jacobian_closure.code_ptr();

        let mut result = Vec::with_capacity(t.len());

        let _lock = FLAG.lock().unwrap();
        for &tout in t.iter() {
            unsafe {
                dlsodes_(
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
