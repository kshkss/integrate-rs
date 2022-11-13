use lazy_static::lazy_static;
use libc::{c_double, c_int};
use libffi::high::Closure4;
use libffi::high::Closure7;
use std::slice;
use std::sync::Mutex;

use super::low;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

pub fn dlsode<'a, 'b>(
    dydt: &(impl 'a + Fn(&'b [f64], f64, &'b mut [f64]) + ?Sized),
    jac: &(impl 'a + Fn(&'b [f64], f64, &'b mut [f64]) + ?Sized),
    y: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let f =
        |n: *const c_int, t_ptr: *const c_double, y_ptr: *mut c_double, dy_ptr: *mut c_double| {
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, *n as usize),
                    slice::from_raw_parts(y_ptr, *n as usize),
                    *t_ptr,
                )
            };
            dydt(y, t, dy);
        };
    let closure = Closure4::new(&f);
    let call = closure.code_ptr();

    let dfdy = |n: *const c_int,
                t_ptr: *const c_double,
                y_ptr: *const c_double,
                _ml_ptr: *const c_int,
                _mu_ptr: *const c_int,
                pd_ptr: *mut c_double,
                nrow_pd: *const c_int| {
        let (pd, y, t) = unsafe {
            let y = slice::from_raw_parts(y_ptr, *n as usize);
            let t = *t_ptr;
            let nrow = *nrow_pd as usize;
            let pd = slice::from_raw_parts_mut(pd_ptr, y.len() * nrow);
            (pd, y, t)
        };
        jac(y, t, pd);
    };
    let closure2 = Closure7::new(&dfdy);
    let call2 = closure2.code_ptr();

    let n = y.len() as i32;

    let itol = 1;
    let itask = 1;
    let iopt = 1;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsode_(
            *call,
            &n,
            y.as_mut_ptr(),
            &mut t0,
            &t1,
            &itol,
            &rtol,
            &atol,
            &itask,
            &mut istate,
            &iopt,
            rwork.as_mut_ptr(),
            &lrw,
            iwork.as_mut_ptr(),
            &liw,
            *call2,
            &mf,
        );
    }
}

pub fn dlsodes<'a, 'b>(
    dydt: &(impl 'a + Fn(&'b [f64], f64, &'b mut [f64]) + ?Sized),
    jac: &(impl 'a + Fn(&'b [f64], f64, usize, &'b mut [f64]) + ?Sized),
    y0: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let f =
        |n: *const c_int, t_ptr: *const c_double, y_ptr: *mut c_double, dy_ptr: *mut c_double| {
            let (dy, y, t) = unsafe {
                (
                    slice::from_raw_parts_mut(dy_ptr, *n as usize),
                    slice::from_raw_parts(y_ptr, *n as usize),
                    *t_ptr,
                )
            };
            dydt(y, t, dy);
        };
    let closure = Closure4::new(&f);
    let call = closure.code_ptr();

    let jac = |neq: *const c_int,
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
        jac(y, t, j - 1, pd);
    };
    let closure2 = Closure7::new(&jac);
    let call2 = closure2.code_ptr();

    let n = y0.len() as i32;

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsodes_(
            *call,
            &n,
            y0.as_mut_ptr(),
            &mut t0,
            &t1,
            &itol,
            &rtol,
            &atol,
            &itask,
            &mut istate,
            &iopt,
            rwork.as_mut_ptr(),
            &lrw,
            iwork.as_mut_ptr(),
            &liw,
            *call2,
            &mf,
        );
    }
}
