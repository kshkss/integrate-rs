use libc::{c_double, c_int};
use once_cell::sync::Lazy;
use std::slice;
use std::sync::Mutex;

use odepack_sys as low;

static FLAG: Lazy<Mutex<()>> = Lazy::new(|| Mutex::new(()));

pub trait LsodeCallback {
    fn f(&self, t: f64, y: &[f64], dy: &mut [f64]);
    fn jac(&self, t: f64, y: &[f64], pd: &mut [f64]);
}

extern "C" fn lsode_f<T: LsodeCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *mut c_double,
    dy_ptr: *mut c_double,
) {
    let (dy, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(dy_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.f(t, y, dy);
}

extern "C" fn lsode_jac<T: LsodeCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    _ml_ptr: *const c_int,
    _mu_ptr: *const c_int,
    pd_ptr: *mut c_double,
    nrow_pd: *const c_int,
) {
    let (pd, y, t) = unsafe {
        let y = slice::from_raw_parts(y_ptr, *n as usize);
        let t = *t_ptr;
        let nrow = *nrow_pd as usize;
        let pd = slice::from_raw_parts_mut(pd_ptr, y.len() * nrow);
        (pd, y, t)
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.jac(t, y, pd);
}

pub fn dlsode<'a, 'b, T: LsodeCallback>(
    callback: &T,
    y: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let n = {
        let mut n = [y.len() as i32; 3];
        let ptr = unsafe { std::mem::transmute::<&T, [i32; 2]>(callback) };
        n[1] = ptr[0];
        n[2] = ptr[1];
        n
    };

    let itol = 1;
    let itask = 1;
    let iopt = 1;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsode_(
            lsode_f::<T>,
            n.as_ptr(),
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
            lsode_jac::<T>,
            &mf,
        );
    }
}

pub trait LsodesCallback {
    fn f(&self, t: f64, y: &[f64], dy: &mut [f64]);
    fn jac(&self, t: f64, y: &[f64], j: usize, pd: &mut [f64]);
}

extern "C" fn lsodes_f<T: LsodesCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *mut c_double,
    dy_ptr: *mut c_double,
) {
    let (dy, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(dy_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.f(t, y, dy);
}

extern "C" fn lsodes_jac<T: LsodesCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    j_ptr: *const c_int,
    _ian: *const c_int,
    _jan: *const c_int,
    pd_ptr: *mut c_double,
) {
    let j = unsafe { *j_ptr as usize };
    let (pd, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(pd_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.jac(t, y, j - 1, pd);
}

pub fn dlsodes<'a, 'b, T: LsodesCallback>(
    callback: &T,
    y0: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let n = {
        let mut n = [y0.len() as i32; 3];
        let ptr = unsafe { std::mem::transmute::<&T, [i32; 2]>(callback) };
        n[1] = ptr[0];
        n[2] = ptr[1];
        n
    };

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsodes_(
            lsodes_f::<T>,
            n.as_ptr(),
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
            lsodes_jac::<T>,
            &mf,
        );
    }
}

pub trait LsodiCallback {
    fn residual(&self, t: f64, y: &[f64], s: &[f64], res: &mut [f64]);
    fn jac(&self, t: f64, y: &[f64], s: &[f64], pd: &mut [f64]);
    fn adda(&self, t: f64, y: &[f64], pd: &mut [f64]);
}

extern "C" fn lsodi_res<T: LsodiCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    s_ptr: *const c_double,
    r_ptr: *mut c_double,
    _ires: *const c_int,
) {
    let (r, y, s, t) = unsafe {
        (
            slice::from_raw_parts_mut(r_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            slice::from_raw_parts(s_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.residual(t, y, s, r);
}

extern "C" fn lsodi_adda<T: LsodiCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    _ml: *const c_int,
    _mu: *const c_int,
    p_ptr: *mut c_double,
    nrow_p: *const c_int,
) {
    let (p, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(p_ptr, *n as usize * *nrow_p as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.adda(t, y, p);
}

extern "C" fn lsodi_jac<T: LsodiCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    s_ptr: *const c_double,
    _ml: *const c_int,
    _mu: *const c_int,
    p_ptr: *mut c_double,
    nrow_p: *const c_int,
) {
    let (p, y, s, t) = unsafe {
        (
            slice::from_raw_parts_mut(p_ptr, *n as usize * *nrow_p as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            slice::from_raw_parts(s_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.jac(t, y, s, p);
}

pub fn dlsodi<'a, 'b, T: LsodiCallback>(
    callback: &T,
    y0: &mut [f64],
    dydt0: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let n = {
        let mut n = [y0.len() as i32; 3];
        let ptr = unsafe { std::mem::transmute::<&T, [i32; 2]>(callback) };
        n[1] = ptr[0];
        n[2] = ptr[1];
        n
    };

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsodi_(
            lsodi_res::<T>,
            lsodi_adda::<T>,
            lsodi_jac::<T>,
            n.as_ptr(),
            y0.as_mut_ptr(),
            dydt0.as_mut_ptr(),
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
            &mf,
        );
    }
}

pub trait LsodisCallback {
    fn residual(&self, t: f64, y: &[f64], s: &[f64], res: &mut [f64]);
    fn jac(&self, t: f64, y: &[f64], s: &[f64], j: usize, pd: &mut [f64]);
    fn adda(&self, t: f64, y: &[f64], j: usize, pd: &mut [f64]);
}

extern "C" fn lsodis_res<T: LsodisCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    s_ptr: *const c_double,
    r_ptr: *mut c_double,
    _ires: *const c_int,
) {
    let (r, y, s, t) = unsafe {
        (
            slice::from_raw_parts_mut(r_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            slice::from_raw_parts(s_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.residual(t, y, s, r);
}

extern "C" fn lsodis_adda<T: LsodisCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    j_ptr: *const c_int,
    _ian: *const c_int,
    _jan: *const c_int,
    pd_ptr: *mut c_double,
) {
    let j = unsafe { *j_ptr as usize };
    let (pd, y, t) = unsafe {
        (
            slice::from_raw_parts_mut(pd_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.adda(t, y, j - 1, pd);
}

extern "C" fn lsodis_jac<T: LsodisCallback>(
    n: *const c_int,
    t_ptr: *const c_double,
    y_ptr: *const c_double,
    s_ptr: *const c_double,
    j_ptr: *const c_int,
    _ian: *const c_int,
    _jan: *const c_int,
    pd_ptr: *mut c_double,
) {
    let j = unsafe { *j_ptr as usize };
    let (pd, y, s, t) = unsafe {
        (
            slice::from_raw_parts_mut(pd_ptr, *n as usize),
            slice::from_raw_parts(y_ptr, *n as usize),
            slice::from_raw_parts(s_ptr, *n as usize),
            *t_ptr,
        )
    };
    let callback = unsafe {
        let ptr = [*n.offset(1), *n.offset(2)];
        std::mem::transmute::<[i32; 2], &T>(ptr)
    };
    callback.jac(t, y, s, j - 1, pd);
}

pub fn dlsodis<'a, 'b, T: LsodisCallback>(
    callback: &T,
    y0: &mut [f64],
    dydt0: &mut [f64],
    mut t0: f64,
    t1: f64,
    rtol: f64,
    atol: f64,
    rwork: &mut [f64],
    iwork: &mut [i32],
    mf: i32,
) {
    let n = {
        let mut n = [y0.len() as i32; 3];
        let ptr = unsafe { std::mem::transmute::<&T, [i32; 2]>(callback) };
        n[1] = ptr[0];
        n[2] = ptr[1];
        n
    };

    let itol = 1;
    let itask = 1;
    let iopt = 0;
    let mut istate = 1;

    let lrw = rwork.len() as i32;
    let liw = iwork.len() as i32;

    let _lock = FLAG.lock().unwrap();
    unsafe {
        low::dlsodis_(
            lsodis_res::<T>,
            lsodis_adda::<T>,
            lsodis_jac::<T>,
            n.as_ptr(),
            y0.as_mut_ptr(),
            dydt0.as_mut_ptr(),
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
            &mf,
        );
    }
}
