//! implicit Runge-Kutta method (Radau IIA) of variable order (switches automatically between orders 5, 9, and 13)
//! for problems of the form My'=f(x,y) with possibly singular matrix M; For the choices IWORK(11)=3 and IWORK(12)=3,
//! the code is mathematically equivalent to RADAU5 (in general a little bit slower than RADAU5).

use libc::{c_double, c_int};

#[link(name = "gfortran")]
extern "C" {
    pub fn radau_(
        N: *const c_int,
        FCN: extern "C" fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *mut c_double,
            *mut c_double,
            *mut c_int,
        ),
        X: *mut c_double,
        Y: *mut c_double,
        XEND: *const c_double,
        H: *const c_double,
        RTOL: *const c_double,
        ATOL: *const c_double,
        ITOL: *const c_int,
        JAC: extern "C" fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *mut c_double,
            *const c_int,
            *mut c_double,
            *mut c_double,
        ),
        IJAC: *const c_int,
        MLJAC: *const c_int,
        MUJAC: *const c_int,
        MAS: extern "C" fn(*const c_int, *mut c_double, *const c_int, *mut c_double, *mut c_int),
        IMAS: *const c_int,
        MLMAS: *const c_int,
        MUMAS: *const c_int,
        SOLOUT: extern "C" fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *const c_double,
            *const c_double,
            *const c_int,
            *const c_int,
            *mut c_double,
            *mut c_int,
            *mut c_int,
        ),
        IOUT: *const c_int,
        WORK: *mut c_double,
        LWORK: *const c_int,
        IWORK: *mut c_int,
        LIWORK: *const c_int,
        RPAR: *mut c_double,
        IPAR: *mut c_int,
        IDID: *mut c_int,
    );

    /// THIS FUNCTION CAN BE USED FOR CONINUOUS OUTPUT. IT PROVIDES AN
    /// APPROXIMATION TO THE I-TH COMPONENT OF THE SOLUTION AT X.
    /// IT GIVES THE VALUE OF THE COLLOCATION POLYNOMIAL, DEFINED FOR
    /// THE LAST SUCCESSFULLY COMPUTED STEP (BY RADAU).
    pub fn contra_(
        I: *const c_int,
        S: *const c_double,
        CONT: *const c_double,
        LRC: *const c_int,
    ) -> c_double;
}
