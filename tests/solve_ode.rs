use integrate::odepack::Lsode;

fn solution_stiff(t: f64) -> [f64; 2] {
    [
        2.0 * (-t).exp() - (-(1000.0 * t)).exp(),
        -((-t).exp()) + (-(1000.0 * t)).exp(),
    ]
}

fn rhs_stiff(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = 998.0 * y[0] + 1998.0 * y[1];
    dy[1] = -999.0 * y[0] - 1999.0 * y[1];
}

#[test]
fn stiff() {
    let y0 = [1.0, 0.0];
    let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();

    let sol = Lsode::new(&rhs_stiff, Default::default()).run(&ts, &y0);

    for (analytical, calculated) in ts.iter().map(|x| solution_stiff(*x)).zip(sol) {
        assert!(
            (analytical[0] - calculated[0]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[0],
            calculated[0]
        );
        assert!(
            (analytical[1] - calculated[1]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[1],
            calculated[1]
        );
    }
}

use ndarray::prelude::*;

#[test]
fn stiff_with_jacobian() {
    let g = |_t: f64, _y: &[f64], mut jac: ArrayViewMut2<f64>| {
        jac.assign(&array![[998., 1998.], [-999., -1999.]])
    };
    let y0 = [1.0, 0.0];
    let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();

    let sol = Lsode::new(&rhs_stiff, Default::default())
        .gen_full_jacobian_by(&g)
        .run(&ts, &y0);

    for (analytical, calculated) in ts.iter().map(|x| solution_stiff(*x)).zip(sol) {
        assert!(
            (analytical[0] - calculated[0]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[0],
            calculated[0]
        );
        assert!(
            (analytical[1] - calculated[1]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[1],
            calculated[1]
        );
    }
}

fn solution_decay(t: f64) -> [f64; 1] {
    [1000.0 * (-t).exp()]
}

fn rhs_decay(_t: f64, y: &[f64], dy: &mut [f64]) {
    dy[0] = -y[0];
}

#[test]
fn decay() {
    let y0 = [1000.0];
    let ts: Vec<f64> = (0..7).map(|i| 1.0 * i as f64).collect();

    let sol = Lsode::new(&rhs_decay, Default::default()).run(&ts, &y0);

    println!("{:?}", sol);

    for (analytical, calculated) in ts.iter().map(|x| solution_decay(*x)).zip(sol) {
        assert!(
            (analytical[0] - calculated[0]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[0],
            calculated[0]
        );
    }
}

#[test]
fn decay_with_jacobian() {
    let y0 = [1000.0];
    let g = |_t: f64, _y: &[f64], mut jac: ArrayViewMut2<f64>| {
        jac[[0, 0]] = -1.;
    };
    let ts: Vec<f64> = (0..7).map(|i| 1.0 * i as f64).collect();

    let sol = Lsode::new(&rhs_decay, Default::default())
        .gen_full_jacobian_by(&g)
        .run(&ts, &y0);

    println!("{:?}", sol);

    for (analytical, calculated) in ts.iter().map(|x| solution_decay(*x)).zip(sol) {
        assert!(
            (analytical[0] - calculated[0]).abs() < 1e-3,
            "|{} - {}| calculated and expected results are suspiciously different",
            analytical[0],
            calculated[0]
        );
    }
}

#[test]
fn closure_rhs() {
    let y0 = [1.0];
    let ts = vec![0.0, 1.0];
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = t * y[0];
    };
    let sol = Lsode::new(&f, Default::default()).run(&ts, &y0);
    println!("{:?}", sol);
    assert!(
        (sol[1][0] - y0[0] * 0.5_f64.exp()).abs() < 1e-3,
        "error too large"
    );
}

#[test]
fn closure_rhs_with_jacobian() {
    let y0 = [1.0];
    let ts = vec![0.0, 1.0];
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = t * y[0];
    };
    let g = |t: f64, _y: &[f64], mut dy: ArrayViewMut2<f64>| dy[[0, 0]] = t;
    let sol = Lsode::new(&f, Default::default())
        .gen_full_jacobian_by(&g)
        .run(&ts, &y0);
    println!("{:?}", sol);
    assert!(
        (sol[1][0] - y0[0] * 0.5_f64.exp()).abs() < 1e-3,
        "error too large"
    );
}

#[test]
fn closure_rhs_estimate_banded_jacobian() {
    let y0 = [1.0, 10.];
    let ts = vec![0.0, 1.0];
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = t * y[0];
        dy[1] = (t + 2.) * y[1];
    };
    let sol = Lsode::new(&f, Default::default())
        .gen_banded_jacobian(0, 0)
        .run(&ts, &y0);
    println!("{:?}", sol);
    assert!(
        (sol[1][0] - y0[0] * 0.5_f64.exp()).abs() < 1e-3,
        "error too large"
    );
}
