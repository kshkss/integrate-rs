fn solution_stiff(t: f64) -> [f64; 2] {
    [
        2.0 * (-t).exp() - (-(1000.0 * t)).exp(),
        -((-t).exp()) + (-(1000.0 * t)).exp(),
    ]
}

fn rhs_stiff(y: &[f64], _t: f64) -> Vec<f64> {
    let mut dy = vec![0.0, 0.0];
    dy[0] = 998.0 * y[0] + 1998.0 * y[1];
    dy[1] = -999.0 * y[0] - 1999.0 * y[1];
    dy
}

#[test]
fn stiff() {
    let y0 = [1.0, 0.0];
    let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();
    let atol = 1e-6;
    let rtol = 1e-8;

    let sol = integrate::BDF::new(rhs_stiff).solve(&y0, &ts, atol, rtol);

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
    let g = |_y: &[f64], _t: f64| -> Array2<f64> { array![[998., 1998.], [-999., -1999.]] };
    let y0 = [1.0, 0.0];
    let ts: Vec<f64> = (0..10).map(|i| 0.1 * i as f64).collect();
    let atol = 1e-6;
    let rtol = 1e-8;

    let sol = integrate::BDF::new(rhs_stiff)
        .gen_full_jacobian_by(g)
        .solve(&y0, &ts, atol, rtol);

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

fn rhs_decay(y: &[f64], _t: f64) -> Vec<f64> {
    let mut dy = vec![0.0];
    dy[0] = -y[0];
    dy
}

#[test]
fn decay() {
    let y0 = [1000.0];
    let ts: Vec<f64> = (0..7).map(|i| 1.0 * i as f64).collect();
    let atol = 1e-6;
    let rtol = 1e-8;

    let sol = integrate::BDF::new(rhs_decay).solve(&y0, &ts, atol, rtol);

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
    let g = |_y: &[f64], _t: f64| -> Array2<f64> { array![[-1.]] };
    let ts: Vec<f64> = (0..7).map(|i| 1.0 * i as f64).collect();
    let atol = 1e-6;
    let rtol = 1e-8;

    let sol = integrate::BDF::new(rhs_decay)
        .gen_full_jacobian_by(g)
        .solve(&y0, &ts, atol, rtol);

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
    let f = |y: &[f64], t: f64| {
        let mut dy = vec![0.0];
        dy[0] = t * y[0];
        dy
    };
    let sol = integrate::BDF::new(f).solve(&y0, &ts, 1e-6, 1e-6);
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
    let f = |y: &[f64], t: f64| {
        let mut dy = vec![0.0];
        dy[0] = t * y[0];
        dy
    };
    let g = |_y: &[f64], t: f64| array![[t]];
    let sol = integrate::BDF::new(f)
        .gen_full_jacobian_by(g)
        .solve(&y0, &ts, 1e-6, 1e-6);
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
    let f = |y: &[f64], t: f64| vec![t * y[0], (t + 2.) * y[1]];
    let sol = integrate::BDF::new(f)
        .gen_banded_jacobian(0, 0)
        .solve(&y0, &ts, 1e-6, 1e-6);
    println!("{:?}", sol);
    assert!(
        (sol[1][0] - y0[0] * 0.5_f64.exp()).abs() < 1e-3,
        "error too large"
    );
}
