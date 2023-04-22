use integrate::radau::Radau;

#[test]
fn test_case01() {
    let f = |t: f64, _y: &[f64], dy: &mut [f64]| {
        dy[0] = 1.;
        dy[1] = t;
    };
    Radau::new(&f, Default::default()).solve(vec![0., 10.], vec![1., 1.]);
}

#[test]
fn test_case02() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -y[0];
    };
    Radau::new(&f, Default::default()).solve(vec![0., 100.], vec![1., 0.]);
}

#[test]
fn test_case03() {
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[0] * 0.5 / (t + 1.) - 2. * t * y[1];
        dy[1] = y[1] * 0.5 / (t + 1.) + 2. * t * y[0];
    };
    Radau::new(&f, Default::default()).solve(vec![0., 10.], vec![1., 0.]);
}

#[test]
fn test_case04() {
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -2. * y[0] + y[1] - t.cos();
        dy[1] = 1998. * y[0] - 1999. * y[1] + 1999. * t.cos() - t.sin();
    };
    Radau::new(&f, Default::default()).solve(vec![0., 10.], vec![1., 2.]);
}

#[test]
fn test_case05() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = 1e-6_f64.powi(2) * ((1. - y[0].powi(2)) * y[1] - y[0]);
    };
    Radau::new(&f, Default::default()).solve(vec![0., 2.], vec![2., -0.66]);
}
