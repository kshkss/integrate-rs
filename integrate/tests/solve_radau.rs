use approx::assert_relative_eq;
use integrate::radau::Radau;
use ndarray::Array1;

#[test]
fn test_case01() {
    let f = |_t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = -0.1 * y[0] - 49.9 * y[1];
        dy[1] = -50. * y[1];
        dy[2] = 70. * y[1] - 120. * y[2];
    };
    let ts = Array1::linspace(0., 10., 100).to_vec();
    let res = Radau::new(&f, Default::default()).solve(ts.clone(), vec![2., 1., 2.]);

    for (&t, y) in ts.iter().zip(&res) {
        assert_relative_eq!(
            y[0],
            (-0.1 * t).exp() + (-50. * t).exp(),
            max_relative = 1e-8,
            epsilon = 1e-8
        );

        assert_relative_eq!(y[1], (-50. * t).exp(), max_relative = 1e-8, epsilon = 1e-8);

        assert_relative_eq!(
            y[2],
            (-50. * t).exp() + (-120. * t).exp(),
            max_relative = 1e-8,
            epsilon = 1e-8
        );
    }
}

#[test]
fn test_case02() {
    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = 100. * (t.sin() - y[0]);
    };
    let ts = Array1::linspace(0., 10., 100).to_vec();
    let res = Radau::new(&f, Default::default()).solve(ts.clone(), vec![0.]);

    for (&t, y) in ts.iter().zip(&res) {
        assert_relative_eq!(
            y[0],
            (t.sin() - 0.01 * t.cos() + 0.01 * (-100. * t).exp()) / 1.0001,
            max_relative = 1e-8,
            epsilon = 1e-8
        );
    }
}
