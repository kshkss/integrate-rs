use integrate::odepack::Lsode;
use ndarray::prelude::*;

/// This example has the setup from:
/// https://docs.juliadiffeq.org/stable/tutorials/ode_example/#Example-3:-Solving-Nonhomogeneous-Equations-using-Parameterized-Functions-1
fn main() {
    let l: f64 = 1.0;
    let m: f64 = 1.0;
    let g: f64 = 9.81;
    let torque = |t: f64| 0.1 * t.sin();

    let f = |t: f64, y: &[f64], dy: &mut [f64]| {
        dy[0] = y[1];
        dy[1] = -3.0 * g / (2.0 * l) * y[0].sin() + 3.0 / (m * l * l) * torque(t);
    };

    let y0 = [0.01, 0.0];
    let ts = Array1::linspace(0.0, 10.0, 1000).to_vec();
    let sol = Lsode::new(&f, Default::default()).run(&y0, &ts);

    for (i, t) in sol.iter().zip(ts) {
        println!("{} {} {}", t, i[0], i[1]);
    }
}
