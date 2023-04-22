pub mod low;
pub mod mid;
use ndarray::ArrayViewMut2;

use mid::CallbackType;
pub use mid::{Callback, Control};

pub struct Radau<'a, T> {
    callback: &'a T,
    control: Control,
}

impl<'a, T: Callback> Radau<'a, T> {
    pub fn new<'b: 'a>(callback: &'b T, control: Control) -> Self {
        Self { callback, control }
    }

    pub fn solve(&self, ts: Vec<f64>, mut y0: Vec<f64>) -> Vec<Vec<f64>> {
        let mut ys = Vec::with_capacity(ts.len());
        let hinit = (ts[1] - ts[0]) * 1e-8;
        let atol = vec![self.control.atol; y0.len()];
        let rtol = vec![self.control.rtol; y0.len()];
        let jac_type = self.callback.jacobian_type();
        let mas_type = self.callback.mass_type();
        let use_solout = self.callback.use_solout();
        let mut rwork = mid::real_work_space(y0.len(), &jac_type, &mas_type, &self.control);
        let mut iwork = mid::integer_work_space(y0.len(), &jac_type, &mas_type, &self.control);

        for (&t0, &t1) in ts.iter().zip(ts.iter().skip(1)) {
            ys.push(y0.clone());
            mid::radau(
                self.callback,
                &mut y0,
                t0,
                t1,
                hinit,
                &rtol,
                &atol,
                &jac_type,
                &mas_type,
                use_solout,
                &mut rwork,
                &mut iwork,
            );
        }
        ys.push(y0);
        ys
    }
}

impl<F: Fn(f64, &[f64], &mut [f64])> Callback for F {
    fn fcn<'b>(&self, t: f64, y: &'b [f64], f: &'b mut [f64]) {
        self(t, y, f);
    }
    fn jac<'b>(&self, _t: f64, _y: &'b [f64], _pd: ArrayViewMut2<'b, f64>) {
        unimplemented!();
    }
    fn mas<'b>(&self, _mas: ArrayViewMut2<'b, f64>) {
        unimplemented!();
    }
    fn solout<'b, CONTRA>(
        &self,
        _nr: usize,
        _t_old: f64,
        _t: f64,
        _y: &'b [f64],
        _contra: &'b CONTRA,
    ) -> i32
    where
        CONTRA: 'b + Fn(f64) -> Vec<f64>,
    {
        unimplemented!();
    }
    fn jacobian_type(&self) -> CallbackType {
        CallbackType::FullMatrix { supplied: false }
    }
    fn mass_type(&self) -> CallbackType {
        CallbackType::FullMatrix { supplied: false }
    }
    fn use_solout(&self) -> bool {
        false
    }
}

impl<F, JAC> Callback for (F, JAC)
where
    F: Fn(f64, &[f64], &mut [f64]),
    JAC: Fn(f64, &[f64], ArrayViewMut2<f64>),
{
    fn fcn<'b>(&self, t: f64, y: &'b [f64], f: &'b mut [f64]) {
        self.0(t, y, f);
    }
    fn jac<'b>(&self, t: f64, y: &'b [f64], pd: ArrayViewMut2<'b, f64>) {
        self.1(t, y, pd);
    }
    fn mas<'b>(&self, _mas: ArrayViewMut2<'b, f64>) {
        unimplemented!();
    }
    fn solout<'b, CONTRA>(
        &self,
        _nr: usize,
        _t_old: f64,
        _t: f64,
        _y: &'b [f64],
        _contra: &'b CONTRA,
    ) -> i32
    where
        CONTRA: 'b + Fn(f64) -> Vec<f64>,
    {
        unimplemented!();
    }
    fn jacobian_type(&self) -> CallbackType {
        CallbackType::FullMatrix { supplied: true }
    }
    fn mass_type(&self) -> CallbackType {
        CallbackType::FullMatrix { supplied: false }
    }
    fn use_solout(&self) -> bool {
        false
    }
}
