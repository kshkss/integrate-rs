#[derive(Debug, Clone)]
pub struct Control {
    pub atol: f64,
    pub rtol: f64,
    pub max_steps: usize,
    pub max_step_size: f64,
}

impl Default for Control {
    fn default() -> Self {
        Self {
            atol: 1e-8,
            rtol: 1e-8,
            max_steps: 500,
            max_step_size: f64::INFINITY,
        }
    }
}

pub mod lsode;
pub mod lsodi;
pub mod mid;

pub use lsode::Lsode;
pub use lsodi::Lsodi;
