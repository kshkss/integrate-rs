use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

#[derive(Debug, Clone)]
pub struct Control {
    pub atol: f64,
    pub rtol: f64,
    pub max_steps: usize,
}

impl Default for Control {
    fn default() -> Self {
        Self {
            atol: 1e-8,
            rtol: 1e-8,
            max_steps: 500,
        }
    }
}

pub mod low;
pub mod lsode;
pub mod lsodi;
pub mod mid;

pub use lsode::Lsode;
pub use lsodi::Lsodi;
