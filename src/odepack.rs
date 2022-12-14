use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

pub mod high;
pub mod low;
pub mod lsodi;
pub mod mid;

pub use high::{Control, Lsode};
pub use lsodi::Lsodi;
