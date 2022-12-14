use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

pub mod low;
pub mod lsode;
pub mod lsodi;
pub mod mid;

pub use lsode::{Control, Lsode};
pub use lsodi::Lsodi;
