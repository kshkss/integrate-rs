use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref FLAG: Mutex<()> = Mutex::<()>::new(());
}

pub mod high;
pub mod mid;
pub mod low;

pub use high::{Adams, BDF, Control};
