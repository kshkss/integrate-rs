use sprs::{CsMat, TriMat};
use ndarray::prelude::*;
use num_traits::Num;

#[derive(Debug, Clone)]
pub struct SparseMatrix<F>(CsMat<F>);

impl<F: Num + Clone> SparseMatrix<F> {
    fn new(shape: (usize, usize), row_indices: Vec<usize>, col_indices: Vec<usize>, values: Vec<F>) -> Self {
        Self(TriMat::from_triplets(shape, row_indices, col_indices, values).to_csc())
    }
}

impl<F: Num + Clone + Default> From<CsMat<F>> for SparseMatrix<F> {
    fn from(mat: CsMat<F>) -> Self {
        if mat.is_csr() {
            Self(mat.to_csc())
        } else {
            Self(mat)
        }
    }
}
