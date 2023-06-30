pub mod matrix;
pub mod galois;

use crate::matrix::Matrix;
use crate::galois::GaloisField;

/// A Struct to represent and store data for Reed Solomon Erasure Coding.
pub struct ReedSolomon {
    data_shard_count: usize,
    parity_shard_count: usize,
    total_shard_count: usize,
    gf: GaloisField,
    matrix: Matrix,
}

impl ReedSolomon {
    /// Create a matrix used for encoding.
    /// Since top square of the matrix is guaranteed to be an identity
    /// matrix, the data shards will remain unchanged after encoding.
    /// # Arguments
    ///
    /// * `data_shards` - No. of Data Shards
    /// * `total_shards` - Total no. of Shards (Data + Parity)
    /// * `gf` - Galois Field where all the arithmetic will take place
    ///
    /// # Example
    /// ```
    /// use crate::ReedSolomon;
    /// use crate::galois::GaloisField;
    /// 
    /// let gf = GaloisField::new();
    /// let matrix = ReedSolomon::build_matrix(4, 6, gf);
    /// ```
    pub(crate) fn build_matrix(data_shards: usize, total_shards: usize, gf: GaloisField) -> Matrix {
        // Start with a Vandermonde matrix but this matrix doesn't have the property 
        // that the data shards are unchanged after encoding.
        let vandermonde = Matrix::new_vandermonde(total_shards, data_shards, gf);

        // Multiply the inverse of the top square of the matrix with matrix.
        // This will make the top square of the matrix be the identity matrix, but
        // will preserve the property that any square subset of rows is invertible.
        let top = vandermonde.new_sub_matrix(0, data_shards, 0, data_shards);
        let top_inv = top.invert(gf);

        vandermonde.mul(top_inv, gf)
    }

    /// Create a new Reed Solomon Erasure Coding to be used to encode data.
    /// # Arguments
    ///
    /// * `data_shards` - No. of Data Shards
    /// * `parity_shards` - No. of Parity Shards i.e. Checksum Shards
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;;
    /// 
    /// let matrix = ReedSolomon::new(4, 2);
    /// ```
    pub fn new(data_shards: usize, parity_shards: usize) -> ReedSolomon {
        if data_shards == 0 {
            panic!("Data Shards can't be zero")
        }
        if parity_shards == 0 {
            panic!("Parity Shards can't be zero")
        }
        // More than 256 will lead to duplicate rows in the Vandermonde matrix,
        // which would then lead to duplicate rows in the built matrix.
        // Any subset of the rows containing the duplicate rows would 
        // be singular and thus non-invertible.
        if data_shards + parity_shards > 256 {
            panic!("too many shards - max is 256")
        }

        let gf = GaloisField::new();
        let total_shards = data_shards + parity_shards;

        let matrix = Self::build_matrix(data_shards, total_shards, gf);
        ReedSolomon {
            data_shard_count: data_shards,
            parity_shard_count: parity_shards,
            total_shard_count: total_shards,
            gf,
            matrix,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let rs = ReedSolomon::new(4, 2);
        let exp_res: [[u8; 4]; 6] = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [27, 28, 18, 20], [28, 27, 20, 18]];

        for (row_index, row) in rs.matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
}