pub mod galois;
pub mod matrix;

use crate::galois::GaloisField;
use crate::matrix::Matrix;

/// A Struct to represent and store data for Reed Solomon Erasure Coding.
pub struct ReedSolomon {
    data_shard_count: usize,
    parity_shard_count: usize,
    total_shard_count: usize,
    parity: Matrix,
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
    /// let rs = ReedSolomon::new(4, 2);
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

        let mut parity = Matrix::new(parity_shards, data_shards);
        for i in 0..parity_shards {
            parity.data[i] = matrix.data[data_shards + i].clone();
        }

        ReedSolomon {
            data_shard_count: data_shards,
            parity_shard_count: parity_shards,
            total_shard_count: total_shards,
            parity,
            gf,
            matrix,
        }
    }

    /// Check the consistency of shards passed to other methods.
    /// # Arguments
    ///
    /// * `shards` - All shards including data and parity shards.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let shards = vec![vec![0, 1, 2], vec![3, 4, 5], vec![200, 201, 203], vec![100, 101, 102]];
    /// rs.check_shard_sizes(shards);
    /// ```
    pub(crate) fn check_shard_sizes(&self, shards: &Vec<Vec<u8>>) {
        if shards.len() != self.total_shard_count {
            panic!("Wrong no. of shards");
        }

        let shard_elem_len = shards[0].len();
        if shard_elem_len == 0 {
            panic!("Empty Shard");
        }
        for elem in shards.iter() {
            if (*elem).len() != shard_elem_len {
                panic!("Length of shards are different");
            }
        }
    }

    /// Encodes checksum shards for a set of data shards.
    /// # Arguments
    ///
    /// * `shards` - All shards including data and parity shards. Parity shards will be overwritten.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let shards = vec![vec![0, 1, 2], vec![3, 4, 5], vec![200, 201, 203], vec![100, 101, 102]];
    /// let encoded_shards = rs.encode(shards);
    /// ```
    pub fn encode(&self, shards: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
        self.check_shard_sizes(&shards);

        let mut inputs = shards[..self.data_shard_count].to_vec();
        let mut outputs = shards[self.data_shard_count..].to_vec();

        self.encode_shards(&inputs, &mut outputs);

        inputs.extend(outputs);

        inputs
    }

    /// Encodes checksum shards for a given input (data shards) and modifies the output.
    /// # Arguments
    ///
    /// * `inputs` - Data shards.
    /// * `outputs` - Parity shards (to be overwritten).
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let inputs = vec![vec![0, 1, 2], vec![3, 4, 5]];
    /// let mut outputs = vec![vec![200, 201, 203], vec![100, 101, 102]];
    /// rs.encode_shards(&inputs, &mut outputs);
    /// ```
    pub(crate) fn encode_shards(&self, inputs: &Vec<Vec<u8>>, outputs: &mut Vec<Vec<u8>>) {
        for inp in 0..self.data_shard_count {
            for out in 0..self.parity_shard_count {
                let parity_byte = self.parity.data[out][inp];
                if inp == 0 {
                    for (i_byte, input) in inputs[inp].iter().enumerate() {
                        outputs[out][i_byte] = self.gf.mul(parity_byte, *input);
                    }
                } else {
                    let mut val: u8;
                    for (i_byte, input) in inputs[inp].iter().enumerate() {
                        val = self.gf.mul(parity_byte, *input);
                        outputs[out][i_byte] = GaloisField::add(outputs[out][i_byte], val);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let rs = ReedSolomon::new(4, 2);
        let exp_res: [[u8; 4]; 6] = [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [27, 28, 18, 20],
            [28, 27, 20, 18],
        ];

        for (row_index, row) in rs.matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
        for (row_index, row) in rs.parity.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[rs.data_shard_count + row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_encode() {
        let rs = ReedSolomon::new(2, 2);
        let shards = vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![200, 201, 203],
            vec![100, 101, 102],
        ];
        let encoded_shard = rs.encode(shards);
        let exp_res: [[u8; 3]; 4] = [[0, 1, 2], [3, 4, 5], [6, 11, 12], [5, 14, 11]];
        for (row_index, row) in encoded_shard.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
}
