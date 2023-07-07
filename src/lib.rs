pub mod error;
pub mod galois;
pub mod matrix;

use crate::error::Error;
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
    pub(crate) fn build_matrix(
        data_shards: usize,
        total_shards: usize,
        gf: GaloisField,
    ) -> Result<Matrix, Error> {
        // Start with a Vandermonde matrix but this matrix doesn't have the property
        // that the data shards are unchanged after encoding.
        let vandermonde = Matrix::new_vandermonde(total_shards, data_shards, gf);

        // Multiply the inverse of the top square of the matrix with matrix.
        // This will make the top square of the matrix be the identity matrix, but
        // will preserve the property that any square subset of rows is invertible.
        let top = vandermonde.new_sub_matrix(0, data_shards, 0, data_shards);
        let top_inv = top.invert(gf)?;

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
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(4, 2);
    /// ```
    pub fn new(data_shards: usize, parity_shards: usize) -> Result<ReedSolomon, Error> {
        if data_shards == 0 {
            return Err(Error::ZeroDataShards);
        }
        if parity_shards == 0 {
            return Err(Error::ZeroParityShards);
        }
        // More than 256 will lead to duplicate rows in the Vandermonde matrix,
        // which would then lead to duplicate rows in the built matrix.
        // Any subset of the rows containing the duplicate rows would
        // be singular and thus non-invertible.
        if data_shards + parity_shards > 256 {
            return Err(Error::ShardsOverflow);
        }

        let gf = GaloisField::new();
        let total_shards = data_shards + parity_shards;

        let matrix = Self::build_matrix(data_shards, total_shards, gf)?;

        let mut parity = Matrix::new(parity_shards, data_shards);
        for i in 0..parity_shards {
            parity.data[i] = matrix.data[data_shards + i].clone();
        }

        Ok(ReedSolomon {
            data_shard_count: data_shards,
            parity_shard_count: parity_shards,
            total_shard_count: total_shards,
            parity,
            gf,
            matrix,
        })
    }

    /// Check the consistency of shards passed to other methods.
    /// # Arguments
    ///
    /// * `shards` - All shards including data and parity shards.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let shards = vec![vec![0, 1, 2], vec![3, 4, 5], vec![200, 201, 203], vec![100, 101, 102]];
    /// rs.check_shard_sizes(shards);
    /// ```
    pub(crate) fn check_shard_sizes(&self, shards: &Vec<Vec<u8>>) -> Result<(), Error> {
        if shards.len() != self.total_shard_count {
            return Err(Error::WrongNoOfShards);
        }

        let shard_elem_len = shards[0].len();
        if shard_elem_len == 0 {
            return Err(Error::EmptyShards);
        }
        for elem in shards.iter() {
            if (*elem).len() != shard_elem_len {
                return Err(Error::InconsistentShards);
            }
        }

        Ok(())
    }

    /// Encodes checksum shards for a set of data shards.
    /// Returns all the shards including all data and parity shards.
    /// # Arguments
    ///
    /// * `shards` - All shards including data and parity shards. Parity shards will be overwritten.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let shards = vec![vec![0, 1, 2], vec![3, 4, 5], vec![200, 201, 203], vec![100, 101, 102]];
    /// let encoded_shards = rs.encode(shards);
    /// ```
    pub fn encode(&self, shards: Vec<Vec<u8>>) -> Result<Vec<Vec<u8>>, Error> {
        self.check_shard_sizes(&shards)?;

        let mut inputs = shards[..self.data_shard_count].to_vec();
        let mut outputs = shards[self.data_shard_count..].to_vec();

        self.encode_shards(&self.parity, &inputs, &mut outputs);

        inputs.extend(outputs);

        Ok(inputs)
    }

    /// Encodes checksum shards for a given input (data shards) and modifies the output.
    /// # Arguments
    ///
    /// * `inputs` - Data shards.
    /// * `outputs` - Parity shards (to be overwritten).
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let inputs = vec![vec![0, 1, 2], vec![3, 4, 5]];
    /// let mut outputs = vec![vec![200, 201, 203], vec![100, 101, 102]];
    /// rs.encode_shards(&inputs, &mut outputs);
    /// ```
    pub(crate) fn encode_shards(
        &self,
        parity: &Matrix,
        inputs: &Vec<Vec<u8>>,
        outputs: &mut Vec<Vec<u8>>,
    ) {
        for inp in 0..self.data_shard_count {
            for out in 0..self.parity_shard_count {
                let parity_byte = (*parity).data[out][inp];
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

    /// Check the no. and consistency of shards passed to decode methods.
    /// # Arguments
    ///
    /// * `shards` - Given shards including data and parity shards. Some shards might be missing.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let shards = vec![vec![0, 1, 2], vec![], vec![6, 11, 12], vec![]];
    /// rs.check_shard_sizes_for_decode(shards);
    /// ```
    pub(crate) fn check_shard_sizes_for_decode(
        &self,
        shards: &Vec<Vec<u8>>,
    ) -> Result<(usize, usize), Error> {
        if shards.len() < self.total_shard_count {
            return Err(Error::TooFewShards);
        }
        if shards.len() > self.total_shard_count {
            return Err(Error::TooManyShards);
        }

        let mut shard_elem_len = 0;
        let mut present: usize = 0;
        for elem in shards.iter() {
            if (*elem).len() != 0 {
                present += 1;
                shard_elem_len = (*elem).len();
            }
        }
        let mut same_size: usize = 0;
        for elem in shards.iter() {
            if (*elem).len() == shard_elem_len {
                same_size += 1;
            }
        }
        if present != same_size {
            return Err(Error::InconsistentShards);
        }
        if present < self.data_shard_count {
            return Err(Error::TooFewShards);
        }

        Ok((present, shard_elem_len))
    }

    /// Takes shards as input and recover any data or parity shards that is missing.
    /// Returns all the shards including all data and parity shards.
    /// # Arguments
    ///
    /// * `shards` - Given shards including data and parity shards. Some shards might be missing.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::ReedSolomon;
    ///
    /// let rs = ReedSolomon::new(2, 2);
    /// let inputs = vec![vec![0, 1, 2], vec![3, 4, 5]];
    /// let mut outputs = vec![vec![200, 201, 203], vec![100, 101, 102]];
    /// rs.encode_shards(&inputs, &mut outputs);
    /// ```
    pub fn decode(&self, shards: Vec<Vec<u8>>) -> Result<Vec<Vec<u8>>, Error> {
        let (present, shard_elem_len) = self.check_shard_sizes_for_decode(&shards)?;

        if present == self.total_shard_count {
            // All of the shards have data so we can return
            return Ok(shards);
        }

        // Pull out the rows of the matrix that correspond
        // to the given shards and build a square matrix.
        // This matrix could be used to generate the shards
        // that we have from the original data.
        //
        // Create an array holding just the shards that
        // correspond to the rows of the submatrix. These
        // shards will be the input to the decoding process
        // that re-creates the missing data shards.
        let mut sub_matrix = Matrix::new(self.data_shard_count, self.data_shard_count);
        let mut sub_shard: Vec<Vec<u8>> = vec![vec![]; self.data_shard_count];
        let mut sub_matrix_row: usize = 0;
        let mut matrix_row: usize = 0;
        while matrix_row < self.total_shard_count && sub_matrix_row < self.data_shard_count {
            if shards[matrix_row].len() != 0 {
                sub_matrix.data[sub_matrix_row] = self.matrix.data[matrix_row].clone();
                sub_shard[sub_matrix_row] = shards[matrix_row].clone();
                sub_matrix_row += 1;
            }
            matrix_row += 1;
        }
        // Invert the matrix, so we can go from the encoded shards
        // back to the original data. Then pull out the row that
        // generates the shard that we want to decode. Since this
        // matrix maps back to the orginal data, it can be used
        // to create a data shard, but not a parity shard.
        let data_decode_matrix = sub_matrix.invert(self.gf)?;

        // Re-create any data shards that were missing.
        //
        // The input to the coding is all of the shards we actually
        // have, and the output is the missing data shards. The computation
        // is done using the special decode matrix we just built.
        let mut matrix_rows = Matrix::new(self.parity_shard_count, self.parity_shard_count);
        let mut outputs: Vec<Vec<u8>> = vec![vec![0; shard_elem_len]; self.parity_shard_count];
        let mut output_count: usize = 0;
        for i in 0..self.data_shard_count {
            if shards[i].len() == 0 {
                matrix_rows.data[output_count] = data_decode_matrix.data[i].clone();
                output_count += 1;
            }
        }
        self.encode_shards(&matrix_rows, &sub_shard, &mut outputs);

        // Filling the missing data shards.
        output_count = 0;
        let mut shards = shards;
        for i in 0..self.data_shard_count {
            if shards[i].len() == 0 {
                shards[i] = outputs[output_count].clone();
                output_count += 1;
            }
        }

        // Filling missing parity shards with placeholder
        for i in self.data_shard_count..self.total_shard_count {
            if shards[i].len() == 0 {
                shards[i] = vec![0; shard_elem_len]
            }
        }
        // Now that we have all of the data shards intact, we can
        // compute any of the parity that is missing.
        //
        // The input to the coding is ALL of the data shards, including
        // any that we just calculated. The output is all parity shards.
        self.encode(shards)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let result = ReedSolomon::new(4, 2);
        let rs = match result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
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
        let result = ReedSolomon::new(2, 2);
        let rs = match result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let shards = vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![200, 201, 203],
            vec![100, 101, 102],
        ];
        let encoded_shard_result = rs.encode(shards);
        let encoded_shard = match encoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let exp_res: [[u8; 3]; 4] = [[0, 1, 2], [3, 4, 5], [6, 11, 12], [5, 14, 11]];
        for (row_index, row) in encoded_shard.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_decode_missing_data() {
        let result = ReedSolomon::new(2, 2);
        let rs = match result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let shards = vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![200, 201, 203],
            vec![100, 101, 102],
        ];
        let encoded_shard_result = rs.encode(shards.clone());
        let encoded_shard = match encoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let broken_shards = vec![vec![0, 1, 2], vec![], vec![6, 11, 12], vec![5, 14, 11]];
        let decoded_shard_result = rs.decode(broken_shards);
        let decoded_shard = match decoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        for (row_index, row) in decoded_shard.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(encoded_shard[row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_decode_missing_parity() {
        let result = ReedSolomon::new(2, 2);
        let rs = match result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let shards = vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![200, 201, 203],
            vec![100, 101, 102],
        ];
        let encoded_shard_result = rs.encode(shards.clone());
        let encoded_shard = match encoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let broken_shards = vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 11, 12], vec![]];
        let decoded_shard_result = rs.decode(broken_shards);
        let decoded_shard = match decoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        for (row_index, row) in decoded_shard.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(encoded_shard[row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_decode_missing_data_and_parity() {
        let result = ReedSolomon::new(2, 2);
        let rs = match result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let shards = vec![
            vec![0, 1, 2],
            vec![3, 4, 5],
            vec![200, 201, 203],
            vec![100, 101, 102],
        ];
        let encoded_shard_result = rs.encode(shards.clone());
        let encoded_shard = match encoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        let broken_shards = vec![vec![0, 1, 2], vec![], vec![], vec![5, 14, 11]];
        let decoded_shard_result = rs.decode(broken_shards);
        let decoded_shard = match decoded_shard_result {
            Ok(x) => x,
            Err(e) => panic!("{}", e),
        };
        for (row_index, row) in decoded_shard.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(encoded_shard[row_index][col_index], elem);
            }
        }
    }
}
