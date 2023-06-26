/// A struct to represent Matrix
pub struct Matrix {
    row: usize,
    col: usize,
    data: Vec<Vec<u8>>,
}

impl Matrix {
    /// Create a new matrix and fill it with 0s.
    /// # Arguments
    ///
    /// * `row_size` - Size of the row of the matrix
    /// * `col_size` - Size of the col of the matrix
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new(3, 3);
    /// ```
    pub fn new(row_size: usize, col_size: usize) -> Matrix {
        let data: Vec<Vec<u8>> = vec![vec![0; col_size]; row_size];

        Matrix {
            row: row_size,
            col: col_size,
            data,
        }
    }

    /// Create a new identity matrix and fill the primary diagonal
    /// with 1 and the rest with 0s.
    /// # Arguments
    ///
    /// * `size` - Size of the identity matrix
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new_identity(3);
    /// ```
    pub fn new_identity(size: usize) -> Matrix {
        let mut data: Vec<Vec<u8>> = vec![vec![0; size]; size];

        for i in 0..size {
            data[i][i] = 1;
        }

        Matrix {
            row: size,
            col: size,
            data,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let matrix = Matrix::new(3, 3);

        assert_eq!(matrix.row, 3);
        assert_eq!(matrix.col, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        for row in matrix.data.iter() {
            for &elem in row.iter() {
                assert_eq!(elem, 0);
            }
        }
    }
    #[test]
    fn test_new_identity() {
        let matrix = Matrix::new_identity(3);

        assert_eq!(matrix.row, 3);
        assert_eq!(matrix.col, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        for (row_index, row) in matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                if row_index == col_index {
                    assert_eq!(elem, 1);
                } else {
                    assert_eq!(elem, 0);
                }
            }
        }
    }
}
