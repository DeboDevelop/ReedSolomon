/// A struct to represent Matrix
pub struct Matrix {
    row: usize,
    col: usize,
    data: Vec<Vec<u8>>
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
            data
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
}
