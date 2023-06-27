use crate::galois::GaloisField;

/// A struct to represent Matrix
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<Vec<u8>>,
}

impl Matrix {
    /// Create a new matrix and fill it with 0s.
    /// # Arguments
    ///
    /// * `rows` - Size of the row of the matrix
    /// * `cols` - Size of the col of the matrix
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new(3, 3);
    /// ```
    pub fn new(rows: usize, cols: usize) -> Matrix {
        let data: Vec<Vec<u8>> = vec![vec![0; cols]; rows];

        Matrix { rows, cols, data }
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
            rows: size,
            cols: size,
            data,
        }
    }

    /// Create a new vandermonde matrix where the given property is guaranteed.
    /// Any subset of rows that forms a square matrix is invertible.
    /// # Arguments
    ///
    /// * `rows` - Size of the row of the matrix
    /// * `cols` - Size of the col of the matrix
    /// * `gf` - Galois Field for the elements of matrix
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    /// use reed_solomon::galois::GaloisField;
    /// 
    /// let gf8 = GaloisField::new();
    /// let matrix = Matrix::new_vandermonde(3, 3, gf8);
    /// ```
    pub fn new_vandermonde(rows: usize, cols: usize, gf: GaloisField) -> Matrix {
        let mut data: Vec<Vec<u8>> = vec![vec![0; cols]; rows];

        for r in 0..rows {
            for c in 0..cols {
                data[r][c] = gf.exp(r as u8, c);
            }
        }

        Matrix { rows, cols, data }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let matrix = Matrix::new(3, 3);

        assert_eq!(matrix.rows, 3);
        assert_eq!(matrix.cols, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        for row in matrix.data.iter() {
            for &elem in row.iter() {
                assert_eq!(0, elem);
            }
        }
    }
    #[test]
    fn test_new_identity() {
        let matrix = Matrix::new_identity(3);

        assert_eq!(matrix.rows, 3);
        assert_eq!(matrix.cols, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        for (row_index, row) in matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                if row_index == col_index {
                    assert_eq!(1, elem);
                } else {
                    assert_eq!(0, elem);
                }
            }
        }
    }
    #[test]
    fn test_new_vandermonde() {
        let gf8 = GaloisField::new();
        let matrix = Matrix::new_vandermonde(3, 3, gf8);
        let exp_res: [[u8; 3]; 3] = [[1, 0, 0], [1, 1, 1], [1, 2, 4]];

        assert_eq!(matrix.rows, 3);
        assert_eq!(matrix.cols, 3);
        assert_eq!(matrix.data.len(), 3);
        assert_eq!(matrix.data[0].len(), 3);
        for (row_index, row) in matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
}
