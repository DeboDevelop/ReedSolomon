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

    /// Create a new sub matrix from the given matrix (self), r_start, 
    /// r_end, c_start, c_end.
    /// # Arguments
    ///
    /// * `r_start` - Starting index of the row in given matrix
    /// * `r_end` - Ending index of the row in given matrix
    /// * `c_start` - Starting index of the col in given matrix
    /// * `c_end` - Ending index of the col in given matrix
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    /// 
    /// let matrix = Matrix::new_identity(3);
    /// let sub_matrix = matrix.new_sub_matrix(1, 3, 1, 3);
    /// ```
    pub fn new_sub_matrix(&self, r_start: usize, r_end: usize, c_start: usize, c_end: usize) -> Matrix {
        let rows = r_end - r_start;
        let cols = c_end - c_start;
        let mut data: Vec<Vec<u8>> = vec![vec![0; cols]; rows];

        for r in r_start..r_end {
            for c in c_start..c_end {
                data[r - r_start][c - c_start] = self.data[r][c];
            }
        }

        Matrix { rows, cols, data }
    }

    /// Create a new augmented matrix from the given Matrices - self, right.
    /// # Arguments
    ///
    /// * `right` - Right side of the augmented matrix.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    /// 
    /// let left = Matrix::new_identity(3);
    /// let right = Matrix::new_identity(3);
    /// let augmented_matrix = left.new_augmented_matrix(right);
    /// ```
    pub fn new_augmented_matrix(&self, right: Matrix) -> Matrix {
        if self.rows != right.rows {
            panic!(
                "Row count of the matrices must match. Current row count, left: {}, right: {}",
                self.rows, right.rows
            )
        }

        let cols = self.cols + right.cols;
        let mut data: Vec<Vec<u8>> = vec![vec![0; cols]; self.rows];
        for r in 0..self.rows {
            for c in 0..self.cols {
                data[r][c] = self.data[r][c];
            }
            for c in 0..right.cols {
                data[r][self.cols + c] = right.data[r][c];
            }
        }

        Matrix { rows: self.rows, cols, data }
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
    #[test]
    fn test_new_sub_matrix() {
        let gf8 = GaloisField::new();
        let matrix = Matrix::new_vandermonde(3, 3, gf8);
        let sub_matrix = matrix.new_sub_matrix(1, matrix.rows, 1, matrix.cols);
        let exp_res: [[u8; 2]; 2] = [[1, 1], [2, 4]];

        assert_eq!(sub_matrix.rows, 2);
        assert_eq!(sub_matrix.cols, 2);
        assert_eq!(sub_matrix.data.len(), 2);
        assert_eq!(sub_matrix.data[0].len(), 2);
        for (row_index, row) in sub_matrix.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_new_augmented_matrix() {
        let gf8 = GaloisField::new();
        let left = Matrix::new_vandermonde(3, 3, gf8);
        let right = Matrix::new_identity(3);
        let res = left.new_augmented_matrix(right);
        let exp_res: [[u8; 6]; 3] = [[1, 0, 0, 1, 0, 0], [1, 1, 1, 0, 1, 0], [1, 2, 4, 0, 0, 1]];

        assert_eq!(res.rows, 3);
        assert_eq!(res.cols, 6);
        assert_eq!(res.data.len(), 3);
        assert_eq!(res.data[0].len(), 6);
        for (row_index, row) in res.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
}
