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

    /// Create a new matrix from the given data.
    /// # Arguments
    ///
    /// * `data` - Matrix data
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new_from_data(vec![vec![1, 2, 3], vec![1, 2, 3]]);
    /// ```
    pub fn new_from_data(data: Vec<Vec<u8>>) -> Matrix {
        Matrix {
            rows: data.len(),
            cols: data[0].len(),
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
    pub fn new_sub_matrix(
        &self,
        r_start: usize,
        r_end: usize,
        c_start: usize,
        c_end: usize,
    ) -> Matrix {
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

        Matrix {
            rows: self.rows,
            cols,
            data,
        }
    }

    /// Multiply given 2 matrices - self, right.
    /// # Arguments
    ///
    /// * `right` - 2nd matrix to be multiplied.
    /// * `gf` - Galois Field where the multiplication will occur.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let left = Matrix::new_identity(3);
    /// let right = Matrix::new_identity(3);
    /// let multiplied_matrix = left.mul(right);
    /// ```
    pub fn mul(&self, right: Matrix, gf: GaloisField) -> Matrix {
        if self.cols != right.rows {
            panic!(
                "Colomn count on left has to be same as row count on right. left column: {}, right row: {}",
                self.cols, right.rows
            )
        }

        let mut res = Matrix::new(self.rows, right.cols);
        for r in 0..self.rows {
            for c in 0..right.cols {
                let mut value: u8 = 0;
                for lc in 0..self.cols {
                    let m = gf.mul(self.data[r][lc], right.data[lc][c]);
                    value = GaloisField::add(value, m);
                }
                res.data[r][c] = value;
            }
        }

        res
    }

    /// Returns the inverted matrix of self.
    /// # Arguments
    ///
    /// * `gf` - Galois Field where the multiplication will occur.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new_from_data(vec![vec![56, 23, 98], vec![3, 100, 200], vec![45, 201, 123]]);
    /// let gf8 = GaloisField::new();
    /// let inv_matrix = matrix.invert(gf8);
    /// ```
    pub fn invert(&self, gf: GaloisField) -> Matrix {
        if self.rows != self.cols {
            panic!("Can't invert a non-square matrix")
        }
        // Create a working matrix by augmenting this one with an identity matrix on the right.
        let mut work = self.new_augmented_matrix(Matrix::new_identity(self.rows));

        // Do Gaussian elimination to transform the left half into an identity matrix.
        work.gauss_elim(gf);

        // The right half is now the inverse.
        work.new_sub_matrix(0, self.rows, self.cols, self.cols * 2)
    }

    /// Swap two given rows of Matrix data.
    /// # Arguments
    ///
    /// * `row1` - 1st row to be swapped in the given matrix.
    /// * `row2` - 2nd row to be swapped in the given matrix.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new_identity(3);
    /// matrix.swap_rows(0, 1);
    /// ```
    fn swap_rows(&mut self, row1: usize, row2: usize) {
        if row1 == row2 {
            return;
        }

        self.data.swap(row1, row2);
    }

    /// Perform Gaussian Elimination on the given matrix (self)
    /// # Arguments
    ///
    /// * `gf` - Galois Field where the multiplication will occur.
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    /// use reed_solomon::matrix::Matrix;
    ///
    /// let matrix = Matrix::new_from_data(vec![vec![56, 23, 98], vec![3, 100, 200], vec![45, 201, 123]]);
    /// let gf8 = GaloisField::new();
    /// matrix.gauss_elim(gf8);
    /// ```
    fn gauss_elim(&mut self, gf: GaloisField) {
        // Clear out the part below the main diagonal and scale the main
        // diagonal to be 1.
        for r in 0..self.rows {
            // If the element on the diagonal is 0, find a row below
            // that has a non-zero and swap them.
            if self.data[r][r] == 0 {
                for r_below in r + 1..self.rows {
                    if self.data[r_below][r] != 0 {
                        self.swap_rows(r_below, r);
                        break;
                    }
                }
            }
            // If we couldn't find one, the matrix is singular.
            if self.data[r][r] == 0 {
                panic!("The given matrix is singular");
            }
            // Scale to 1.
            if self.data[r][r] != 1 {
                let scale = gf.div(1, self.data[r][r]);
                for c in 0..self.cols {
                    self.data[r][c] = gf.mul(self.data[r][c], scale)
                }
            }
            // Make everything below the 1 be a 0 by subtracting
            // a multiple of it.  (Subtraction and addition are
            // both exclusive or in the Galois field.)
            for r_below in r + 1..self.rows {
                if self.data[r_below][r] != 0 {
                    let scale = self.data[r_below][r];
                    for c in 0..self.cols {
                        let m = gf.mul(scale, self.data[r][c]);
                        self.data[r_below][c] = GaloisField::add(self.data[r_below][c], m);
                    }
                }
            }
        }
        // Now clear the part above the main diagonal.
        for d in 0..self.rows {
            for r_above in 0..d {
                if self.data[r_above][d] != 0 {
                    let scale = self.data[r_above][d];
                    for c in 0..self.cols {
                        let m = gf.mul(scale, self.data[d][c]);
                        self.data[r_above][c] = GaloisField::add(self.data[r_above][c], m);
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
    fn test_new_from_data() {
        let matrix = Matrix::new_from_data(vec![vec![1, 1, 1], vec![1, 1, 1]]);

        assert_eq!(matrix.rows, 2);
        assert_eq!(matrix.cols, 3);
        assert_eq!(matrix.data.len(), 2);
        assert_eq!(matrix.data[0].len(), 3);
        for row in matrix.data.iter() {
            for &elem in row.iter() {
                assert_eq!(1, elem);
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
    #[test]
    fn test_mul() {
        let gf8 = GaloisField::new();
        let left = Matrix::new_from_data(vec![vec![1, 2], vec![3, 4]]);
        let right = Matrix::new_from_data(vec![vec![5, 6], vec![7, 8]]);
        let res = left.mul(right, gf8);
        let exp_res: [[u8; 2]; 2] = [[11, 22], [19, 42]];

        assert_eq!(res.rows, 2);
        assert_eq!(res.cols, 2);
        assert_eq!(res.data.len(), 2);
        assert_eq!(res.data[0].len(), 2);
        for (row_index, row) in res.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }
    }
    #[test]
    fn test_new_swap() {
        let gf8 = GaloisField::new();
        let mut matrix = Matrix::new_vandermonde(3, 3, gf8);
        matrix.swap_rows(0, 1);
        let exp_res: [[u8; 3]; 3] = [[1, 1, 1], [1, 0, 0], [1, 2, 4]];

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
    fn test_invert() {
        let gf8 = GaloisField::new();
        let matrix = Matrix::new_from_data(vec![
            vec![56, 23, 98],
            vec![3, 100, 200],
            vec![45, 201, 123],
        ]);
        let res = matrix.invert(gf8);
        let exp_res: [[u8; 3]; 3] = [[175, 133, 33], [130, 13, 245], [112, 35, 126]];
        let iden = Matrix::new_identity(matrix.rows);

        for (row_index, row) in res.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(exp_res[row_index][col_index], elem);
            }
        }

        let mul = matrix.mul(res, gf8);

        for (row_index, row) in iden.data.iter().enumerate() {
            for (col_index, &elem) in row.iter().enumerate() {
                assert_eq!(mul.data[row_index][col_index], elem);
            }
        }
    }
}
