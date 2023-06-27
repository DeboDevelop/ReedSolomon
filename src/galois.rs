/// This size of the field i.e. 2^8.
const FIELD_SIZE: usize = 256;

/// Size of the exponent table. The highest log value is FIELD_SIZE - 2
/// so we decreased 2 from FIELD_SIZE. The table is repeated a 2nd time
/// so fn mul() doesn't have to check bounds.
const EXP_TABLE_SIZE: usize = FIELD_SIZE * 2 - 2;

/// The irreducible polynomial which is used to generate the log and
/// exp table. The possibilities are 29, 43, 45, 77, 95, 99, 101, 105,
/// 113, 135, 141, 169, 195, 207, 231, and 245.
const IRREDUCIBLE_POLYNOMIAL: usize = 29;

pub struct GaloisField {
    field_size: usize,
    irre_poly: usize,
    exp_table_size: usize,
    log_table: [u8; FIELD_SIZE],
    exp_table: [u8; EXP_TABLE_SIZE],
}

/// Generate the log table given an irreducible polynomial which maps
/// the elements of the Galois field to their discrete logarithm. Since there is
/// no log for 0 so the entry in 0th index can be ignored.
/// # Arguments
///
/// * `irre_poly` - An irreducible polynomial for GF(2^8)
///
/// # Example
/// ```
/// use reed_solomon::galois::gen_log_table;
///
/// let log_table = gen_log_table(29);
/// ```
pub fn gen_log_table(irre_poly: usize) -> [u8; FIELD_SIZE] {
    let mut res = [0 as u8; FIELD_SIZE];
    // Primitive element
    let mut b: usize = 1;

    for log in 0..FIELD_SIZE - 1 {
        res[b] = log as u8;

        // raising power of the element
        b = b << 1;

        // modulo the element so that it remain inside the field
        if FIELD_SIZE <= b {
            b = (b - FIELD_SIZE) ^ irre_poly;
        }
    }

    res
}
/// Generate the exp table given the log table. The exp table maps logarithms
/// to elements of the Galois field.
/// # Arguments
///
/// * `log_table` - The log table for GF(2^8)
///
/// # Example
/// ```
/// use reed_solomon::galois::gen_exp_table;
/// use reed_solomon::galois::gen_log_table;
///
/// let log_table = gen_log_table(29);
/// let exp_table = gen_exp_table(&log_table);
/// ```
pub fn gen_exp_table(log_table: &[u8; FIELD_SIZE]) -> [u8; EXP_TABLE_SIZE] {
    let mut res = [0 as u8; EXP_TABLE_SIZE];

    for i in 1..FIELD_SIZE {
        let log = log_table[i] as usize;
        res[log] = i as u8;
        // Populating the repeated table
        res[log + FIELD_SIZE - 1] = i as u8;
    }

    res
}

impl GaloisField {
    /// Create a new GaloisField(2^8)
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    ///
    /// let gf8 = GaloisField::new();
    /// ```
    pub fn new() -> GaloisField {
        let log_table = gen_log_table(IRREDUCIBLE_POLYNOMIAL);
        let exp_table = gen_exp_table(&log_table);

        GaloisField {
            field_size: FIELD_SIZE,
            irre_poly: IRREDUCIBLE_POLYNOMIAL,
            exp_table_size: EXP_TABLE_SIZE,
            log_table,
            exp_table,
        }
    }

    /// Adds 2 elements in the field.
    /// # Arguments
    ///
    /// * `a` - First element to be added
    /// * `b` - Second element to be added
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    ///
    /// let res = GaloisField::add(1, 1);
    /// ```
    pub fn add(a: u8, b: u8) -> u8 {
        a ^ b
    }

    /// Subtract 1 element from another in the field.
    /// # Arguments
    ///
    /// * `a` - Minuend
    /// * `b` - Subtrahend
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    ///
    /// let res = GaloisField::sub(1, 1);
    /// ```
    pub fn sub(a: u8, b: u8) -> u8 {
        a ^ b
    }

    /// Multiplies 2 elements in the field.
    /// # Arguments
    ///
    /// * `log_table` - The log table for GF(2^8)
    /// * `exp_table` - The exp table for GF(2^8)
    /// * `a` - First element to be multiplied
    /// * `b` - Second element to be multiplied
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    /// 
    /// let gf8 = GaloisField::new();
    /// let res = gf8.mul(1, 1);
    /// ```
    pub fn mul(&self, a: u8, b: u8) -> u8 {
        if a == 0 || b == 0 {
            0
        } else {
            let log_a = self.log_table[a as usize];
            let log_b = self.log_table[b as usize];
            let log_res = log_a as usize + log_b as usize;
            self.exp_table[log_res]
        }
    }

    /// Computes a^n in Galois field
    /// # Arguments
    ///
    /// * `a` - Base element
    /// * `b` - Exponent element
    ///
    /// # Example
    /// ```
    /// use reed_solomon::galois::GaloisField;
    /// 
    /// let gf8 = GaloisField::new();
    /// let res = gf8.exp(2, 2);
    /// ```
    pub fn exp(
        &self,
        a: u8,
        n: usize,
    ) -> u8 {
        if n == 0 {
            1
        } else if a == 0 {
            0
        } else {
            let log_a = self.log_table[a as usize];
            let mut log_res = log_a as usize * n;
            while 255 <= log_res {
                log_res -= 255;
            }
            self.exp_table[log_res]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    static EXPECTED_LOG_RES: [u8; FIELD_SIZE] = [
        0, 0, 1, 25, 2, 50, 26, 198, 3, 223, 51, 238, 27, 104, 199, 75, 4, 100, 224, 14, 52, 141,
        239, 129, 28, 193, 105, 248, 200, 8, 76, 113, 5, 138, 101, 47, 225, 36, 15, 33, 53, 147,
        142, 218, 240, 18, 130, 69, 29, 181, 194, 125, 106, 39, 249, 185, 201, 154, 9, 120, 77,
        228, 114, 166, 6, 191, 139, 98, 102, 221, 48, 253, 226, 152, 37, 179, 16, 145, 34, 136, 54,
        208, 148, 206, 143, 150, 219, 189, 241, 210, 19, 92, 131, 56, 70, 64, 30, 66, 182, 163,
        195, 72, 126, 110, 107, 58, 40, 84, 250, 133, 186, 61, 202, 94, 155, 159, 10, 21, 121, 43,
        78, 212, 229, 172, 115, 243, 167, 87, 7, 112, 192, 247, 140, 128, 99, 13, 103, 74, 222,
        237, 49, 197, 254, 24, 227, 165, 153, 119, 38, 184, 180, 124, 17, 68, 146, 217, 35, 32,
        137, 46, 55, 63, 209, 91, 149, 188, 207, 205, 144, 135, 151, 178, 220, 252, 190, 97, 242,
        86, 211, 171, 20, 42, 93, 158, 132, 60, 57, 83, 71, 109, 65, 162, 31, 45, 67, 216, 183,
        123, 164, 118, 196, 23, 73, 236, 127, 12, 111, 246, 108, 161, 59, 82, 41, 157, 85, 170,
        251, 96, 134, 177, 187, 204, 62, 90, 203, 89, 95, 176, 156, 169, 160, 81, 11, 245, 22, 235,
        122, 117, 44, 215, 79, 174, 213, 233, 230, 231, 173, 232, 116, 214, 244, 234, 168, 80, 88,
        175,
    ];

    static EXPECTED_EXP_RES: [u8; EXP_TABLE_SIZE] = [
        1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117,
        234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181,
        119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222,
        161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231,
        211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175, 67, 134, 17,
        34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59,
        118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66,
        132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115,
        230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219,
        171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141,
        7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155,
        43, 86, 172, 69, 138, 9, 18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139,
        11, 22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142,
        // Repeat the table a second time
        1, 2, 4, 8, 16, 32, 64, 128, 29, 58, 116, 232, 205, 135, 19, 38, 76, 152, 45, 90, 180, 117,
        234, 201, 143, 3, 6, 12, 24, 48, 96, 192, 157, 39, 78, 156, 37, 74, 148, 53, 106, 212, 181,
        119, 238, 193, 159, 35, 70, 140, 5, 10, 20, 40, 80, 160, 93, 186, 105, 210, 185, 111, 222,
        161, 95, 190, 97, 194, 153, 47, 94, 188, 101, 202, 137, 15, 30, 60, 120, 240, 253, 231,
        211, 187, 107, 214, 177, 127, 254, 225, 223, 163, 91, 182, 113, 226, 217, 175, 67, 134, 17,
        34, 68, 136, 13, 26, 52, 104, 208, 189, 103, 206, 129, 31, 62, 124, 248, 237, 199, 147, 59,
        118, 236, 197, 151, 51, 102, 204, 133, 23, 46, 92, 184, 109, 218, 169, 79, 158, 33, 66,
        132, 21, 42, 84, 168, 77, 154, 41, 82, 164, 85, 170, 73, 146, 57, 114, 228, 213, 183, 115,
        230, 209, 191, 99, 198, 145, 63, 126, 252, 229, 215, 179, 123, 246, 241, 255, 227, 219,
        171, 75, 150, 49, 98, 196, 149, 55, 110, 220, 165, 87, 174, 65, 130, 25, 50, 100, 200, 141,
        7, 14, 28, 56, 112, 224, 221, 167, 83, 166, 81, 162, 89, 178, 121, 242, 249, 239, 195, 155,
        43, 86, 172, 69, 138, 9, 18, 36, 72, 144, 61, 122, 244, 245, 247, 243, 251, 235, 203, 139,
        11, 22, 44, 88, 176, 125, 250, 233, 207, 131, 27, 54, 108, 216, 173, 71, 142,
    ];

    #[test]
    fn test_gen_log_table() {
        let res = gen_log_table(IRREDUCIBLE_POLYNOMIAL);
        for i in 0..FIELD_SIZE {
            assert_eq!(EXPECTED_LOG_RES[i], res[i]);
        }
    }
    #[test]
    fn test_gen_exp_table() {
        let log_table = gen_log_table(IRREDUCIBLE_POLYNOMIAL);
        let res = gen_exp_table(&log_table);
        for i in 0..EXP_TABLE_SIZE {
            assert_eq!(EXPECTED_EXP_RES[i], res[i]);
        }
    }
    #[test]
    fn test_gf_new() {
        let gf8 = GaloisField::new();
        assert_eq!(FIELD_SIZE, gf8.field_size);
        assert_eq!(IRREDUCIBLE_POLYNOMIAL, gf8.irre_poly);
        assert_eq!(EXP_TABLE_SIZE, gf8.exp_table_size);
        for i in 0..FIELD_SIZE {
            assert_eq!(EXPECTED_LOG_RES[i], gf8.log_table[i]);
        }
        for i in 0..EXP_TABLE_SIZE {
            assert_eq!(EXPECTED_EXP_RES[i], gf8.exp_table[i]);
        }
    }
    #[test]
    fn test_add() {
        assert_eq!(0, GaloisField::add(1, 1));
        assert_eq!(26, GaloisField::add(21, 15));
        assert_eq!(68, GaloisField::add(120, 60));
    }
    #[test]
    fn test_sub() {
        assert_eq!(0, GaloisField::sub(1, 1));
        assert_eq!(26, GaloisField::sub(21, 15));
        assert_eq!(68, GaloisField::sub(120, 60));
    }
    #[test]
    fn test_mul() {
        let gf8 = GaloisField::new();
        assert_eq!(12, gf8.mul(3, 4));
        assert_eq!(21, gf8.mul(7, 7));
        assert_eq!(41, gf8.mul(23, 45));
    }
    #[test]
    fn test_exp() {
        let gf8 = GaloisField::new();
        assert_eq!(4, gf8.exp(2, 2));
        assert_eq!(235, gf8.exp(5, 20));
        assert_eq!(43, gf8.exp(13, 7));
    }
}
