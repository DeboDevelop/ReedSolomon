use std::fmt;

#[derive(Debug)]
pub enum Error {
    RowsMustMatch(usize, usize),
    RowColMustMatch(usize, usize),
    NonSquareMatrix,
    SingularMatrix,
    ZeroDataShards,
    ZeroParityShards,
    ShardsOverflow,
    WrongNoOfShards,
    EmptyShards,
    InconsistentShards,
    TooFewShards,
    TooManyShards,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::RowsMustMatch(left_rows, right_rows) => write!(
                f,
                "Row count of the matrices must match. Current row count, left: {}, right: {}",
                *left_rows, *right_rows
            ),
            Error::RowColMustMatch(left_cols, right_rows) => write!(
                f,
                "Column count on left has to be same as row count on right. left column: {}, right row: {}",
                *left_cols, *right_rows
            ),
            Error::NonSquareMatrix =>  write!(f, "The given matrix is non-square matrix and they are not invertible"),
            Error::SingularMatrix =>  write!(f, "The given matrix is singular"),
            Error::ZeroDataShards =>  write!(f, "Data Shards can't be zero"),
            Error::ZeroParityShards =>  write!(f, "Parity Shards can't be zero"),
            Error::ShardsOverflow =>  write!(f, "More than 256 shards are nt allowed"),
            Error::WrongNoOfShards =>  write!(f, "Wrong no. of shards"),
            Error::EmptyShards =>  write!(f, "There is a empty shards"),
            Error::InconsistentShards =>  write!(f, "Length of the given shards are different"),
            Error::TooFewShards =>  write!(f, "Too few no. of shards"),
            Error::TooManyShards =>  write!(f, "Too many no. of shards"),
        }
    }
}
