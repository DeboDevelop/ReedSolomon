## Reed Solomon

This is an implementation of Reed Solomon Erasure Coding based on [Backblaze's Java implementation](https://github.com/Backblaze/JavaReedSolomon) and loosely on [Klaus Post's Go implementation](https://github.com/klauspost/reedsolomon). Reed-Solomon erasure coding is an error correction technique used to recover lost or corrupted data in storage or transmission systems. It works by performing matrix arithmetic over Galois Field. You can get an overview of Reed Solomon in [this blog post](https://www.backblaze.com/blog/reed-solomon/). You can get more details in [this paper](http://web.eecs.utk.edu/~jplank/plank/papers/SPE-9-97.html) by James S. Plank as well as [the correction paper](http://web.eecs.utk.edu/~jplank/plank/papers/CS-03-504.html) by James S. Plank and Ying Ding.

### Usage Example

```
extern crate reed_solomon;
use reed_solomon::ReedSolomon;

// Creating a new ReedSolomon Struct. First parameter is 
// no. of data shard and 2nd parameter is no. of parity shard.
let result = ReedSolomon::new(2, 2);
let rs = match result {
    Ok(x) => x,
    Err(e) => panic!("{}", e),
};
// We can also read files and put them as shard. Size of all shards should be same.
let shards = vec![
    // data shards
    vec![0, 1, 2],
    vec![3, 4, 5],
    // parity shards (to be overwritten)
    vec![0, 0, 0],
    vec![0, 0, 0],
];
// Encoding shards
let encoded_shard_result = rs.encode(shards.clone());
let encoded_shard = match encoded_shard_result {
    Ok(x) => x,
    Err(e) => panic!("{}", e),
};
// Reconstructing missing shards
let broken_shards = vec![vec![0, 1, 2], vec![], vec![], vec![5, 14, 11]];
let decoded_shard_result = rs.decode(broken_shards);
let decoded_shard = match decoded_shard_result {
    Ok(x) => x,
    Err(e) => panic!("{}", e),
};
```

### Special Thanks To

Backblaze

Klaus Post

James S. Plank

Ying Ding

### License

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details.

### Author

[Debajyoti Dutta](https://github.com/DeboDevelop)