//! Invertible integer hash function described by Uriel Elias Wiebelitz in his masters thesis.
//! Works for u64 keys and produces an (address, fingerprint) tuple.
//! The size of the fingerprint depends on the size of the desired hash table:
//! |fp| = 64 - \floor(log2(|HT|)) bits
//!
//!
//! ```
//!
//! use bpht_eval::hash_function::*;
//! let key = 42;
//! let uneven_multiplier = 7;
//! let ht_size = 512;
//!
//! let (pos, fp) = inv_mult_hash(key, uneven_multiplier, ht_size);
//! let multiplicative_inverse = get_inverse_constant(uneven_multiplier);
//! let restored_key = invert_mult_hash(
//!     pos, fp, ht_size, multiplicative_inverse
//! );
//! assert_eq!(key, restored_key);
//! ```
//!
//! To generate a hash function family / universal hash function use
//! different uneven values for the uneven multiplier.
use mod_exp;
use rand::Rng;
use serde::{Deserialize, Serialize};

#[derive(Clone, Serialize, Deserialize)]
pub struct InvMultParams {
    pub a: u32,
}

impl InvMultParams {
    pub fn new() -> Self {
        // restrict the random value to be uneven and "large" so that the hash function scatters well
        let mut rng = rand::thread_rng();
        let mut mult = rng.gen_range(1_000_000, 4_294_967_295);
        if mult % 2 == 0 {
            mult += 1;
        }
        InvMultParams { a: mult }
    }

    pub fn with_params(a: u32) -> Self {
        InvMultParams { a }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct HlinParams {
    pub a: u64,
    pub b: u64,
}

impl HlinParams {
    pub fn new() -> Self {
        HlinParams {
            a: rand::random::<u64>(),
            b: rand::random::<u64>(),
        }
    }

    pub fn with_params(a: u64, b: u64) -> Self {
        HlinParams {
            a,
            b,
        }
    }
}

/// Compute an (address, fingerprint) tuple from a u64 key
///
/// m: (uneven) multiplier (for universal hashing)
/// n: hash table size (# slots)
///
/// NOTE: size of key universe u needs to be a power of 2
/// m needs to be uneven
pub fn inv_mult_hash(key: u64, m: u64, n: u64) -> (u64, u64) {
    let hv = key.wrapping_mul(m);
    let p = hv % n;
    let f = hv / n;
    (p, f)
}

/// Compute an inverse constant required to reverse the hash function
/// This needs to be computed once per multiplier.
pub fn get_inverse_constant(m: u64) -> u64 {
    let u = 2_u128.pow(64);
    let exp = (u / 2) - 1;
    (mod_exp::mod_exp(m as u128, exp, u)) as u64
}

/// Restore a key from an (address, fingerprint) tuple
/// using the hash table size and a multiplicative inverse
/// computed with the `get_inverse_constant` function.
pub fn invert_mult_hash(p: u64, f: u64, n: u64, inverse_constant: u64) -> u64 {
    (n * f + p).wrapping_mul(inverse_constant)
}

/// Compute an (address, fingerprint) tuple from a u32 key
///
/// NOTE: size of key universe u needs to be a power of 2
/// m needs to be uneven
pub fn inv_mult_hash_32(key: u32, m: u32, n: u32) -> (u32, u32) {
    let hv = key.wrapping_mul(m);
    let p = hv % n;
    let f = hv / n;
    (p, f)
}

pub fn simplified_inv_mult_hash_32(key: u32, a: u32) -> u32 {
    key.wrapping_mul(a)
}

/// Compute an inverse constant required to reverse the hash function
/// This needs to be computed once per multiplier.
pub fn get_inverse_constant_32(m: u32) -> u32 {
    let u = 2_u64.pow(32);
    let exp = (u / 2) - 1;
    (mod_exp::mod_exp(m as u64, exp, u)) as u32
}

/// Restore a key from an (address, fingerprint) tuple
/// using the hash table size and a multiplicative inverse
/// computed with the `get_inverse_constant` function.
pub fn invert_mult_hash_32(p: u32, f: u32, n: u32, inverse_constant: u32) -> u32 {
    (n * f + p).wrapping_mul(inverse_constant)
}

/// Compute the fingerprint length for a given parameter set
pub fn fp_length(n: u64) -> u32 {
    (2_u128.pow(64) as f64 / n as f64).log2().floor() as u32
}

/// Compute the fingerprint length for a given parameter set
pub fn fp_length_32(n: u32) -> u32 {
    (2_u64.pow(32) as f64 / n as f64).log2().floor() as u32
}

// pub fn get_hlin32_inverse_constant(m: u64) -> u64 {
//     let u = 2_u128.pow(64);
//     let exp = (u / 2) - 1;
//     (mod_exp::mod_exp(m as u128, exp, u)) as u64
// }

/// Compute an (address, fingerprint) tuple from a u32 key
///
/// NOTE: size of key universe u needs to be a power of 2
/// m needs to be uneven
pub fn hlin_32(key: u32, a: u64, b: u64) -> u32 {
    ((key as u64).wrapping_mul(a).wrapping_add(b) >> 32) as u32
}

// pub fn invert_hlin_32(key: u32, b: u64, inverse_constant: u64) -> u32 {
//     (key as u64).wrapping_sub(b).wrapping_mul(inverse_constant) as u32
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    /// Is hash funciton inversion working properly for a small example?
    fn hf_inversion_crafted() {
        let m = 3;
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let keys: Vec<u64> = (0..1000)
            .map(|_| rng.gen_range(0, (2_u128.pow(64) - 1) as u64))
            .collect();
        let n = 512;
        let inv_const = get_inverse_constant(m);
        for key in keys {
            let (address, fingerprint) = inv_mult_hash(key, m, n);
            let restored_key = invert_mult_hash(address, fingerprint, n, inv_const);
            assert_eq!(key, restored_key);
        }
    }

    #[test]
    /// Is hash function inversion working properly for random keys and sizes?
    fn hf_inversion_randomized() {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            let ht_size = rng.gen_range(8, 2_u64.pow(63)); // N
                                                           // generate an uneven multiplier
            let mut multiplier = rng.gen_range(1, (2_u128.pow(64) - 1) as u64);
            if multiplier % 2 == 0 {
                multiplier += 1;
            }

            // generate random keys
            let keys: Vec<u64> = (0..1000)
                .map(|_| rng.gen_range(0, (2_u128.pow(64) - 1) as u64))
                .collect();

            // Compute inverse constant for the uneven multiplier
            let inv_const = get_inverse_constant(multiplier);
            for key in keys {
                // Test restoring the keys
                let (address, fingerprint) = inv_mult_hash(key, multiplier, ht_size);
                let restored_key = invert_mult_hash(address, fingerprint, ht_size, inv_const);
                assert_eq!(key, restored_key);
            }
        }
    }

    #[test]
    /// Is hash function inversion working properly for random keys and sizes with 32 bit?
    fn hf_inversion_randomized_32() {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            let ht_size = rng.gen_range(8, 2_u32.pow(31)); // N
                                                           // generate an uneven multiplier
            let mut multiplier = rng.gen_range(1, (2_u64.pow(32) - 1) as u32);
            if multiplier % 2 == 0 {
                multiplier += 1;
            }

            // generate random keys
            let keys: Vec<u32> = (0..1000)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            // Compute inverse constant for the uneven multiplier
            let inv_const = get_inverse_constant_32(multiplier);
            for key in keys {
                // Test restoring the keys
                let (address, fingerprint) = inv_mult_hash_32(key, multiplier, ht_size);
                let restored_key =
                    invert_mult_hash_32(address, fingerprint, ht_size, inv_const);
                assert_eq!(key, restored_key);
            }
        }
    }

    #[test]
    /// Test if the fingerprint length holds up for 64 bit values.
    fn fingerprint_length() {
        use rand::Rng;
        let mut rng = rand::thread_rng();

        for _ in 0..1000 {
            let mut longest_fingerprint = 0;
            let ht_power = rng.gen_range(24, 32);
            let ht_size = 2_u64.pow(ht_power); // N

            let computed_length = fp_length(ht_size as u64);

            // generate an uneven multiplier
            let mut multiplier = rng.gen_range(1, (2_u128.pow(64) - 1) as u64);
            if multiplier % 2 == 0 {
                multiplier += 1;
            }

            // generate random keys
            let keys: Vec<u64> = (0..1000)
                .map(|_| rng.gen_range(0, (2_u128.pow(64) - 1) as u64))
                .collect();

            // Compute inverse constant for the uneven multiplier
            let _inv_const = get_inverse_constant(multiplier);
            for key in keys {
                // Test restoring the keys
                let (_address, fingerprint) = inv_mult_hash(key, multiplier, ht_size);
                let fp_length = 64 - fingerprint.leading_zeros();
                if fp_length > longest_fingerprint {
                    longest_fingerprint = fp_length;
                }
            }
            assert!(longest_fingerprint <= computed_length);
        }
    }

    #[test]
    /// Test if the fingerprint length holds up for 32 bit values.
    fn fingerprint_length_32() {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            let mut longest_fingerprint = 0;
            let ht_power = rng.gen_range(5, 32);
            let ht_size = 2_u32.pow(ht_power); // N

            let computed_length = fp_length_32(ht_size);

            // generate an uneven multiplier
            let mut multiplier = rng.gen_range(1, (2_u64.pow(32) - 1) as u32);
            if multiplier % 2 == 0 {
                multiplier += 1;
            }

            // generate random keys
            let keys: Vec<u32> = (0..1000)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            // Compute inverse constant for the uneven multiplier
            let inv_const = get_inverse_constant_32(multiplier);
            for key in keys {
                // Test restoring the keys
                let (address, fingerprint) = inv_mult_hash_32(key, multiplier, ht_size);
                let _restored_key =
                    invert_mult_hash_32(address, fingerprint, ht_size, inv_const);
                let fp_length = 32 - fingerprint.leading_zeros();
                if fp_length > longest_fingerprint {
                    longest_fingerprint = fp_length;
                }
                assert!(longest_fingerprint <= computed_length);
            }
        }
    }

    // hlin cannot be (easily) inverted, since bits are lost during the >> 32 step

    // #[test]
    // /// Is hash function inversion working properly for random keys and sizes with 32 bit?
    // fn hlin_inversion_randomized_32() {
    //     use rand::Rng;
    //     let mut rng = rand::thread_rng();
    //     for _ in 0..1000 {

    //         let mut multiplier = rng.gen_range(1, (2_u128.pow(64) - 1) as u64);
    //         let offset = rng.gen_range(1, (2_u128.pow(74) - 1) as u64);
    //         if multiplier % 2 == 0 {
    //             multiplier += 1;
    //         }

    //         // generate random keys
    //         let keys: Vec<u32> = (0..1000)
    //             .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
    //             .collect();

    //         // Compute inverse constant for the uneven multiplier
    //         let inv_const = get_hlin32_inverse_constant(multiplier);
    //         for key in keys {
    //             // Test restoring the keys
    //             let hash = hlin_32(key, multiplier, offset);
    //             let restored_key = invert_hlin_32(hash, offset, inv_const);
    //             assert_eq!(key, restored_key);
    //         }
    //     }
    // }
}
