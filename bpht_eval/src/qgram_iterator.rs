//! This module contains a 2-bit q-gram encoder and enums to modulate its behaviour.
//!
//! The encoder itself is an iterator generated from:
//!
//!   * A reference to a text
//!   * A value for q (1 <= q <= 32)
//!   * A canonicity function (selected using calues of the Canonical enum)
//!   * A hash mixing function (selected using calues of the HashMixing enum)
//!
//! A q-gram is encoded in up to 64 bits and returned as Option<u64>.
//! The 2-bit encoding can only encode A, C, G, and T, hence all other
//! characters cause an invalid q-gram. This is signified by a returned None calue.
//! Note that this makes the item type of the iterator Option<u64> and consequently
//! the iterator returns items of the type Option<Option<u64>>, where:
//!
//!   * None signifies the end of the sequence and terminates the iterator
//!   * Some(None) signifies the end of a valid subsequence of q-grams in the text
//!   * Some(Some(g)) contains the valid q-gram g.

use crate::hash_function;
use crate::hash_function::{InvMultParams, HlinParams};
use bincode::{deserialize_from, serialize_into};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};

/// Return true if the given base is in ACGT or in acgt.
fn is_acgt(base: &u8) -> bool {
    [65, 97, 67, 99, 71, 103, 84, 116].contains(base)
}

/// Used to keep tabs on the state of the q-gram encoder.
/// It enteres an invalid state if a non-ACGT character is read.
pub enum State {
    Valid,
    Invalid,
}

/// Used to pick a canonicity function. The effect of the function
/// is implemented in the result function of the QGrams struct.
/// No means no canonical q-gram is computed, i.e. the q-gram seen
/// in the text is returned directly. Min/max mean the minimum/ maximum
///  of the q-gram and its reverse complement are returned.
#[derive(Clone)]
pub enum Canonical {
    No,
    Min,
    Max,
}

/// Hash mixing functions used to increase the entropy in the q-grams
/// before computing caonical q-grams. All hash mixing functions are implemented in
/// the hash_function module.
#[derive(Clone, Serialize, Deserialize)]
pub enum HashFunction {
    No,
    Hlin(HlinParams),
    InvMult(InvMultParams),
    Tab32Simple(tab_hash::Tab32Simple),
    Tab32Twisted(tab_hash::Tab32Twisted),
}

impl HashFunction {
    pub fn load(path: &str) -> HashFunction {
        let loaded_hf: HashFunction;
        {
            let mut f =
                BufReader::new(File::open(path).unwrap_or_else(|_| {
                    panic!("Opening the saved hf {} did not work", path)
                }));
            loaded_hf = deserialize_from(&mut f).unwrap_or_else(|_| {
                panic!(
                    "Deserializing the hf file {} did not work",
                    path
                )
            });
        }
        loaded_hf
    }

    pub fn save(&mut self, path: &str) {
        eprintln!("Saving hash function to {}", path);
        let mut f = BufWriter::new(File::create(path).unwrap_or_else(|_| {
            panic!(
                "Opening file to save hf at {} did not work. Check that the path exists.",
                path
            )
        }));
        serialize_into(&mut f, self).expect("Serialization of the hf for saving did not work.");
    }
}

/// Iterator over 2-bit encoded q-grams
pub struct QGrams<'a> {
    seq: &'a [u8],
    q: usize,
    canonical: Canonical,
    mask: u32,
    pos: usize,
    len: usize,
    last_qgram_fw: Option<u32>,
    last_qgram_rc: Option<u32>,
    state: State,
    hash_function: HashFunction,
    canonisize_qgrams: bool,
}

impl<'a> QGrams<'a> {
    /// Create a new iterator over 2-bit encoded q-grams with a canonicity function
    /// and a hash mixing function. Note that q has to be 1 <= q <= 32.
    pub fn new(
        seq: &'a [u8],
        q: usize,
        canonical: Canonical,
        hash_function: HashFunction,
        canonisize_qgrams: bool,
    ) -> Self {
        if q == 0 || q > 16 || seq.len() < q {
            panic!(
                "q has to be at most 16, at least 1, and smaller than the sequence length ({}).",
                seq.len()
            );
        }
        QGrams {
            seq,
            q,
            canonical,
            mask: (2_u64.pow(2 * q as u32) - 1) as u32,
            pos: 0,
            len: seq.len(),
            last_qgram_fw: None,
            last_qgram_rc: None,
            state: State::Valid,
            hash_function,
            canonisize_qgrams,
        }
    }

    /// Encode the qgram starting at pos.
    /// Return the qgram and its rc, if possible, return Err if non ACGT is present.
    /// This is used to initially encode a q-gram either at the start of a
    /// sequence or after an invalid subsequence has been skipped.
    ///
    /// Note: This does not modify the internal state of the iterator, but only
    /// returns the encoded q-grams for forward and reverse complement.
    fn encode(&self, pos: usize) -> Result<(u32, u32), ()> {
        let mut qgram_fw: u32 = 0;
        let mut qgram_rc: u32 = 0;

        // Encode the forward q-gram
        for base in self.seq[pos..pos + self.q].iter() {
            qgram_fw = qgram_fw.overflowing_shl(2).0; // Shift left (possibly shift out the two bits no longer needed)
            qgram_fw |= match base {
                65 | 97 => 0b00,
                67 | 99 => 0b01,
                71 | 103 => 0b10,
                84 | 116 => 0b11,
                _ => return Err(()),
            }
        }

        // Encode the reverse complement q-gram
        for base in self.seq[pos..pos + self.q].iter().rev() {
            qgram_rc = qgram_rc.overflowing_shl(2).0;
            qgram_rc |= match base {
                65 | 97 => 0b11,  // Complement
                67 | 99 => 0b10,  // Complement
                71 | 103 => 0b01, // Complement
                84 | 116 => 0b00, // Complement
                _ => panic!(
                    "Encountered non ACGT char. This shouldn't happen here, but in the fw case!"
                ),
            }
        }
        Ok((qgram_fw & self.mask, qgram_rc & self.mask))
    }

    /// Keep both forward and reverse complement kmer up to date after
    /// advancing one position/ base.
    ///
    /// Seq: ATGCA
    /// Pos   î
    /// q: 4
    ///
    /// fw:
    ///     [ A  T  G  C] A
    ///      00 11 10 01
    ///     << shift out A
    ///     [ T  G  C  -] A
    ///      11 10 01 ??
    ///     add in A
    ///     [ T  G  C  A]
    ///      11 10 01 00
    ///
    /// rc:
    ///    T [ G  C  A  T]
    ///       10 01 00 11
    ///    >> shift out T to the right
    ///    T [ -  G  C  A]
    ///       ?? 10 01 00
    ///    << shift 11 for new T to the left to the fourth base position
    ///    T [ -  G  C  A]
    ///       ?? 10 01 00
    ///       11
    ///      add in T
    ///      [ T  G  C  A]
    ///       11 10 01 00
    ///
    /// Note that this does not modify the internal state
    fn update(&self) -> Result<(u32, u32), ()> {
        let mut qgram_fw;
        let qgram_rc;
        if let Some(old_qgram_fw) = self.last_qgram_fw {
            // shift out leftmost base
            qgram_fw = old_qgram_fw.overflowing_shl(2).0 & self.mask;
            // add in new base
            qgram_fw |= match self.seq[self.pos + self.q - 1] {
                65 | 97 => 0b00,
                67 | 99 => 0b01,
                71 | 103 => 0b10,
                84 | 116 => 0b11,
                _ => return Err(()), // Shifted in charcter is non ACGT
            }
        } else {
            panic!("Ran update on uninitialized fw qgram");
        }

        if let Some(mut old_qgram_rc) = self.last_qgram_rc {
            let offset = 2 * self.q as u32 - 2;
            // shift out old base
            old_qgram_rc >>= 2;
            // add in the new base, move to the leftmost side since it is RC
            old_qgram_rc |= match self.seq[self.pos + self.q - 1] {
                65 | 97 => 0b11 << offset,  // Complement
                67 | 99 => 0b10 << offset,  // Complement
                71 | 103 => 0b01 << offset, // Complement
                84 | 116 => 0b00 << offset, // Complement
                _ => panic!(
                    "Encountered non ACGT char. This shouldn't happen here, but in the fw case!"
                ),
            };
            qgram_rc = old_qgram_rc;
        } else {
            panic!("Ran update on uninitialized rc qgram");
        }
        Ok((qgram_fw, qgram_rc))
    }

    /// Find the next valid q-gram, i.e. the next full q-gram that does not
    /// contain any non-ACGT characters. If none exists, return None.
    fn find_next_valid_qgram(&self) -> Option<usize> {
        let mut last_invalid = None;

        // Find the last/ rightmost invalid position in the offending q-gram
        for (i, base) in self.seq[self.pos..self.pos + self.q].iter().enumerate() {
            if !is_acgt(base) {
                last_invalid = Some(self.pos + i);
            }
        }

        if last_invalid.is_some() {
            // Advance the q-gram through all remaining valid starting positions
            // and monitor if an invalid character is still present
            for qgram_start_pos in self.pos..self.seq.len() - self.q + 1 {
                let new_base_pos = qgram_start_pos + self.q - 1;
                // Look at the last (rightmost) base in the new q-gram
                let base = self.seq[new_base_pos];

                // If the new base is non-ACGT (i.e. invalid), update last invalid position.
                // This means longer searching is required.
                if !is_acgt(&base) {
                    last_invalid = Some(new_base_pos);
                } else {
                    // Check if the last invalid position has been pushed out by this step
                    if let Some(last_invalid) = last_invalid {
                        if last_invalid < qgram_start_pos && new_base_pos < self.len {
                            // found a new staring position that is completely valid
                            // and contained within the sequence
                            return Some(qgram_start_pos);
                        }
                    } else {
                        panic!("Last invalid was None. Something is flawed.");
                    }
                }
            }
            // If the above loop completely passes without returning,
            // then no new valid q-gram can be found.
            // This None is interpreted in the next method to terminate the iterator
            None
        } else {
            panic!("Did not find invalid character in first window. Something is flawed.");
        }
    }

    /// Compute the returned q-gram using the specified canonicity function
    /// and the hash mixing function.
    fn result(&self) -> Option<Option<u32>> {
        if let (Some(fw), Some(rc)) = (self.last_qgram_fw, self.last_qgram_rc) {
            let mixed_fw;
            let mixed_rc;

            // Apply hash mixing to both forward and reverse q-gram
            match self.hash_function {
                HashFunction::No => {
                    mixed_fw = fw;
                    mixed_rc = rc;
                }
                HashFunction::Hlin(ref params) => {
                    mixed_fw = hash_function::hlin_32(fw, params.a, params.b);
                    mixed_rc = hash_function::hlin_32(rc, params.a, params.b);
                }
                HashFunction::Tab32Simple(ref tab) => {
                    mixed_fw = tab.hash(fw);
                    mixed_rc = tab.hash(rc);
                }
                HashFunction::Tab32Twisted(ref tab) => {
                    mixed_fw = tab.hash(fw);
                    mixed_rc = tab.hash(rc);
                }
                HashFunction::InvMult(ref params) => {
                    mixed_fw = hash_function::simplified_inv_mult_hash_32(fw, params.a);
                    mixed_rc = hash_function::simplified_inv_mult_hash_32(rc, params.a);
                }
            }

            if self.canonisize_qgrams {
                // canonical => smallest of two q-rgams
                match self.canonical {
                    Canonical::No => Some(Some(mixed_fw)),
                    Canonical::Min => {
                        if fw <= rc {
                            Some(Some(mixed_fw))
                        } else {
                            Some(Some(mixed_rc))
                        }
                    }
                    Canonical::Max => {
                        if fw <= rc {
                            Some(Some(mixed_rc))
                        } else {
                            Some(Some(mixed_fw))
                        }
                    }
                }
            } else {
                // canonical => smallest of two hashes
                // Pick and return canonical q-gram based on the hashmixed q-grams
                match self.canonical {
                    Canonical::No => Some(Some(mixed_fw)),
                    Canonical::Min => Some(Some(mixed_fw.min(mixed_rc))),
                    Canonical::Max => Some(Some(mixed_fw.max(mixed_rc))),
                }
            }
        } else {
            None
        }
    }
}

/// Implements and handles all the plumbing for the iterator trait.
/// The item type is Option<u32> to encode invalid subsequences without
/// terminating the iterator.
impl<'a> Iterator for QGrams<'a> {
    type Item = Option<u32>;

    /// Call the encode/ update update methods of
    /// the iterator. The the result is valid, return a q-gram.
    /// If not, manage invalid sequences until either a valid
    /// q-gram turns up, or the sequence ends.
    fn next(&mut self) -> Option<Self::Item> {
        // If the state is valid, i.e. the last read base was an ACGT character,
        // try to update. If that fails, return None, do not advance the position.
        // During the following call to next the Invalid case will then try to find
        // a valid q-gram or terminate the iteration.
        match self.state {
            State::Valid => {
                // handle first qgram
                if self.pos == 0 {
                    match self.encode(0) {
                        Ok((fw, rc)) => {
                            // successfully encoded
                            self.last_qgram_fw = Some(fw);
                            self.last_qgram_rc = Some(rc);
                            self.pos += 1;
                            self.result()
                        }
                        Err(()) => {
                            // Return None and search for the next valid q-gram in the next step
                            self.state = State::Invalid;
                            Some(None)
                        }
                    }
                } else if self.pos <= self.len - self.q {
                    // handle all other valid q-grams
                    match self.update() {
                        Ok((fw, rc)) => {
                            // successfully updated
                            self.last_qgram_fw = Some(fw);
                            self.last_qgram_rc = Some(rc);
                            self.pos += 1;
                            self.result()
                        }
                        Err(()) => {
                            // new base was non ACGT, return None and search for next valid q-gram in next step
                            self.state = State::Invalid;
                            Some(None)
                        }
                    }
                } else {
                    None // Sequence consumed. Finish iteration.
                }
            }

            // Handle the invalid state. Note that the position has not advanced
            // from the invalid position yet.
            State::Invalid => {
                if self.pos <= self.len - self.q {
                    // Recover from invalid state, i.e. a non ACGT character
                    match self.find_next_valid_qgram() {
                        Some(pos) => {
                            // Found a new valid q-gram starting position
                            self.pos = pos;
                            match self.encode(pos) {
                                Ok((fw, rc)) => {
                                    // successfully encoded, 
                                    self.last_qgram_fw = Some(fw);
                                    self.last_qgram_rc = Some(rc);
                                    self.pos += 1;
                                    self.state = State::Valid;
                                    self.result()
                                },
                                Err(()) => panic!("After finding a valid qg, encoding should never fail. Error at position {}", self.pos),

                            }
                        }
                        None => None, // cannot find any more valid q-grams, sequence consumed
                    }
                } else {
                    None // Sequence consumed. Finish iteration.
                }
            }
        }
    }
}
