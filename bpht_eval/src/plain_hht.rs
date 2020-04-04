//! This is not a doctest
//!
//! Computing address and fp from keys:
//!
//! key: 32-bit
//! | 32 - fp_len address-bits         | floor(log2(|ht|)) fingerprint-bits  |
//! | power of 2 of the address space  |                                     |
//!
//! Bitpacking hash table entries:
//!
//! |             32 bit payload     |     ~24 bit fingerprint    |~ 8 hop bits|
//! |pppppppppppppppppppppppppppppppp|ffffffffffffffffffffffff...<|>...hhhhhhhh|
//!
//! Hop bits: 0 means free, 1 means filled
//! Hop bits are read from left to right
//!      1011 means, this position is filled (___1),
//!      as is the next (__1_) and the one with offset 3 (1___)
//!
//! The number of fingerprint bits depends on the size of the hash table.
//!
//! Multiple entries possible:
//! ht[42] = 23
//! ht[42] = 17
//! should r [23, 17]?
use bincode::{deserialize_from, serialize_into, serialized_size};
use bincode;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};


#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PlainHHT {
    pub h: usize,
    pub u: usize,
    pub table: Vec<u64>,
    pub hop_array: Vec<u32>,
    pub hop_bits_mask: u64,
    pub fingerprint_mask: u64,
    pub fingerprint_shift: usize,
    pub admin_bits: u64,
    pub fp_bits_in_key: usize,
    pub key_fp_mask: u32,
    pub in_resize: bool,
    pub allow_resize: bool,
}

impl PlainHHT {
    /// Create a new hopscotch hash table using
    /// a page size of h and u + h - 1 slots in total.
    /// the last h-1 slots are used to keep overflowing entries
    /// from the last position without wraparound
    ///
    /// Note that u > h is required, to cram all bits used for
    /// administration into 32 bits.
    pub fn new(h: usize, u: usize, allow_resize: bool) -> Result<Self, &'static str> {
        // Make sure that we can work with 32-bit values here.
        // Not sure if this is helpful as it restricts the HT size.
        //
        // TODO validate if fp_length_32 and fp_length_64 are a sensible differentiation
        assert!(u <= 2_u64.pow(32) as usize);
        // assure that u is a power of 2
        if u.count_ones() != 1 {
            return Err("The parameter u is not a power of 2");
        }

        // Get the number of fingerprint bits required for a hash table
        // of size u
        let required_fp_bits = PlainHHT::fp_length_for(u as u32);

        // make sure that the sum of fingerprint bits and hop bits
        // does not exceed the 32 bits allocated for them.
        // To circumvent this, either use less hop bits or a larger
        // hash table, i.e. a larger u
        //
        // NOTE: The number of hop bits could be automatically calculated.
        // however, 8 is a reasonable value, since it is equal to a 64 byte
        // cache line
        // NOTE: This is for evaluation purposes only
        // this implementation could work just as well without the h in this term
        // since the hop bits are not part of the entry
        let total_admin_bits = required_fp_bits as usize + h;
        if total_admin_bits > 32 {
            return Err("Total sum of admin bits is >32");
        }

        let key_fp_mask = 2_u32.pow(required_fp_bits as u32) - 1;
        // To reach the fp bits, we need to shift out hop bits
        let fingerprint_shift = 0;
        // The remaining 32 bits not taken up by hop bits store fingerprint info
        // these might not all contain fingerprint bits in the current setup
        // depending on the size of u
        // NOTE: This might be counterproductive later
        let entry_fp_bits = 32;
        // a contiguous mask of fp_bits 1-bits, shifted to the right position
        let fingerprint_mask = (2_u64.pow(entry_fp_bits as u32) - 1) << fingerprint_shift;
        // Create the actual hash table
        Ok(PlainHHT {
            h,
            u,
            table: vec![0_u64; u + h - 1], // u hash values, plus h-1 shifting positions for the last hv
            hop_array: vec![0_u32; u + h - 1], // u hash values, plus h-1 shifting positions for the last hv
            hop_bits_mask: 2_u64.pow(h as u32) - 1,
            fingerprint_mask,
            fingerprint_shift,
            admin_bits: 2_u64.pow(32) - 1,
            fp_bits_in_key: required_fp_bits,
            key_fp_mask,
            in_resize: false,
            allow_resize,
        })
    }

    pub fn load(path: &str) -> PlainHHT {
        let loaded_hht: PlainHHT;
        {
            let mut f =
                BufReader::new(File::open(path).unwrap_or_else(|_| {
                    panic!("Opening the saved hash table {} did not work", path)
                }));
            // loaded_hht = deserialize_from(&mut f).unwrap_or_else(|_| {
            //     panic!(
            //         "Deserializing the index from saved index file {} did not work",
            //         path
            //     )
            // }
            loaded_hht = deserialize_from(&mut f).unwrap();
        }
        loaded_hht
    }

    pub fn save(&mut self, path: &str) {
        eprintln!("Saving hash table to {}", path);
        let mut f = BufWriter::new(File::create(path).unwrap_or_else(|_| {
            panic!(
                "Opening file to save index at {} did not work. Check that the path exists.",
                path
            )
        }));
        eprintln!("This would take {:?}", serialized_size(self));
        serialize_into(&mut f, self).expect("Serialization of the index for saving did not work.");
    }

    pub fn get_h(&self) -> usize{
        self.h
    }
    
    /// Compute the fingerprint length for a given size u
    fn fp_length_for(u: u32) -> usize {
        (2_u64.pow(32) as f64 / f64::from(u)).log2().floor() as usize
    }

    #[inline]
    pub fn split_key(&self, key: u32) -> (usize, u32) {
        // get lowest bit from key to use as fingerprint
        // shift rest and cast into usize to use as address

        // this is the | fp | address| version
        // let address = (key as usize) % self.u;
        // let fp = (key as usize) / self.u;

        // this is the | address | fp | version
        let fp = key & self.key_fp_mask;
        let (address, _) = key.overflowing_shr(self.fp_bits_in_key as u32);

        (address as usize, fp as u32)
    }

    fn _restore_key(&self, address: usize, fp: u32) -> u32 {
        (address << self.fp_bits_in_key) as u32 | (fp & self.key_fp_mask)
    }

    fn restore_key_with(address: usize, fp: u32, fp_bits_in_key: usize, key_fp_mask: u32) -> u32 {
        (address << fp_bits_in_key) as u32 | (fp & key_fp_mask)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Hop bit alteration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Retrieve the hop bit vector for the given address
    pub fn get_hop_bits(&self, address: usize) -> u64 {
        self.hop_array[address] as u64 & self.hop_bits_mask
    }

    /// Get initial hop bits for an address, taking into account the h-1
    /// slots before the target address
    fn initialize_insert_hop_bits(&self, address: usize) -> u64 {
        let start = if (self.h - 1) > address {
            // lower edge of table
            0
        } else {
            // start h-1 position before the address
            // to pass all positions that can influence
            // the address slot
            address - (self.h - 1)
        };
        // extract bits for the first
        let mut shifting_hop_bits = self.get_hop_bits(start);
        // NOTE this used to be (start+1)..=address, however, the tests didn't change
        // with this! Is this correct? Is starting from address instead of address+1
        // a computation that is overwritten by later steps?
        for i in start..=address {
            shifting_hop_bits = (shifting_hop_bits >> 1) | self.get_hop_bits(i);
        }
        shifting_hop_bits
    }

    /// Set a specific hop bit of the address to 1
    ///
    /// Example: hop_bits of address 42: 10001
    /// set_hop_bit_in_table(42, 1)
    /// new hop bits of address 42: 10011
    #[inline]
    fn set_hop_bit_in_table(&mut self, address: usize, offset: usize) {
        self.hop_array[address] |= 0b_1 << offset;
    }

    /// Replace the hop bits of the address with the given hop_bits vector
    #[inline]
    fn replace_hop_bits(&mut self, address: usize, hop_bits: u64) {
        self.hop_array[address] = (((self.hop_array[address] as u64) & (!self.hop_bits_mask)) | hop_bits) as u32;
    }

    /// For a given hop bit vector, set a specific position to 0
    #[inline]
    pub fn unset_hop_bit(&self, hop_bits: u64, pos: usize) -> u64 {
        let inverted_mask = self.hop_bits_mask ^ (0b_1 << pos);
        hop_bits & inverted_mask
    }

    /// For a given hop bit vector, set a specific position to 1
    #[inline]
    pub fn set_hop_bit(&self, hop_bits: u64, pos: usize) -> u64 {
        hop_bits | (0b_1 << pos)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Packing, unpacking, and value alteration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Pack a value, fingerprint and hop bits into one 64-bit integer.
    /// The amount of hop bits and fingerprint bits used depends on the size of
    /// the hash table. Some bits might remain 'empty'.
    fn pack(&self, value: u32, fp: u32) -> u64 {
        (u64::from(value) << 32) | u64::from(fp)
    }

    /// Change the value of an entry without changing the hop bits
    fn repack_value(&self, value: u32, fp: u32) -> u64 {
        (u64::from(value) << 32)
            | u64::from(fp) << self.fingerprint_shift
    }

    /// unpack an entry into value, fingerprint, hop_bits
    #[inline]
    fn unpack(&self, entry: u64) -> (u32, u32) {
        (
            (entry >> 32) as u32,
            (entry & self.fingerprint_mask) as u32,
        )
    }

    /// Extract a payload value from an entry by shifting
    /// out the 32 fingerprint and hop bits
    #[inline]
    fn extract_value(&self, entry: u64) -> u32 {
        (entry >> 32) as u32
    }

    /// Store a value at a position in the table
    #[inline]
    fn set_value(&mut self, value: u32, fp: u32, address: usize) {
        self.table[address] = self.pack(value, fp);
    }

    /// For a given bit packed hash table entry, check,
    /// if it has the same fingerprint as the fp provided.
    #[inline]
    fn has_fp(&self, entry: u64, fp: u64) -> bool {
        let target_fp = (entry & self.fingerprint_mask) >> self.fingerprint_shift;
        target_fp == fp
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Helper functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Check if the given address can store into the bin with given
    /// target_offset by evaluating the hop bits. Returns the offset
    /// relative to address of the item that can be moved into target address.
    ///
    /// Example: address 42, with h=4 and hop bits 1101, is queried:
    /// find_offset_to_replace(42, 0) -> None    (offset of 0 in the hop bits is not free)
    /// find_offset_to_replace(42, 1) -> Some(0) (offset of 1 in the hop bits is 0, offset 0 can be moved there)
    /// find_offset_to_replace(42, 2) -> None    (offset of 2 in the hop bits is not free)
    /// find_offset_to_replace(42, 3) -> None    (offset of 3 in the hop bits is not free)
    ///
    /// Helper function for the `free_up_positions` stage of `insert`.
    fn find_offset_to_replace(&self, address: usize, target_offset: usize) -> Option<usize> {
        assert!(
            target_offset < self.h,
            "target_offset < h\ntarget offset: {}\nself.h: {}",
            target_offset,
            self.h
        );
        let hb = self.get_hop_bits(address);
        // make sure the target position is empty and there exists a valid shifting candidate
        if ((hb >> target_offset) & 1 == 0) && (hb != 0) {
            // NOTE: this does not need to take slots before into account.
            // The target slot is guaranteed to be empty.

            // start from the smallest, go to the largest
            // to minimize the number of replacements needed by
            // choosing the biggest possible step size
            for i in 0..target_offset {
                // find the 1-bit making the most way towards the insert point
                if ((hb >> i) & 0b_1) == 1 {
                    return Some(i);
                }
            }
            // no valid offset was found
            None
        } else {
            // either the target position is not free or
            // the hop bits are empty
            None
        }
    }

    /// Get the offsets for filled position for the given address.
    /// Associated positions are extracted from the hop bits
    /// and returned as a vector of offset from the address.
    ///
    /// Example:
    /// hop bits: 0111 -> offsets: [0, 1, 2]
    /// so that
    /// for o in offsets {
    ///     table[address+offset]
    /// }
    /// Yields an entry associated with the address.
    /// Note that these can be soft collisions.
    ///
    /// Helper function for `get` and `delete`.
    #[inline]
    fn occupied_positions(&self, address: usize) -> Vec<usize> {
        let mut occupied_positions = Vec::with_capacity(self.h);
        let positions = self.get_hop_bits(address);

        let mut offset = 0;
        loop {
            if (positions >> offset & 1) == 1 {
                occupied_positions.push(offset);
            }

            offset += 1;
            if offset >= self.h {
                return occupied_positions;
            }
        }
    }


    // /// Shift entries towards the address within one neighbourhood
    // /// This should only be possible when deletions occurred.
    // /// 
    // pub fn compact(&mut self, address: usize) -> Option<usize>{
    //     // get hop bit mask showing free (0) and filled (1) positions for the current slot
    //     let mut shifting_hop_bits = self.initialize_insert_hop_bits(address);
    //     let mut occupied_positions = self.occupied_positions(address);
    //     // let mut target_offset = 0;
    //     let mut highest_occupied = *occupied_positions.iter().max().unwrap();
    //     let mut moved = 0;

    //     let mut target_offset = 0;
    //     loop {
    //         if target_offset >= highest_occupied {
    //             break;
    //         }
    //         // eprintln!("Checking offset {} with hop bits {:b}", target_offset, shifting_hop_bits);
    //         // look for the first genuinely empty slot
    //         if (shifting_hop_bits & 0b_1) == 0 {
    //             // eprintln!("Could move {} + {} to {} + {}", address, highest_occupied, address, target_offset);
    //             // remove highest occupied, move into address + offset
    //             moved += 1;

    //             let address_from = address + highest_occupied;
    //             // let address_to = address + target_offset;
                
    //             let entry = self.table[address_from];
    //             let (value, fp, tmp_hop_bits) = self.unpack(entry);

    //             let hop_bits = self.get_hop_bits(address);

    //             let mut new_hop_bits_for_mca = self.unset_hop_bit(hop_bits, highest_occupied);
    //             new_hop_bits_for_mca =
    //                 self.set_hop_bit(new_hop_bits_for_mca, target_offset);

    //             self.table[address_from] = self.pack(0, 0, tmp_hop_bits);
    //             self.replace_hop_bits(address, new_hop_bits_for_mca);
    //             self.set_value(value, fp, address+target_offset);

    //             // TODO This changes the shifting hop bits etc.
    //             // reset them by calling self.initialize_insert_hop_bits(address);
    //             // again.
                
    //             occupied_positions.retain(|x| x != &highest_occupied);
    //             occupied_positions.push(target_offset);
    //             highest_occupied = *occupied_positions.iter().max().unwrap();

    //         } 
    //         // the current position is full
    //         // shift to look at a farther offset
    //         target_offset += 1;
    //         let new_addr = address + target_offset;
    //         shifting_hop_bits = (shifting_hop_bits >> 1) | self.get_hop_bits(new_addr);

    //     }
    //     if moved > 0 {
    //         Some(moved)
    //     } else {
    //         None
    //     }
    // }
    
    /// Create a free address for insertion by shifting items to higher
    /// addresses in their respective pages.
    /// In other words: Try to shift an empty slot towards the target address
    /// so that a new key can be inserted.
    ///
    /// The `address` parameter is the address for the newly inserted key
    /// `free_offset` is the distance from said address to the next free slot
    /// in the hash table.
    ///
    /// Helper function for `insert`.
    pub fn free_up_positions(
        &mut self,
        address: usize,
        free_offset: usize,
    ) -> Result<usize, &'static str> {
        // http://people.csail.mit.edu/shanir/publications/disc2008_submission_98.pdf

        // address is the POSITION in the HT at which the new item should be inserted
        // free_offset is the DISTANCE/ OFFSET from address to the next free slot
        //
        // starting from (address + free_offset) work backwards to move an empty slot
        // into range h of address.
        //
        // for j = (address + free_offset), j > (address + h - 1)
        // check positions (j - h + 1)..j
        // if one of these contains an item that can be moved into j
        // do it
        //
        // repeat until j is closer than h - 1 positions to the initial address
        let mut j = address + free_offset;
        // eprintln!("j = {} = addr {} + free_offset {}", j, address, free_offset);
        // move backwards from the first free slot by h positions
        // and try to find an entry that can be moved into the free slot.
        loop {
            // sub-loop to check H-1 slots below address + active_free_offset
            // i.e. the addresses that proncipially can move items into the free slot.
            let mut successfully_shifted = false;
            for move_candidate_address in (j - (self.h - 1))..j {
                // if current address has items that can be moved into
                // the current free spot (which is guaranteed to be free
                // due to previous steps of this loop or the initial free slot)
                // move it.
                // eprintln!("Evaluating move candidate address {}", move_candidate_address);
                if let Some(moveable_offset) =
                    self.find_offset_to_replace(move_candidate_address, j - move_candidate_address)
                {
                    // move the identified offset into the current free slot (j)
                    // update the current free slot and move on

                    // The address from which a value can be moved to free up space
                    // Note that this is the address computed for said key
                    // plus an offset at which is was inserted relative to the initial address
                    let address_from = move_candidate_address + moveable_offset;

                    // Extract the entry that is moved
                    let entry = self.table[address_from];
                    let tmp_hop_bits = self.hop_array[address_from];
                    let (value, fp) = self.unpack(entry);

                    // compute the new hop bits for the original address () of the moved value
                    // by removing the old offset of said entru and adding the offset to the new position
                    let hop_bits = self.get_hop_bits(move_candidate_address);
                    let address_offset_to_j = j - move_candidate_address;

                    // Assemble new hop bits for the move_candidate_address.
                    // These will be stored at the real address of the key that is moved
                    let mut new_hop_bits_for_mca = self.unset_hop_bit(hop_bits, moveable_offset);
                    new_hop_bits_for_mca =
                        self.set_hop_bit(new_hop_bits_for_mca, address_offset_to_j);
                    // eprintln!("Setting new hop bits {:b} for position {}", new_hop_bits_for_mca, move_candidate_address);

                    // Clear the slot that the item was moved out of, only adding its hop bits
                    // back in. These are not changed, unless moveable offset = 0
                    if moveable_offset == 0 {
                        // this is the case, where addres_from and movable offset are the same slot
                        assert_eq!(address_from, move_candidate_address);
                        self.table[address_from] = self.pack(0, 0);
                        self.hop_array[address_from] = new_hop_bits_for_mca as u32;
                    } else {
                        // add the hop bits that were present at the slot from which the key
                        // was extracted back into the table
                        self.table[address_from] = self.pack(0, 0);
                        self.hop_array[address_from] = tmp_hop_bits;
                        self.replace_hop_bits(move_candidate_address, new_hop_bits_for_mca);
                    }

                    // enter the shifted value at the target position j
                    self.set_value(value, fp, j);

                    // eprintln!("Shifted from {}+{} to {}", move_candidate_address, moveable_offset, j);
                    // make sure a new free slot is set and terminate the sub-loop
                    j = address_from;
                    successfully_shifted = true;
                    break;
                }
            }

            // check if the sub-loop above (h-1 slots before the current free position)
            // could shift, if not, stop here and trigger a resize in insert.
            if !successfully_shifted || j < address {
                return Err("No freeable slot for address. Needs a resize.");
            }

            // stop, when a freed up slot is close enough to the target address
            if j < (address + self.h) {
                return Ok(j - address); // the offset from the keys target address to the freed up slot
            }
        }
    }

    /// Double the size of the HT to accomodate more entries
    ///
    /// Iterate through all slots from u down to h
    /// for each slot:
    ///   extract all entries (key-value pairs) stored at this address (up to h)
    ///   restore their keys from address and remainder
    ///   reinsert the key-value pairs using the new table parameters
    /// for the last h slots (0 .. h-1)
    /// extract all values into a vector
    /// reinsert the key value pairs in the vec
    ///
    /// Note that the order of keys within a sequence of slots is stable.
    /// They are extracted in a certain order and reinserted in the same order.
    fn resize(&mut self) {
        // NOTE: Allocate the new size. Start from the largest hash value, pull a new bit out of the FP
        // and reenter the key. This potentially allows to resize without allocating |old HT| + |new HT| and recomputing the hashes
        // but do with |new HT| and only repacking.

        // Proof that in place shifting works for addresses larger than h:
        //
        // If a given address a receives an additional (least significant) bit,
        // the new address a' is either a' = 2a (0-bit) or a' = 2a+1 (1-bit).
        // Unless a <= h-1, a new address can always be inserted without touching old values.
        // Since:
        //  [new address]        a'  >  (a-1) + h - 1   [highest non shifted entry; rightmost entry ((h-1)-th soft collision) in the bucket of (a-1)]
        //                      2a   >  a + h - 2
        //                       a   >  h - 2
        //

        // TODO: Can resizing create more collisions?
        // In other words: Can a valid table that is being resized become invalid?
        // I don't think so: Neighbouring address in the old table are spaced apart in the new one.
        //
        // Assume a >> h and both fp bits are 0.
        // a and a+1 turn into 2a and 2(a+1) = 2a+2
        //
        // |   a  |  a+1 |  ->  |  2a  | 2a+1 | 2a+2 |
        // +------+------+      +------+------+------+
        // |  e1  |  e2  |      |  e1  |      |  e2  |
        //
        // Therefore, there exist \floor(h/2) additional and free slots
        // that can take up entries.
        //
        // For fp(e2) = 0 and fp(1) = 1:
        //
        // |   a  |  a+1 |  ->  |  2a  | 2a+1 | 2a+2 |
        // +------+------+      +------+------+------+
        // |  e1  |  e2  |      |      |  e1  |  e2  |
        // This results in the same behaviour as before and cannot become worse
        //
        // for a <= h:
        // This is hard and (hopefully) an unlikely edge case
        eprintln!("RESIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%  {} -> {} slots  (~{}MB) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", self.u, 2* self.u, (2* self.u * 8) / 1000_usize.pow(2));
        // set flag that a resize is in progress.
        // this flag is used to prevent a call to insert made during the resize
        // process that triggers another resize. This should not happen and is
        // always an error.
        self.in_resize = true;

        // Update all administrative parameters, keeping a backup of the old ones
        // needed to unpack and restore old entries.
        let old_len = self.table.len();
        self.table.reserve_exact((2 * self.u) + self.h - 1);
        self.table.extend(vec![0; self.u]);
        self.hop_array.reserve_exact((2 * self.u) + self.h - 1);
        self.hop_array.extend(vec![0; self.u]);

        let old_fp_bits_in_key = self.fp_bits_in_key;
        let old_key_fp_mask = self.key_fp_mask;

        self.fp_bits_in_key -= 1;
        self.key_fp_mask >>= 1;

        self.u *= 2;

        // NOTE this could theoretically go wrong, if during a resize, another resize is triggered.
        // this will mess up the unpacking with old values.
        // A solution for this would be to make the resize function recursive
        // however, there is still no way to track, in which iteration, which key was
        // (re)inserted.
        //
        // However, this should not arise in the first place.
        // Resizing cannot invalidate a table and the given table is valid
        // before the item that triggered the resize is added.
        // Hence, this should be a valid table that is moved to a larger space.
        //
        // To assert this, the variable in_resize is checked, before a resize is triggered.

        // only run until h, put the rest into a vec and reinsert them piecewise
        // this is to prevent reinserted keys touching positions still
        // occupied by not updated entries.
        for old_address in (self.h..old_len).rev() {
            let hb = self.get_hop_bits(old_address);
            if hb > 0 {
                // unpack
                // shift one bit from fingerprint to address
                // insert back into table
                for offset in self.occupied_positions(old_address) {
                    // Current problem:
                    // if a cluster of keys all has trailing ones as fingerprints,
                    // resizing does not solve the resize issue.
                    // Also there are incosistencies concerning the order of fingerprints and
                    // addressbits in the hash value.

                    let extracted_entry = self.table[old_address + offset];
                    self.set_value(0, 0, old_address + offset);

                    // unpack, with old params
                    let (value, fp) = self.unpack(extracted_entry);
                    // restore original key used for insertion
                    let key = PlainHHT::restore_key_with(
                        old_address,
                        fp,
                        old_fp_bits_in_key,
                        old_key_fp_mask,
                    );
                    // dbg!(key, old_address, fp);
                    // eprintln!("Restored: {:>32b} = {:>32b} | {:>32b}", key, old_address, fp);

                    self.insert(key, value).unwrap();
                }
                // after this, all entries stored for the addres were removed
                // and its hop bits are 0
                self.replace_hop_bits(old_address, 0);
            }
        }

        // print!("Resized\n\n");

        // the last h addresses as well as their soft collisions cannot be
        // reinserted directly, but need to
        let mut kv_pairs_left = Vec::with_capacity(self.h.pow(2)); // 2h-1 should sufffice
        for old_address in (0..self.h).rev() {
            let hb = self.get_hop_bits(old_address);
            if hb > 0 {
                for offset in self.occupied_positions(old_address) {
                    let extracted_entry = self.table[old_address + offset];
                    self.set_value(0, 0, old_address + offset);
                    let (value, fp) = self.unpack(extracted_entry);
                    let old_key = PlainHHT::restore_key_with(
                        old_address,
                        fp,
                        old_fp_bits_in_key,
                        old_key_fp_mask,
                    );
                    // dbg!(old_key);

                    // instead of re-inserting directly, put the kv pair into a vector that
                    // will be drained later for insertion
                    kv_pairs_left.push((old_key, value));
                }
                self.replace_hop_bits(old_address, 0);
            }
        }

        // insert all leftover key value pairs back into the table.
        for (key, value) in kv_pairs_left.drain(..) {
            self.insert(key, value).unwrap();
        }

        // Resize has finished. Set the flag accordingly before continueing
        self.in_resize = false;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Basic operation
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// Insert a new (key, value) pair into the HT
    pub fn insert(&mut self, key: u32, value: u32) -> Result<(), &'static str> {
        let (address, fp) = self.split_key(key);

        // eprintln!("Inserted: {:>32b} = {:>32b} | {:>32b}", key, address, fp);
        // eprintln!("Inserted:   key:  {:} = new address: {} | new fp: {}", key, address, fp);

        // Check if address already has all hop bits set.
        // If that is the case, immediately resize, there is
        // no way we can make room for the new key in the current setup.
        if self.get_hop_bits(address) == self.hop_bits_mask {
            if self.allow_resize {
                self.resize();
                self.insert(key, value).unwrap();
                return Ok(());
            } else {
                return Err("Resizes not allowed; full slot");
            }
        }

        // Use linear probing to find the first empty slot.
        // Extract hop bits, shift them after each iteration
        // to maintain a list of occupied positions
        // start with the accumulated hop bits from the H-1 positions
        // left from the address, which could have filled the slot
        let mut shifting_hop_bits = self.initialize_insert_hop_bits(address);
        let mut probe_offset = 0;
        loop {
            // look for the first genuinely empty slot
            if (shifting_hop_bits & 0b_1) == 0 {
                // found empty position
                if probe_offset >= self.h {
                    // start shifting process
                    if let Ok(freed_offset) = self.free_up_positions(address, probe_offset) {
                        // if a valid free position was found, fill it with the given value and fp
                        self.set_value(value, fp, address + freed_offset);
                        self.set_hop_bit_in_table(address, freed_offset);
                        return Ok(());
                    } else {
                        // prevent triggering a resize during an active resize
                        if self.in_resize {
                            panic!("Double resize");
                        }
                        if self.allow_resize {
                            self.resize();
                        } else {
                            return Err("Resizes not allowed, couldn't move free slot");
                        }
                        // eprintln!("After resize");
                        // self.print_ht_fw();

                        self.insert(key, value).unwrap();
                        // eprintln!("After insert");
                        // self.print_ht_fw();

                        return Ok(());
                    }
                } else {
                    // all is fine. insert at address + offset
                    // set probe_offset bit in hop_bits(address)
                    self.set_value(value, fp, address + probe_offset);
                    self.set_hop_bit_in_table(address, probe_offset);
                    return Ok(());
                }
            } else {
                // the current position is full
                // shift to look at a farther offset
                probe_offset += 1;

                let new_addr = address + probe_offset;
                // Check if the end of the table is reached
                if new_addr >= self.table.len() {
                    // prevent triggering a resize during an active resize
                    if self.in_resize {
                        panic!("Double resize");
                    }
                    if self.allow_resize {
                        self.resize();
                    } else {
                        return Err("Resizes not allowed, ran over last slot");
                    }
                    self.insert(key, value).unwrap();
                    return Ok(());
                } else {
                    // shift already collected hop bits one position. This
                    // shifts out the last active position. By OR-ing in the
                    // hop bits for the new active position all filled position
                    // bits are combined.
                    shifting_hop_bits = (shifting_hop_bits >> 1) | self.get_hop_bits(new_addr);
                }
            }
        }
    }

    /// Increment the count for the supplied key.
    /// This is only a valid operation when the HT is used as a counter.
    pub fn increment_count(&mut self, key: u32) -> Result<(), &'static str> {
        let mut hit_addresses = Vec::with_capacity(self.h);

        // Get address and fingerprint
        let (address, query_fp) = self.split_key(key);

        // identify positions with target address.
        // these can contain soft collisions with different
        // fingerprint
        for offset in self.occupied_positions(address) {
            // Check if the entry shares its fingerprint with
            // the query to weed out soft collisions
            if self.has_fp(self.table[address + offset], u64::from(query_fp)) {
                hit_addresses.push(address + offset);
            }
        }

        match hit_addresses.len() {
            0 => {
                // this key was not yet present
                // insert it with a count of 1
                return self.insert(key, 1);
            }
            1 => {
                // This key was already present
                // increment its count by 1
                let hit_address = hit_addresses[0];
                let (value, fingerprint) = self.unpack(self.table[hit_address]);
                self.table[hit_address] = self.pack(value + 1, fingerprint);
                return Ok(());
            }
            x => {
                // This should not arise when using the table for q-gram counting
                panic!(
                    "More than one hit ({}). This should not be possible with a counting table.",
                    x
                );
            }
        }
    }

    /// Get the count for the supplied key.
    /// This is only a valid operation when the HT is used as a counter.
    pub fn get_count(&self, key: u32) -> Option<u32> {
        let mut hit_address = None;

        // Get address and fingerprint
        let (address, query_fp) = self.split_key(key);
        let query_fp = u64::from(query_fp);
        
        // identify positions with target address.
        // these can contain soft collisions with different
        // fingerprint
        for offset in self.occupied_positions(address) {
            // Check if the entry shares its fingerprint with
            // the query to weed out soft collisions
            if self.has_fp(self.table[address + offset], query_fp) {
                hit_address = Some(address + offset);
            }
        }

        match hit_address {
            None => None,
            Some(hit_address) => {
                let (value, _) = self.unpack(self.table[hit_address]);
                Some(value)
            }
        }

    }

    // /// Mybe call this replace?
    // ///
    // pub fn update(&self, key: u32) -> Option<Vec<u32>> {

    // }

    // increment and decrement?

    /// Get all entries for the given key in a Option<Vector>.
    /// If no entry is found, return None.
    pub fn get(&self, key: u32) -> Option<Vec<u32>> {
        // Initialize output. At most h hits can be found.
        let mut hits = Vec::with_capacity(self.h);
        
        // Get address and fingerprint
        let (address, query_fp) = self.split_key(key);
        // prefetch
        // if self.h > 8 {
        //     unsafe {
        //         _mm_prefetch(self.table[address + 8] as *const i8, _MM_HINT_T0);
        //     }
        // }
        // identify positions with target address.
        // these can contain soft collisions with different
        // fingerprint
        for offset in self.occupied_positions(address) {
            let candidate = self.table[address + offset];

            // Check if the entry shares its fingerprint with
            // the query to weed out soft collisions
            if self.has_fp(candidate, u64::from(query_fp)) {
                hits.push(self.extract_value(candidate))
            }
        }

        if !hits.is_empty() {
            Some(hits)
        } else {
            None
        }
    }

    /// Remove the occurrences of this key from the hash table
    ///
    /// What should the signature be? (key) or (key, value)?
    /// Currently: (key) removes all occurences of key
    pub fn delete(&mut self, key: u32) -> Result<(), &'static str> {
        let (address, query_fp) = self.split_key(key);

        let mut updated_hop_bits = self.get_hop_bits(address);

        for offset in self.occupied_positions(address) {
            let candidate = self.table[address + offset];

            // Check if the entry shares its fingerprint with
            // the query to weed out soft collisions
            if self.has_fp(candidate, u64::from(query_fp)) {
                // println!("Deleting offset {}", offset);
                // take a one bit, shift it by the current offset.
                // invert an h-bit bitvector using the mask XOR the set one-bit
                // AND it to the current hop bits to set the current offset to 0
                updated_hop_bits &= self.hop_bits_mask ^ (1 << offset);
                self.set_value(0, 0, address + offset);
            }
        }

        self.replace_hop_bits(address, updated_hop_bits);

        Ok(())
    }

    //////////////////////////////////////////////////////////////
    // DEBUG. Remove this after debugging the insert bug is done
    pub fn _print_ht_fw(&self) {
        let mut last_empty = false;
        for (i, entry) in self.table.iter().enumerate() {
            match (*entry == 0, last_empty) {
                (false, true) => {
                    println!("{:3}  {:>64b}", i, entry);
                    last_empty = false;
                }
                (false, false) => {
                    println!("{:3}  {:>64b}", i, entry);
                }
                (true, true) => (),
                (true, false) => {
                    println!("{:3}  {:>64b}\n[...]", i, entry);
                    last_empty = true;
                }
            }
        }
    }

    // telemetry
    pub fn get_size(&self) -> usize {
        self.u
    }

    pub fn fill_rate(&self) -> f64 {
        let mut nonzero = 0_usize;
        for (_addr, value) in self.table.iter().enumerate() {
            if (value & (!self.hop_bits_mask)) != 0 {
                nonzero += 1;
            }
        }
        nonzero as f64 / (self.table.len() as f64)
    }
}
