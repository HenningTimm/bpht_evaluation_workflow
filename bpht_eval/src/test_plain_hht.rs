#[cfg(test)]
mod tests {
    use crate::plain_hht::PlainHHT;
    use itertools::iproduct;
    use rand::seq::SliceRandom;
    use std::collections::HashSet;
    use std::iter::FromIterator;
    use bpht::BPHT;

    /// Test implementing debug methods for a hopscotch HT
    trait HopscotchDebug {
        fn count_total_hop_bits(&self) -> u32;
        fn print_ht_fw(&self);
        fn print_ht(&self);
        fn is_valid(&self) -> bool;
        fn nonzero_entries(&self) -> usize;
        fn key_from_parts(&self, address: u32, remainder: u32) -> u32;
    }
    
    impl HopscotchDebug for PlainHHT {
        fn count_total_hop_bits(&self) -> u32 {
            let mut total_hop_bits = 0;
            for (addr, _value) in self.table.iter().enumerate() {
                let hop_bits = self.get_hop_bits(addr);
                total_hop_bits += hop_bits.count_ones();
            }
            total_hop_bits
        }

        /// Print the hash table, showing only filled buckets
        /// formatted into 64 bits
        fn print_ht_fw(&self) {
            let mut last_empty = false;
            for (i, entry) in self.table.iter().enumerate() {
                match (*entry == 0, last_empty) {
                    (false, true) => {
                        println!("{:3}  {:>64b} {:>32b}", i, entry, self.get_hop_bits(i));
                        last_empty = false;
                    }
                    (false, false) => {
                        println!("{:3}  {:>64b} {:>32b}", i, entry, self.get_hop_bits(i));
                    }
                    (true, true) => (),
                    (true, false) => {
                        println!("{:3}  {:>64b} {:>32b}\n[...]", i, entry, self.get_hop_bits(i));
                        last_empty = true;
                    }
                }
            }
        }

        /// Print a full hash table, formatted into 42 bits.
        fn print_ht(&self) {
            for (i, entry) in self.table.iter().enumerate() {
                println!("{:3}  {:>42b} {:>32b}", i, entry, self.get_hop_bits(i));
            }
        }

        /// Check if the following things hold for the hash table:
        ///
        /// - There are no invalid hop bits, i.e. no two hop bits point to
        ///   the same bucket.
        ///
        fn is_valid(&self) -> bool {
            // check that no conflicting hop-bits are present
            let mut shifting_hop_bits = 0;

            for (addr, _value) in self.table.iter().enumerate() {
                shifting_hop_bits = shifting_hop_bits >> 1;
                let hop_bits = self.get_hop_bits(addr);
                if (shifting_hop_bits & hop_bits) != 0 {
                    panic!("Invalid hop bits!");
                }
            }

            // if this is reached, no conflicting hop bits were found
            true

            // when done, add this as last step to all tests.
        }

        /// count the non-zero entries in the HT by counting
        /// all slots that contain a 1-bit that is not part of the hop
        /// bits.
        /// Note that this can miss counting the value 0 with remainder 0.
        fn nonzero_entries(&self) -> usize {
            let mut nonzero = 0;
            for (_addr, value) in self.table.iter().enumerate() {
                if *value != 0 {
                    nonzero += 1;
                }
            }
            nonzero
        }

        /// assemble a key that will be split into the given address and
        /// remainder for this HT configuration
        fn key_from_parts(&self, address: u32, remainder: u32) -> u32 {
            (address << self.fp_bits_in_key) | remainder
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Test cases
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #[test]
    fn resize_automated() {
        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        for (size_power, h) in iproduct!(5..16, 3..16) {
            eprintln!("Parameter Set: u = 2^{}  h = {}", size_power, h);
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, true) {
                Ok(ht) => ht,
                Err(_) => continue, // skip all invalid configurations
            };

            // generate somewhat random keys
            // that force the HT into resizing
            // by adding >= h items with identical hash value

            // add as many items as possible without resize
            let mut values: Vec<u32> = (0..=h)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            // extract the value that will trigger the resize (the last one)
            let overflow_value = values.pop().unwrap();

            // Assemble test keys. Target setup:
            // |0..0aaa|fff0...0| so that the next resize
            // disperses the clumped keys.
            // Shift left until the msb of the largest
            // fp is next to the current  addr|fp interface
            // so that at least on key will be redistributed
            // to a new address after resize
            let address = rng.gen_range(0, ht.u);
            let addr_shift = 32 - size_power;

            let largest_fp = h - 1;
            let largest_fp_bit_length = (largest_fp as f64).log2().floor() as u32;
            let fp_shift = addr_shift - (1 + largest_fp_bit_length);

            // insert keys
            for (fp, _) in values.iter().enumerate() {
                let key = (address << addr_shift) | (fp << fp_shift);

                ht.insert(key as u32, (fp + 1) as u32).unwrap();
            }

            // ht before resize
            eprintln!("right before resize");
            ht.print_ht_fw();
            // assemble a key that will trigger a resize
            let key = (address << (32 - size_power)) | (2_u32.pow(32 - size_power) - 1) as usize;
            ht.insert(key as u32, overflow_value).unwrap();
            eprintln!("after resize");
            ht.print_ht_fw();

            // test that the table is twice as big (plus the h overflow slots)
            eprintln!(
                "\nht.u: {}\ninitial_u: {}\n2*initial_u: {}\n",
                ht.u,
                initial_u,
                initial_u * 2
            );
            assert_eq!(ht.u, initial_u * 2);
            assert_eq!(ht.table.len(), (initial_u * 2) + h - 1);

            ht.is_valid();
        }
    }

    #[test]
    fn insert_get_identity_unique_keys_automated() {
        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, true) {
                Ok(ht) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    ht
                }
                Err(_) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };

            // generate somewhat random keys
            // that force the HT into resizing
            // by adding >= h items with identical hash value

            // aim for a number of keys that would fill the initial hash table to 90%
            let n = (0.9 * 2_u32.pow(size_power) as f64) as usize;
            // add as many items as possible without resize
            let keys: HashSet<u32> = (0..n)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let values: Vec<u32> = (0..keys.len())
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            for (_step, (key, value)) in keys.iter().zip(values.iter()).enumerate() {
                ht.insert(*key, *value).unwrap();
                // if step+1 != ht.count_total_hop_bits() as usize {
                //     let (address, fp) = ht.split_key(*key);
                //     eprintln!("Hop bits in disarray after inserting key {} = addr {} | fp: {:32b}; value = {:32b}", key, address, fp, value);
                //     ht.print_ht_fw();
                //     eprintln!("Instance for reproduction:");
                //     eprintln!("HopscotchHT::new({}, {});", ht.h, ht.u);
                //     eprintln!("let keys = vec!{:?}", keys);
                //     eprintln!("let values = vec!{:?};", values);
                // }
                // assert_eq!(step+1, ht.count_total_hop_bits() as usize);
                // assert_eq!(step+1, ht.nonzero_entries());
                // eprintln!("State after inserting {}\n___", key);
            }

            for (key, exp_value) in keys.iter().zip(values.iter()) {
                if let Some(value) = ht.get(*key) {
                    assert_eq!(value, vec![*exp_value]);
                } else {
                    eprintln!("\n\n==========================================================================================\n\n");
                    eprintln!(
                        "Error! Could not find a value for key {} with expected value {:b}",
                        key, exp_value
                    );
                    eprintln!("Parameter Set: u = 2^{}  h = {}", size_power, h);
                    let (addr, fp) = ht.split_key(*key);
                    eprintln!("addr: {} fp: {:b}", addr, fp);
                    eprintln!("\n");
                    ht.print_ht_fw();
                    panic!("COULD NOT FIND PREVIOUSLY INSERTED KEY.")
                }
            }
            ht.is_valid();
            assert_eq!(keys.len(), ht.count_total_hop_bits() as usize);
            assert_eq!(keys.len(), ht.nonzero_entries());

            eprintln!("\n##########################################################################################");
            eprintln!("##### EVALUATION FINISHED SUCCESSFULLY  #######");
            eprintln!("##########################################################################################");
        }
        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }

    /// Test if multiplicities are counted correctly.
    #[test]
    fn counting_automated() {
        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, true) {
                Ok(ht) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    ht
                }
                Err(_) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };

            // aim for a number of keys that would fill the initial hash table to 90%
            let n = (0.5 * 2_u32.pow(size_power) as f64) as usize;
            // add as many items as possible without resize
            let keys: HashSet<u32> = (0..n)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let multiplicities: Vec<usize> = (0..n).map(|_| rng.gen_range(1, 100)).collect();

            for (_step, (key, multiplicity)) in keys.iter().zip(multiplicities.iter()).enumerate() {
                for _ in 0..*multiplicity {
                    ht.increment_count(*key).unwrap();
                }
            }

            for (_step, (key, multiplicity)) in keys.iter().zip(multiplicities.iter()).enumerate() {
                assert_eq!(ht.get_count(*key), Some(*multiplicity as u32));
            }
        }
        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }

    #[test]
    fn insert_delete_empty() {
        // inserts certain number of keys
        // delete them all
        // check that all slots are absolutely empty

        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, true) {
                Ok(ht) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    ht
                }
                Err(_) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };

            // aim for a number of keys that would fill the initial hash table to 90%
            let n = (0.9 * 2_u32.pow(size_power) as f64) as usize;
            // add as many items as possible without resize
            let keys: HashSet<u32> = (0..n)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let values: Vec<u32> = (0..keys.len())
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            // insert a certain number of key-value pairs
            for (_step, (key, value)) in keys.iter().zip(values.iter()).enumerate() {
                ht.insert(*key, *value).unwrap();
            }

            // remove all key value pairs
            for key in keys.iter() {
                ht.delete(*key).unwrap();
            }

            // assert that the table is completely empty
            assert_eq!(ht.table, vec![0; ht.u + h - 1]);
            assert_eq!(ht.hop_array, vec![0; ht.u + h - 1]);
        }

        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }

    #[test]
    fn insert_delete_validity() {
        // inserts certain number of keys
        // delete some of them
        // check if the expected keys are contained and/ or removed

        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, true) {
                Ok(ht) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    ht
                }
                Err(_) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };

            // aim for a number of keys that would fill the initial hash table to 90%
            let n = (0.9 * 2_u32.pow(size_power) as f64) as usize;
            // add as many items as possible without resize
            let keys: HashSet<u32> = (0..n)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let mut candidate_keys_to_delete: Vec<&u32> = Vec::from_iter(keys.iter());
            candidate_keys_to_delete.shuffle(&mut rng);

            let values: Vec<u32> = (0..keys.len())
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let m = rng.gen_range(0, n);

            let keys_to_delete: HashSet<u32> = candidate_keys_to_delete
                .iter()
                .take(m)
                .map(|x| **x)
                .collect();

            let remaining_keys: HashSet<u32> = candidate_keys_to_delete
                .iter()
                .skip(m)
                .map(|x| **x)
                .collect();

            // insert a certain number of key-value pairs
            for (_step, (key, value)) in keys.iter().zip(values.iter()).enumerate() {
                ht.insert(*key, *value).unwrap();
            }

            // remove all key value pairs
            for key in keys_to_delete.iter() {
                ht.delete(*key).unwrap();
            }

            for (key, value) in keys.iter().zip(values.iter()) {
                // assert that the two hash sets are a partition of the keys
                match (keys_to_delete.contains(key), remaining_keys.contains(key)) {
                    (true, true) => {
                        panic!("A key cannot be contained in both deleted and not deleted keys!");
                    }
                    (true, false) => {
                        assert!(ht.get(*key).is_none());
                    }
                    (false, true) => {
                        assert_eq!(ht.get(*key).unwrap(), vec![*value]);
                    }
                    (false, false) => {
                        panic!("A key has to be contained in either deleted or not deleted keys!");
                    }
                }
            }
        }

        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }


    #[test]
    fn complete_saturation() {
        color_backtrace::install();
        let mut rng = rand::thread_rng();


        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let mut ht = match PlainHHT::new(h, initial_u, false) {
                Ok(ht) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    ht
                }
                Err(_) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };
        
            let mut keys: Vec<u32> = (0..2_u32.pow(size_power))
                .map(|x| x << (32-size_power))
                .collect();
            keys.shuffle(&mut rng);

            
            for key in keys.iter() {
                if let Ok(()) = ht.insert(*key, 42) {
                } else {
                    // ht._print_ht_fw();
                    eprintln!("Crashed with fill rate: {}", ht.fill_rate());
                    
                    ht.insert(*key, 42).unwrap();
                    panic!("Overflow")
                };
            }
            eprintln!("Expected fill rate: {}\nGot: {}", ht.fill_rate(), (initial_u as f64) / ((initial_u as u64 + h as u64 - 1) as f64));
            assert_eq!(ht.fill_rate(), (initial_u as f64) / (initial_u as u64 + h as u64 - 1) as f64);
        }
        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }


    #[test]
    fn insert_get_identity_plain_vs_bitpacked() {
        color_backtrace::install();
        use rand::Rng;
        let mut rng = rand::thread_rng();

        let mut nr_evaluated = 0;
        for (size_power, h) in iproduct!(5..16, 3..16) {
            let initial_u = 2_usize.pow(size_power);

            let (mut bpht, mut plht) = match (BPHT::new(h, initial_u, true), PlainHHT::new(h, initial_u, true)) {
                (Ok(bpht), Ok(plht)) => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} is valid",
                        size_power, h
                    );
                    nr_evaluated += 1;
                    (bpht, plht)
                }
                _ => {
                    eprintln!(
                        "\n\n\nParameter Set: u = 2^{}  h = {} SKIPPED",
                        size_power, h
                    );
                    continue; // skip all invalid configurations
                }
            };

            // generate somewhat random keys
            // that force the HT into resizing
            // by adding >= h items with identical hash value

            // aim for a number of keys that would fill the initial hash table to 90%
            let n = (0.9 * 2_u32.pow(size_power) as f64) as usize;
            // add as many items as possible without resize
            let keys: HashSet<u32> = (0..n)
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            let values: Vec<u32> = (0..keys.len())
                .map(|_| rng.gen_range(0, (2_u64.pow(32) - 1) as u32))
                .collect();

            for (_step, (key, value)) in keys.iter().zip(values.iter()).enumerate() {
                bpht.insert(*key, *value).unwrap();
                plht.insert(*key, *value).unwrap();
                // if step+1 != ht.count_total_hop_bits() as usize {
                //     let (address, fp) = ht.split_key(*key);
                //     eprintln!("Hop bits in disarray after inserting key {} = addr {} | fp: {:32b}; value = {:32b}", key, address, fp, value);
                //     ht.print_ht_fw();
                //     eprintln!("Instance for reproduction:");
                //     eprintln!("HopscotchHT::new({}, {});", ht.h, ht.u);
                //     eprintln!("let keys = vec!{:?}", keys);
                //     eprintln!("let values = vec!{:?};", values);
                // }
                // assert_eq!(step+1, ht.count_total_hop_bits() as usize);
                // assert_eq!(step+1, ht.nonzero_entries());
                // eprintln!("State after inserting {}\n___", key);
            }

            for (key, exp_value) in keys.iter().zip(values.iter()) {
                if let Some(value) = bpht.get(*key) {
                    assert_eq!(value, vec![*exp_value]);
                } else {
                    eprintln!("\n\n==========================================================================================\n\n");
                    eprintln!(
                        "Error! Could not find a value for key {} with expected value {:b}",
                        key, exp_value
                    );
                    eprintln!("Parameter Set: u = 2^{}  h = {}", size_power, h);
                    let (addr, fp) = bpht.split_key(*key);
                    eprintln!("addr: {} fp: {:b}", addr, fp);
                    eprintln!("\n");
                    panic!("COULD NOT FIND PREVIOUSLY INSERTED KEY.")
                }

                if let Some(value) = plht.get(*key) {
                    assert_eq!(value, vec![*exp_value]);
                } else {
                    eprintln!("\n\n==========================================================================================\n\n");
                    eprintln!(
                        "Error! Could not find a value for key {} with expected value {:b}",
                        key, exp_value
                    );
                    eprintln!("Parameter Set: u = 2^{}  h = {}", size_power, h);
                    let (addr, fp) = plht.split_key(*key);
                    eprintln!("addr: {} fp: {:b}", addr, fp);
                    eprintln!("\n");
                    plht.print_ht_fw();
                    panic!("COULD NOT FIND PREVIOUSLY INSERTED KEY.")
                }
            }

            plht.is_valid();
            assert_eq!(keys.len(), plht.count_total_hop_bits() as usize);
            assert_eq!(keys.len(), plht.nonzero_entries());

            assert_eq!(bpht.fill_rate(), plht.fill_rate());

            eprintln!("\n##########################################################################################");
            eprintln!("##### EVALUATION FINISHED SUCCESSFULLY  #######");
            eprintln!("##########################################################################################");
        }
        // there are 88 valid parameter combinations that can be evaluated
        // for the parameters iproduct!(5..16, 3..16)
        // assert that they are all visited
        assert_eq!(nr_evaluated, 88);
    }
}
