use rand;
use rand::seq::IteratorRandom;

use std::fs::File;
use std::io::prelude::*;
use std::io::{BufReader,BufWriter};

use bincode::{deserialize_from, serialize_into};
use std::time::Instant;
use bpht::BPHT;

pub enum Hashtable {
    BPHT(BPHT),
    PLHT(crate::plain_hht::PlainHHT),
}


pub fn evaluate_ht(
    bpht_path:Option<String>,
    plht_path:Option<String>,
    keys_path: &str,
    _stats_path: &str,
) {
    let mut counted = 0_usize;
    let mut stashed = 0_usize;
    let mut sum = 0;

    let ht = match (bpht_path, plht_path) {
        (Some(s), None) => {Hashtable::BPHT(BPHT::load(&s))},
        (None, Some(s)) => {Hashtable::PLHT(crate::plain_hht::PlainHHT::load(&s))},
        _ => panic!("Invalid combination. Can only evaluate either BPHT or PLHR at a time"),
    };
    let keys:Vec<u32> = deserialize_from(BufReader::new(File::open(keys_path).unwrap())).unwrap();

    eprintln!("Read in {} keys", keys.len());

    let total_time = match ht {
        Hashtable::BPHT(bpht) => {
            let start = Instant::now();
            for key in keys {
                // counted += 1;
                // stashed += 1;
                if let Some(count) = bpht.get_count(key) {
                    counted += 1;
                    sum += count;
                } else {
                    stashed += 1;
                }
            }
            start.elapsed()
        }
        Hashtable::PLHT(plht) => {
            let start = Instant::now();
            for key in keys {
                // counted += 1;
                // stashed += 1;
                if let Some(count) = plht.get_count(key) {
                    counted += 1;
                    sum += count;
                } else {
                    stashed += 1;
                }
            }
            start.elapsed()
        }
    };

    eprintln!("Counterd: {}     stashed: {}   sum: {}    time: {}", counted, stashed, sum, total_time.as_nanos());
}


pub fn prepare_indices(
    size_power: usize,
    h: usize,
    stats_path: &str,
    keys_path: &str,
    bpht_path: &str,
    plht_path: &str,
) {
    let mut rng = rand::thread_rng();
    let size = 2_usize.pow(size_power as u32);
    let mut results = Vec::new();
    results.push(format!("{},{},{},{},{},{},{}\n", "size_power", "h", "nr_keys", "bpht_fill_rate", "bpht_fill_rate", "bpht_keys_stashed", "plht_keys_stashed"));
    

    let mut bpht = BPHT::new(h, size, false).unwrap();
    let mut plht = crate::plain_hht::PlainHHT::new(h, size, false).unwrap();
    let mut bpht_stashed = 0;
    let mut plht_stashed = 0;

    // aim for a fill rate of 1
    let keys: Vec<u32> = (0..(2_u64.pow(32)-1) as u32)
        .choose_multiple(&mut rng, 2_u64.pow(size_power as u32) as usize);
    
    for key in keys.iter() {
        if let Ok(()) = bpht.increment_count(*key) {
        } else {
            bpht_stashed += 1;
        }
        if let Ok(()) = plht.increment_count(*key) {
        } else {
            plht_stashed += 1;
        }
    }
    
    let bpht_fill_rate = bpht.fill_rate();
    let plht_fill_rate = plht.fill_rate();
    results.push(format!("{},{},{},{},{},{},{}\n", size_power, h, 2_u64.pow(size_power as u32), bpht_fill_rate, plht_fill_rate, bpht_stashed, plht_stashed));

    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(stats_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(stats_file)
        }
        Err(e) => panic!("Could not open output file for stats: {}", e),
    };

    for result in results {
        stats_file
            .write_all(result.as_bytes())
            .unwrap();
    }

    // hash_function.save(hf_path);
    plht.save(plht_path);
    bpht.save(bpht_path);

    let mut f = BufWriter::new(File::create(keys_path).unwrap());
    serialize_into(&mut f, &keys).unwrap();
}


pub fn compare_in_memory(
    size_power: usize,
    h: usize,
    stats_path: &str,
) {
 let mut rng = rand::thread_rng();
    let size = 2_usize.pow(size_power as u32);

    let mut tidy_results = Vec::new();
    tidy_results.push(
        format!(
            "{},{},{},{},{},{},{},{},{}\n",
            "size_power",
            "h",
            "nr_keys",
            "table",
            "fill_rate",
            "keys_stashed",
            "time (ns)",
            "counted",
            "total_sum",
        )
    );
    
    let mut results = Vec::new();
    results.push(format!("{},{},{},{},{},{},{}\n", "size_power", "h", "nr_keys", "bpht_fill_rate", "bpht_fill_rate", "bpht_keys_stashed", "plht_keys_stashed"));

    let mut bpht = BPHT::new(h, size, false).unwrap();
    let mut plht = crate::plain_hht::PlainHHT::new(h, size, false).unwrap();
    let mut bpht_stashed = 0;
    let mut plht_stashed = 0;

    // aim for a fill rate of 1
    let keys: Vec<u32> = (0..(2_u64.pow(32)-1) as u32)
        .choose_multiple(&mut rng, 2_u64.pow(size_power as u32) as usize);
    
    for key in keys.iter() {
        if let Ok(()) = bpht.increment_count(*key) {
        } else {
            bpht_stashed += 1;
        }
        if let Ok(()) = plht.increment_count(*key) {
        } else {
            plht_stashed += 1;
        }
    }
    
    let bpht_fill_rate = bpht.fill_rate();
    let plht_fill_rate = plht.fill_rate();
    results.push(format!("{},{},{},{},{},{},{}\n", size_power, h, 2_u64.pow(size_power as u32), bpht_fill_rate, plht_fill_rate, bpht_stashed, plht_stashed));


    let mut bpht_counted = 0_usize;
    let mut bpht_sum = 0;
    let mut bpht_stashed = 0_usize;
    let start_bpht = Instant::now();
    for key in keys.iter() {
        // counted += 1;
        // stashed += 1;
        if let Some(count) = bpht.get_count(*key) {
            bpht_counted += 1;
            bpht_sum += count;
        } else {
            bpht_stashed += 1;
        }
    }
    let time_bpht = start_bpht.elapsed();

    let mut plht_counted = 0_usize;
    let mut plht_sum = 0;
    let mut plht_stashed = 0_usize;
    let start_plht = Instant::now();
    for key in keys.iter() {
        // counted += 1;
        // stashed += 1;
        if let Some(count) = plht.get_count(*key) {
            plht_counted += 1;
            plht_sum += count;
        } else {
            plht_stashed += 1;
        }
    }
    let time_plht = start_plht.elapsed();
    results.push(format!("BPHT Counted {} keys with total counts of {}.\n", bpht_counted, bpht_stashed));
    results.push(format!("BPHT time (ns): {}.\n", time_bpht.as_nanos()));
    results.push(format!("PLHT Counted {} keys with total counts of {}.\n", plht_counted, plht_stashed));
    results.push(format!("PLHT time (ns): {}.\n", time_plht.as_nanos()));


    // Summarize results in long form data
    tidy_results.push(
        format!(
            "{},{},{},{},{},{},{},{},{}\n",
            size_power,
            h,
            keys.len(),
            "BPHT",
            bpht.fill_rate(),
            bpht_stashed,
            time_bpht.as_nanos(),
            bpht_counted,
            bpht_sum,
        )
    );

    tidy_results.push(
        format!(
            "{},{},{},{},{},{},{},{},{}\n",
            size_power,
            h,
            keys.len(),
            "PLHT",
            plht.fill_rate(),
            plht_stashed,
            time_plht.as_nanos(),
            plht_counted,
            plht_sum,
        )
    );
    
    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(stats_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(stats_file)
        }
        Err(e) => panic!("Could not open output file for stats: {}", e),
    };

    for result in tidy_results {
        stats_file
            .write_all(result.as_bytes())
            .unwrap();
    }
}
