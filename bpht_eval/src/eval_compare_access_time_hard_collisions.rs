//! Compare access times for BPHTs and PLHTs while keeping them in memory (no creation->saving->reloading).
//! Hard collisions are introduced by sampling from the multiset where
//! 0..2^{size_power} all have a multiplicity of 2, i.e. {{0, 0, 1, 1, ...,
//! 2^{size_power}-1, 2^{size_power}-1 }}.
//! This mimics the behavior of soft collisions due to quotienting to maintain
//! comparability with other tests.
//!
use rand;
use rand::seq::IteratorRandom;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use std::time::Instant;

use bpht::BPHT;

pub enum Hashtable {
    BPHT(BPHT),
    PLHT(crate::plain_hht::PlainHHT),
}

/// Create a BPHT and a PLHT, insert the same keys.
/// Subseqeuntly, query with the whole key set, first for BPHT then for PLHT.
pub fn compare_hard_collision_in_memory(size_power: usize, h: usize, stats_path: &str) {
    let mut rng = rand::thread_rng();
    let size = 2_usize.pow(size_power as u32);

    let mut tidy_results = Vec::new();
    tidy_results.push(format!(
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
    ));

    let mut results = Vec::new();
    results.push(format!(
        "{},{},{},{},{},{},{}\n",
        "size_power",
        "h",
        "nr_keys",
        "bpht_fill_rate",
        "bpht_fill_rate",
        "bpht_keys_stashed",
        "plht_keys_stashed"
    ));

    let mut bpht = BPHT::new(h, size, false).unwrap();
    let mut plht = crate::plain_hht::PlainHHT::new(h, size, false).unwrap();
    let mut bpht_stashed = 0;
    let mut plht_stashed = 0;

    eprintln!("Simulating keys");

    let keys: Vec<u32> = (0..(2_u64.pow(32) - 1) as u32)
        .chain(0..(2_u64.pow(32) - 1) as u32)
        .choose_multiple(&mut rng, 2_u64.pow(size_power as u32) as usize);

    for key in keys.iter() {
        if let Ok(()) = bpht.insert(*key, 1) {
        } else {
            bpht_stashed += 1;
        }
        if let Ok(()) = plht.insert(*key, 1) {
        } else {
            plht_stashed += 1;
        }
    }

    let bpht_fill_rate = bpht.fill_rate();
    let plht_fill_rate = plht.fill_rate();
    results.push(format!(
        "{},{},{},{},{},{},{}\n",
        size_power,
        h,
        2_u64.pow(size_power as u32),
        bpht_fill_rate,
        plht_fill_rate,
        bpht_stashed,
        plht_stashed
    ));

    let mut bpht_counted = 0_usize;
    let mut bpht_sum = 0;
    let mut bpht_stashed = 0_usize;
    let start_bpht = Instant::now();
    for key in keys.iter() {
        if let Some(hits) = bpht.get(*key) {
            bpht_counted += 1;
            bpht_sum += hits.len();
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
        if let Some(hits) = plht.get(*key) {
            plht_counted += 1;
            plht_sum += hits.len();
        } else {
            plht_stashed += 1;
        }
    }
    let time_plht = start_plht.elapsed();
    results.push(format!(
        "BPHT Counted {} keys with total counts of {}.\n",
        bpht_counted, bpht_stashed
    ));
    results.push(format!("BPHT time (ns): {}.\n", time_bpht.as_nanos()));
    results.push(format!(
        "PLHT Counted {} keys with total counts of {}.\n",
        plht_counted, plht_stashed
    ));
    results.push(format!("PLHT time (ns): {}.\n", time_plht.as_nanos()));

    // Summarize results in long form data
    tidy_results.push(format!(
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
    ));

    tidy_results.push(format!(
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
    ));

    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(stats_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(stats_file)
        }
        Err(e) => panic!("Could not open output file for stats: {}", e),
    };

    for result in tidy_results {
        stats_file.write_all(result.as_bytes()).unwrap();
    }
}
