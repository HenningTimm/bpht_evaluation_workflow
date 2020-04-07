use rand;
use rand::seq::{SliceRandom, IteratorRandom};
// use rand_distr::{Poisson, Distribution};
    
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;

use std::time::Instant;
use bpht::BPHT;

pub enum Hashtable {
    BPHT(BPHT),
    PLHT(crate::plain_hht::PlainHHT),
}


pub fn compare_hard_collision_in_memory(
    size_power: usize,
    h: usize,
    stats_path: &str,
) {
    let mut rng = rand::thread_rng();
    // let pd = Poisson::new(0.01).unwrap();
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

    eprintln!("Simulating keys");



    // let mut key_sets = Vec::new();
    // let initial_keys: Vec<u32> = (0..(2_u64.pow(32)-1) as u32)
    //     .choose_multiple(&mut rng, 2_u64.pow(size_power as u32) as usize);
    // key_sets.push(initial_keys);

    // inverse collision rate describes the number of collisions added, like:
    // number of colliding keys = total number of keys * (1 / inv_col_rate)^step
    // i.e. 100 -> 1 percent (0.01) of keys receive a collision
    // let inverse_collision_rate = 4_u64;
    // let collision_depth = 6;
    // for i in 0..collision_depth {
    //     let sample_size = 2_u64.pow(size_power as u32) / inverse_collision_rate.pow(i as u32 + 1);
    //     key_sets.push(
    //         key_sets[i]
    //             .choose_multiple(&mut rng, sample_size as usize)
    //             .map(|x| *x)
    //             .collect()
    //     );
    // }

    // let keys: Vec<u32> = key_sets.into_iter().rev().flatten().collect();




    let keys: Vec<u32> = (0..(2_u64.pow(32)-1) as u32)
        .chain(0..(2_u64.pow(32)-1) as u32)
        .choose_multiple(&mut rng, 2_u64.pow(size_power as u32) as usize);
    
    // for key in keys.iter().take(100) {
    //     eprintln!("Inserting: {}", key);
    // }
    // eprintln!("Inserting keys");
    let mut foo = 0;
    // for key in keys.iter() {
    for key in keys.iter() {
        foo += 1;
        if let Ok(()) = bpht.insert(*key, 1) {
        } else {
            bpht_stashed += 1;
        }
        if let Ok(()) = plht.insert(*key, 1) {
        } else {
            plht_stashed += 1;
        }
    }
    
    eprintln!("foo: {}", foo);
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
        // counted += 1;
        // stashed += 1;
        if let Some(hits) = plht.get(*key) {
            plht_counted += 1;
            plht_sum += hits.len();
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
