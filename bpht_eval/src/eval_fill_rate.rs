use rand;
use rand::seq::IteratorRandom;

use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;


/// Generate random, unique keys, insert them putting overflows into a stash
pub fn evaluate_fill_rate(size_power:usize, h: usize, steps: usize, stats_path: &str) {
    let mut rng = rand::thread_rng();
    let size = 2_usize.pow(size_power as u32);
    let mut results = Vec::new();
    results.push(format!("{},{},{},{},{}\n", "size_power", "h", "nr_keys", "fill_rate", "keys_stashed"));
    
    // compute key fractions
    let step_size = 1.0 / (steps as f64);
    let nr_keys_iter = (1..=steps)
        .map(
            |x| {
                (x as f64 * step_size * (size as f64)) as usize
            }
        );
    for nr_keys in nr_keys_iter {
        let mut ht = bpht::BPHT::new(h, size, false).unwrap();
        let mut stashed = 0;
        let key_iter = (0..(2_u64.pow(32)-1) as u32)
            .choose_multiple(&mut rng, nr_keys);
        for key in key_iter {
            if let Ok(()) = ht.increment_count(key) {
            } else {
                stashed += 1;
            }
        }
        let fill_rate = ht.fill_rate();
        results.push(format!("{},{},{},{},{}\n", size_power, h, nr_keys, fill_rate, stashed));
    }


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

}
