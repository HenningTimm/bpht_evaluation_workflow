use crate::qgram_iterator::{Canonical, HashFunction};
use crate::hash_function::{InvMultParams, HlinParams};
use crate::qgram_iterator;
use bpht::BPHT;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use tab_hash::{Tab32Simple,Tab32Twisted};
use std::time::Instant;
use needletail::parse_sequence_path;
use std::ops::Sub;
    
pub fn evaluate_genome_fill_time(
    genome_path: &str,
    size_power: usize,
    h: usize,
    q: usize,
    hf: &str,
    // out_path: &str,
    stats_path: &str,
    hf_path: &str,
) {
    let ht_size = 2_usize.pow(size_power as u32);
    eprintln!("Analyzing genome {} to build a HHT of size {} using H={} containing {}-grams", genome_path, ht_size, h, q);
    let mut results = vec![format!("{},{},{},{},{},{},{},{},{},{}\n", "size_power", "h", "q", "hf", "fill_rate", "nr_qgrams", "qgrams_stashed", "total_time", "fill_time", "init_time")];
    
    let mut hash_function = match hf {
        "no" => HashFunction::No,
        "mult" => HashFunction::InvMult(InvMultParams::new()),
        "hlin" => HashFunction::Hlin(HlinParams::new()),
        "tab_simple" => HashFunction::Tab32Simple(Tab32Simple::new()),
        "tab_twisted" => HashFunction::Tab32Twisted(Tab32Twisted::new()),
        _ => panic!("Invalid hash function"),
    };

    let canonical = Canonical::No;
    let canonisize_qgrams = false;
    let mut nr_of_qgrams = 0_usize;
    let mut qgrams_stashed = 0_usize;

    let start_init = Instant::now();
    let mut hht = BPHT::new(h, ht_size, false).unwrap();
    let start_insertion = Instant::now();
    parse_sequence_path(
        genome_path,
        |_| {},
        |seq| {
            if seq.seq.len() > q {
                let q_grams = qgram_iterator::QGrams::new(
                    &seq.seq,
                    q,
                    canonical.clone(),
                    hash_function.clone(),
                    canonisize_qgrams,
                );

                for q_gram in q_grams {
                    if let Some(q_gram) = q_gram {
                        nr_of_qgrams += 1;
                        if let Ok(()) = hht.increment_count(q_gram) {
                        } else {
                            qgrams_stashed += 1;
                        }
                    }
                }
            }
        },
    )
    .expect("failed to iterate through FASTA file.");
    let total_time = start_init.elapsed();
    let filling_time = start_insertion.elapsed();
    let init_time = start_insertion.sub(start_init);

    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(distr_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(distr_file)
        }
        Err(e) => panic!("Could not open output file for distribution: {}", e),
    };

    results.push(
        format!("{},{},{},{},{},{},{},{},{},{}\n", size_power, h, q, hf, hht.fill_rate(), nr_of_qgrams, qgrams_stashed, total_time.as_secs(), filling_time.as_secs(), init_time.as_secs())
    );

    for result in results {
        stats_file
            .write_all(result.as_bytes())
            .unwrap();
    }
    hash_function.save(hf_path);
    // hht.save(out_path);
}
