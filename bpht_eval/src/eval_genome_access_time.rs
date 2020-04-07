use crate::qgram_iterator::{Canonical, HashFunction};
use crate::hash_function::{InvMultParams, HlinParams};
use crate::qgram_iterator;
use bpht::BPHT;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use tab_hash::{Tab32Simple,Tab32Twisted};
use std::time::{Instant, Duration};
use needletail::parse_sequence_path;
use std::ops::Sub;




pub fn evaluate_genome_access_time(
    genome_path: &str,
    bpht_path: &str,
    hf_path: &str,
    q: usize,
    stats_path: &str,
) {

    eprintln!(
        "Querying index {} with {}-grams from genome {}. Results go to {}",
        bpht_path, genome_path, q, stats_path
    );
    let bpht = BPHT::load(bpht_path);
    let canonical = Canonical::No;
    let hash_function = HashFunction::load(hf_path);
    let canonisize_qgrams = false;

    let mut nr_of_qgrams = 0_usize;
    let mut distinct_counted = 0_usize;
    let mut total_counted = 0_usize;
    let mut nr_stashed = 0_usize;
    let mut results = vec![
        format!(
            "{},{},{},{},{},{},{},{},{},{},{}\n",
            "size_power",
            "h",
            "q",
            "hf",
            "fill_rate",
            "nr_qgrams",
            "total_time (ns)",
            "access_time (ns)",
            "distinct_qgrams_counted",
            "total_counted",
            "nr_stashed",
        )
    ];
    
    let total_start = Instant::now();
    let mut chromosome_times = Vec::with_capacity(50);
    parse_sequence_path(
        genome_path,
        |_| {},
        |seq| {
            if seq.seq.len() > q {
                let q_grams: Vec<Option<u32>> = qgram_iterator::QGrams::new(
                    &seq.seq,
                    q,
                    canonical.clone(),
                    hash_function.clone(),
                    canonisize_qgrams,
                ).collect();

                let chromosome_start = Instant::now();
                for q_gram in q_grams {
                    if let Some(q_gram) = q_gram {
                        nr_of_qgrams += 1;

                        if let Some(count) = bpht.get_count(q_gram) {
                            distinct_counted += 1;
                            total_counted += count as usize;
                        } else {
                            nr_stashed += 1;
                        };
                    }
                }
                chromosome_times.push(chromosome_start.elapsed());
            }
        },
    )
    .expect("failed to iterate through FASTA file.");
    let total_time = total_start.elapsed();
    let accumulated_insert_times: Duration = chromosome_times.iter().sum();

    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(distr_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(distr_file)
        }
        Err(e) => panic!("Could not open output file for distribution: {}", e),
    };

    let hf = match hash_function {
        HashFunction::No => "no",
        HashFunction::InvMult(_) => "mult",
        HashFunction::Hlin(_) => "hlin",
        HashFunction::Tab32Simple(_) => "tab_simple",
        HashFunction::Tab32Twisted(_) => "tab_twisted",
    };
    
    results.push(
        format!(
            "{},{},{},{},{},{},{},{},{},{},{}\n",
            bpht.get_size(),
            bpht.get_h(),
            q,
            hf,
            bpht.fill_rate(),
            nr_of_qgrams,
            total_time.as_nanos(),
            accumulated_insert_times.as_nanos(),
            distinct_counted,
            total_counted,
            nr_stashed,
        )
    );
    
    for result in results {
        stats_file
            .write_all(result.as_bytes())
            .unwrap();
    }
}


pub fn prepare_index(
    genome_path: &str,
    size_power: usize,
    h: usize,
    q: usize,
    hf: &str,
    stats_path: &str,
    hf_path: &str,
    bpht_path: &str,
) {
    let ht_size = 2_usize.pow(size_power as u32);
    eprintln!("Analyzing genome {} to build a BPHT of size {} using H={} containing {}-grams", genome_path, ht_size, h, q);
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
    let mut bpht = BPHT::new(h, ht_size, false).unwrap();
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
                        if let Ok(()) = bpht.increment_count(q_gram) {
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
        format!("{},{},{},{},{},{},{},{},{},{}\n", size_power, h, q, hf, bpht.fill_rate(), nr_of_qgrams, qgrams_stashed, total_time.as_secs(), filling_time.as_secs(), init_time.as_secs())
    );

    for result in results {
        stats_file
            .write_all(result.as_bytes())
            .unwrap();
    }
    hash_function.save(hf_path);
    bpht.save(bpht_path);
}
