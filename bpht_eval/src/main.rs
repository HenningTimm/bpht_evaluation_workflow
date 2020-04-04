///
/// https://doc.rust-lang.org/core/arch/x86_64/fn._mm_prefetch.html
use bpht_eval;
use bpht_eval::hash_function::InvMultParams;
use bpht_eval::qgram_iterator;
use bpht_eval::qgram_iterator::{Canonical, HashFunction};
use bpht_eval::eval_fill_rate::evaluate_fill_rate;
use bpht_eval::eval_genome_fill_time::evaluate_genome_fill_time;
use bpht_eval::eval_genome_access_time::{evaluate_genome_access_time, prepare_index};
use bpht_eval::eval_cache_misses::{prepare_caching_index, evaluate_caching};
use bpht_eval::eval_compare_access_time::{prepare_indices, evaluate_ht, compare_in_memory};
use bpht_eval::eval_compare_access_time_hard_collisions::compare_hard_collision_in_memory;
use needletail::parse_sequence_path;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufWriter;
use std::time::Instant;
use structopt;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "HHT evaluation")]
enum HhtEval {
    CreateTable {
        /// Input file
        #[structopt(name = "genome-path")]
        genome_path: String,

        /// Size of slots in the HT. Has to be a power of 2
        #[structopt(short = "s", long = "ht-size", required = true)]
        ht_size: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// q-gram size
        #[structopt(short = "q", required = true)]
        q: usize,

        /// Path to save the built index
        #[structopt(short = "o", long = "out", required = true)]
        out_path: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,

        /// Path to save the used hash function
        #[structopt(short = "f", long = "hf", required = true)]
        hf_path: String,
    },

    QueryTable {
        /// Input file
        #[structopt(name = "genome-path")]
        genome_path: String,

        /// HHT file
        #[structopt(name = "hht-path")]
        hht_path: String,

        /// q-gram size
        #[structopt(short = "q", required = true)]
        q: usize,

        /// Path to results
        #[structopt(short = "o", long = "out", required = true)]
        out_path: String,

        /// Path to save the used hash function
        #[structopt(short = "f", long = "hf", required = true)]
        hf_path: String,
    },


    EvalFillRate {
        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// nr of steps evaluated
        #[structopt(short = "s", long = "steps", required = true)]
        steps: usize,
        
        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },


    GenomeFillTime {
        /// Input file
        #[structopt(name = "genome-path")]
        genome_path: String,

        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// q-gram size
        #[structopt(short = "q", required = true)]
        q: usize,

        /// q-gram size
        #[structopt(long = "hf", required = true)]
        hf: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,

        /// Path to save the used hash function
        #[structopt(short = "f", long = "hf-path", required = true)]
        hf_path: String,
    },


    PrepareQueryIndex {
        /// Input file
        #[structopt(name = "genome-path")]
        genome_path: String,

        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// q-gram size
        #[structopt(short = "q", required = true)]
        q: usize,

        /// q-gram size
        #[structopt(long = "hf", required = true)]
        hf: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,

        /// Path to save the used hash function
        #[structopt(short = "f", long = "hf-path", required = true)]
        hf_path: String,

        /// BPHT file
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: String,
    },
    

    EvalQueryTime {
        /// Input file
        #[structopt(name = "genome-path")]
        genome_path: String,

        /// BPHT file
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: String,

        /// q-gram size
        #[structopt(long = "hf-path", required = true)]
        hf_path: String,
        
        /// q-gram size
        #[structopt(short = "q", required = true)]
        q: usize,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },

    PrepareCacheMisses {
        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// q-gram size
        #[structopt(long = "hf", required = true)]
        hf: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,

        /// Path to save the used hash function
        #[structopt(short = "f", long = "hf-path", required = true)]
        hf_path: String,


        #[structopt(short = "k", long = "keys-path", required = true)]
        keys_path: String,
        
        /// BPHT file
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: String,
    },

    EvalCacheMisses {
        /// BPHT file
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: String,

        #[structopt(short = "k", long = "keys-path", required = true)]
        keys_path: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },

    PrepareComparisonIndices {
        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,

        #[structopt(short = "k", long = "keys-path", required = true)]
        keys_path: String,
        
        /// BPHT file
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: String,

        /// BPHT file
        #[structopt(long = "plht-path", name = "plht-path")]
        plht_path: String,
    },

    AccessTimeComparison {
        #[structopt(long = "bpht-path", name = "bpht-path")]
        bpht_path: Option<String>,

        #[structopt(long = "plht-path", name = "plht-path")]
        plht_path: Option<String>,

        #[structopt(short = "k", long = "keys-path", required = true)]
        keys_path: String,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },

    CompareInMemory {
        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },

    CompareHardColl {
        /// Size of slots in the HT. As a power of 2
        #[structopt(short = "p", long = "ht-power", required = true)]
        ht_pow: usize,

        /// H-parameter for hopscotch hashing
        #[structopt(short = "H", long = "hop-positions", required = true)]
        h: usize,

        /// Path to save the stats file
        #[structopt(short = "t", long = "stats", required = true)]
        stats_path: String,
    },
    
}

fn main() {
    match HhtEval::from_args() {
        HhtEval::CreateTable {
            genome_path,
            ht_size,
            h,
            q,
            out_path,
            stats_path,
            hf_path,
        } => {
            create_index(
                &genome_path,
                ht_size,
                h,
                q,
                &out_path,
                &stats_path,
                &hf_path,
            );
        }
        HhtEval::QueryTable {
            genome_path,
            hht_path,
            q,
            out_path,
            hf_path,
        } => {
            query_index(&genome_path, &hht_path, q, &out_path, &hf_path);
        }
        HhtEval::EvalFillRate {
            ht_pow,
            h,
            steps,
            stats_path,
        } => {
            evaluate_fill_rate(ht_pow, h, steps, &stats_path);
        }
        HhtEval::GenomeFillTime {
            genome_path,
            ht_pow,
            h,
            q,
            hf,
            // out_path,
            stats_path,
            hf_path,
        } => {
            evaluate_genome_fill_time(
                &genome_path,
                ht_pow,
                h,
                q,
                &hf,
                // &out_path,
                &stats_path,
                &hf_path,
            );
        }

        HhtEval::PrepareQueryIndex {
            genome_path,
            ht_pow,
            h,
            q,
            hf,
            stats_path,
            hf_path,
            bpht_path,
        } => {
            prepare_index(
                &genome_path,
                ht_pow,
                h,
                q,
                &hf,
                &stats_path,
                &hf_path,
                &bpht_path,
            );
        }
        
        HhtEval::EvalQueryTime {
            genome_path,
            bpht_path,
            hf_path,
            q,
            stats_path,
        } => {
            evaluate_genome_access_time(&genome_path, &bpht_path, &hf_path, q, &stats_path);
        }

        HhtEval::PrepareCacheMisses {
            ht_pow,
            h,
            hf,
            stats_path,
            hf_path,
            keys_path,
            bpht_path,
        } => {
            prepare_caching_index(
                ht_pow,
                h,
                &hf,
                &stats_path,
                &hf_path,
                &keys_path,
                &bpht_path,
            );
        }

        HhtEval::EvalCacheMisses {
            bpht_path,
            keys_path,
            stats_path,
        } => {
            evaluate_caching(&bpht_path, &keys_path, &stats_path);
        }

        HhtEval::PrepareComparisonIndices {
            ht_pow,
            h,
            stats_path,
            keys_path,
            bpht_path,
            plht_path,
        } => {
            prepare_indices(
                ht_pow,
                h,
                &stats_path,
                &keys_path,
                &bpht_path,
                &plht_path,
            );
        }

        HhtEval::AccessTimeComparison {
            bpht_path,
            plht_path,
            keys_path,
            stats_path,
        } => {
            evaluate_ht(bpht_path, plht_path, &keys_path, &stats_path);
        }
        
        HhtEval::CompareInMemory {
            ht_pow,
            h,
            stats_path,
        } => {
            compare_in_memory(ht_pow, h, &stats_path);
        }

        HhtEval::CompareHardColl {
            ht_pow,
            h,
            stats_path,
        } => {
            compare_hard_collision_in_memory(ht_pow, h, &stats_path);
        }
        
    }
}

fn create_index(
    genome_path: &str,
    ht_size: usize,
    h: usize,
    q: usize,
    out_path: &str,
    stats_path: &str,
    hf_path: &str,
) {
    eprintln!("Analyzing genome {} to build a HHT of size {} using H={} containing {}-grams which is written to {}", genome_path, ht_size, h, q, out_path);
    let mut hht = bpht::BPHT::new(h, ht_size, true).unwrap();

    let canonical = Canonical::No;
    let mut hash_function = HashFunction::InvMult(InvMultParams::new());
    let canonisize_qgrams = false;
    let mut nr_of_qgrams = 0;
    let start = Instant::now();
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
                        hht.increment_count(q_gram).unwrap();
                    }
                }
            }
        },
    )
    .expect("failed to iterate through FASTA file.");
    let stop = start.elapsed();

    // Write statistics_file
    let mut stats_file = match File::create(stats_path) {
        Ok(distr_file) => {
            eprintln!("Writing statistics file to {}", stats_path);
            BufWriter::new(distr_file)
        }
        Err(e) => panic!("Could not open output file for distribution: {}", e),
    };

    stats_file
        .write_all(
            format!(
                "key: value\nSize: {}\nFill rate: {}\nq-grams: {}\nTime (s): {}",
                hht.get_size(),
                hht.fill_rate(),
                nr_of_qgrams,
                stop.as_secs()
            )
            .as_bytes(),
        )
        .unwrap();
    hash_function.save(hf_path);
    hht.save(out_path);
}

fn query_index(genome_path: &str, hht_path: &str, q: usize, out_path: &str, hf_path: &str) {
    eprintln!(
        "Querying index {} with {}-grams from genome {}. Results go to {}",
        hht_path, genome_path, q, out_path
    );
    let hht = bpht::BPHT::load(hht_path);
    // print!("{:?}", hht);
    let canonical = Canonical::No;
    let hash_function = HashFunction::load(hf_path);
    let canonisize_qgrams = false;

    let start = Instant::now();
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
                        hht.get_count(q_gram).unwrap();
                    }
                }
            }
        },
    )
    .expect("failed to iterate through FASTA file.");
    let stop = start.elapsed();
    // Write output_file
    let mut stats_file = match File::create(out_path) {
        Ok(distr_file) => {
            eprintln!("Writing statistics file to {}", out_path);
            BufWriter::new(distr_file)
        }
        Err(e) => panic!("Could not open output file for distribution: {}", e),
    };

    stats_file
        .write_all(
            format!(
                "Size: {}\nFill rate: {}\nq-grams: {}\nTime (s): {}",
                hht.get_size(),
                hht.fill_rate(),
                "placeholder",
                stop.as_secs()
            )
            .as_bytes(),
        )
        .unwrap();
}
