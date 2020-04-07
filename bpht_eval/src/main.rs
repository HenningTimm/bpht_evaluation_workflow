use bpht_eval::eval_compare_access_time::{compare_in_memory, evaluate_ht, prepare_indices};
use bpht_eval::eval_compare_access_time_hard_collisions::compare_hard_collision_in_memory;
use bpht_eval::eval_fill_rate::evaluate_fill_rate;
use bpht_eval::eval_genome_access_time::{evaluate_genome_access_time, prepare_index};
use bpht_eval::eval_genome_fill_time::evaluate_genome_fill_time;
use structopt;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "BPHT evaluation")]
enum BPHTEval {
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
    match BPHTEval::from_args() {
        BPHTEval::EvalFillRate {
            ht_pow,
            h,
            steps,
            stats_path,
        } => {
            evaluate_fill_rate(ht_pow, h, steps, &stats_path);
        }
        BPHTEval::GenomeFillTime {
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

        BPHTEval::PrepareQueryIndex {
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

        BPHTEval::EvalQueryTime {
            genome_path,
            bpht_path,
            hf_path,
            q,
            stats_path,
        } => {
            evaluate_genome_access_time(&genome_path, &bpht_path, &hf_path, q, &stats_path);
        }

        BPHTEval::CompareInMemory {
            ht_pow,
            h,
            stats_path,
        } => {
            compare_in_memory(ht_pow, h, &stats_path);
        }

        BPHTEval::CompareHardColl {
            ht_pow,
            h,
            stats_path,
        } => {
            compare_hard_collision_in_memory(ht_pow, h, &stats_path);
        }
    }
}
