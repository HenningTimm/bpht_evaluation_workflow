rule compile_bpht_eval:
    input:
        f"{config['bpht_eval_path']}/src/lib.rs",
        f"{config['bpht_eval_path']}/src/main.rs",
    output:
        bin="bpht_eval",
    conda:
        "../envs/rust.yaml"
    params:
        target_path=f"{config['bpht_eval_path']}/target/release/bitpacked_hopscotch",
        manifest_path=f"{config['bpht_eval_path']}/Cargo.toml",
    shell:
        "cargo build --release --manifest-path {params.manifest_path} &&"
        "cp {params.target_path} {output.bin}"

rule create_index:
    input:
        bin="bpht_eval",
        genome_path=lambda w: g[w.genome],
    output:
        index=f"{config['tmp_path']}/genome={{genome}}_s={{s}}_H={{H}}_q={{q}}.hht",
        stats="results/index_creation/genome={genome}_s={s}_H={H}_q={q}.tsv",
        hf="hfs/genome={genome}_s={s}_H={H}_q={q}.hf",
    shell:
        "./{input.bin} create-table {input.genome_path} -s {wildcards.s} -H {wildcards.H} -q {wildcards.q} "
        "-o {output.index} -t {output.stats} -f {output.hf}"

rule evaluate_index:
    input:
        bin="bpht_eval",
        genome_path=lambda w: g[w.genome],
        index=f"{config['tmp_path']}/genome={{genome}}_s={{s}}_H={{H}}_q={{q}}.hht",
        hf="hfs/genome={genome}_s={s}_H={H}_q={q}.hf",
    output:
        tsv="results/genome={genome}_s={s}_H={H}_q={q}.tsv"
    shell:
        "./{input.bin} query-table {input.genome_path} {input.index} -q {wildcards.q} -o {output.tsv} -f {input.hf}"

rule eval_fill_rate:
    input:
        bin="bpht_eval"
    output:
        csv="results/random_fill_rate/p={p}_H={H}_run={run}.tsv",
    params:
        steps=config["fill_rate"]["steps"]
    shell:
        "./{input.bin} eval-fill-rate -H {wildcards.H} -p {wildcards.p} -t {output.csv} --steps {params.steps}"

rule eval_fill_time_genome:
    input:
        bin="bpht_eval",
        genome_path=lambda w: g[w.genome],
    output:
        stats="results/genome_fill_time/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.csv",
        hf="hfs/genome_fill_time/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.hf",
    shell:
        "./{input.bin} genome-fill-time {input.genome_path} -p 31 -H {wildcards.H} "
        "--hf {wildcards.hf} -q {wildcards.q} -t {output.stats} -f {output.hf}"

rule prepare_bpht_access_time:
    input:
        bin="bpht_eval",
        genome_path=lambda w: g[w.genome],
    output:
        stats="results/genome_access_time/index/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.csv",
        hf="hfs/genome_access_time/index/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.hf",
        bpht="hfs/genome_access_time/index/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.bpht",
    shell:
        "./{input.bin} prepare-query-index {input.genome_path} -p 31 -H {wildcards.H} "
        "--hf {wildcards.hf} -q {wildcards.q} -t {output.stats} -f {output.hf} --bpht-path {output.bpht} "

rule evaluate_bpht_access_time:
    input:
        bin="bpht_eval",
        genome_path=lambda w: g[w.genome],
        hf="hfs/genome_access_time/index/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.hf",
        bpht="hfs/genome_access_time/index/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.bpht",
    output:
        stats="results/genome_access_time/eval/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.csv",
    shell:
        "./{input.bin} eval-query-time {input.genome_path} --bpht-path {input.bpht} "
        "--hf-path {input.hf} -q {wildcards.q} -t {output.stats} "

rule compare_access_times_numa0:
    input:
        bin="bpht_eval",
    output:
        stats="results/access_time_comparison/p={p}_H={H}_run={run}_numa=0.csv",
    resources:
        numa_0=1
    shell:
        "numactl -N 0 -m 0 ./{input.bin} compare-in-memory -H {wildcards.H} -p {wildcards.p} -t {output.stats} "

rule compare_access_times_numa1:
    input:
        bin="bpht_eval",
    output:
        stats="results/access_time_comparison/p={p}_H={H}_run={run}_numa=1.csv",
    resources:
        numa_1=1
    shell:
        "numactl -N 1 -m 1 ./{input.bin} compare-in-memory -H {wildcards.H} -p {wildcards.p} -t {output.stats} "

rule new_compare_access_times_numa0:
    input:
        bin="bpht_eval",
    output:
        stats="results/new_access_time_comparison/p={p}_H={H}_run={run}_numa=0.csv",
    resources:
        numa_0=1
    shell:
        "numactl -N 0 -m 0 ./{input.bin} compare-hard-coll -H {wildcards.H} -p {wildcards.p} -t {output.stats} "

rule new_compare_access_times_numa1:
    input:
        bin="bpht_eval",
    output:
        stats="results/new_access_time_comparison/p={p}_H={H}_run={run}_numa=1.csv",
    resources:
        numa_1=1
    shell:
        "numactl -N 1 -m 1 ./{input.bin} compare-hard-coll -H {wildcards.H} -p {wildcards.p} -t {output.stats} "

