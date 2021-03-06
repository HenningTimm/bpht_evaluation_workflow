all_h, all_p = zip(*rfr_valid_h_p_combinations())

rule plot_random_fill_rate:
    input:
        stats=expand(
        expand(
            "results/random_fill_rate/p={p}_H={H}_run={{run}}.tsv",
            zip,
            p=all_p,
            H=all_h,
            ),
                run=range(config["fill_rate"]["runs"]),
            )
    output:
        fill_rate_pdf="plots/fill_rate_random/fill_rate_statistics.pdf",
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plot_fill_rate_random.py"


rule plot_all_genome_fill_times:
    input:
        stats=expand(
            "results/genome_fill_time/genome={genome}_H={H}_hf={hf}_q={q}_run={run}.csv",
            genome=config["target_genomes"],
            H=config["genomes_fill_time"]["H"],
            hf=config["genomes_fill_time"]["hfs"],
            q=config["genomes_fill_time"]["q"],
            run=range(config["genomes_fill_time"]["runs"]),
            )
    output:
        fill_time_pdf="plots/genome_fill_time/genome_fill_times.pdf",
        time_per_insert="plots/genome_fill_time/time_per_insert.pdf"
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plot_genome_fill_time.py"


def comparison_valid_h_p_combinations():
    return [(h, p) for h, p in itertools.product(config["comparison"]["H"], config["comparison"]["p"]) if h <= p]

all_comparison_h, all_comparison_p = zip(*comparison_valid_h_p_combinations())

rule plot_access_time_comparison:
    input:
        stats=expand(expand(
            "results/access_time_comparison/p={p}_H={H}_run={{run}}_numa={{numa}}.csv",
            zip,
            p=all_comparison_p,
            H=all_comparison_h,
            ),
                run=range(config["comparison"]["runs"]),
                numa=config["numa_blocks"],
            )
    output:
        atc_pdf="plots/counting_access_time_comparison.pdf",
        speedup_pdf="plots/counting_speedup.pdf"
    params:
        title="Access Time Comparison for Counting Mode"
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plot_access_time_comparison.py"


def new_comparison_valid_h_p_combinations():
    return [(h, p) for h, p in itertools.product(config["new_comparison"]["H"], config["new_comparison"]["p"]) if h <= p]

all_new_comparison_h, all_new_comparison_p = zip(*new_comparison_valid_h_p_combinations())


rule plot_new_access_time_comparison:
    input:
        stats=expand(expand(
            "results/new_access_time_comparison/p={p}_H={H}_run={{run}}_numa={{numa}}.csv",
            zip,
            p=all_new_comparison_p,
            H=all_new_comparison_h,
            ),
                run=range(config["new_comparison"]["runs"]),
                numa=config["numa_blocks"],
            )
    output:
        atc_pdf="plots/normal_access_time_comparison.pdf",
        speedup_pdf="plots/normal_speedup.pdf"
    params:
        title="Access Time Comparison for Normal Mode"
    conda:
        "../envs/plotting.yaml"
    script:
        "../scripts/plot_access_time_comparison.py"
