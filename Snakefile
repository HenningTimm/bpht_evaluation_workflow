"""Evaluation workflow for BPHT.

Requires numctl to be installed.

Note that the accesss time evaluation is tailored
to a specific server with two numa blocks.
Performance on other devices might differ.

The hg38 and hop genomes evaluated are too big for this repository.
They need to be downloaded seperately, placed in the genomes folder
and added to genomes.tsv.
"""
import pandas as pd
import socket
import itertools

# Change this file to the server config when working on a suitably powerful machine
configfile: "config_desktop.yaml"
# Server config requires: Processors with two seperate NUMA blocks, about 2 TB of free disk space and 128GB of RAM
# configfile: "config_server.yaml"

genomes = pd.read_csv("genomes.tsv", sep="\t", comment="#", dtype=str)
g = {name: {"path": path, "download_link": download_link} for name, path, download_link in zip(genomes["name"], genomes["path"], genomes["download_link"])}
print(g)

def rfr_valid_h_p_combinations():
    return [(h, p) for h, p in itertools.product(config["fill_rate"]["H"], config["fill_rate"]["p"]) if h <= p]

rule all:
    input:
        "plots/fill_rate_random/fill_rate_statistics.pdf",
        "plots/genome_fill_time/genome_fill_times.pdf",
        "plots/genome_fill_time/time_per_insert.pdf",
        "plots/normal_access_time_comparison.pdf",
        "plots/counting_access_time_comparison.pdf",
        "plots/counting_speedup.pdf",

include: "rules/analysis.smk"
include: "rules/plotting.smk"
include: "rules/files.smk"
