import glob
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()


rule get_mxanthus_genome:
    input:
        FTP.remote(
            g["mxanthus"]["download_link"],
            keep_local=True,
            )
    output:
        g["mxanthus"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"


rule get_pfalciparum_genome:
    input:
        FTP.remote(
            g["pfalciparum"]["download_link"],
            keep_local=True,
            )
    output:
        g["pfalciparum"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

rule get_hops_genome:
    input:
        HTTP.remote(
            g["hops"]["download_link"],
            keep_local=True,
            insecure=True,  # Hopbase download is HTTP, not HTTPS
        )
    output:
        g["hops"]["path"],
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

# TODO refactor untared path
rule get_human_genome:
    input:
        FTP.remote(
            g["hg38"]["download_link"],
            keep_local=True,
        )
    output:
        g["hg38"]["path"],
    params:
        untared_path=glob.glob(g["hg38"]["download_link"][:-3] + "/*/GRCh38.p13_genomic.fna.gz")
    shell:
        "mkdir -p genomes && tar -xf {input} && zcat {params.untared_path} > {output}"
