import glob
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

HTTP = HTTPRemoteProvider()
FTP = FTPRemoteProvider()


rule get_mxanthus_genome:
    input:
        FTP.remote(
            g["mxanthus"]["download_link"],
            )
    output:
        g["mxanthus"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"


rule get_pfalciparum_genome:
    input:
        FTP.remote(
            g["pfalciparum"]["download_link"],
            )
    output:
        g["pfalciparum"]["path"]
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

rule get_hops_genome:
    input:
        HTTP.remote(
            g["hops"]["download_link"],
            insecure=True,  # Hopbase download is HTTP, not HTTPS
        )
    output:
        g["hops"]["path"],
    shell:
        "mkdir -p genomes && zcat {input} > {output}"

rule get_human_genome:
    input:
        FTP.remote(
            g["hg38"]["download_link"],
        )
    output:
        g["hg38"]["path"],
    shell:
        "mkdir -p genomes && zcat {input} > {output}"
