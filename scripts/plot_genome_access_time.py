import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import glob


def file_to_df(file_path):
    pattern = "genome=([a-zA-Z0-9]+)_.*_run=(\d+)"
    genome, run = re.search(pattern, file_path).groups()
    df = pd.read_csv(file_path)
    df["fraction_stashed"] = df["nr_stashed"] / max(df["nr_qgrams"])
    df["Time per access (μs)"] = (df["access_time (ns)"] / df["nr_qgrams"]) / (10**3)
    df["Total Runtime (s)"] = df["total_time (ns)"] / (10**9)
    df["Total Access Time (s)"] = df["access_time (ns)"] / (10**9)
    df["h"] = [f"H={h}" for h in df.h]
    df["run"] = int(run)
    df["genome"] = genome
    return df


def get_df():
    dfs = []

    for index_file in snakemake.input.stats:
        dfs.append(file_to_df(index_file))
    df = pd.concat(dfs)
    df["genome"] = pd.Categorical(df["genome"], ["mxanthus", "pfalciparum", "hops", "hg38"])
    df["h"] = pd.Categorical(df["h"], ["H=8", "H=16", "H=24"])
    df.sort_values("genome")
    return df


def plot():
    df = get_df()
    df = df.rename(
        columns={
            "h": "Hopscotch Neighborhood",
            "total_time (ns)": "Total Runtime (ns)",
            "genome": "Genome",
        }
    )

    sns.set(
        style="whitegrid",
        font_scale=1.2,
    )
    sns.despine()

    g = sns.catplot(
        x="Hopscotch Neighborhood",
        y="Time per access (μs)",
        row="q",
        col="Genome",
        hue="hf",
        kind="point",
        size=5,
        aspect=0.8,
        data=df,
        legend_out=False,
        margin_titles=True,
        dodge=True,
    )

    plt.savefig(snakemake.output.access_time_pdf)
    

if __name__ == "__main__":
    plot()
