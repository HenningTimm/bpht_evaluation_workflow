import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import matplotlib
matplotlib.use('Agg')



def file_to_df(file_path):
    pattern = "genome=([a-zA-Z0-9]+)_.*_run=(\d+)"
    genome, run = re.search(pattern, file_path).groups()
    df = pd.read_csv(file_path)
    df["fraction_stashed"] = df["qgrams_stashed"] / max(df["nr_qgrams"])
    df["Time per insert (μs)"] = (df["fill_time"] / df["nr_qgrams"])*(10**6)
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
    df.sort_values("genome")
    return df


def plot():
    df = get_df()

    df = df.rename(
        columns={
            "h": "Hopscotch Neighborhood",
            "total_time": "Total Runtime (s)",
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
        y="Total Runtime (s)",
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

    plt.savefig(snakemake.output.fill_time_pdf)

    plt.clf()
    sns.set(
        style="whitegrid",
        font_scale=1.2,
    )
    sns.despine()
    g = sns.catplot(
        x="Hopscotch Neighborhood",
        y="Time per insert (μs)",
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
    plt.savefig(snakemake.output.time_per_insert)


if __name__ == "__main__":
    plot()
