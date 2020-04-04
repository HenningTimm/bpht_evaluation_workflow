import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re

def file_to_df(file_path):
    pattern = "genome=(.+)_s=(\d+)_H=(\d+)_q=(\d+)"
    genome, size, h, q = re.search(pattern, file_path).groups()

    keys = ["genome", "size", "h", "q"]
    values = [genome, int(size), int(h), int(q)]
    with open(file_path, "r") as index_file:
        for line in index_file:
            key, value = line.rsplit(" ", 1)
            if key == "key:" or not key:
                continue
            keys.append(key)
            values.append(float(value.strip()))

    df = pd.DataFrame(
        {key: [value] for key, value in zip(keys, values)}
    )
    return df


def plot():
    dfs = []
    for index_file in snakemake.input.index_tsv:
        dfs.append(file_to_df(index_file))
    df = pd.concat(dfs)

    sns.set(style="whitegrid")
    p = sns.pointplot(
        x="q",
        y="Fill rate:",
        style="h",
        alpha=0.5,
        data=df,
    )
    sns.despine()
    plt.savefig(snakemake.output.fill_rate_pdf)

    sns.pairplot(data=df)
    sns.despine()
    plt.savefig(snakemake.output.runtime_pdf)

if __name__ == "__main__":
    plot()
