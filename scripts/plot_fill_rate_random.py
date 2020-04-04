import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import matplotlib
matplotlib.use('Agg')



def file_to_df(file_path):
    pattern = "p=(\d+)_H=(\d+)_run=(\d+)"
    size_power, h, run = (int(x) for x in re.search(pattern, file_path).groups())
    df = pd.read_csv(file_path)
    df["fraction_keys"] = df["nr_keys"] / max(df["nr_keys"])
    df["fraction_stashed"] = df["keys_stashed"] / max(df["nr_keys"])
    df["h"] = [f"H={h}" for h in df.h]
    df["run"] = run
    return df


def get_df():
    dfs = []
    for index_file in snakemake.input.stats:
        dfs.append(file_to_df(index_file))
    df = pd.concat(dfs)
    print(df)
    return df


def plot():
    df = get_df()

    sns.set(style="whitegrid")

    df = df.rename(
        columns={
            "size_power": r"$\frak{a}$",
            "fraction_keys": "Fraction of inserted keys",
            "fraction_stashed": "Fraction of keys stashed",
            "h": "H",
            "fill_rate": "Fill rate",
        }
    )
    df["X"] = "Measurement"
    print(df)
    f = sns.relplot(
        row=r"$\frak{a}$",
        x="Fraction of inserted keys",
        y="Fill rate",
        hue="H",
        height=5,
        aspect=2,
        alpha=0.7,
        data=df,
        kind="line",
        legend=False,
        style="X",
        markers=True,
        dashes=False,
    )

    sp_groups = [pd.DataFrame(y) for x, y in df.groupby(r'$\frak{a}$', as_index=False)]

    for row, data in zip(f.axes, sp_groups):
        ax_stash = row[0].twinx()

        sns.lineplot(
            x="Fraction of inserted keys",
            y="Fraction of keys stashed",
            hue="H",
            data=data,
            ax=ax_stash,
            alpha=0.7,
            style="X",
            markers=True,
            dashes=False,
            legend="brief",
        )
        ax_stash.legend(loc='upper left')
        row[0].set_ylim(0, 1.0)
        for t in ax_stash.get_legend().texts:
            if t.get_text() == "X":
                t.set_text("")
        ax_stash.set_ylim(0, 1.0)
        for line in ax_stash.lines:
            line.set_linestyle("--")
    sns.despine()
    plt.savefig(snakemake.output.fill_rate_pdf, bbox_inches='tight')


if __name__ == "__main__":
    plot()
