import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re


def file_to_df(file_path):
    pattern = "p=(\d+)_H=(\d+)_run=(\d+)"
    size_power, h, run = (int(x) for x in re.search(pattern, file_path).groups())
    df = pd.read_csv(file_path)
    df["fraction_keys"] = df["nr_keys"] / max(df["nr_keys"])
    df["fraction_stashed"] = df["keys_stashed"] / max(df["nr_keys"])
    df["Time (s)"] = df["time (ns)"] / (10**9)
    df["h"] = [f"H={h}" for h in df.h]
    df["run"] = run
    return df


def get_df():
    dfs = []
    for index_file in snakemake.input.stats:
        dfs.append(file_to_df(index_file))
    df = pd.concat(dfs)
    return df


def create_speedup_df(df):
    speedup_df = pd.DataFrame(columns=[r"$\frak{a}$", "Hopscotch neighborhood", "BPHT", "PLHT", "Speedup"])
    for (a, H), d in df.groupby([r"$\frak{a}$", "H"]):
        (table_1, df_1), (table_2, df_2) = d.groupby(["Table"])
    
        m_1 = df_1["time (ns)"].mean()
        m_2 = df_2["time (ns)"].mean()
        if table_1 == "BPHT":
            time_bpht = m_1
            time_plht = m_2
        else:
            time_bpht = m_2
            time_plht = m_1
               
        speedup_df = speedup_df.append(
        pd.DataFrame(
            {r"$\frak{a}$": [a], "Hopscotch neighborhood": H, "BPHT": time_bpht, "PLHT": time_plht, "Speedup": time_plht/time_bpht})
        )
    speedup_df["Hopscotch neighborhood"] = pd.Categorical(speedup_df["Hopscotch neighborhood"], categories=sorted(list(set(speedup_df["Hopscotch neighborhood"])), key=lambda x: int(x.split("=")[1])), ordered=True)
    print(f" min speedup: {min(speedup_df['Speedup'])}  max speedup: {max(speedup_df['Speedup'])}")
    return speedup_df


def plot():
    df = get_df()
    print(df)
    sns.set(style="whitegrid")

    df = df.rename(
        columns={
            "size_power": r"$\frak{a}$",
            "h": "H",
            "table": "Table",
            }
    )
    
    g = sns.catplot(
        row=r"$\frak{a}$",
        x="H",
        y="Time (s)",
        hue="Table",
        kind="swarm",
        data=df,
        height=3,
        aspect=2.5
    )
    plt.gcf().suptitle(snakemake.params.title, y=1.01)
    
    plt.savefig(snakemake.output.atc_pdf, bbox_inches='tight')

    plt.clf()
    speedup_df = create_speedup_df(df)
    custom_palette = sns.color_palette("GnBu", 5)
    sns.set_palette(custom_palette)
    g = sns.catplot(
        hue="$\frak{a}$",
        x="Hopscotch neighborhood",
        y="Speedup",
        data=speedup_df,
        kind="point",
    )
    g.set(ylim=(1.0, 1.7))
    plt.savefig(snakemake.output.speedup_pdf, bbox_inches='tight')


if __name__ == "__main__":
    plot()
