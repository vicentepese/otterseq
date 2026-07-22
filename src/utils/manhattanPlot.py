#!/usr/bin/env python
"""Manhattan plot for GWAS results."""

import math
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
import yaml

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
)

with open("settings.yaml") as f:
    settings = yaml.safe_load(f)

assoc_path = (
    settings["plinkFiles"]["GWASQC"]
    + settings["plinkFiles"]["prefix"]
    + ".assoc.logistic"
)
df = pd.read_csv(assoc_path, delim_whitespace=True)
df = df[df["TEST"] == "ADD"].dropna(subset=["P"]).reset_index(drop=True)

# Alternating colors per chromosome
colors_cycle = ["#F8766D", "#00BFC4"]
color_map = {
    chrom: colors_cycle[i % 2] for i, chrom in enumerate(df["CHR"].unique())
}
df["color"] = df["CHR"].map(color_map)
df["pos"] = range(len(df))
df["log_p"] = df["P"].apply(lambda p: -math.log10(p))

# Chromosome label positions (center of each chromosome's points)
axis_df = df.groupby("CHR")["pos"].agg(
    center=lambda x: (x.max() + x.min()) / 2
)

plt.figure(figsize=(19.2, 10.8), dpi=100)
plt.scatter(df["pos"], df["log_p"], c=df["color"], alpha=1, s=1)
plt.axhline(y=-math.log10(5e-8), color="black", linewidth=0.8)
plt.xticks(axis_df["center"], axis_df.index)
plt.xlabel("Chromosome")
plt.ylabel("-log10(pval)")
plt.title("GWAS")
plt.legend().remove()
plt.tight_layout()
plt.savefig("GWAS_py.png")
print("Saved GWAS_py.png")  # noqa: T201
