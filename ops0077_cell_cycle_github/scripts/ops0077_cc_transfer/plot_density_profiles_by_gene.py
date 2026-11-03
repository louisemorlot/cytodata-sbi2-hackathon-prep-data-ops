"""Cell density and marker profiles across the cell cycle for each gene target, in the style of marker_profiles.png.

Per panel (NTC and each knockout gene):
  top    : median (and IQR) of phosphorylated Rb, p21, pRb and DNA content (robust z) along the phase, plus log nuclear area (dashed)
  bottom : density of cells along the phase (histogram), the NTC density as a grey line for reference, NTC-defined sector shading

Two versions: the stain score (4i markers + DNA, frozen NTC curve) and the phase-contrast score (all cells).

Inputs : reports/ops0077_cc_transfer/reference/, reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv.gz
Outputs: reports/ops0077_cc_transfer/ko3_density_profiles/

Usage: python plot_density_profiles_by_gene.py [project_root] [--scores CSV] [--reference DIR] [--out DIR]
"""
import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cc_transfer_lib as L  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--scores", default="reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv.gz")
ap.add_argument("--reference", default="reports/ops0077_cc_transfer/reference")
ap.add_argument("--out", default="reports/ops0077_cc_transfer/ko3_density_profiles")
args = ap.parse_args()
ROOT = Path(args.root)
OUT = ROOT / args.out
OUT.mkdir(parents=True, exist_ok=True)
SECCOL = {"G1": "#4c78a8", "S": "#f58518", "G2/M": "#e45756", "M-to-G1 exit": "#72b7b2"}
AXES = [("z_phos", "C3"), ("z_p21", "C0"), ("z_pRb", "C1"), ("z_dna", "k")]
GENES = ["CCNB1", "MCM3", "CENPC"]
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})


def profile(ax, d, phase_col, nb, min_cells, title):
    edges = np.linspace(0, 1, nb + 1)
    mid = (edges[:-1] + edges[1:]) / 2
    b = pd.cut(d[phase_col], edges, include_lowest=True)
    size = d.groupby(b, observed=False).size().values
    keep = size >= min_cells
    for f, c in AXES:
        g = d.groupby(b, observed=False)[f]
        q50, q25, q75 = g.median().values, g.quantile(.25).values, g.quantile(.75).values
        ax.plot(mid[keep], q50[keep], c=c, label=f, lw=1.8)
        ax.fill_between(mid[keep], q25[keep], q75[keep], color=c, alpha=.10)
    a2 = ax.twinx()
    la = np.log(d["nuc_area_px"]).groupby(b, observed=False).median().values
    a2.plot(mid[keep], la[keep], "--", c="gray", lw=1.2)
    a2.set_ylim(7.3, 8.3)
    a2.spines["top"].set_visible(False)
    a2.set_yticks([])
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(f"{phase_col} (0 = G1 start, wraps)")


def density(ax, ph, ref_ph, ref, bins, color, title, n):
    ph = ph[~np.isnan(ph)]
    ax.hist(ph, bins=bins, range=(0, 1), density=True, color=color, alpha=.85)
    if ref_ph is not None:
        h, e = np.histogram(ref_ph[~np.isnan(ref_ph)], bins=bins, range=(0, 1), density=True)
        ax.step(np.r_[e[:-1], 1], np.r_[h, h[-1]], where="post", color="0.25", lw=1.6, label="NTC")
        ax.legend(fontsize=8, loc="upper center")
    for s, x0, x1 in L.sector_intervals(ref):
        ax.axvspan(x0, x1, alpha=.07, color=SECCOL[s], lw=0)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("phase")
    ax.text(.98, .95, f"n = {n}", transform=ax.transAxes, ha="right", va="top", fontsize=8)


def main():
    ref = L.load_reference(ROOT / args.reference / "ntc_cc_reference.json")
    ntc = pd.read_csv(ROOT / args.reference / "ntc_reference_scores.csv.gz", index_col=0)
    ko = pd.read_csv(ROOT / args.scores, index_col=0)
    genes = [g for g in GENES if g in set(ko["gene"])]
    ntc_ph = ntc["cc_phase_stain"].values

    # ---- stain score: NTC + genes
    fig, ax = plt.subplots(2, 1 + len(genes), figsize=(4.6 * (1 + len(genes)), 7.4), sharey="row")
    d0 = ntc.dropna(subset=["cc_phase_stain"])
    profile(ax[0, 0], d0, "cc_phase_stain", 30, 30, f"NTC (n={len(d0)})")
    density(ax[1, 0], ntc_ph, None, ref, 40, "0.5", "NTC: density of cells", len(d0))
    for k, g in enumerate(genes, start=1):
        d = ko[(ko["gene"] == g) & ko["cc_phase_stain"].notna()]
        profile(ax[0, k], d, "cc_phase_stain", 15, 8, f"{g} (n={len(d)})")
        density(ax[1, k], d["cc_phase_stain"].values, ntc_ph, ref, 40, {"CCNB1": "#1f77b4", "MCM3": "#d62728", "CENPC": "#2ca02c"}[g],
                f"{g}: density of cells", len(d))
    ax[0, 0].legend(fontsize=8, loc="upper right")
    ax[0, 0].set_ylabel("median robust z (IQR band); dashed: log nuclear area")
    ax[1, 0].set_ylabel("cell density")
    fig.suptitle("Stain score (4i markers + DNA, frozen NTC curve). Knockout profiles use bins with >= 8 cells; grey step = NTC density", fontsize=10)
    plt.tight_layout()
    plt.savefig(OUT / "density_profiles_by_gene_stain_score.png", dpi=110, bbox_inches="tight")
    plt.close()

    # ---- phase-contrast score: all cells (markers shown for the cells that have them)
    fig, ax = plt.subplots(2, len(genes), figsize=(4.8 * len(genes), 7.4), sharey="row")
    for k, g in enumerate(genes):
        d = ko[ko["gene"] == g]
        dm = d[d["has_markers"]]
        profile(ax[0, k], dm, "cc_phase_phase", 15, 8, f"{g}: markers (n={len(dm)} with markers)")
        density(ax[1, k], d["cc_phase_phase"].values, ntc_ph, ref, 40, {"CCNB1": "#1f77b4", "MCM3": "#d62728", "CENPC": "#2ca02c"}[g],
                f"{g}: density of cells (all cells)", int(d["cc_phase_phase"].notna().sum()))
    ax[0, 0].legend(fontsize=8, loc="upper right")
    ax[0, 0].set_ylabel("median robust z (IQR band); dashed: log nuclear area")
    ax[1, 0].set_ylabel("cell density")
    fig.suptitle("Phase-contrast score (NTC-trained model, all cells). Grey step = NTC density of the true stain phase, for reference", fontsize=10)
    plt.tight_layout()
    plt.savefig(OUT / "density_profiles_by_gene_phase_contrast_score.png", dpi=110, bbox_inches="tight")
    plt.close()

    # ---- fraction of cells per sector per gene, as a table
    rows = []
    for name, d, col in [("NTC", ntc.assign(sector_stain=ntc["sector"]), "sector_stain")] + [(g, ko[ko["gene"] == g], "sector_stain") for g in genes]:
        s = d[col].dropna()
        s = s[s != ""]
        rows.append(dict(group=name, score="stain", n=len(s), **s.value_counts(normalize=True).reindex(L.SECTORS, fill_value=0).round(4).to_dict()))
    for g in genes:
        s = ko.loc[ko["gene"] == g, "sector_phase"].dropna()
        rows.append(dict(group=g, score="phase_contrast", n=len(s), **s.value_counts(normalize=True).reindex(L.SECTORS, fill_value=0).round(4).to_dict()))
    pd.DataFrame(rows).to_csv(OUT / "sector_fractions_by_gene.csv", index=False)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
