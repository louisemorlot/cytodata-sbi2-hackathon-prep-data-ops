"""Visual evaluation of the transferred cell-cycle phase scores in gene-knockout cells against the NTC reference.

Compares, per knockout gene, the phase distribution of the knockout cells (scored with the frozen NTC reference, see
apply_cc_transfer.py) with the NTC cells, using the phase sectors defined from the NTC curve only.
Expectation from the experiment design (stated by the user): CCNB1 -> G1, MCM3 -> S, CENPC -> G2/M. These are drawn on the
figures but nothing is tuned to them and there is no pass/fail: effect sizes come with bootstrap intervals.

Inputs : reports/ops0077_cc_transfer/reference/, reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv
Outputs: reports/ops0077_cc_transfer/ko3_evaluation/ (figures + tables)

Usage: python evaluate_cc_transfer.py [project_root] [--scores CSV] [--reference DIR] [--out DIR] [--n-boot N] [--n-perm N]
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cc_transfer_lib as L  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--scores", default="reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv")
ap.add_argument("--reference", default="reports/ops0077_cc_transfer/reference")
ap.add_argument("--out", default="reports/ops0077_cc_transfer/ko3_evaluation")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--n-perm", type=int, default=5000)
args = ap.parse_args()
ROOT = Path(args.root)
OUT = ROOT / args.out
OUT.mkdir(parents=True, exist_ok=True)
REF_DIR = ROOT / args.reference
RNG = np.random.default_rng(0)
EXPECTED = {"CCNB1": "G1", "MCM3": "S", "CENPC": "G2/M"}  # stated by the user (CCNB1 -> G1 is unusual: see README notes)
COL = {"NTC": "0.35", "CCNB1": "#1f77b4", "MCM3": "#d62728", "CENPC": "#2ca02c"}
SECT = L.SECTORS
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})


def valid(sectors):
    """Sector labels without missing values (cells with no score are '' in memory and NaN after a CSV round trip)."""
    s = pd.Series(sectors, dtype=object)
    return s[s.notna() & (s != "")]


def fractions(sectors):
    s = valid(sectors)
    return s.value_counts(normalize=True).reindex(SECT, fill_value=0.0).values, len(s)


def boot_ci(sectors, n_boot):
    s = valid(sectors).values.astype(object)
    idx = RNG.integers(0, len(s), (n_boot, len(s)))
    fr = np.array([[np.mean(s[i] == k) for k in SECT] for i in idx])
    return np.percentile(fr, [2.5, 97.5], axis=0)


def perm_test_shift(gene_phase, ntc_phase, n_perm):
    """Permutation test on the distance between the mean unit vectors of the two phase samples."""
    g, n = gene_phase[~np.isnan(gene_phase)], ntc_phase[~np.isnan(ntc_phase)]
    pooled = np.exp(2j * np.pi * np.r_[g, n])
    obs = abs(pooled[:len(g)].mean() - pooled[len(g):].mean())
    cnt = 0
    for _ in range(n_perm):
        perm = RNG.permutation(len(pooled))
        cnt += abs(pooled[perm[:len(g)]].mean() - pooled[perm[len(g):]].mean()) >= obs
    return float(obs), float((cnt + 1) / (n_perm + 1))


def polar_hist(ax, phase, color, label, bins=36, lw=2):
    h, e = np.histogram(2 * np.pi * phase[~np.isnan(phase)], bins=bins, range=(0, 2 * np.pi), density=True)
    th = np.r_[(e[:-1] + e[1:]) / 2, (e[0] + e[1]) / 2 + 2 * np.pi]
    ax.plot(th, np.r_[h, h[0]], color=color, lw=lw, label=label)


def sector_shading(ax, ref, polar=False, alpha=.07):
    cols = {"G1": "#4c78a8", "S": "#f58518", "G2/M": "#e45756", "M-to-G1 exit": "#72b7b2"}
    for s, a, b in L.sector_intervals(ref):
        if polar:
            ax.fill_between(np.linspace(2 * np.pi * a, 2 * np.pi * b, 30), 0, ax.get_ylim()[1], color=cols[s], alpha=alpha, lw=0)
        else:
            ax.axvspan(a, b, color=cols[s], alpha=alpha, lw=0)


def main():
    ref = L.load_reference(REF_DIR / "ntc_cc_reference.json")
    ntc = pd.read_csv(REF_DIR / "ntc_reference_scores.csv", index_col=0)
    prof = pd.read_csv(REF_DIR / "ntc_profiles_by_phase.csv")
    ko = pd.read_csv(ROOT / args.scores, index_col=0)
    genes = [g for g in EXPECTED if g in set(ko["gene"])] + [g for g in sorted(set(ko["gene"])) if g not in EXPECTED]
    print("genes:", genes, "| NTC cells scored:", int(ntc["cc_phase_stain"].notna().sum()))
    ntc_stain = ntc["cc_phase_stain"].values
    ntc_sec = ntc["sector"].values  # NaN where the NTC cell has no complete score

    # ---------------------------------------------------------------- 0. coverage
    cov = ko.groupby("gene").agg(cells=("cc_phase", "size"), has_markers=("has_markers", "sum"), dna_reliable=("dna_reliable", "sum"),
                                 stain_score=("cc_phase_stain", lambda s: int(s.notna().sum())),
                                 phase_score=("cc_phase_phase", lambda s: int(s.notna().sum())),
                                 out_of_ntc_support=("out_of_ntc_support", "sum"))
    cov["fraction_with_stain"] = (cov["stain_score"] / cov["cells"]).round(3)
    cov.to_csv(OUT / "coverage_by_gene.csv")
    print(cov)

    # ---------------------------------------------------------------- 1. phase distributions vs NTC
    fig = plt.figure(figsize=(19, 5.2))
    for k, (score, ttl) in enumerate([("cc_phase_stain", "stain score (4i markers + DNA, frozen NTC curve)"),
                                      ("cc_phase_phase", "phase-contrast score (NTC-trained model)")]):
        ax = fig.add_subplot(1, 3, 1 + k, projection="polar")
        polar_hist(ax, (ntc_stain if score == "cc_phase_stain" else ntc["cc_phase_stain"].values), COL["NTC"], "NTC", lw=2.5) if score == "cc_phase_stain" else None
        for g in genes:
            polar_hist(ax, ko.loc[ko["gene"] == g, score].values, COL.get(g, "k"), f"{g} (n={int(ko.loc[ko['gene'] == g, score].notna().sum())})")
        sector_shading(ax, ref, polar=True)
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_title(ttl, fontsize=10, pad=14)
        ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.12), fontsize=8)
    ax = fig.add_subplot(1, 3, 3)
    xs = np.linspace(0, 1, 400)
    for name, ph in [("NTC", ntc_stain)] + [(g, ko.loc[ko["gene"] == g, "cc_phase_stain"].values) for g in genes]:
        ph = ph[~np.isnan(ph)]
        dens = gaussian_kde(np.r_[ph - 1, ph, ph + 1], bw_method=0.06)(xs) * 3  # wrapped KDE
        ax.plot(xs, dens, color=COL.get(name, "k"), lw=2.5 if name == "NTC" else 1.8, label=name)
    sector_shading(ax, ref); ax.set_xlabel("stain phase (0 = G1 start; 0 and 1 are the same point)"); ax.set_ylabel("density")
    ax.set_title("stain score, linear view (sector shading from NTC)"); ax.legend(fontsize=8)
    plt.tight_layout(); plt.savefig(OUT / "fig1_phase_distributions_vs_ntc.png", dpi=110, bbox_inches="tight"); plt.close()

    # ---------------------------------------------------------------- 2. sector fractions with bootstrap CIs + shift tests
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.6), sharey=True)
    for ax, (score, secol, ttl) in zip(axes, [("cc_phase_stain", "sector_stain", "stain score"), ("cc_phase_phase", "sector_phase", "phase-contrast score")]):
        groups = [("NTC", ntc_sec if score == "cc_phase_stain" else None)] + [(g, ko.loc[ko["gene"] == g, secol].values) for g in genes]
        w = 0.8 / len(groups)
        for gi, (name, sec) in enumerate(groups):
            if sec is None:
                continue
            fr, n = fractions(sec)
            ci = boot_ci(sec, args.n_boot)
            ax.bar(np.arange(len(SECT)) + gi * w, fr, w, yerr=[np.clip(fr - ci[0], 0, None), np.clip(ci[1] - fr, 0, None)], color=COL.get(name, "k"), label=f"{name} (n={n})", capsize=2)
            for si, s in enumerate(SECT):
                rows.append(dict(score=score, group=name, sector=s, n=n, fraction=fr[si], ci_low=ci[0][si], ci_high=ci[1][si]))
        ax.set_xticks(np.arange(len(SECT)) + 0.4 - w / 2); ax.set_xticklabels(SECT)
        for g, s in EXPECTED.items():
            if g in genes and s in SECT:
                ax.annotate(f"expected\n{g}", (SECT.index(s) + (genes.index(g) + 1) * w, 0.0), xytext=(0, -42), textcoords="offset points",
                            ha="center", fontsize=7, color=COL[g], annotation_clip=False)
        ax.set_title(ttl + " (error bars: 95% bootstrap CI)"); ax.legend(fontsize=8)
    axes[0].set_ylabel("fraction of cells")
    plt.tight_layout(); plt.savefig(OUT / "fig2_sector_fractions.png", dpi=110, bbox_inches="tight"); plt.close()
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "sector_fractions_with_ci.csv", index=False)

    shift = []
    for g in genes:
        s = ko.loc[ko["gene"] == g, "cc_phase_stain"].values
        for score, nt in [("cc_phase_stain", ntc_stain)]:
            d, p = perm_test_shift(s, nt, args.n_perm)
            fr_g, n_g = fractions(ko.loc[ko["gene"] == g, "sector_stain"].values)
            fr_n, _ = fractions(ntc_sec)
            ci = boot_ci(ko.loc[ko["gene"] == g, "sector_stain"].values, args.n_boot)
            for si, sct in enumerate(SECT):
                shift.append(dict(gene=g, sector=sct, ntc_fraction=fr_n[si], gene_fraction=fr_g[si], delta=fr_g[si] - fr_n[si],
                                  delta_ci_low=ci[0][si] - fr_n[si], delta_ci_high=ci[1][si] - fr_n[si], n_gene=n_g,
                                  expected_sector=EXPECTED.get(g, ""), mean_vector_shift=d, permutation_p=p))
    shift = pd.DataFrame(shift)
    shift.to_csv(OUT / "sector_shift_vs_ntc.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    hm = shift.pivot(index="gene", columns="sector", values="delta").reindex(index=genes, columns=SECT)
    lo = shift.pivot(index="gene", columns="sector", values="delta_ci_low").reindex(index=genes, columns=SECT)
    hi = shift.pivot(index="gene", columns="sector", values="delta_ci_high").reindex(index=genes, columns=SECT)
    v = np.nanmax(np.abs(hm.values))
    im = ax.imshow(hm.values, cmap="RdBu_r", vmin=-v, vmax=v)
    for i, g in enumerate(genes):
        for j, s in enumerate(SECT):
            ax.text(j, i, f"{hm.values[i, j]:+.2f}\n[{lo.values[i, j]:+.2f}, {hi.values[i, j]:+.2f}]", ha="center", va="center", fontsize=7.5,
                    fontweight="bold" if EXPECTED.get(g) == s else "normal")
    ax.set_xticks(range(len(SECT))); ax.set_xticklabels(SECT); ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"{g}\n(p={shift.loc[shift.gene == g, 'permutation_p'].iloc[0]:.3g})" for g in genes])
    ax.set_title("change in fraction of cells vs NTC (95% CI); bold = stated expectation", fontsize=9); plt.colorbar(im, ax=ax)
    plt.tight_layout(); plt.savefig(OUT / "fig3_sector_shift_heatmap.png", dpi=110); plt.close()

    # ---------------------------------------------------------------- 4. replicate consistency: per guide and per well
    fig, ax = plt.subplots(1, 2, figsize=(15, 4.6))
    # guide = the 20-nt protospacer (the barcode column is truncated to 9 or 10 nt, which would split one guide in two)
    ko["guide"] = ko["gene"] + "|" + ko["sgRNA"].astype(str)
    gt = ko.groupby(["guide", "gene"])["sector_stain"].apply(lambda s: pd.Series(fractions(s)[0], index=SECT)).unstack().reset_index()
    gt["n"] = ko.groupby("guide")["cc_phase_stain"].apply(lambda s: s.notna().sum()).reindex(gt["guide"]).values
    gt = gt[gt["n"] >= 10].sort_values(["gene", "guide"])
    gt.to_csv(OUT / "fractions_per_guide.csv", index=False)
    bottom = np.zeros(len(gt))
    for s, c in zip(SECT, ["#4c78a8", "#f58518", "#e45756", "#72b7b2"]):
        ax[0].bar(range(len(gt)), gt[s], bottom=bottom, color=c, label=s); bottom += gt[s].values
    ax[0].set_xticks(range(len(gt))); ax[0].set_xticklabels([f"{g.split('|')[0]}\n{g.split('|')[1][:8]}\nn={n}" for g, n in zip(gt["guide"], gt["n"])], fontsize=7)
    ax[0].set_title("sector fractions per sgRNA (stain score, guides with >= 10 scored cells)"); ax[0].legend(fontsize=7, ncol=2)
    wt = ko.groupby(["gene", "well"])["sector_stain"].apply(lambda s: pd.Series(fractions(s)[0], index=SECT)).unstack().reset_index()
    wt["n"] = ko.groupby(["gene", "well"])["cc_phase_stain"].apply(lambda s: s.notna().sum()).values
    nw = ntc.groupby("well")["sector"].apply(lambda s: pd.Series(fractions(s)[0], index=SECT)).unstack().reset_index()
    nw["gene"], nw["n"] = "NTC", ntc.groupby("well")["cc_phase_stain"].apply(lambda s: s.notna().sum()).values
    wt = pd.concat([nw, wt]).sort_values(["gene", "well"], key=lambda c: c.map({"NTC": "0"}).fillna(c) if c.name == "gene" else c)
    wt.to_csv(OUT / "fractions_per_well.csv", index=False)
    bottom = np.zeros(len(wt))
    for s, c in zip(SECT, ["#4c78a8", "#f58518", "#e45756", "#72b7b2"]):
        ax[1].bar(range(len(wt)), wt[s], bottom=bottom, color=c); bottom += wt[s].values
    ax[1].set_xticks(range(len(wt))); ax[1].set_xticklabels([f"{g}\n{w}\nn={int(n)}" for g, w, n in zip(wt["gene"], wt["well"], wt["n"])], fontsize=6)
    ax[1].set_title("sector fractions per well (NTC and knockouts)")
    plt.tight_layout(); plt.savefig(OUT / "fig4_replicates_guide_and_well.png", dpi=110); plt.close()

    # ---------------------------------------------------------------- 5. do knockout cells follow the NTC cycle? profiles + support
    fig, ax = plt.subplots(1, 4, figsize=(20, 4))
    for a_, (f, ttl) in zip(ax, [("z_dna", "DNA content (z)"), ("z_phos", "phosphorylated Rb fraction (z)"), ("z_p21", "p21 (z)"), ("z_pRb", "pRb (z)")]):
        a_.fill_between(prof["phase_mid"], prof[f + "_q25"], prof[f + "_q75"], color="0.7", alpha=.4)
        a_.plot(prof["phase_mid"], prof[f + "_q50"], color="0.2", lw=2.5, label="NTC (median, IQR)")
        for g in genes:
            d = ko[(ko["gene"] == g) & ko["cc_phase_stain"].notna()]
            b = pd.cut(d["cc_phase_stain"], np.linspace(0, 1, 16), include_lowest=True)
            m = d.groupby(b, observed=True)[f].median()
            n = d.groupby(b, observed=True).size()
            mid = (np.linspace(0, 1, 16)[:-1] + 1 / 30)[:len(m)]
            a_.plot(mid[n.values >= 5], m.values[n.values >= 5], "o-", ms=3, color=COL[g], label=g)
        a_.set_title(ttl); a_.set_xlabel("stain phase")
    ax[0].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(OUT / "fig5_profiles_knockout_vs_ntc.png", dpi=110); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.2))
    q99 = ref["ntc_curve_dist_q99"]
    xs = np.linspace(0, max(q99 * 2.2, ko["cc_curve_dist"].quantile(.99)), 200)
    for name, dd in [("NTC", ntc["cc_curve_dist"].dropna().values)] + [(g, ko.loc[ko["gene"] == g, "cc_curve_dist"].dropna().values) for g in genes]:
        ax[0].plot(xs, gaussian_kde(dd)(xs), color=COL.get(name, "k"), lw=2.5 if name == "NTC" else 1.8, label=name)
    ax[0].axvline(q99, c="k", ls="--", lw=.8); ax[0].set_xlabel("distance to the NTC cycle curve (dashed: NTC 99th pct)"); ax[0].legend(fontsize=8)
    fr = ko.groupby("gene")["out_of_ntc_support"].mean().reindex(genes)
    ax[1].bar(genes, fr.values, color=[COL[g] for g in genes]); ax[1].axhline(.01, c="k", ls="--", lw=.8)
    ax[1].set_ylabel("fraction beyond NTC 99th pct"); ax[1].set_title("cells off the NTC cycle (dashed: 1% expected)")
    for g in genes:
        d = ko[ko["gene"] == g]
        ax[2].scatter(d["z_dna"], d["z_p21"], s=4, alpha=.45, color=COL[g], label=g)
    sub = ntc.dropna(subset=["cc_phase_stain"]).sample(4000, random_state=0)
    ax[2].scatter(sub["z_dna"], sub["z_p21"], s=2, color="0.75", zorder=0, label="NTC (sample)")
    ax[2].set_xlabel("z_dna"); ax[2].set_ylabel("z_p21"); ax[2].legend(fontsize=7, markerscale=2)
    plt.tight_layout(); plt.savefig(OUT / "fig6_support_and_distance_to_ntc_curve.png", dpi=110); plt.close()

    # ---------------------------------------------------------------- 6. stain vs phase-contrast score agreement
    both = ko.dropna(subset=["cc_phase_stain", "cc_phase_phase"]).copy()
    both["err"] = L.circ_diff(both["cc_phase_phase"], both["cc_phase_stain"])
    agree = both.groupby("gene").apply(lambda d: pd.Series(dict(n=len(d), mae_turns=d["err"].abs().mean(), within_0p10=(d["err"].abs() <= .1).mean(),
                                                                circ_agreement=abs(np.exp(2j * np.pi * d["err"]).mean())))).round(3)
    agree.to_csv(OUT / "stain_vs_phase_agreement.csv")
    fig, ax = plt.subplots(1, len(genes) + 1, figsize=(5 * (len(genes) + 1), 4.3))
    for a_, g in zip(ax, genes):
        d = both[both["gene"] == g]
        a_.scatter(d["cc_phase_stain"], d["cc_phase_phase"], s=5, alpha=.5, color=COL[g])
        for off in (-1, 0, 1):
            a_.plot([0, 1], [off, 1 + off], "k--", lw=.6)
        a_.set_xlim(0, 1); a_.set_ylim(0, 1); a_.set_xlabel("stain score"); a_.set_ylabel("phase-contrast score")
        a_.set_title(f"{g}: MAE {agree.loc[g, 'mae_turns']:.3f}, n={len(d)}")
    for g in genes:
        ax[-1].hist(both.loc[both["gene"] == g, "err"], bins=40, histtype="step", color=COL[g], lw=2, label=g, density=True)
    ax[-1].set_xlabel("circular error: phase-contrast minus stain (turns)"); ax[-1].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(OUT / "fig7_stain_vs_phase_contrast_scores.png", dpi=110); plt.close()

    # ---------------------------------------------------------------- 7. comparison with the reference columns shipped in the input file
    fig, ax = plt.subplots(1, 3, figsize=(19, 4.8))
    for a_, (rc, ttl) in zip(ax[:2], [("ref_cycle_stage", "ref_cycle_stage (DNA x pRb gating)"), ("ref_phase_predicted_stage", "ref_phase_predicted_stage")]):
        if rc not in ko:
            a_.axis("off"); continue
        d = ko[(ko["sector_stain"] != "")]
        ct = pd.crosstab(d[rc].astype(str), d["sector_stain"], normalize="index").reindex(columns=SECT, fill_value=0)
        n = d[rc].astype(str).value_counts().reindex(ct.index)
        im = a_.imshow(ct.values, cmap="viridis", vmin=0, vmax=1)
        for i in range(ct.shape[0]):
            for j in range(ct.shape[1]):
                a_.text(j, i, f"{ct.values[i, j]:.2f}", ha="center", va="center", color="w", fontsize=8)
        a_.set_xticks(range(len(SECT))); a_.set_xticklabels(SECT, fontsize=8); a_.set_yticks(range(ct.shape[0]))
        a_.set_yticklabels([f"{i} (n={int(n[i])})" for i in ct.index], fontsize=8)
        a_.set_title(f"our stain sector (columns) within {ttl}", fontsize=9); a_.set_xlabel("row-normalised")
        ct.to_csv(OUT / f"crosstab_{rc}_vs_our_sector.csv")
    if "ref_phase_predicted_cc_phase_ridge" in ko:
        d = ko.dropna(subset=["cc_phase_phase", "ref_phase_predicted_cc_phase_ridge"])
        err = L.circ_diff(d["cc_phase_phase"], d["ref_phase_predicted_cc_phase_ridge"])
        ax[2].hexbin(d["ref_phase_predicted_cc_phase_ridge"], d["cc_phase_phase"], gridsize=35, cmap="viridis", mincnt=1, bins="log")
        ax[2].set_xlabel("ref_phase_predicted_cc_phase_ridge (in file)"); ax[2].set_ylabel("our cc_phase_phase")
        ax[2].set_title(f"two phase-contrast predictors (different NTC clocks): MAE {np.abs(err).mean():.3f} turns")
    plt.tight_layout(); plt.savefig(OUT / "fig8_comparison_with_file_reference_columns.png", dpi=110); plt.close()

    # ---------------------------------------------------------------- 8. DNA content reference check (QC only; not used for any score)
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.2), sharey=True)
    xs = np.linspace(-1, 2, 300)
    peak_rows = []
    for a_, w in zip(ax, sorted(ko["well"].unique())):
        for name, dna in [("NTC", ntc.loc[ntc["well"] == w, "dna_rel"].dropna().values)] + [(g, ko.loc[(ko["gene"] == g) & (ko["well"] == w), "dna_content"].dropna().values) for g in genes]:
            if len(dna) < 20:
                continue
            kd = gaussian_kde(dna, bw_method=0.15)(xs)
            low = (xs > -0.6) & (xs < 0.35)
            peak_rows.append(dict(well=w, group=name, n=len(dna), g1_peak_mode=float(xs[low][np.argmax(kd[low])]), median=float(np.median(dna))))
            a_.plot(xs, kd, color=COL.get(name, "k"), lw=2.5 if name == "NTC" else 1.5, label=f"{name} (n={len(dna)})")
        a_.set_title(f"well {w}"); a_.set_xlabel("relative DNA content (log2)"); a_.legend(fontsize=7)
    ax[0].set_ylabel("density")
    plt.suptitle("DNA content reference check: the lowest-DNA (G1, 2N) peak should sit at the same place in NTC and knockout cells", fontsize=10)
    plt.tight_layout(); plt.savefig(OUT / "fig9_dna_reference_check.png", dpi=110); plt.close()
    pk = pd.DataFrame(peak_rows)
    pk.to_csv(OUT / "dna_g1_peak_check.csv", index=False)
    print(pk.round(3).to_string())

    summ = dict(genes=genes, sectors=[list(x) for x in L.sector_intervals(ref)],
                shift_vs_ntc=shift[["gene", "sector", "ntc_fraction", "gene_fraction", "delta", "delta_ci_low", "delta_ci_high", "permutation_p"]].round(4).to_dict("records"),
                stain_vs_phase=agree.reset_index().to_dict("records"))
    json.dump(summ, open(OUT / "evaluation_summary.json", "w"), indent=2, default=str)
    print(shift[["gene", "sector", "ntc_fraction", "gene_fraction", "delta", "delta_ci_low", "delta_ci_high", "permutation_p"]].round(3).to_string())
    print(agree)


if __name__ == "__main__":
    main()
