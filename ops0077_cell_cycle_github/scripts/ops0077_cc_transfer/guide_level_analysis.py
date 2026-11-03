"""Guide-level analysis of the transferred cell-cycle phase in gene-knockout cells.

Why: the gene-level average of the knockout cells shows no clear arrest, but guides may differ (some effective, some not).
A guide is the 20-nt protospacer (obs['sgRNA']); the barcode column is truncated to 9 or 10 nt and would split one guide in two.

Questions answered
  1. Per guide: sector fractions and the change against NTC, standardised to the guide's own well mix (NTC cells re-weighted
     to the wells in which the guide's cells sit), with bootstrap intervals; pooled circular-shift test with FDR.
  2. Do guides of the same gene differ more than chance allows? (chi-square, null from permuting guide labels within wells)
  3. Is a guide effect real? Compare the same guide's phase composition between two DISJOINT sets of cells measured by
     different instruments: cells with a stain score (stain sectors) vs cells without one (phase-contrast sectors).
  4. Do other phenotypes (nuclear size, p21, DNA beyond 4N, off-cycle cells) rank the guides the same way?
  5. Statistical power: how large a change can a guide of n cells reveal?

Inputs : reports/ops0077_cc_transfer/reference/, reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv.gz
Outputs: reports/ops0077_cc_transfer/ko3_guide_analysis/ (figures g1..g6 + tables)

Usage: python guide_level_analysis.py [project_root] [--scores CSV] [--reference DIR] [--out DIR] [--min-cells N] [--n-boot N] [--n-perm N]
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
from scipy.stats import gaussian_kde, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cc_transfer_lib as L  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--scores", default="reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv.gz")
ap.add_argument("--reference", default="reports/ops0077_cc_transfer/reference")
ap.add_argument("--out", default="reports/ops0077_cc_transfer/ko3_guide_analysis")
ap.add_argument("--min-cells", type=int, default=10, help="guides with fewer scored cells are not analysed")
ap.add_argument("--n-boot", type=int, default=2000)
ap.add_argument("--n-perm", type=int, default=5000)
args = ap.parse_args()
ROOT = Path(args.root)
OUT = ROOT / args.out
OUT.mkdir(parents=True, exist_ok=True)
REF_DIR = ROOT / args.reference
RNG = np.random.default_rng(0)
SECT = L.SECTORS
EXPECTED = {"CCNB1": "G1", "MCM3": "S", "CENPC": "G2/M"}  # stated by the user
COL = {"NTC": "0.35", "CCNB1": "#1f77b4", "MCM3": "#d62728", "CENPC": "#2ca02c"}
SECCOL = {"G1": "#4c78a8", "S": "#f58518", "G2/M": "#e45756", "M-to-G1 exit": "#72b7b2"}
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 110})


def valid(x):
    s = pd.Series(x, dtype=object)
    return s[s.notna() & (s != "")]


def frac(sectors):
    s = valid(sectors)
    return s.value_counts(normalize=True).reindex(SECT, fill_value=0.0).values, len(s)


def short(guide):
    g, seq = guide.split("|")
    return f"{g}|{seq[:8]}"


def benjamini_hochberg(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def standardised_delta(sec, wells, ntc_w):
    """Guide sector fractions minus the NTC fractions re-weighted to the guide's own well mix. ntc_w: well -> fraction vector."""
    s = np.asarray(sec, object)
    w = np.asarray(wells, object)
    ok = pd.notna(s) & (s != "")
    s, w = s[ok], w[ok]
    f = np.array([np.mean(s == k) for k in SECT])
    comp = pd.Series(w).value_counts(normalize=True)
    exp = sum(comp[k] * ntc_w[k] for k in comp.index if k in ntc_w)
    return f, f - exp, exp


def boot_delta(sec, wells, ntc_w, n_boot):
    s = np.asarray(sec, object)
    w = np.asarray(wells, object)
    ok = pd.notna(s) & (s != "")
    s, w = s[ok], w[ok]
    out = np.empty((n_boot, len(SECT)))
    for b in range(n_boot):
        i = RNG.integers(0, len(s), len(s))
        out[b] = standardised_delta(s[i], w[i], ntc_w)[1]
    return np.percentile(out, [2.5, 97.5], axis=0)


def mean_vec_test(g_phase, ntc_phase, n_perm):
    g, n = g_phase[~np.isnan(g_phase)], ntc_phase[~np.isnan(ntc_phase)]
    pooled = np.exp(2j * np.pi * np.r_[g, n])
    obs = abs(pooled[:len(g)].mean() - pooled[len(g):].mean())
    cnt = sum(abs(pooled[p[:len(g)]].mean() - pooled[p[len(g):]].mean()) >= obs for p in (RNG.permutation(len(pooled)) for _ in range(n_perm)))
    return float(obs), float((cnt + 1) / (n_perm + 1))


def chi2_stat(table):
    t = np.asarray(table, float)
    e = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.nansum((t - e) ** 2 / e))


def heterogeneity(df, sec_col, n_perm):
    """Chi-square across the guides of one gene; null from permuting guide labels within wells."""
    d = df[df[sec_col].notna() & (df[sec_col] != "")]
    obs_t = pd.crosstab(d["guide"], d[sec_col]).reindex(columns=SECT, fill_value=0)
    obs = chi2_stat(obs_t.values)
    guides = d["guide"].values
    codes = pd.Categorical(d[sec_col], categories=SECT).codes
    wells = d["well"].values
    idx_by_well = [np.where(wells == w)[0] for w in np.unique(wells)]
    ng = len(obs_t)
    gcode = pd.Categorical(guides, categories=obs_t.index).codes
    cnt = 0
    for _ in range(n_perm):
        pg = gcode.copy()
        for ix in idx_by_well:
            pg[ix] = gcode[ix][RNG.permutation(len(ix))]
        t = np.zeros((ng, len(SECT)))
        np.add.at(t, (pg, codes), 1)
        cnt += chi2_stat(t) >= obs
    return obs, float((cnt + 1) / (n_perm + 1)), obs_t


def main():
    ref = L.load_reference(REF_DIR / "ntc_cc_reference.json")
    ntc = pd.read_csv(REF_DIR / "ntc_reference_scores.csv.gz", index_col=0)
    ko = pd.read_csv(ROOT / args.scores, index_col=0)
    ko["guide"] = ko["gene"] + "|" + ko["sgRNA"].astype(str)
    genes = [g for g in EXPECTED if g in set(ko["gene"])]
    ko["has_stain"] = ko["cc_phase_stain"].notna()
    n_scored = ko.groupby("guide")["cc_phase"].size()
    guides = [g for g in sorted(ko["guide"].unique()) if n_scored[g] >= args.min_cells]
    dropped = [g for g in ko["guide"].unique() if g not in guides]
    print(f"guides analysed: {len(guides)} (dropped {dropped})")
    ko = ko[ko["guide"].isin(guides)].copy()

    ntc_sc = ntc[ntc["sector"].notna() & (ntc["sector"] != "")]
    ntc_w = {w: frac(d["sector"])[0] for w, d in ntc_sc.groupby("well")}
    ntc_all = frac(ntc_sc["sector"])[0]

    # ------------------------------------------------------------------ 1. per-guide sector shifts (stain score), well-standardised
    rows = []
    for g in guides:
        d = ko[ko["guide"] == g]
        ds = d[d["has_stain"]]
        f, delta, exp = standardised_delta(ds["sector_stain"].values, ds["well"].values, ntc_w)
        ci = boot_delta(ds["sector_stain"].values, ds["well"].values, ntc_w, args.n_boot)
        obs_shift, p = mean_vec_test(ds["cc_phase_stain"].values, ntc["cc_phase_stain"].values, args.n_perm)
        # raw fractions with all cells (phase-contrast sectors for the cells without a stain score)
        f_all, n_all = frac(d["sector_final"])
        f_ph, n_ph = frac(d["sector_phase"])
        for i, s in enumerate(SECT):
            rows.append(dict(guide=g, gene=d["gene"].iloc[0], sector=s, n_stain=len(ds), n_all=n_all, frac_stain=f[i], ntc_expected=exp[i],
                             delta=delta[i], ci_low=ci[0][i], ci_high=ci[1][i], frac_final_all_cells=f_all[i], frac_phase_score_all_cells=f_ph[i],
                             mean_vector_shift=obs_shift, perm_p=p))
    gs = pd.DataFrame(rows)
    first = gs.groupby("guide").first()
    q = pd.Series(benjamini_hochberg(first.loc[guides, "perm_p"].values), index=guides)
    gs["fdr_q"] = gs["guide"].map(q)
    gs["short"] = gs["guide"].map(short)
    gs.to_csv(OUT / "guide_sector_shifts.csv", index=False)
    wide = gs.pivot(index="guide", columns="sector", values="delta")[SECT].round(3)
    wide.insert(0, "n_stain", first.loc[wide.index, "n_stain"].astype(int))
    wide["perm_p"], wide["fdr_q"] = first.loc[wide.index, "perm_p"].round(4), q.round(4)
    print(wide.to_string())

    # ------------------------------------------------------------------ 2. heterogeneity across guides within a gene
    het = []
    for g in genes:
        for score, col in [("stain", "sector_stain"), ("all_cells_final", "sector_final")]:
            d = ko[(ko["gene"] == g)]
            chi, p, t = heterogeneity(d, col, args.n_perm)
            het.append(dict(gene=g, score=score, n_guides=len(t), n_cells=int(t.values.sum()), chi2=chi, perm_p_within_wells=p))
    het = pd.DataFrame(het)
    het.to_csv(OUT / "heterogeneity_across_guides.csv", index=False)
    print(het.round(4).to_string())

    # ------------------------------------------------------------------ 3. is a guide effect real? disjoint cell sets, different instruments
    a = ko[ko["has_stain"]].groupby("guide")["sector_stain"].apply(lambda s: pd.Series(frac(s)[0], index=SECT)).unstack()
    b = ko[~ko["has_stain"]].groupby("guide")["sector_phase"].apply(lambda s: pd.Series(frac(s)[0], index=SECT)).unstack()
    nb = ko[~ko["has_stain"]].groupby("guide").size().reindex(guides).fillna(0)
    na = ko[ko["has_stain"]].groupby("guide").size().reindex(guides)
    c = ko[ko["has_stain"]].groupby("guide")["sector_phase"].apply(lambda s: pd.Series(frac(s)[0], index=SECT)).unstack()
    conc = []
    for s in SECT:
        ok = nb[guides] >= args.min_cells
        r1, p1 = spearmanr(a.loc[ok[ok].index, s], b.loc[ok[ok].index, s])
        r2, p2 = spearmanr(a[s], c.loc[a.index, s])
        conc.append(dict(sector=s, n_guides_disjoint=int(ok.sum()), spearman_stain_cells_vs_phase_score_other_cells=r1, p_disjoint=p1,
                         spearman_stain_vs_phase_score_same_cells=r2, p_same=p2))
    conc = pd.DataFrame(conc)
    conc.to_csv(OUT / "guide_effect_concordance.csv", index=False)
    print(conc.round(3).to_string())

    # ------------------------------------------------------------------ 4. other phenotypes per guide
    def phen(d):
        return pd.Series(dict(n=len(d), median_log_nuclear_area=np.log(d["nuc_area_px"]).median(), median_z_p21=d["z_p21"].median(),
                              pct_p21_high_z_gt_2=100 * (d["z_p21"] > 2).mean(), pct_dna_beyond_4N_z_ge_3p5=100 * (d["z_dna"] >= 3.5).mean(),
                              pct_off_ntc_cycle=100 * d["out_of_ntc_support"].mean(), median_z_phos=d["z_phos"].median()))
    ph = ko[ko["has_markers"]].groupby("guide").apply(phen)
    ntc_row = phen(ntc.assign(out_of_ntc_support=ntc["cc_curve_dist"] > ref["ntc_curve_dist_q99"]))
    ph.loc["NTC"] = ntc_row
    ph.to_csv(OUT / "phenotypes_by_guide.csv")
    g2m = wide["G2/M"]
    pcorr = pd.DataFrame({m: [spearmanr(g2m.loc[ph.index.intersection(g2m.index)], ph.loc[ph.index.intersection(g2m.index), m])[0]]
                          for m in ph.columns if m != "n"}, index=["spearman_with_delta_G2/M_across_guides"]).T
    pcorr.to_csv(OUT / "phenotype_vs_g2m_shift_across_guides.csv")
    print(pcorr.round(3).to_string())

    # ------------------------------------------------------------------ 5. power
    pw = pd.DataFrame({"n_stain_cells": na.reindex(guides)})
    pw["halfwidth_95CI_at_fraction_0.25"] = 1.96 * np.sqrt(0.25 * 0.75 / pw["n_stain_cells"])
    pw["halfwidth_95CI_at_fraction_0.63"] = 1.96 * np.sqrt(0.63 * 0.37 / pw["n_stain_cells"])
    pw.round(3).to_csv(OUT / "power_by_guide.csv")

    plot_all(ko, gs, wide, first, ntc, ref, genes, guides, a, b, c, na, nb, ph, het)
    summ = dict(guides=guides, dropped=dropped, heterogeneity=het.round(4).to_dict("records"), concordance=conc.round(3).to_dict("records"),
                guide_shift=wide.reset_index().to_dict("records"))
    json.dump(summ, open(OUT / "guide_analysis_summary.json", "w"), indent=2, default=str)


def plot_all(ko, gs, wide, first, ntc, ref, genes, guides, a, b, c, na, nb, ph, het):
    # ---- g1: forest plots of the standardised change per guide, rows = sectors, columns = genes
    sects = ["G1", "S", "G2/M"]
    fig, ax = plt.subplots(len(sects), len(genes), figsize=(5.2 * len(genes), 3.1 * len(sects)), sharex="row")
    for si, s in enumerate(sects):
        for gi, g in enumerate(genes):
            a_ = ax[si, gi]
            d = gs[(gs["gene"] == g) & (gs["sector"] == s)].reset_index(drop=True)
            y = np.arange(len(d))
            a_.errorbar(d["delta"], y, xerr=[d["delta"] - d["ci_low"], d["ci_high"] - d["delta"]], fmt="o", color=COL[g], capsize=3, label="stain-scored cells")
            a_.plot(d["frac_final_all_cells"] - d["ntc_expected"], y + .22, "s", mfc="none", color=COL[g], ms=5, label="all cells (phase contrast where no stain)")
            a_.axvline(0, c="k", lw=.8)
            a_.set_yticks(y); a_.set_yticklabels([f"{x}  n={n}" for x, n in zip(d["short"].str.split("|").str[1], d["n_stain"])], fontsize=7); a_.invert_yaxis()
            a_.set_title(f"{g}: {s}" + ("  (stated expectation)" if EXPECTED[g] == s else ""), fontsize=9, fontweight="bold" if EXPECTED[g] == s else "normal")
            if si == len(sects) - 1:
                a_.set_xlabel("fraction minus NTC (well-matched), 95% bootstrap CI")
    ax[0, 0].legend(fontsize=6.5, loc="lower right")
    plt.tight_layout(); plt.savefig(OUT / "g1_guide_shift_forest.png", dpi=110); plt.close()

    # ---- g2: density along the cycle per guide (stain score on top, phase-contrast score for all cells below)
    xs = np.linspace(0, 1, 300)
    fig, ax = plt.subplots(2, len(genes), figsize=(5.4 * len(genes), 7), sharey="row")
    nt = ntc["cc_phase_stain"].dropna().values
    kd0 = gaussian_kde(np.r_[nt - 1, nt, nt + 1], bw_method=.06)(xs) * 3
    for gi, g in enumerate(genes):
        for ri, (score, ttl) in enumerate([("cc_phase_stain", "stain score"), ("cc_phase_phase", "phase-contrast score (all cells)")]):
            a_ = ax[ri, gi]
            if ri == 0:
                a_.plot(xs, kd0, color="0.35", lw=3.5, label="NTC")
            for k, gd in enumerate([x for x in guides if x.startswith(g + "|")]):
                ph_ = ko.loc[ko["guide"] == gd, score].dropna().values
                if len(ph_) < 10:
                    continue
                a_.plot(xs, gaussian_kde(np.r_[ph_ - 1, ph_, ph_ + 1], bw_method=.09)(xs) * 3, lw=1.6, color=plt.cm.tab10(k), label=f"{gd.split('|')[1][:8]} (n={len(ph_)})")
            for s_, x0, x1 in ref["sector_intervals"]:
                a_.axvspan(x0, x1, alpha=.06, color=SECCOL[s_])
            a_.set_title(f"{g}: {ttl}", fontsize=9); a_.set_xlabel("phase (0 and 1 are the same point)"); a_.legend(fontsize=6.5)
    ax[0, 0].set_ylabel("density"); ax[1, 0].set_ylabel("density")
    plt.tight_layout(); plt.savefig(OUT / "g2_guide_density_along_cycle.png", dpi=110); plt.close()

    # ---- g3: is the guide effect real? same guide, two disjoint cell sets
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
    for a_, s in zip(ax, ["G1", "G2/M", "S"]):
        okg = [x for x in guides if nb[x] >= 10]
        for x in okg:
            gene = x.split("|")[0]
            a_.scatter(a.loc[x, s], b.loc[x, s], s=18 + nb[x] / 2, color=COL[gene], alpha=.8)
            a_.annotate(x.split("|")[1][:5], (a.loc[x, s], b.loc[x, s]), fontsize=6, xytext=(3, 3), textcoords="offset points")
        r, p = spearmanr(a.loc[okg, s], b.loc[okg, s])
        lim = [0, max(a.loc[okg, s].max(), b.loc[okg, s].max()) * 1.1]
        a_.plot(lim, lim, "k:", lw=.8); a_.set_xlabel(f"fraction in {s}: cells WITH a stain score (stain sectors)")
        a_.set_ylabel(f"fraction in {s}: other cells (phase-contrast sectors)"); a_.set_title(f"{s}: Spearman {r:.2f} (p={p:.2g}, {len(okg)} guides)", fontsize=9)
    for gname in genes:
        ax[0].scatter([], [], color=COL[gname], label=gname)
    ax[0].legend(fontsize=8)
    plt.tight_layout(); plt.savefig(OUT / "g3_guide_effect_concordance_disjoint_cells.png", dpi=110); plt.close()

    # ---- g4: other phenotypes per guide, ranked, against NTC
    metrics = [("median_log_nuclear_area", "median log nuclear area"), ("median_z_p21", "median p21 (z)"),
               ("pct_p21_high_z_gt_2", "% cells p21 z > 2"), ("pct_dna_beyond_4N_z_ge_3p5", "% DNA content beyond 4N (z >= 3.5)"),
               ("pct_off_ntc_cycle", "% off the NTC cycle")]
    fig, ax = plt.subplots(1, len(metrics) + 1, figsize=(4.4 * (len(metrics) + 1), 5.4), sharey=True)
    order = list(wide.sort_values("G2/M").index)
    y = np.arange(len(order))
    ax[0].barh(y, wide.loc[order, "G2/M"], color=[COL[o.split("|")[0]] for o in order]); ax[0].axvline(0, c="k", lw=.8)
    ax[0].set_yticks(y); ax[0].set_yticklabels([f"{short(o).split('|')[0]} {short(o).split('|')[1]} (n={int(wide.loc[o, 'n_stain'])})" for o in order], fontsize=7)
    ax[0].set_title("change in G2/M fraction vs NTC", fontsize=9)
    for a_, (m, t) in zip(ax[1:], metrics):
        a_.barh(y, ph.loc[order, m], color=[COL[o.split("|")[0]] for o in order])
        a_.axvline(ph.loc["NTC", m], c="k", ls="--", lw=1); a_.set_title(t + "\n(dashed: NTC)", fontsize=9)
        lo, hi = ph.loc[order, m].min(), ph.loc[order, m].max(); a_.set_xlim(min(lo, ph.loc["NTC", m]) - .1 * (hi - lo + 1e-9), hi + .1 * (hi - lo + 1e-9))
    plt.tight_layout(); plt.savefig(OUT / "g4_guide_phenotypes_ranked_by_g2m_shift.png", dpi=110); plt.close()

    # ---- g5: heterogeneity tests
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    het.pivot(index="gene", columns="score", values="perm_p_within_wells").reindex(genes).plot.bar(ax=ax[0])
    ax[0].axhline(.05, c="k", ls="--", lw=.8); ax[0].set_ylabel("permutation p (guides differ?)"); ax[0].set_yscale("log"); ax[0].set_title("do guides of one gene differ more than chance?")
    ax[0].tick_params(axis="x", rotation=0)
    n_tab = ko[ko["has_stain"]].groupby("gene")["guide"].apply(lambda s: s.value_counts().sort_index().values)
    ax[1].axis("off")
    ax[1].text(0, .9, "stratified permutation: guide labels shuffled within wells\n(removes well composition as an explanation)\n\n" +
               "\n".join(f"{g}: stain-scored cells per guide {list(n_tab[g])}" for g in genes), fontsize=9, va="top")
    plt.tight_layout(); plt.savefig(OUT / "g5_heterogeneity_tests.png", dpi=110); plt.close()

    # ---- g6: guide x well composition (confound check)
    fig, ax = plt.subplots(figsize=(10, 3.6))
    ct = pd.crosstab(ko["guide"], ko["well"], normalize="index").reindex(guides)
    ct.plot.bar(stacked=True, ax=ax, width=.8); ax.set_xticklabels([short(x) for x in ct.index], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("fraction of the guide's cells"); ax.set_title("well composition of each guide's cells (NTC fractions are re-weighted to this mix)")
    plt.tight_layout(); plt.savefig(OUT / "g6_guide_well_composition.png", dpi=110); plt.close()


if __name__ == "__main__":
    main()
