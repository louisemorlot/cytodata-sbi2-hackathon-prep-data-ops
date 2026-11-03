"""Continuous (open-ended) cell-cycle stage score from p21 / Rb / pRb intensities only, NO DAPI / DNA content.

Reference analysis to show later what DNA content adds (see ops0077_cell_cycle_score.py).

Stage 1: intensity-only score (open principal curve through biology-oriented axes).
Stage 2: add per-channel / pairwise PCs that track the score without being batch or size driven.
All outputs carry the suffix `_only_markers` and go to reports/ops0077_cell_cycle_score_only_markers/.

Usage: python ops0077_cell_cycle_score_only_markers.py [project_root]
"""
import json
import sys
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from scipy.stats import kruskal, spearmanr
from sklearn.mixture import GaussianMixture

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
H5AD = ROOT / "data/altair/ops0077_4i/OPS0077_NTC_OrganelleProfiler_pRb_Rb_p21.h5ad"
PCS = ROOT / "reports/ops0077_cell_cycle_clock/per_cell_pca_embeddings.csv"
OUT = ROOT / "reports/ops0077_cell_cycle_score_only_markers"
OUT.mkdir(parents=True, exist_ok=True)
SUF = "_only_markers"
RNG = np.random.default_rng(0)
MARKERS = ["pRb", "Rb", "p21"]


# ---------------------------------------------------------------- features
def load_intensities():
    a = ad.read_h5ad(H5AD)
    v = list(a.var.index)
    g = lambda n: np.asarray(a.X[:, v.index(n)], dtype=float)
    df = pd.DataFrame(index=a.obs_names)
    df["well"] = a.obs["well_canonical"].values
    for m in MARKERS:
        p = f"cp_single_object_{m}"
        cell_int, cell_mean = g(f"{p}_cell_Intensity_IntegratedIntensity"), g(f"{p}_cell_Intensity_MeanIntensity")
        cyt_int, cyt_mean = g(f"{p}_cytoplasm_Intensity_IntegratedIntensity"), g(f"{p}_cytoplasm_Intensity_MeanIntensity")
        nuc_area = cell_int / cell_mean - cyt_int / cyt_mean  # px
        df[f"{m}_nuc"] = (cell_int - cyt_int) / nuc_area  # nuclear mean intensity
        df[f"{m}_cyt"] = cyt_mean
        if m == "pRb":
            df["nuc_area_px"] = nuc_area
    df.loc[~(df["nuc_area_px"] > 50), [f"{m}_nuc" for m in MARKERS] + ["nuc_area_px"]] = np.nan
    return df


def preprocess(df):
    """Per-well background subtraction, log, robust z. Returns log-intensity and z axes."""
    logi = pd.DataFrame(index=df.index)
    for m in MARKERS:
        col = df[f"{m}_nuc"]
        # background = per-well 0.5th percentile of the nuclear mean; floor at the per-well 1st percentile of the
        # positive remainder so a few dim cells cannot pile up at log(0)
        bg = col.groupby(df["well"]).transform(lambda x: np.nanpercentile(x, 0.5))
        sig = col - bg
        floor = sig.where(sig > 0).groupby(df["well"]).transform(lambda x: np.nanpercentile(x, 1))
        sig = sig.clip(lower=floor)
        logi[m] = np.log(sig)
    z = pd.DataFrame(index=df.index)
    z["phos"] = logi["pRb"] - logi["Rb"]  # phosphorylated fraction of Rb
    z["p21"] = logi["p21"]
    z["pRb"] = logi["pRb"]
    for c in z:
        grp = z[c].groupby(df["well"])
        med = grp.transform("median")
        mad = (z[c] - med).abs().groupby(df["well"]).transform("median") * 1.4826
        z[c] = ((z[c] - med) / mad).clip(-4, 4)
    return logi, z


# ---------------------------------------------------------------- open principal curve
def principal_curve(X, init_dir, weights=None, n_nodes=12, n_iter=30, smooth=4):
    """Open principal curve; returns nodes and weights."""
    w = np.ones(X.shape[1]) if weights is None else np.asarray(weights)
    Xw = X * w
    mu = Xw.mean(0)
    t = (Xw - mu) @ (init_dir * w) / np.linalg.norm(init_dir * w)
    for _ in range(n_iter):
        t = (t - t.min()) / (t.max() - t.min())
        # equal-count bins by rank (robust to ties at the curve ends)
        rank = np.argsort(np.argsort(t, kind="stable"), kind="stable")
        idx = np.minimum(rank * n_nodes // len(t), n_nodes - 1)
        nodes = np.array([Xw[idx == k].mean(0) for k in range(n_nodes)])
        for _ in range(smooth):  # light smoothing, keep endpoints
            nodes[1:-1] = 0.25 * nodes[:-2] + 0.5 * nodes[1:-1] + 0.25 * nodes[2:]
        t = project(Xw, nodes)
    return nodes, w


def project(Xw, nodes):
    seg = np.diff(nodes, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    cum = np.r_[0, np.cumsum(seglen)]
    best_d = np.full(len(Xw), np.inf)
    best_t = np.zeros(len(Xw))
    for i in range(len(seg)):
        s = seg[i]
        u = np.clip(((Xw - nodes[i]) @ s) / (s @ s + 1e-12), 0, 1)
        d = np.linalg.norm(Xw - (nodes[i] + u[:, None] * s), axis=1)
        m = d < best_d
        best_d[m] = d[m]
        best_t[m] = cum[i] + u[m] * seglen[i]
    return best_t / cum[-1]


def fit_score(X, orient_dir, weights=None, fit_n=20000):
    sub = RNG.choice(len(X), min(fit_n, len(X)), replace=False)
    nodes, w = principal_curve(X[sub], orient_dir, weights)
    t = project(X * w, nodes)
    # percentile-rank so the score is uniform-ish on [0, 1] (comparable across fits)
    return pd.Series(t).rank(pct=True).values, t, nodes


# ---------------------------------------------------------------- soft states (biology prior)
STATES = ["G0/early G1", "late G1", "S", "G2/M"]
# prior means in z space of [phos, p21, pRb]
PRIOR = np.array([[-1.0, 1.0, -1.0], [0.0, 0.0, 0.0], [0.8, -0.8, 0.8], [1.5, 0.3, 1.5]])


def soft_states(Z):
    gm = GaussianMixture(4, covariance_type="full", means_init=PRIOR, random_state=0, max_iter=300).fit(Z)
    return gm.predict_proba(Z), gm.means_


def eta2(x, groups):
    ok = ~np.isnan(x)
    g = [x[ok & (groups == k)] for k in np.unique(groups)]
    h = kruskal(*g).statistic
    return (h - len(g) + 1) / (ok.sum() - len(g))


def main():
    df = load_intensities()
    logi, z = preprocess(df)
    ok = z.notna().all(axis=1) & df["nuc_area_px"].notna()
    print(f"cells total {len(df)}, complete intensity {ok.sum()}")
    Z = z[ok].values
    well = df.loc[ok, "well"].values
    area = np.log(df.loc[ok, "nuc_area_px"].values)

    # ---- stage 1: intensity-only
    orient = np.array([1.0, 0.0, 1.0])  # pRb/phos increase along the cycle; p21 left free (non-monotonic)
    s1, t1, nodes1 = fit_score(Z, orient)
    if spearmanr(s1, Z[:, 0]).statistic < 0:
        s1 = 1 - s1
    post, gm_means = soft_states(Z)
    stage_idx = post @ np.arange(4) / 3
    print("GMM means (phos,p21,pRb):\n", np.round(gm_means, 2))
    print("score vs soft-state expected stage rho:", round(spearmanr(s1, stage_idx).statistic, 3))
    print("score vs nuclear area rho:", round(spearmanr(s1, area).statistic, 3))

    res = pd.DataFrame(index=z.index)
    res["cc_score_intensity"] = np.nan
    res.loc[ok, "cc_score_intensity"] = s1
    for k, s in enumerate(STATES):
        res[f"p_{s}"] = np.nan
        res.loc[ok, f"p_{s}"] = post[:, k]
    res["soft_stage"] = np.nan
    res.loc[ok, "soft_stage"] = stage_idx
    res = res.join(z.add_prefix("z_")).join(df[["well", "nuc_area_px"]])

    # ---- stage 2: PC layer
    pcs = pd.read_csv(PCS, index_col=0).reindex(df.index)
    rows = []
    for c in pcs.columns:
        x = pcs.loc[ok, c].values
        m = ~np.isnan(x)
        if m.sum() < 5000:
            rows.append(dict(pc=c, n=m.sum()))
            continue
        rows.append(dict(
            pc=c, n=int(m.sum()), coverage=m.mean(),
            rho_score=spearmanr(x[m], s1[m]).statistic,
            rho_area=spearmanr(x[m], area[m]).statistic,
            eta2_well=eta2(x, well),
        ))
    scr = pd.DataFrame(rows)
    scr["keep"] = (scr["coverage"] > 0.5) & (scr["rho_score"].abs() >= 0.15) & (scr["eta2_well"] < 0.10) & (scr["rho_area"].abs() < 0.6)
    scr.to_csv(OUT / f"pc_screening{SUF}.csv", index=False)
    print(scr.round(3).to_string())
    keep = scr.loc[scr["keep"], "pc"].tolist()
    print("kept PCs:", keep)

    res["cc_score_combined"] = np.nan
    if keep:
        P = pcs.loc[ok, keep]
        full = P.notna().all(axis=1).values
        Pz = (P - P.median()) / (P.quantile(0.75) - P.quantile(0.25)) * 0.7413
        Pz = Pz.clip(-4, 4)
        Xc = np.hstack([Z[full], Pz.values[full]])
        # intensity axes weight 1; PCs share a total weight of 1 so they refine rather than dominate
        w = np.r_[np.ones(3), np.full(len(keep), 1 / np.sqrt(len(keep)))]
        init = np.r_[orient, np.zeros(len(keep))]
        s2, _, nodes2 = fit_score(Xc, init, w)
        if spearmanr(s2, Z[full, 0]).statistic < 0:
            s2 = 1 - s2
        idx = np.where(ok)[0][full]
        res.iloc[idx, res.columns.get_loc("cc_score_combined")] = s2
        print(f"combined score cells {full.sum()}; rho with intensity-only: {spearmanr(s2, s1[full]).statistic:.3f}")
        print("rho combined vs area:", round(spearmanr(s2, area[full]).statistic, 3))
        # independent check: PC-only trajectory
        s3, _, _ = fit_score(Pz.values[full], np.ones(len(keep)))
        rho_pc = spearmanr(s3, s1[full]).statistic
        if rho_pc < 0:
            s3 = 1 - s3
        res["cc_score_pc_only"] = np.nan
        res.iloc[idx, res.columns.get_loc("cc_score_pc_only")] = s3
        print("PC-only vs intensity-only |rho|:", round(abs(rho_pc), 3))
    res["cc_score"] = res["cc_score_combined"].fillna(res["cc_score_intensity"])
    res.to_csv(OUT / f"cell_cycle_scores{SUF}.csv")

    plot_all(res, ok, Z, logi, keep)
    summ = dict(
        cells_total=int(len(df)), cells_intensity_score=int(ok.sum()),
        cells_combined_score=int(res["cc_score_combined"].notna().sum()), kept_pcs=keep,
        rho_intensity_vs_softstage=float(spearmanr(s1, stage_idx).statistic),
        rho_intensity_vs_area=float(spearmanr(s1, area).statistic),
        well_eta2_intensity_score=float(eta2(s1, well)),
    )
    json.dump(summ, open(OUT / f"summary{SUF}.json", "w"), indent=2)
    print(summ)


def plot_all(res, ok, Z, logi, keep):
    r = res[ok].copy()
    sub = RNG.choice(len(r), min(25000, len(r)), replace=False)
    emb = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(Z[sub])
    rs = r.iloc[sub]
    fig, ax = plt.subplots(2, 4, figsize=(20, 9))
    cols = [("cc_score_intensity", "viridis"), ("z_phos", "magma"), ("z_p21", "magma"), ("z_pRb", "magma")]
    for a_, (c, cm) in zip(ax[0], cols):
        h = a_.scatter(emb[:, 0], emb[:, 1], c=rs[c], s=2, cmap=cm, rasterized=True)
        a_.set_title(c); plt.colorbar(h, ax=a_)
    h = ax[1, 0].scatter(emb[:, 0], emb[:, 1], c=np.log(rs["nuc_area_px"]), s=2, cmap="cividis", rasterized=True)
    ax[1, 0].set_title("log nuclear area"); plt.colorbar(h, ax=ax[1, 0])
    for w_, c_ in zip(sorted(rs["well"].unique()), ["C0", "C1", "C2"]):
        m = (rs["well"] == w_).values
        ax[1, 1].scatter(emb[m, 0], emb[m, 1], s=1, c=c_, label=w_, alpha=.5, rasterized=True)
    ax[1, 1].legend(markerscale=6); ax[1, 1].set_title("well")
    h = ax[1, 2].scatter(emb[:, 0], emb[:, 1], c=rs["soft_stage"], s=2, cmap="viridis", rasterized=True)
    ax[1, 2].set_title("soft-state expected stage (prior GMM)"); plt.colorbar(h, ax=ax[1, 2])
    c2 = "cc_score_combined" if r["cc_score_combined"].notna().any() else "cc_score_intensity"
    h = ax[1, 3].scatter(emb[:, 0], emb[:, 1], c=rs[c2], s=2, cmap="viridis", rasterized=True)
    ax[1, 3].set_title(c2); plt.colorbar(h, ax=ax[1, 3])
    for a_ in ax.ravel(): a_.set_xticks([]); a_.set_yticks([])
    plt.tight_layout(); plt.savefig(OUT / f"umap_overview{SUF}.png", dpi=110); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    L = logi[ok].iloc[sub]
    for a_, (x, y, nm) in zip(ax, [(L["Rb"], L["pRb"], "log Rb vs log pRb"), (rs["z_phos"], rs["z_p21"], "phos fraction vs p21 (z)"), (rs["z_pRb"], rs["z_p21"], "pRb vs p21 (z)")]):
        h = a_.scatter(x, y, c=rs["cc_score_intensity"], s=2, cmap="viridis", rasterized=True); a_.set_title(nm)
    plt.colorbar(h, ax=ax); plt.savefig(OUT / f"phase_planes{SUF}.png", dpi=110, bbox_inches="tight"); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))
    mid = np.linspace(0, 1, 26)[:-1] + 0.02
    for a_, c_ in zip(ax[:2], ["cc_score_intensity", c2]):
        bins = pd.cut(r[c_], np.linspace(0, 1, 26))
        for f, col in [("z_phos", "C3"), ("z_p21", "C0"), ("z_pRb", "C1")]:
            g = r.groupby(bins, observed=True)[f]
            a_.plot(mid[: len(g.median())], g.median().values, label=f, c=col)
            a_.fill_between(mid[: len(g.median())], g.quantile(.25).values, g.quantile(.75).values, color=col, alpha=.15)
        g = r.groupby(bins, observed=True)["nuc_area_px"].median()
        a2 = a_.twinx(); a2.plot(mid[: len(g)], g.values, "k--", label="nuc area"); a2.set_ylabel("nuclear area px")
        a_.set_title(f"marker profiles along {c_}"); a_.set_xlabel("score"); a_.legend()
    bins = pd.cut(r["cc_score_intensity"], np.linspace(0, 1, 26))
    for s in STATES:
        ax[2].plot(mid, r.groupby(bins, observed=True)[f"p_{s}"].mean().values, label=s)
    ax[2].set_title("soft-state probability along intensity score"); ax[2].legend()
    plt.tight_layout(); plt.savefig(OUT / f"marker_profiles{SUF}.png", dpi=110); plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for w_ in sorted(r["well"].unique()):
        ax[0].hist(r.loc[r["well"] == w_, "cc_score_intensity"], bins=40, alpha=.5, label=w_, density=True)
    ax[0].legend(); ax[0].set_title("intensity score by well")
    both = r.dropna(subset=["cc_score_combined"])
    if len(both):
        ax[1].hexbin(both["cc_score_intensity"], both["cc_score_combined"], gridsize=40, cmap="viridis", mincnt=1)
        ax[1].set_xlabel("intensity-only"); ax[1].set_ylabel("intensity + PCs"); ax[1].set_title("score comparison")
    plt.tight_layout(); plt.savefig(OUT / f"score_distributions{SUF}.png", dpi=110); plt.close()


if __name__ == "__main__":
    main()
