"""Circular cell-cycle phase score for OPS0077 NTC cells.

Inputs per cell: nuclear p21 / Rb / pRb intensities (4i), relative DNA content (DAPI cpmeasure),
optionally per-channel / pairwise PCs (reports/ops0077_cell_cycle_clock/per_cell_pca_embeddings.csv.gz; skipped if absent).

Stages
  1. markers   : closed principal curve through [phos fraction, p21, pRb]          (QC only, weak without DNA)
  2. dna       : closed principal curve through [phos, p21, pRb, DNA content]       (core score)
  3. dna + PCs : stage 2 plus PCs that track the phase without being batch/size driven
cc_phase is in [0, 1): 0 = G1 start (lowest DNA content), increasing through S (p21 minimum) to G2/M,
and wraps back to G1. The cycle is circular: 0 and 1 are the same point.

Usage: python ops0077_cell_cycle_score.py [project_root]
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

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
H5AD = ROOT / "data/altair/ops0077_4i/OPS0077_NTC_OrganelleProfiler_pRb_Rb_p21.h5ad"
DAPI = ROOT / "data/altair/OPS0077_NTC_DAPI_cpmeasure.h5ad"
PCS = ROOT / "reports/ops0077_cell_cycle_clock/per_cell_pca_embeddings.csv.gz"
OUT = ROOT / "reports/ops0077_cell_cycle_score"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(0)
MARKERS = ["pRb", "Rb", "p21"]
AXES = ["phos", "p21", "pRb", "dna"]  # column order of the core feature matrix
DNA_COL, P21_COL = 3, 1


# ---------------------------------------------------------------- features
def load_intensities():
    a = ad.read_h5ad(H5AD)
    v = list(a.var.index)
    g = lambda n: np.asarray(a.X[:, v.index(n)], dtype=float)
    df = pd.DataFrame(index=a.obs_names)
    df["cell_key"] = a.obs["cell_key"].values
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
    # 31 cell_keys are duplicated in the 4i file: ambiguous match to the DAPI table, drop them
    df.loc[df["cell_key"].duplicated(keep=False), [f"{m}_nuc" for m in MARKERS]] = np.nan
    return df


def attach_dapi(df):
    d = ad.read_h5ad(DAPI)
    o = d.obs
    v = list(d.var.index)
    area = np.nanmean(np.log(np.column_stack([np.asarray(d.X[:, v.index(f"nucleus_R{k}_Area")], float) for k in (2, 3, 4)])), 1)
    keep = ["dna_rel_R234", "dna_reliable", "dense_field", "nucleus_found_R234", "dna_R1_any_saturated", "pRb_mean", "x_cp", "y_cp"]
    t = o[keep].copy()
    t["dapi_log_area"] = area
    t.index.name = "cell_key"
    df = df.join(t, on="cell_key")
    # DNA content is only trusted where the DAPI table flags it reliable
    df["dna"] = df["dna_rel_R234"].where(df["dna_reliable"].fillna(False).astype(bool))
    return df


def preprocess(df):
    """Per-well background subtraction, log, robust z. Returns log-intensities and z axes."""
    logi = pd.DataFrame(index=df.index)
    for m in MARKERS:
        col = df[f"{m}_nuc"]
        # background = per-well 0.5th percentile of the nuclear mean; floor at the per-well 1st percentile of the
        # positive remainder so dim cells cannot pile up at log(0)
        bg = col.groupby(df["well"]).transform(lambda x: np.nanpercentile(x, 0.5))
        sig = col - bg
        floor = sig.where(sig > 0).groupby(df["well"]).transform(lambda x: np.nanpercentile(x, 1))
        logi[m] = np.log(sig.clip(lower=floor))
    z = pd.DataFrame(index=df.index)
    z["phos"] = logi["pRb"] - logi["Rb"]  # phosphorylated fraction of Rb
    z["p21"] = logi["p21"]
    z["pRb"] = logi["pRb"]
    z["dna"] = df["dna"]
    for c in z:
        med = z[c].groupby(df["well"]).transform("median")
        mad = (z[c] - med).abs().groupby(df["well"]).transform("median") * 1.4826
        z[c] = ((z[c] - med) / mad).clip(-4, 4)
    return logi, z


# ---------------------------------------------------------------- closed principal curve
def project_closed(Xw, nodes):
    """Project onto the closed polyline through nodes; returns phase in [0,1) and distance."""
    P = np.vstack([nodes, nodes[:1]])
    seg = np.diff(P, axis=0)
    seglen = np.linalg.norm(seg, axis=1)
    cum = np.r_[0, np.cumsum(seglen)]
    best_d = np.full(len(Xw), np.inf)
    best_t = np.zeros(len(Xw))
    for i in range(len(seg)):
        s = seg[i]
        u = np.clip(((Xw - P[i]) @ s) / (s @ s + 1e-12), 0, 1)
        d = np.linalg.norm(Xw - (P[i] + u[:, None] * s), axis=1)
        m = d < best_d
        best_d[m], best_t[m] = d[m], cum[i] + u[m] * seglen[i]
    return (best_t / cum[-1]) % 1.0, best_d


def node_positions(nodes):
    P = np.vstack([nodes, nodes[:1]])
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    return cum[:-1] / cum[-1]


def fit_closed_curve(Xw, init_phase, n_nodes=16, n_iter=30, smooth=3):
    phase = init_phase.copy()
    for _ in range(n_iter):
        rank = np.argsort(np.argsort(phase, kind="stable"), kind="stable")
        idx = np.minimum(rank * n_nodes // len(phase), n_nodes - 1)  # equal-count bins along the phase
        nodes = np.array([Xw[idx == k].mean(0) for k in range(n_nodes)])
        for _ in range(smooth):  # periodic smoothing
            nodes = 0.25 * np.roll(nodes, 1, 0) + 0.5 * nodes + 0.25 * np.roll(nodes, -1, 0)
        phase, dist = project_closed(Xw, nodes)
    return nodes, phase, dist


def anchor(nodes, w, phase):
    """Rotate/flip phase: 0 = node with lowest DNA content; forward direction passes the p21 minimum (S) before the DNA maximum."""
    pos = node_positions(nodes)
    dna_n, p21_n = nodes[:, DNA_COL] / w[DNA_COL], nodes[:, P21_COL] / w[P21_COL]
    p0, pmax, pp21 = pos[np.argmin(dna_n)], pos[np.argmax(dna_n)], pos[np.argmin(p21_n)]
    rel = lambda x: (x - p0) % 1.0
    flip = not (rel(pp21) < rel(pmax))  # p21 minimum must lie on the arc between G1 start and DNA maximum
    out = rel(phase)
    if flip:
        out = (-out) % 1.0
    return out, flip


def init_phase_from_pca(Xw):
    Xc = Xw - Xw.mean(0)
    u, s, vt = np.linalg.svd(Xc[RNG.choice(len(Xc), min(len(Xc), 20000), replace=False)], full_matrices=False)
    p = Xc @ vt[:2].T
    return (np.arctan2(p[:, 1], p[:, 0]) / (2 * np.pi)) % 1.0


def fit_score(X, w, init=None, fit_n=20000):
    """Fit on a subsample, project all cells. Returns anchored phase, distance to curve, nodes."""
    Xw = X * w
    init = init_phase_from_pca(Xw) if init is None else init
    sub = RNG.choice(len(X), min(fit_n, len(X)), replace=False)
    nodes, _, _ = fit_closed_curve(Xw[sub], init[sub])
    phase, dist = project_closed(Xw, nodes)
    if len(w) > DNA_COL:  # without a DNA axis the phase origin is arbitrary; aligned afterwards
        phase, _ = anchor(nodes, w, phase)
    return phase, dist, nodes


def circ_diff(a, b):
    return (a - b + 0.5) % 1.0 - 0.5


def circ_corr(a, b):
    """Resultant length of the phase difference (1 = identical, 0 = unrelated); rotation-sensitive."""
    return float(np.abs(np.mean(np.exp(2j * np.pi * circ_diff(a, b)))))


def eta2(x, groups):
    ok = ~np.isnan(x)
    g = [x[ok & (groups == k)] for k in np.unique(groups)]
    h = kruskal(*g).statistic
    return (h - len(g) + 1) / (ok.sum() - len(g))


def harmonic_r2(x, phase):
    """Variance in x explained by the first two harmonics of the phase (circular association)."""
    m = ~np.isnan(x)
    t = 2 * np.pi * phase[m]
    D = np.column_stack([np.ones(m.sum()), np.sin(t), np.cos(t), np.sin(2 * t), np.cos(2 * t)])
    beta, *_ = np.linalg.lstsq(D, x[m], rcond=None)
    return 1 - np.var(x[m] - D @ beta) / np.var(x[m])


def main():
    df = attach_dapi(load_intensities())
    logi, z = preprocess(df)
    z_mk = z[["phos", "p21", "pRb"]]
    ok_mk = z_mk.notna().all(axis=1)
    ok_dna = ok_mk & z["dna"].notna()
    print(f"cells {len(df)} | complete markers {ok_mk.sum()} | + reliable DNA {ok_dna.sum()}")

    res = pd.DataFrame(index=z.index)
    res = res.join(z.add_prefix("z_")).join(
        df[["cell_key", "well", "nuc_area_px", "dapi_log_area", "dna_rel_R234", "dna_reliable", "dense_field",
            "nucleus_found_R234", "dna_R1_any_saturated", "x_cp", "y_cp"]])

    # ---- stage 2 (core): markers + DNA content
    Zd = z.loc[ok_dna, AXES].values
    w_d = np.array([1.0, 1.0, 1.0, 1.5])
    ph_d, dist_d, nodes_d = fit_score(Zd, w_d)
    res.loc[ok_dna, "cc_phase_dna"] = ph_d
    res.loc[ok_dna, "cc_curve_dist_dna"] = dist_d
    well_d = df.loc[ok_dna, "well"].values
    print("node means (phos,p21,pRb,dna) around the loop, from G1 start:")
    pos = node_positions(nodes_d)
    order = np.argsort((pos - pos[np.argmin(nodes_d[:, DNA_COL])]) % 1.0)
    print(np.round((nodes_d / w_d)[order][::2], 2))
    print("phase vs DNA: harmonic R2", round(harmonic_r2(Zd[:, 3], ph_d), 3), "| p21", round(harmonic_r2(Zd[:, 1], ph_d), 3),
          "| phos", round(harmonic_r2(Zd[:, 0], ph_d), 3))
    print("phase vs area: harmonic R2", round(harmonic_r2(df.loc[ok_dna, "dapi_log_area"].values, ph_d), 3),
          "| well eta2", round(eta2(ph_d, well_d), 4))

    # ---- stage 1: markers only (QC / fallback), aligned to the DNA-anchored phase
    Zm = z_mk[ok_mk].values
    ph_m, dist_m, _ = fit_score(Zm, np.ones(3))
    both = ok_dna[ok_mk].values
    ref = res.loc[ok_mk, "cc_phase_dna"].values
    # align the (rotation-free) markers-only loop to the DNA-anchored phase using the shared cells
    best = (-1, 0, 1)
    for f in (1, -1):
        for sh in np.linspace(0, 1, 200, endpoint=False):
            r = circ_corr((f * ph_m[both] + sh) % 1.0, ref[both])
            if r > best[0]:
                best = (r, sh, f)
    ph_m_al = (best[2] * ph_m + best[1]) % 1.0
    res.loc[ok_mk, "cc_phase_markers"] = ph_m_al
    res.loc[ok_mk, "cc_curve_dist_markers"] = dist_m
    print("markers-only vs DNA-anchored phase: circular agreement", round(best[0], 3),
          "| mean abs phase diff", round(float(np.mean(np.abs(circ_diff(ph_m_al[both], ref[both])))), 3))

    # ---- stage 3: add PCs
    keep = []
    if PCS.exists():
        pcs = pd.read_csv(PCS, index_col=0).reindex(df.index)
        rows = []
        for c in pcs.columns:
            x = pcs.loc[ok_dna, c].values
            m = ~np.isnan(x)
            if m.sum() < 5000:
                rows.append(dict(pc=c, n=int(m.sum())))
                continue
            rows.append(dict(pc=c, n=int(m.sum()), coverage=m.mean(), harmonic_r2_phase=harmonic_r2(x, ph_d),
                             rho_dna=spearmanr(x[m], Zd[m, 3]).statistic,
                             rho_area=spearmanr(x[m], df.loc[ok_dna, "dapi_log_area"].values[m], nan_policy="omit").statistic,
                             eta2_well=eta2(x, well_d)))
        scr = pd.DataFrame(rows)
        scr["keep"] = ((scr["coverage"] > 0.5) & (scr["harmonic_r2_phase"] >= 0.03) & (scr["eta2_well"] < 0.10)
                       & (scr["rho_area"].abs() < 0.6))
        scr.to_csv(OUT / "pc_screening.csv", index=False)
        print(scr.round(3).to_string())
        keep = scr.loc[scr["keep"], "pc"].tolist()
        print("kept PCs:", keep)
    else:
        print(f"PC file not found ({PCS}); skipping the PC layer")
    res["cc_phase_dna_pcs"] = np.nan
    if keep:
        P = pcs.loc[ok_dna, keep]
        full = P.notna().all(axis=1).values
        Pz = ((P - P.median()) / (P.quantile(0.75) - P.quantile(0.25)) * 0.7413).clip(-4, 4)
        Xc = np.hstack([Zd[full], Pz.values[full]])
        w_c = np.r_[w_d, np.full(len(keep), 1 / np.sqrt(len(keep)))]
        ph_c, dist_c, _ = fit_score(Xc, w_c, init=ph_d[full])
        idx = np.where(ok_dna)[0][full]
        res.iloc[idx, res.columns.get_loc("cc_phase_dna_pcs")] = ph_c
        print(f"dna+PCs cells {full.sum()} | agreement with dna-only: {circ_corr(ph_c, ph_d[full]):.3f} "
              f"| mean abs diff {np.mean(np.abs(circ_diff(ph_c, ph_d[full]))):.3f}")

    # ---- final score: best available
    res["cc_phase"] = res["cc_phase_dna_pcs"]
    res["cc_phase_source"] = np.where(res["cc_phase"].notna(), "dna+pcs", "")
    # markers-only phase is kept as a separate QC column: it agrees poorly with the DNA-anchored phase, so it is not a fallback
    for col, name in [("cc_phase_dna", "dna")]:
        miss = res["cc_phase"].isna() & res[col].notna()
        res.loc[miss, "cc_phase"] = res.loc[miss, col]
        res.loc[miss, "cc_phase_source"] = name
    res["cc_phase_angle_rad"] = 2 * np.pi * res["cc_phase"]
    res.to_csv(OUT / "cell_cycle_scores.csv.gz")
    print(res["cc_phase_source"].value_counts().to_dict())

    plot_all(res, z, logi, ok_dna)
    summ = dict(
        cells_total=int(len(df)), cells_markers=int(ok_mk.sum()), cells_dna=int(ok_dna.sum()),
        cells_dna_pcs=int(res["cc_phase_dna_pcs"].notna().sum()), kept_pcs=keep,
        source_counts=res["cc_phase_source"].value_counts().to_dict(),
        markers_vs_dna_circular_agreement=float(best[0]),
        well_eta2_phase=float(eta2(ph_d, well_d)),
    )
    json.dump(summ, open(OUT / "summary.json", "w"), indent=2)
    print(summ)


def plot_all(res, z, logi, ok_dna):
    r = res[ok_dna].copy()
    sub = RNG.choice(len(r), min(25000, len(r)), replace=False)
    Z = z.loc[ok_dna, AXES].values
    emb = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=0).fit_transform(Z[sub] * np.array([1, 1, 1, 1.5]))
    rs = r.iloc[sub]
    fig, ax = plt.subplots(2, 4, figsize=(20, 9))
    spec = [("cc_phase_dna", "twilight"), ("z_dna", "magma"), ("z_phos", "magma"), ("z_p21", "magma"),
            ("z_pRb", "magma"), ("dapi_log_area", "cividis"), ("cc_phase_markers", "twilight"), ("cc_phase", "twilight")]
    for a_, (c, cm) in zip(ax.ravel(), spec):
        h = a_.scatter(emb[:, 0], emb[:, 1], c=rs[c], s=2, cmap=cm, rasterized=True)
        a_.set_title(c); a_.set_xticks([]); a_.set_yticks([]); plt.colorbar(h, ax=a_)
    plt.tight_layout(); plt.savefig(OUT / "umap_overview.png", dpi=110); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(17, 5))
    for a_, (x, y) in zip(ax, [("z_dna", "z_phos"), ("z_dna", "z_p21"), ("z_pRb", "z_p21")]):
        h = a_.scatter(rs[x], rs[y], c=rs["cc_phase_dna"], s=2, cmap="twilight", rasterized=True)
        a_.set_xlabel(x); a_.set_ylabel(y)
    plt.colorbar(h, ax=ax, label="cc_phase_dna"); plt.savefig(OUT / "phase_planes.png", dpi=110, bbox_inches="tight"); plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(15, 4.5))
    bins = pd.cut(r["cc_phase_dna"], np.linspace(0, 1, 31))
    mid = np.linspace(0, 1, 31)[:-1] + 1 / 60
    for f, col in [("z_phos", "C3"), ("z_p21", "C0"), ("z_pRb", "C1"), ("z_dna", "k")]:
        g = r.groupby(bins, observed=True)[f]
        ax[0].plot(mid, g.median().values, c=col, label=f)
        ax[0].fill_between(mid, g.quantile(.25).values, g.quantile(.75).values, color=col, alpha=.12)
    g = r.groupby(bins, observed=True)["dapi_log_area"].median()
    a2 = ax[0].twinx(); a2.plot(mid, g.values, "--", c="gray"); a2.set_ylabel("log nuclear area (dashed)")
    ax[0].set_xlabel("cc_phase_dna (0 = G1 start, wraps)"); ax[0].legend(); ax[0].set_title("marker profiles around the cycle")
    ax[1].hist(r["cc_phase_dna"], bins=60, color="C4"); ax[1].set_title("cell density around the cycle")
    plt.tight_layout(); plt.savefig(OUT / "marker_profiles.png", dpi=110); plt.close()

    fig, ax = plt.subplots(1, 3, figsize=(17, 4.5))
    for w_ in sorted(r["well"].unique()):
        ax[0].hist(r.loc[r["well"] == w_, "cc_phase_dna"], bins=40, alpha=.5, label=w_, density=True)
    ax[0].legend(); ax[0].set_title("phase by well")
    b = r.dropna(subset=["cc_phase_markers"])
    ax[1].hexbin(b["cc_phase_dna"], b["cc_phase_markers"], gridsize=40, cmap="viridis", mincnt=1)
    ax[1].set_xlabel("dna-anchored"); ax[1].set_ylabel("markers only (aligned)"); ax[1].set_title("markers-only vs DNA-anchored")
    b = r.dropna(subset=["cc_phase_dna_pcs"])
    if len(b):
        ax[2].hexbin(b["cc_phase_dna"], b["cc_phase_dna_pcs"], gridsize=40, cmap="viridis", mincnt=1)
        ax[2].set_xlabel("dna only"); ax[2].set_ylabel("dna + PCs"); ax[2].set_title("effect of PC layer")
    else:
        ax[2].set_title("no PC layer"); ax[2].axis("off")
    plt.tight_layout(); plt.savefig(OUT / "score_distributions.png", dpi=110); plt.close()


if __name__ == "__main__":
    main()
