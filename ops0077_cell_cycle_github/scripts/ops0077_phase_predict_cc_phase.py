"""Predict the cell-cycle phase score (cc_phase) from phase-contrast features only, by circular regression.

Label  : cc_phase from ops0077_cell_cycle_score.py (markers + DNA content). Phase contrast is NOT used to build it.
Target : (sin 2*pi*phase, cos 2*pi*phase); predicted phase = angle of the predicted pair, so 0 and 1 are the same point.
Inputs : OPS0077_NTC_phase_OrganelleProfiler.h5ad (single-cell phase-contrast features).

Features are restricted to those derived from the phase / focus images. Fluorescence-derived features (gfp, mcherry)
are always excluded and mask-geometry features (cell / nucleus shape, taken from the 4i segmentation) are excluded by
default (--include-mask-morphology to add them), so the predictor is label-free phase contrast.

Validation is out-of-fold, grouped so that neighbouring cells never straddle train and test:
  * primary  : GroupKFold over tiles (tile_pheno)
  * stringent: leave-one-well-out
Outputs go to reports/ops0077_phase_predicts_cc_phase/.

Usage: python ops0077_phase_predict_cc_phase.py [project_root] [--quick] [--include-mask-morphology] [--max-reg-dist PX]
"""
import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import QuantileTransformer

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--quick", action="store_true", help="subsample cells and skip the slow parts (testing)")
ap.add_argument("--include-mask-morphology", action="store_true")
ap.add_argument("--max-reg-dist", type=float, default=None, help="train/evaluate only on cells with phase-to-fixed distance <= this (px)")
args = ap.parse_args()

ROOT = Path(args.root)
PHASE = ROOT / "data/altair/OPS0077_NTC_phase_OrganelleProfiler.h5ad"
SCORES = ROOT / "reports/ops0077_cell_cycle_score/cell_cycle_scores.csv"
OUT = ROOT / "reports/ops0077_phase_predicts_cc_phase"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(0)
MAX_MISSING = 0.05  # keep features that are present in >= 95% of cells


# ---------------------------------------------------------------- feature groups
def feature_group(name):
    """Assign a feature to a group; None = never used."""
    if name.startswith(("op_contact_gfp", "op_contact_mcherry")):
        return None  # fluorescence-derived, not phase contrast
    if name.startswith("cp_single_object_Phase"):
        return "phase_texture_intensity"  # CellProfiler on the phase image
    if name.startswith(("op_phase2d", "op_nucleoli_phase2d", "op_network_phase2d", "op_contact_phase2d", "op_contact_nucleoli")):
        return "organelle_phase2d"
    if name.startswith(("op_focus3d", "op_nucleoli_focus3d", "op_network_focus3d", "op_contact_focus3d")):
        return "organelle_focus3d"
    if name.startswith(("cp_cell_", "cp_cytoplasm_", "cp_nucleus_", "op_cell_", "op_cp_", "op_nuclei_", "op_contact_nuclei")):
        return "mask_morphology"  # geometry of the 4i cell / nucleus masks (not phase-image derived)
    return None


# ---------------------------------------------------------------- helpers
def circ_diff(a, b):
    return (a - b + 0.5) % 1.0 - 0.5


def to_phase(sc):
    return (np.arctan2(sc[:, 0], sc[:, 1]) / (2 * np.pi)) % 1.0  # sin first, cos second


def metrics(true, pred, conf=None):
    err = circ_diff(pred, true)
    return dict(
        n=int(len(true)), mae_turns=float(np.abs(err).mean()), median_abs_err=float(np.median(np.abs(err))),
        within_0p05=float((np.abs(err) <= 0.05).mean()), within_0p10=float((np.abs(err) <= 0.10).mean()),
        circ_agreement=float(np.abs(np.exp(2j * np.pi * err).mean())),  # 1 = perfect, 0 = unrelated
        acc_4bins=float((np.floor(pred * 4) % 4 == np.floor(true * 4) % 4).mean()),
        acc_8bins=float((np.floor(pred * 8) % 8 == np.floor(true * 8) % 8).mean()),
    )  # a random guess has mae 0.25, within_0p10 = 0.20, acc_4bins = 0.25, acc_8bins = 0.125


def harmonic_r2(x, phase):
    m = ~np.isnan(x)
    t = 2 * np.pi * phase[m]
    D = np.column_stack([np.ones(m.sum()), np.sin(t), np.cos(t), np.sin(2 * t), np.cos(2 * t)])
    beta, *_ = np.linalg.lstsq(D, x[m], rcond=None)
    return 1 - np.var(x[m] - D @ beta) / np.var(x[m])


def ridge_model():
    return make_pipeline(SimpleImputer(strategy="median"), QuantileTransformer(n_quantiles=500, output_distribution="normal", subsample=20000, random_state=0),
                         RidgeCV(alphas=[1e1, 1e2, 1e3, 3e3, 1e4, 3e4, 1e5]))


class HgbCircular:
    """PCA(100) of the quantile-normalised features, then one gradient-boosting regressor for sin and one for cos."""

    def __init__(self, n_pca=100):
        self.n_pca = n_pca

    def fit(self, X, Y):
        self.pre = make_pipeline(SimpleImputer(strategy="median"),
                                 QuantileTransformer(n_quantiles=500, output_distribution="normal", subsample=20000, random_state=0),
                                 PCA(min(self.n_pca, X.shape[1]), random_state=0))
        P = self.pre.fit_transform(X)
        self.m = [HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, max_leaf_nodes=24, l2_regularization=1.0,
                                                random_state=0).fit(P, Y[:, k]) for k in range(2)]
        return self

    def predict(self, X):
        P = self.pre.transform(X)
        return np.column_stack([m.predict(P) for m in self.m])


def cv_predict(make, X, Y, groups, splitter):
    oof = np.full(Y.shape, np.nan)
    for tr, te in splitter.split(X, Y, groups):
        m = make().fit(X[tr], Y[tr])
        oof[te] = m.predict(X[te])
    return oof


# ---------------------------------------------------------------- data
def load():
    a = ad.read_h5ad(PHASE)
    sc = pd.read_csv(SCORES, index_col=0)
    sc = sc.drop_duplicates("cell_key", keep=False).set_index("cell_key")  # duplicated keys are ambiguous
    obs = a.obs.copy()
    obs.index = obs["cell_key"].values
    names = np.array(a.var_names)
    grp = np.array([feature_group(n) for n in names], dtype=object)
    X = np.asarray(a.X, dtype=np.float32)
    X[~np.isfinite(X)] = np.nan
    miss = np.isnan(X).mean(0)
    sd = np.nanstd(X, 0)
    use = (miss <= MAX_MISSING) & (sd > 0) & (grp != None)  # noqa: E711
    if not args.include_mask_morphology:
        use &= grp != "mask_morphology"
    feats = pd.DataFrame({"feature": names, "group": grp, "missing_frac": miss, "used": use})
    # context / size covariates for baselines (not phase-image features)
    obs["cells_in_tile"] = obs.groupby("tile_pheno")["cell_key"].transform("size")
    size_cols = [n for n in ("op_cell_area", "op_nuclei_area_mean") if n in names]
    size = pd.DataFrame(np.asarray(a.X[:, [list(names).index(c) for c in size_cols]], dtype=float), index=obs.index, columns=size_cols)
    obs = obs.join(sc[["cc_phase", "cc_phase_source", "cc_curve_dist_dna", "z_dna", "z_p21", "z_phos", "z_pRb"]])
    return obs, X, names, feats, size


def main():
    obs, X, names, feats, size = load()
    print(f"cells {len(obs)} | features used {int(feats['used'].sum())} of {len(feats)}")
    print(feats[feats["used"]].groupby("group").size().to_dict())
    feats.to_csv(OUT / "feature_table.csv", index=False)

    lab = obs["cc_phase"].notna().values
    if args.max_reg_dist is not None:
        lab &= (obs["pheno_to_fixed_distance_px"] <= args.max_reg_dist).values
    if args.quick:
        lab &= RNG.random(len(obs)) < 0.25
    print("labelled cells for training/validation:", lab.sum())
    cols = np.where(feats["used"].values)[0]
    XL = X[lab][:, cols]
    phase = obs.loc[lab, "cc_phase"].values
    Y = np.column_stack([np.sin(2 * np.pi * phase), np.cos(2 * np.pi * phase)])
    tile = obs.loc[lab, "tile_pheno"].values
    well = obs.loc[lab, "well"].values

    tile_cv = GroupKFold(n_splits=5)
    results, oof = {}, {}
    models = {"ridge": ridge_model, "hgb_pca": lambda: HgbCircular()}
    for name, make in models.items():
        for cvname, splitter, g in [("tile_cv", tile_cv, tile), ("leave_one_well_out", LeaveOneGroupOut(), well)]:
            if args.quick and cvname == "leave_one_well_out":
                continue
            o = cv_predict(make, XL, Y, g, splitter)
            oof[(name, cvname)] = o
            results[f"{name}|{cvname}"] = metrics(phase, to_phase(o))
            print(f"{name:8s} {cvname:20s} MAE {results[f'{name}|{cvname}']['mae_turns']:.3f} "
                  f"within0.1 {results[f'{name}|{cvname}']['within_0p10']:.3f} 8bin {results[f'{name}|{cvname}']['acc_8bins']:.3f}")

    # ---- baselines (confound checks), gradient boosting on a handful of non-phase-image covariates
    def small_hgb():
        return make_pipeline(SimpleImputer(strategy="median"), _MultiHgb())

    base_sets = {"baseline_size_only": size.loc[obs.index[lab]].values,
                 "baseline_context_only": obs.loc[lab, ["cells_in_tile", "pheno_to_fixed_distance_px"]].values}
    for bname, Xb in base_sets.items():
        o = cv_predict(small_hgb, Xb, Y, tile, tile_cv)
        oof[(bname, "tile_cv")] = o
        results[f"{bname}|tile_cv"] = metrics(phase, to_phase(o))
        print(f"{bname:24s} MAE {results[f'{bname}|tile_cv']['mae_turns']:.3f}")
    results["random_guess"] = dict(mae_turns=0.25, within_0p10=0.20, acc_4bins=0.25, acc_8bins=0.125)

    best = min(models, key=lambda k: results[f"{k}|tile_cv"]["mae_turns"])
    print("best model (tile CV):", best)

    # ---- final fit for cells without a label (and out-of-fold values for labelled ones)
    final = models[best]().fit(XL, Y)
    allpred = final.predict(X[:, cols])
    pred = pd.DataFrame(index=obs.index)
    for name in models:
        pred[f"pred_sin_{name}"] = np.nan
        pred[f"pred_cos_{name}"] = np.nan
    idx_lab = np.where(lab)[0]
    for name in models:
        o = oof[(name, "tile_cv")]
        pred.iloc[idx_lab, pred.columns.get_loc(f"pred_sin_{name}")] = o[:, 0]
        pred.iloc[idx_lab, pred.columns.get_loc(f"pred_cos_{name}")] = o[:, 1]
    pred["pred_sin"] = allpred[:, 0]
    pred["pred_cos"] = allpred[:, 1]
    pred.iloc[idx_lab, pred.columns.get_loc("pred_sin")] = oof[(best, "tile_cv")][:, 0]
    pred.iloc[idx_lab, pred.columns.get_loc("pred_cos")] = oof[(best, "tile_cv")][:, 1]
    pred["pred_phase"] = to_phase(pred[["pred_sin", "pred_cos"]].values)
    pred["pred_confidence"] = np.hypot(pred["pred_sin"], pred["pred_cos"])  # length of the predicted vector
    pred["pred_is_out_of_fold"] = lab
    out = obs[["well", "tile_pheno", "pheno_to_fixed_distance_px", "cells_in_tile", "cc_phase", "cc_phase_source", "cc_curve_dist_dna",
               "z_dna", "z_p21", "z_phos", "z_pRb"]].join(pred)
    out["error_turns"] = circ_diff(out["pred_phase"], out["cc_phase"])
    out.to_csv(OUT / "phase_predicted_cc_phase.csv")

    # ---- feature importance: univariate circular association + ridge coefficients on quantile-normalised features
    rid = ridge_model().fit(XL, Y)
    coef = rid[-1].coef_  # (2, n_features)
    fi = feats[feats["used"]].copy().reset_index(drop=True)
    fi["ridge_importance"] = np.hypot(coef[0], coef[1])
    fi["harmonic_r2_phase"] = [harmonic_r2(XL[:, j].astype(float), phase) for j in range(XL.shape[1])]
    fi.sort_values("harmonic_r2_phase", ascending=False).to_csv(OUT / "feature_importance.csv", index=False)

    # ---- ablation by feature group (ridge, 3-fold tile CV on a subsample)
    abl = []
    if not args.quick:
        gl = feats[feats["used"]]["group"].values
        sub = RNG.choice(len(phase), min(15000, len(phase)), replace=False)
        cv3 = GroupKFold(n_splits=3)
        for g in sorted(set(gl)):
            for mode in ("only", "without"):
                keepc = (gl == g) if mode == "only" else (gl != g)
                if keepc.sum() == 0:
                    continue
                o = cv_predict(ridge_model, XL[sub][:, keepc], Y[sub], tile[sub], cv3)
                m = metrics(phase[sub], to_phase(o))
                abl.append(dict(group=g, mode=mode, n_features=int(keepc.sum()), **m))
                print(f"ablation {g:26s} {mode:8s} MAE {m['mae_turns']:.3f}")
        pd.DataFrame(abl).to_csv(OUT / "ablation_by_group.csv", index=False)

    # ---- CellProfiler phase-texture features exist for only ~22% of cells (the better-registered ones).
    # Compare, on exactly those cells, the OrganelleProfiler features alone vs. plus the CellProfiler ones.
    gall = feats["group"].values
    cp_cols = np.where((gall == "phase_texture_intensity") & (np.nanstd(X, 0) > 0))[0]
    if len(cp_cols):
        ref = cp_cols[0]
        sel = lab & ~np.isnan(X[:, ref])
        if sel.sum() > 1000:
            Ys = np.column_stack([np.sin(2 * np.pi * obs.loc[sel, "cc_phase"].values), np.cos(2 * np.pi * obs.loc[sel, "cc_phase"].values)])
            ph_s, tile_s = obs.loc[sel, "cc_phase"].values, obs.loc[sel, "tile_pheno"].values
            cp_ok = cp_cols[np.isnan(X[sel][:, cp_cols]).mean(0) <= MAX_MISSING]
            for tag, cc in [("ridge_op_only", cols), ("ridge_op_plus_cp", np.r_[cols, cp_ok])]:
                o = cv_predict(ridge_model, X[sel][:, cc], Ys, tile_s, tile_cv)
                results[f"{tag}|cp_subset"] = metrics(ph_s, to_phase(o))
                print(f"cp subset ({sel.sum()} cells) {tag:18s} {len(cc)} features  MAE {results[f'{tag}|cp_subset']['mae_turns']:.3f}")

    # ---- QC breakdown tables
    d = out[out["pred_is_out_of_fold"]].copy()
    d["abs_err"] = d["error_turns"].abs()
    d["true_bin"] = pd.cut(d["cc_phase"], np.linspace(0, 1, 9), include_lowest=True)
    d.groupby("true_bin", observed=True)["abs_err"].agg(["mean", "median", "count"]).to_csv(OUT / "error_by_true_phase.csv")
    d["reg_bin"] = pd.qcut(d["pheno_to_fixed_distance_px"], 5, duplicates="drop")
    d.groupby("reg_bin", observed=True)["abs_err"].agg(["mean", "median", "count"]).to_csv(OUT / "error_by_registration_distance.csv")
    d.groupby("well")["abs_err"].agg(["mean", "median", "count"]).to_csv(OUT / "error_by_well.csv")

    summ = dict(best_model=best, n_features=int(feats["used"].sum()), n_labelled=int(lab.sum()), n_cells=int(len(obs)),
                feature_groups={k: int(v) for k, v in feats[feats["used"]].groupby("group").size().items()},
                results=results, options=dict(quick=args.quick, include_mask_morphology=args.include_mask_morphology,
                                              max_reg_dist=args.max_reg_dist))
    json.dump(summ, open(OUT / "summary.json", "w"), indent=2)
    plot_all(out, results, oof, phase, pd.DataFrame(abl), fi, best)
    print(json.dumps({k: {m: round(v, 3) for m, v in r.items() if m in ("mae_turns", "within_0p10", "acc_8bins")}
                      for k, r in results.items()}, indent=1))


class _MultiHgb:
    """Small gradient boosting for the baselines (sin and cos)."""

    def fit(self, X, Y):
        self.m = [HistGradientBoostingRegressor(max_iter=150, learning_rate=0.06, max_leaf_nodes=16, random_state=0).fit(X, Y[:, k]) for k in range(2)]
        return self

    def predict(self, X):
        return np.column_stack([m.predict(X) for m in self.m])

    def get_params(self, deep=True):
        return {}

    def set_params(self, **p):
        return self


# ---------------------------------------------------------------- plots
def plot_all(out, results, oof, phase, abl, fi, best):
    d = out[out["pred_is_out_of_fold"]].copy()
    d["abs_err"] = d["error_turns"].abs()
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False})

    # 1. predicted vs true (circular), error distribution, polar error
    fig = plt.figure(figsize=(18, 4.8))
    ax = fig.add_subplot(1, 3, 1)
    h = ax.hexbin(d["cc_phase"], d["pred_phase"], gridsize=50, cmap="viridis", mincnt=1, bins="log")
    for off in (-1, 0, 1):
        ax.plot([0, 1], [off, 1 + off], "w--", lw=.8)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("true cc_phase (markers + DNA)"); ax.set_ylabel("predicted from phase contrast")
    ax.set_title(f"predicted vs true ({best}, out-of-fold)"); plt.colorbar(h, ax=ax, label="log10 cells")
    ax = fig.add_subplot(1, 3, 2)
    ax.hist(d["error_turns"], bins=80, color="C0"); ax.axvline(0, c="k", lw=.8)
    ax.set_xlabel("circular error (turns; +-0.5 = opposite side of the cycle)"); ax.set_title(f"MAE {d['abs_err'].mean():.3f} turns (random 0.25)")
    axp = fig.add_subplot(1, 3, 3, projection="polar")
    cnt, edges = np.histogram(2 * np.pi * d["error_turns"], bins=48, range=(-np.pi, np.pi))
    axp.bar(edges[:-1], cnt, width=np.diff(edges), align="edge", color="C0", alpha=.85)
    axp.set_title("polar error (0 at top = no error)"); axp.set_theta_zero_location("N")
    plt.tight_layout(); plt.savefig(OUT / "qc_pred_vs_true.png", dpi=110); plt.close()

    # 2. error by true phase, confusion matrix (8 bins)
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    d["true_bin"] = pd.cut(d["cc_phase"], np.linspace(0, 1, 9), include_lowest=True)
    d.boxplot(column="abs_err", by="true_bin", ax=ax[0], showfliers=False, grid=False)
    ax[0].set_title("absolute error by true phase bin"); ax[0].set_xlabel("true phase bin (8)"); ax[0].set_ylabel("abs error (turns)"); ax[0].tick_params(axis="x", rotation=45)
    fig.suptitle("")
    tb = np.floor(d["cc_phase"] * 8).astype(int) % 8
    pb = np.floor(d["pred_phase"] * 8).astype(int) % 8
    cm = pd.crosstab(tb, pb, normalize="index").reindex(index=range(8), columns=range(8), fill_value=0)
    im = ax[1].imshow(cm.values, cmap="viridis", vmin=0, vmax=1)
    ax[1].set_xlabel("predicted bin"); ax[1].set_ylabel("true bin"); ax[1].set_title("row-normalised confusion (8 phase bins)"); plt.colorbar(im, ax=ax[1])
    for i in range(8):
        for j in range(8):
            ax[1].text(j, i, f"{cm.values[i, j]:.2f}", ha="center", va="center", color="w", fontsize=7)
    plt.tight_layout(); plt.savefig(OUT / "qc_error_by_phase.png", dpi=110); plt.close()

    # 3. confidence, registration, well, label quality
    fig, ax = plt.subplots(2, 3, figsize=(17, 8.5))
    q = pd.qcut(d["pred_confidence"], 10, duplicates="drop")
    g = d.groupby(q, observed=True)["abs_err"].mean()
    ax[0, 0].plot(range(len(g)), g.values, "o-"); ax[0, 0].axhline(.25, c="gray", ls="--"); ax[0, 0].set_xlabel("prediction confidence decile (length of predicted (sin, cos) vector)")
    ax[0, 0].set_ylabel("mean abs error"); ax[0, 0].set_title("is confidence informative?")
    order = np.argsort(-d["pred_confidence"].values)
    cum = np.cumsum(d["abs_err"].values[order]) / np.arange(1, len(order) + 1)
    ax[0, 1].plot(np.arange(1, len(order) + 1) / len(order), cum); ax[0, 1].axhline(.25, c="gray", ls="--")
    ax[0, 1].set_xlabel("fraction of cells kept (most confident first)"); ax[0, 1].set_ylabel("mean abs error"); ax[0, 1].set_title("error vs coverage")
    q = pd.qcut(d["pheno_to_fixed_distance_px"], 8, duplicates="drop")
    g = d.groupby(q, observed=True)["abs_err"].agg(["mean", "sem"])
    ax[0, 2].errorbar(range(len(g)), g["mean"], g["sem"] * 1.96, fmt="o-"); ax[0, 2].set_xticks(range(len(g)))
    ax[0, 2].set_xticklabels([f"{i.right:.0f}" for i in g.index], rotation=45)
    ax[0, 2].set_xlabel("phase-to-fixed distance (px, upper edge of octile)"); ax[0, 2].set_ylabel("mean abs error"); ax[0, 2].set_title("error vs registration distance")
    ax[1, 0].boxplot([d.loc[d["well"] == w, "abs_err"] for w in sorted(d["well"].unique())], showfliers=False)
    ax[1, 0].set_xticklabels(sorted(d["well"].unique())); ax[1, 0].set_ylabel("abs error"); ax[1, 0].set_title("error by well")
    dd = d.dropna(subset=["cc_curve_dist_dna"])
    q = pd.qcut(dd["cc_curve_dist_dna"], 8, duplicates="drop")
    g = dd.groupby(q, observed=True)["abs_err"].mean()
    ax[1, 1].plot(range(len(g)), g.values, "o-"); ax[1, 1].set_xlabel("distance of the label cell to the fitted cycle (octile)")
    ax[1, 1].set_ylabel("mean abs error"); ax[1, 1].set_title("error vs label quality")
    q = pd.qcut(d["cells_in_tile"], 6, duplicates="drop")
    g = d.groupby(q, observed=True)["abs_err"].mean()
    ax[1, 2].plot(range(len(g)), g.values, "o-"); ax[1, 2].set_xlabel("cells in the tile (sextile)"); ax[1, 2].set_ylabel("mean abs error"); ax[1, 2].set_title("error vs local cell number")
    plt.tight_layout(); plt.savefig(OUT / "qc_error_drivers.png", dpi=110); plt.close()

    # 4. model comparison incl. baselines and stringent CV
    keys = [k for k in results if "|" in k]
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.5))
    ax[0].bar(range(len(keys)), [results[k]["mae_turns"] for k in keys], color=["C3" if "baseline" in k else "C2" if "cp_subset" in k else "C0" if "tile" in k else "C1" for k in keys])
    ax[0].axhline(.25, c="gray", ls="--", label="random guess"); ax[0].set_xticks(range(len(keys))); ax[0].set_xticklabels(keys, rotation=60, ha="right", fontsize=8)
    ax[0].set_ylabel("MAE (turns)"); ax[0].set_title("lower is better"); ax[0].legend()
    ax[1].bar(range(len(keys)), [results[k]["acc_8bins"] for k in keys], color=["C3" if "baseline" in k else "C2" if "cp_subset" in k else "C0" if "tile" in k else "C1" for k in keys])
    ax[1].axhline(.125, c="gray", ls="--"); ax[1].set_xticks(range(len(keys))); ax[1].set_xticklabels(keys, rotation=60, ha="right", fontsize=8)
    ax[1].set_ylabel("accuracy, 8 phase bins"); ax[1].set_title("higher is better")
    plt.tight_layout(); plt.savefig(OUT / "qc_model_comparison.png", dpi=110); plt.close()

    # 5. do the predicted phases recover the marker / DNA profiles? (true vs predicted phase)
    fig, ax = plt.subplots(1, 2, figsize=(14, 4.2), sharey=True)
    for a, col, ttl in [(ax[0], "cc_phase", "along the TRUE phase"), (ax[1], "pred_phase", "along the PREDICTED phase")]:
        b = pd.cut(d[col], np.linspace(0, 1, 26)); mid = np.linspace(0, 1, 26)[:-1] + 0.02
        for f, c in [("z_dna", "k"), ("z_phos", "C3"), ("z_p21", "C0"), ("z_pRb", "C1")]:
            gm = d.groupby(b, observed=True)[f].median()
            a.plot(mid[:len(gm)], gm.values, c=c, label=f)
        a.set_title(ttl); a.set_xlabel("phase")
    ax[0].legend(); ax[0].set_ylabel("median robust z")
    plt.tight_layout(); plt.savefig(OUT / "qc_profiles_true_vs_predicted.png", dpi=110); plt.close()

    # 6. ablation and top features
    nrow = 2 if len(abl) else 1
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    if len(abl):
        p = abl.pivot(index="group", columns="mode", values="mae_turns")
        p.plot(kind="bar", ax=ax[0]); ax[0].axhline(.25, c="gray", ls="--"); ax[0].set_ylabel("MAE (turns)"); ax[0].set_title("feature groups: alone vs left out (ridge, 3-fold tile CV)")
        ax[0].tick_params(axis="x", rotation=30)
    else:
        ax[0].axis("off")
    top = fi.sort_values("harmonic_r2_phase", ascending=False).head(20)[::-1]
    ax[1].barh(top["feature"].str.slice(-60), top["harmonic_r2_phase"], color=top["group"].map({g: f"C{i}" for i, g in enumerate(sorted(fi["group"].unique()))}))
    ax[1].set_xlabel("circular R2 with the true phase (univariate)"); ax[1].set_title("top single features"); ax[1].tick_params(axis="y", labelsize=7)
    plt.tight_layout(); plt.savefig(OUT / "qc_feature_groups_and_top_features.png", dpi=110); plt.close()


if __name__ == "__main__":
    main()
