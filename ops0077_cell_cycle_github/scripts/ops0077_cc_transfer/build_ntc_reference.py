"""Freeze the NTC cell-cycle phase score into a reusable reference and train the phase-contrast predictor.

Reads the NTC data only. Writes (never overwrites the existing score / phase-prediction results):
  reports/ops0077_cc_transfer/reference/
    ntc_cc_reference.json          frozen per-well constants, closed-curve nodes, phase anchoring, sector definitions
    ntc_reference_scores.csv.gz       NTC cells re-scored with the frozen reference (axes, phase, sector)
    ntc_profiles_by_phase.csv      NTC marker / DNA profiles along the phase (for comparison with other cells)
    phase_ridge_model.joblib       ridge circular regression: phase-contrast features -> (sin, cos) of the NTC phase
    phase_model_features.json      names of the phase features the model expects
    reference_summary.json         checks (reproduction of the existing NTC scores) and model cross-validation

Usage: python build_ntc_reference.py [project_root] [--skip-phase-model]
"""
import argparse
import json
import sys
from datetime import date
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import QuantileTransformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cc_transfer_lib as L  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--skip-phase-model", action="store_true")
args = ap.parse_args()
ROOT = Path(args.root)
NTC_4I = ROOT / "data/altair/ops0077_4i/OPS0077_NTC_OrganelleProfiler_pRb_Rb_p21.h5ad"
NTC_DAPI = ROOT / "data/altair/OPS0077_NTC_DAPI_cpmeasure.h5ad"
NTC_PHASE = ROOT / "data/altair/OPS0077_NTC_phase_OrganelleProfiler.h5ad"
EXISTING = ROOT / "reports/ops0077_cell_cycle_score/cell_cycle_scores.csv.gz"
OUT = ROOT / "reports/ops0077_cc_transfer/reference"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(0)
MAX_MISSING = 0.05


# ---------------------------------------------------------------- NTC inputs
def load_ntc_4i():
    a = ad.read_h5ad(NTC_4I)
    v = list(a.var.index)
    get = lambda n: np.asarray(a.X[:, v.index(n)], dtype=float)
    df = L.nuclear_table(get, a.obs["well_canonical"].values, a.obs_names)
    df["cell_key"] = a.obs["cell_key"].values
    # duplicated cell_keys are ambiguous in the DAPI table: not scored (as in the existing NTC analysis)
    df.loc[df["cell_key"].duplicated(keep=False), [f"{m}_nuc" for m in L.MARKERS]] = np.nan
    return df


def attach_dna(df):
    o = ad.read_h5ad(NTC_DAPI).obs[["dna_rel_R234", "dna_reliable"]]
    o.index.name = "cell_key"
    t = df[["cell_key"]].join(o, on="cell_key")
    return t["dna_rel_R234"].where(t["dna_reliable"].fillna(False).astype(bool)).values


def phase_feature_group(name):
    """Phase-image-derived feature groups (same rule as ops0077_phase_predict_cc_phase.py, default settings)."""
    if name.startswith(("op_contact_gfp", "op_contact_mcherry")):
        return None  # fluorescence
    if name.startswith(("op_phase2d", "op_nucleoli_phase2d", "op_network_phase2d", "op_contact_phase2d", "op_contact_nucleoli")):
        return "organelle_phase2d"
    if name.startswith(("op_focus3d", "op_nucleoli_focus3d", "op_network_focus3d", "op_contact_focus3d")):
        return "organelle_focus3d"
    return None  # CellProfiler phase texture (only ~22% of cells) and mask geometry are not used


def ridge_model():
    return make_pipeline(SimpleImputer(strategy="median"),
                         QuantileTransformer(n_quantiles=500, output_distribution="normal", subsample=20000, random_state=0),
                         RidgeCV(alphas=[1e1, 1e2, 1e3, 3e3, 1e4, 3e4, 1e5]))


def main():
    summary = {"created": str(date.today())}
    df = load_ntc_4i()
    dna = attach_dna(df)

    # ---- 1. frozen per-well constants and axes on the NTC cells
    const = L.fit_frozen_constants(df, dna)
    z, _ = L.frozen_axes(df, dna, const)
    ok = z[L.AXES].notna().all(axis=1).values
    print(f"NTC cells {len(df)} | complete axes (markers + reliable DNA) {ok.sum()}")

    # ---- 2. closed curve (same fit, same RNG order as the existing NTC analysis)
    nodes, phase, dist, p0, flip = L.fit_curve_on_ntc(z.loc[ok, L.AXES].values, L.WEIGHTS, RNG)
    ref = {
        "version": 1, "created": str(date.today()), "axes": L.AXES, "weights": L.WEIGHTS.tolist(), "z_clip": L.Z_CLIP,
        "min_nuc_area_px": L.MIN_NUC_AREA_PX, "constants": const, "curve_nodes": nodes.tolist(),
        "anchor": {"p0": p0, "flip": flip}, "n_ntc_cells_fit": int(ok.sum()),
        "ntc_curve_dist_q95": float(np.percentile(dist, 95)), "ntc_curve_dist_q99": float(np.percentile(dist, 99)),
    }
    labels, frac = L.define_sectors(nodes, p0, flip, L.WEIGHTS)
    ref["sector_grid"] = {"n": L.GRID_N, "labels": labels, "dna_fraction_along_curve": frac}
    ref["sector_intervals"] = [[s, round(a, 4), round(b, 4)] for s, a, b in L.sector_intervals(ref)]
    print("sectors (phase intervals):", ref["sector_intervals"])

    # ---- 3. verify against the existing NTC result (cc_phase_dna)
    ex = pd.read_csv(EXISTING, index_col=0)["cc_phase_dna"].reindex(df.index)
    d = np.abs(L.circ_diff(phase, ex.values[ok]))
    have = ~np.isnan(d)
    summary["reproduction_vs_existing_cc_phase_dna"] = dict(
        n_compared=int(have.sum()), max_abs_circ_diff=float(d[have].max()), mean_abs_circ_diff=float(d[have].mean()),
        existing_cells_with_score=int(ex.notna().sum()), reproduced_cells_with_score=int(ok.sum()))
    print("reproduction vs existing cc_phase_dna:", summary["reproduction_vs_existing_cc_phase_dna"])

    # ---- 4. NTC reference scores, sectors and profiles
    sc = pd.DataFrame(index=df.index)
    sc["cell_key"], sc["well"] = df["cell_key"].values, df["well"].values
    sc = sc.join(z.add_prefix("z_"))
    sc["dna_rel"] = dna
    sc["nuc_area_px"] = df["nuc_area_px"].values
    sc["cc_phase_stain"], sc["cc_curve_dist"] = np.nan, np.nan
    sc.loc[ok, "cc_phase_stain"], sc.loc[ok, "cc_curve_dist"] = phase, dist
    sc["sector"] = L.sector_of(sc["cc_phase_stain"].values, ref)
    sc.to_csv(OUT / "ntc_reference_scores.csv.gz")
    b = pd.cut(sc["cc_phase_stain"], np.linspace(0, 1, 31), include_lowest=True)
    prof = sc.groupby(b, observed=True)[["z_dna", "z_phos", "z_p21", "z_pRb", "dna_rel"]].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
    prof.columns = [f"{a}_{ {'median': 'q50', '<lambda_0>': 'q25', '<lambda_1>': 'q75'}[c] }" for a, c in prof.columns]
    prof.insert(0, "phase_mid", np.linspace(0, 1, 31)[:-1] + 1 / 60)
    prof.to_csv(OUT / "ntc_profiles_by_phase.csv", index=False)
    ref["ntc_sector_fractions"] = sc.loc[ok, "sector"].value_counts(normalize=True).round(4).to_dict()
    L.save_reference(ref, OUT / "ntc_cc_reference.json")
    print("NTC sector fractions:", ref["ntc_sector_fractions"])

    # ---- 5. phase-contrast predictor trained on the NTC phase (label = cc_phase_dna, the quantity the stain projection reproduces)
    if not args.skip_phase_model:
        p = ad.read_h5ad(NTC_PHASE)
        names = np.array(p.var_names)
        X = np.asarray(p.X, dtype=np.float32)
        X[~np.isfinite(X)] = np.nan
        grp = np.array([phase_feature_group(n) for n in names], dtype=object)
        use = (np.isnan(X).mean(0) <= MAX_MISSING) & (np.nanstd(X, 0) > 0) & (grp != None)  # noqa: E711
        feats = names[use].tolist()
        obs = p.obs.copy()
        obs.index = obs["cell_key"].values
        lab = pd.Series(sc["cc_phase_stain"].values, index=sc["cell_key"].values)
        lab = lab[~lab.index.duplicated(keep=False)].reindex(obs.index)
        has = lab.notna().values
        Xl = X[has][:, use]
        phase_l = lab.values[has]
        Y = np.column_stack([np.sin(2 * np.pi * phase_l), np.cos(2 * np.pi * phase_l)])
        print(f"phase model: {len(feats)} features, {has.sum()} labelled NTC cells")

        # tile-grouped cross-validation (honest estimate for cells the model never saw)
        oof = np.zeros_like(Y)
        for tr, te in GroupKFold(n_splits=5).split(Xl, Y, obs.loc[has, "tile_pheno"].values):
            oof[te] = ridge_model().fit(Xl[tr], Y[tr]).predict(Xl[te])
        err = np.abs(L.circ_diff(L.to_phase(oof), phase_l))
        summary["phase_model_tile_cv"] = dict(n=int(has.sum()), mae_turns=float(err.mean()), within_0p10=float((err <= 0.10).mean()),
                                              random_guess_mae=0.25)
        print("phase model tile-grouped CV:", summary["phase_model_tile_cv"])
        model = ridge_model().fit(Xl, Y)
        joblib.dump(model, OUT / "phase_ridge_model.joblib")
        json.dump({"features": feats, "target": "(sin, cos) of 2*pi*cc_phase_stain",
                   "note": "trained on NTC cells only; apply to raw (unscaled) phase features"}, open(OUT / "phase_model_features.json", "w"))
    json.dump(summary, open(OUT / "reference_summary.json", "w"), indent=2)
    print("saved to", OUT)


if __name__ == "__main__":
    main()
