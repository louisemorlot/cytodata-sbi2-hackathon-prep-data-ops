"""Predict the cell-cycle phase score in NEW cells (e.g. gene knockouts) from the frozen NTC reference.

Nothing is estimated from the new cells: they are normalised with the frozen NTC per-well constants, projected onto the
frozen NTC cycle curve (stain score) and passed through the phase-contrast model trained on NTC cells (phase score).

Input h5ad (same layout as OPS0077_KO3_CCNB1_MCM3_CENPC_Multimodal_4i_Phase_DNA.h5ad):
  layers['raw'] (or X)       4i features, same columns as the NTC 4i file (raw, NOT pre-scaled)
  obsm['X_phase_raw']        phase features in the order of uns['phase_var_names'] (raw, NOT pre-scaled)
  obs['well']                'A/1/0'-style well matching the NTC wells
  obs['dna_content_log2'], obs['dna_reliable']   relative DNA content and its reliability flag
Existing derived columns of the input (cycle_stage, phase_predicted_*) are only carried as ref_* for comparison.

Outputs in --out (default reports/ops0077_cc_transfer/ko3_scores/):
  cc_phase_new_cells.csv.gz   per cell: cc_phase_stain, cc_phase_phase, cc_phase (final), QC flags, sectors, ref_* columns
  <input name>_cc_phase.h5ad   copy of the input with these columns added to obs (input untouched)
  summary.json

Usage: python apply_cc_transfer.py [project_root] [--input H5AD] [--reference DIR] [--out DIR] [--dna-col COL] [--dna-reliable-col COL]
"""
import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cc_transfer_lib as L  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("root", nargs="?", default=".")
ap.add_argument("--input", default="data/altair/OPS0077_KO3_CCNB1_MCM3_CENPC_Multimodal_4i_Phase_DNA.h5ad")
ap.add_argument("--reference", default="reports/ops0077_cc_transfer/reference")
ap.add_argument("--out", default="reports/ops0077_cc_transfer/ko3_scores")
ap.add_argument("--dna-col", default="dna_content_log2")
ap.add_argument("--dna-reliable-col", default="dna_reliable")
args = ap.parse_args()
ROOT = Path(args.root)
INPUT, REF_DIR, OUT = ROOT / args.input, ROOT / args.reference, ROOT / args.out
OUT.mkdir(parents=True, exist_ok=True)

# columns of the input that were derived by other pipelines: kept for comparison only, never used as inputs
REF_COLS = ["cycle_stage", "phase_predicted_cc_phase_ridge", "phase_predicted_cc_phase_hgb", "phase_predicted_dna",
            "phase_predicted_4N_prob", "phase_predicted_stage"]
CARRY = ["gene", "perturbation", "well", "tile_pheno", "sgRNA", "barcode", "nuclear_area_4i", "cell_area_4i"]


def main():
    ref = L.load_reference(REF_DIR / "ntc_cc_reference.json")
    a = ad.read_h5ad(INPUT)
    print(f"input: {INPUT.name} | {a.n_obs} cells")
    summary = {"input": str(INPUT.name), "n_cells": int(a.n_obs)}

    # ---- stain score: markers + DNA, frozen NTC normalisation and curve
    raw4i = a.layers["raw"] if "raw" in a.layers else a.X
    if "raw" not in a.layers:
        print("WARNING: no layers['raw']; using X (must be raw, unscaled 4i features)")
    v = list(a.var_names)
    get = lambda n: np.asarray(raw4i[:, v.index(n)], dtype=float)
    df = L.nuclear_table(get, a.obs["well"].astype(str).values, a.obs_names)
    reliable = a.obs[args.dna_reliable_col].fillna(False).astype(bool).values
    dna = np.where(reliable, a.obs[args.dna_col].astype(float).values, np.nan)
    z, unseen_well = L.frozen_axes(df, dna, ref["constants"])
    phase_stain, dist = L.project_to_reference(z, ref)
    if unseen_well.any():
        print(f"WARNING: {unseen_well.sum()} cells in wells not seen in the NTC reference (pooled NTC constants used)")

    out = pd.DataFrame(index=a.obs_names)
    for c in CARRY:
        if c in a.obs:
            out[c] = a.obs[c].values
    out = out.join(z.add_prefix("z_"))
    out["dna_content"] = dna
    out["nuc_area_px"] = df["nuc_area_px"].values
    out["has_markers"] = df[[f"{m}_nuc" for m in L.MARKERS]].notna().all(axis=1).values
    out["dna_reliable"] = reliable
    out["cc_phase_stain"] = phase_stain
    out["cc_curve_dist"] = dist
    out["out_of_ntc_support"] = dist > ref["ntc_curve_dist_q99"]  # farther from the NTC cycle than 99% of NTC cells
    out["unseen_well"] = unseen_well

    # ---- phase score: phase-contrast ridge model trained on NTC
    feats = json.load(open(REF_DIR / "phase_model_features.json"))["features"]
    pnames = list(a.uns["phase_var_names"])
    idx = [pnames.index(f) for f in feats]
    Xp = np.asarray(a.obsm["X_phase_raw"], dtype=np.float32)[:, idx]
    Xp[~np.isfinite(Xp)] = np.nan
    model = joblib.load(REF_DIR / "phase_ridge_model.joblib")
    sc = model.predict(Xp)
    out["cc_phase_phase"] = L.to_phase(sc)
    out["phase_confidence"] = np.hypot(sc[:, 0], sc[:, 1])
    out["phase_feature_nan_frac"] = np.isnan(Xp).mean(1)

    # ---- final score: stain where available, otherwise the phase-contrast prediction
    out["cc_phase"] = out["cc_phase_stain"].fillna(out["cc_phase_phase"])
    out["cc_phase_source"] = np.where(out["cc_phase_stain"].notna(), "stain", "phase")
    for col in ("cc_phase_stain", "cc_phase_phase", "cc_phase"):
        out[f"sector_{col.replace('cc_phase_', '').replace('cc_phase', 'final')}"] = L.sector_of(out[col].values, ref)
    for c in REF_COLS:
        if c in a.obs:
            out[f"ref_{c}"] = a.obs[c].values
    out.to_csv(OUT / "cc_phase_new_cells.csv.gz")

    # ---- copy of the input with the new columns
    add = out[[c for c in out.columns if c.startswith(("cc_", "sector_", "phase_conf")) or c in ("out_of_ntc_support",)]].copy()
    for c in add.columns:
        a.obs[f"tr_{c}"] = add[c].values
    for c in a.obs.columns:
        if a.obs[c].dtype == object:
            a.obs[c] = a.obs[c].fillna("none").astype(str)
    a.write_h5ad(OUT / f"{INPUT.stem}_cc_phase.h5ad", compression="gzip")

    summary.update(
        n_stain=int(out["cc_phase_stain"].notna().sum()), n_phase_only=int((out["cc_phase_source"] == "phase").sum()),
        n_out_of_ntc_support=int(out["out_of_ntc_support"].sum()), ntc_curve_dist_q99=ref["ntc_curve_dist_q99"],
        per_gene=out.groupby("gene").agg(n=("cc_phase", "size"), with_stain=("cc_phase_stain", lambda s: int(s.notna().sum())),
                                         out_of_support=("out_of_ntc_support", "sum")).astype(int).to_dict("index") if "gene" in out else {},
    )
    json.dump(summary, open(OUT / "summary.json", "w"), indent=2, default=str)
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
