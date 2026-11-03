#!/usr/bin/env python3
"""Compute an exploratory per-cell G1-to-S progression proxy from pRb, Rb, and p21.

The output is a hypothesis-based marker score, not a validated cell-cycle phase
classifier. DNA content or trusted phase labels are required for validation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

CHANNELS = {
    "pRb": "4i_r3_prb",
    "Rb": "4i_r3_rb",
    "p21": "4i_r4_p21",
}


def feature_index(var: pd.DataFrame, organelle: str) -> int:
    selected = (
        var["source"].eq("organelle_profiler")
        & var["organelle"].eq(organelle)
        & var["category"].eq("intensity")
        & var["metric"].eq("intensity_mean")
        & var["aggregation"].eq("mean")
    )
    matches = np.flatnonzero(selected.to_numpy())
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one intensity-mean feature for {organelle}; found {len(matches)}")
    return int(matches[0])


def robust_z_by_group(values: pd.Series, groups: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    for _, indices in groups.groupby(groups, dropna=False).groups.items():
        group_values = values.loc[indices]
        center = group_values.median()
        iqr = group_values.quantile(0.75) - group_values.quantile(0.25)
        scale = iqr if np.isfinite(iqr) and iqr > 0 else 1.0
        result.loc[indices] = ((group_values - center) / scale).clip(-5, 5)
    return result


def calculate_score(values: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    complete = np.isfinite(values[["pRb", "Rb", "p21"]].to_numpy()).all(axis=1)
    result = pd.DataFrame(index=values.index)
    result["pRb_intensity"] = values["pRb"]
    result["Rb_intensity"] = values["Rb"]
    result["p21_intensity"] = values["p21"]
    result["pRb_to_Rb_ratio"] = values["pRb"] / values["Rb"].where(values["Rb"] > 0)
    result["log_pRb_to_Rb_proxy"] = np.log1p(values["pRb"]) - np.log1p(values["Rb"])
    result["complete_marker_measurements"] = complete
    result["pRb_Rb_evidence_z"] = (
        robust_z_by_group(np.log1p(values["pRb"]), groups)
        - robust_z_by_group(np.log1p(values["Rb"]), groups)
    )
    result["low_p21_evidence_z"] = -robust_z_by_group(np.log1p(values["p21"]), groups)
    result["marker_score_raw"] = 0.65 * result["pRb_Rb_evidence_z"] + 0.35 * result["low_p21_evidence_z"]
    result["cell_cycle_stage_score"] = np.nan
    result["stage_proxy"] = "not_scored_missing_marker"
    result["marker_evidence_discordant"] = False

    for _, indices in groups.groupby(groups, dropna=False).groups.items():
        valid_indices = [index for index in indices if complete[values.index.get_loc(index)]]
        if not valid_indices:
            continue
        score_values = result.loc[valid_indices, "marker_score_raw"]
        result.loc[valid_indices, "cell_cycle_stage_score"] = score_values.rank(method="average", pct=True)

    scored = result["cell_cycle_stage_score"].notna()
    discordant = (
        result["pRb_Rb_evidence_z"].mul(result["low_p21_evidence_z"]).lt(0)
        & result["pRb_Rb_evidence_z"].abs().ge(0.5)
        & result["low_p21_evidence_z"].abs().ge(0.5)
    )
    result.loc[scored, "marker_evidence_discordant"] = discordant.loc[scored]
    result.loc[scored & result["cell_cycle_stage_score"].lt(1 / 3), "stage_proxy"] = "G1_like_low_progression_proxy"
    result.loc[scored & result["cell_cycle_stage_score"].between(1 / 3, 2 / 3, inclusive="both"), "stage_proxy"] = "intermediate_uncertain"
    result.loc[scored & result["cell_cycle_stage_score"].gt(2 / 3), "stage_proxy"] = "G1_S_entry_like_high_progression_proxy"
    result.loc[scored & discordant, "stage_proxy"] = "ambiguous_discordant_markers"
    result.loc[~scored, ["cell_cycle_stage_score", "marker_score_raw", "pRb_Rb_evidence_z", "low_p21_evidence_z"]] = np.nan
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/altair/ops0077_4i/OPS0077_NTC_OrganelleProfiler_pRb_Rb_p21.h5ad"),
    )
    parser.add_argument("--outdir", type=Path, default=Path("reports/ops0077_cell_cycle_score"))
    parser.add_argument("--group-column", default="well_canonical", help="Column used for within-group robust scaling and percentile scores.")
    args = parser.parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(args.input)
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input, backed="r")
    obs = adata.obs.copy()
    var = adata.var.copy()
    if args.group_column not in obs:
        raise ValueError(f"Grouping column {args.group_column!r} is not present in obs")
    indices = [feature_index(var, organelle) for organelle in CHANNELS.values()]
    intensities = np.asarray(adata.X[:, indices], dtype=np.float64)
    values = pd.DataFrame(intensities, index=obs.index, columns=list(CHANNELS))
    adata.file.close()

    groups = obs[args.group_column].astype("string").fillna("missing_group").astype(str)
    scores = calculate_score(values, groups)
    output = obs[[column for column in ("well_canonical", "well", "tile_pheno", "sgRNA", "gene_name", "control_type") if column in obs]].copy()
    output.insert(0, "cell_id", obs.index.astype(str))
    output = output.join(scores)
    output.to_csv(args.outdir / "per_cell_marker_stage_scores.csv.gz", index=False)

    summary = {
        "input": str(args.input),
        "cells": int(len(obs)),
        "scored_cells": int(scores["cell_cycle_stage_score"].notna().sum()),
        "unscored_cells": int(scores["cell_cycle_stage_score"].isna().sum()),
        "group_column": args.group_column,
        "group_counts": groups.value_counts().to_dict(),
        "feature_names": {marker: str(adata_feature) for marker, adata_feature in zip(CHANNELS, [str(var.index[index]) for index in indices])},
        "formula": "within-group percentile of 0.65*(robust_z(log1p(pRb))-robust_z(log1p(Rb))) + 0.35*robust_z(-log1p(p21)); each robust z uses median/IQR and is clipped to +/-5",
        "stage_proxy_thresholds": {
            "below_1_3": "G1_like_low_progression_proxy",
            "1_3_to_2_3": "intermediate_uncertain",
            "above_2_3": "G1_S_entry_like_high_progression_proxy",
            "discordance_rule": "If both marker evidence terms have magnitude >=0.5 and opposite signs, label ambiguous_discordant_markers.",
        },
        "stage_counts": scores["stage_proxy"].value_counts().to_dict(),
        "median_raw_pRb_to_Rb_ratio": float(scores["pRb_to_Rb_ratio"].median()),
        "limitations": [
            "Heuristic G1-to-S progression proxy; not a validated G1/S/G2/M classifier.",
            "pRb-to-Rb intensity ratio is not a calibrated phosphorylation fraction because antibody, exposure, and channel scales may differ.",
            "Within-group scaling and percentiles make scores relative to the selected well/group, not an absolute cross-batch scale.",
            "Total Rb abundance is not assumed to change monotonically across the cycle; Rb enters the score through the pRb/Rb proxy.",
            "High score should not be interpreted as S, G2, or M phase without independent DNA-content or phase-marker validation.",
            "Confirm the pRb antibody epitope and acquisition linearity before biological interpretation.",
        ],
        "pca_role": "PCA is not calculated by this scoring script and must not contribute to the score; it is for visualization only.",
    }
    (args.outdir / "score_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
