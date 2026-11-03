# How the predictions are made

Two kinds of prediction are built on the NTC cell cycle phase (`cc_phase`, see [SCORE_CALCULATION.md](SCORE_CALCULATION.md)):

1. **Phase contrast → phase (NTC cells).** A model predicts the stain-derived phase from label-free phase-contrast features only. Phase contrast is never used to build the score it predicts.
2. **NTC score → new cells (gene knockouts).** The NTC score is frozen into a reference and applied to cells that were not part of it, **without using any signal from the new cells to define it**.

Contents: [Circular regression](#circular-regression-used-in-both) · [Part 3: phase contrast in NTC](#part-3-phase-contrast--cc_phase-in-ntc-cells) · [Part 4: transfer to new cells](#part-4-transfer-to-new-cells-the-gene-knockouts) · [Output columns](#output-columns)

## Circular regression (used in both)

The phase lives on a circle, so a model must not predict the raw 0-to-1 number (cells at 0.99 and 0.01 are neighbours). Instead the model predicts the pair

> `(sin 2π·phase, cos 2π·phase)`

and the **predicted phase is the angle of the predicted pair**. The length of the predicted pair, `sqrt(sin² + cos²)`, is a **confidence**: a model that cannot decide between two phases predicts a short vector. Errors are measured on the circle (the signed circular difference `((pred − true + 0.5) mod 1) − 0.5`, in turns; ±0.5 is the opposite side of the cycle).

## Part 3: phase contrast → `cc_phase` in NTC cells

Script: `scripts/ops0077_phase_predict_cc_phase.py`. Runner: `notebooks/16_run_ops0077_phase_predict_cc_phase.ipynb`. Explore: `notebooks/17_ops0077_phase_prediction_exploration.ipynb`. Results: `reports/ops0077_phase_predicts_cc_phase/`.

**Label.** `cc_phase` of the 35,633 scored NTC cells.

**Features.** Only features derived from the phase images, from the OrganelleProfiler phase file:

| group | used | number |
|---|---|---|
| structures in the 2-D phase image (`op_phase2d_*`, nucleoli, network, contact) | yes | 765 |
| structures in the focus-stack image (`op_focus3d_*`, nucleoli, network, contact) | yes | 759 |
| fluorescence-derived (`op_contact_gfp*`, `op_contact_mcherry*`) | **never** | |
| geometry of the cell / nucleus masks (taken from the 4i segmentation, including DAPI) | **off by default** (`--include-mask-morphology`) | |
| CellProfiler phase texture (`cp_single_object_Phase_*`) | only in a dedicated comparison | 513 |

Only features present in at least 95% of cells are used (1,524 in total). The CellProfiler phase texture exists for only 22% of the cells (the better-registered ones), so it is compared separately on exactly those cells.

**Models** (both predict `(sin, cos)`):
- **Ridge**: median imputation → quantile transform to a normal distribution → ridge regression (penalty chosen by efficient leave-one-out cross-validation inside each training fold). This is the main model.
- **Gradient boosting**: the same transformation, 100 principal components, then one gradient-boosting regressor for sin and one for cos.

**Validation (always out-of-fold).** Neighbouring cells share a field and are correlated, so cells are never split randomly:
- primary: 5-fold cross-validation grouped by **tile** (`tile_pheno`);
- stringent: **leave one well out**.

**Baselines** that show the model learns more than a confound:
- *size only*: nuclear / cell area (gradient boosting);
- *context only*: registration distance and number of cells per tile (gradient boosting);
- a random guess (mean circular error 0.25 turns).

Predictions for labelled cells in `phase_predicted_cc_phase.csv` are the out-of-fold values; cells without a label get the value of a model trained on all labelled cells (`pred_is_out_of_fold` says which).

## Part 4: transfer to new cells (the gene knockouts)

Scripts: `scripts/ops0077_cc_transfer/` (`cc_transfer_lib.py`, `build_ntc_reference.py`, `apply_cc_transfer.py`). Notebooks: `18_...` (build the reference), `19_...` (score new cells). Results: `reports/ops0077_cc_transfer/reference/` and `ko3_scores/`.

### The principle
A score defined on the NTC cells must be applied to new cells **using NTC statistics only**. In particular the new cells are **not** re-normalised on themselves: if a knockout population were re-centred per well on its own medians, an arrest would be partly removed by the normalisation and the shift being measured would disappear. Everything that was estimated from data is therefore frozen from the NTC cells.

### Step 1 – freeze the NTC reference (`build_ntc_reference.py`)
Saved to `ntc_cc_reference.json` (plain JSON, safe to share):
- the per-well, per-marker background and floor, and the per-well median and MAD of each axis (normalisation constants of section 1 of [SCORE_CALCULATION.md](SCORE_CALCULATION.md));
- the nodes of the closed curve, the axis weights, and the anchoring (phase 0 and direction);
- the **phase sectors**, defined from the NTC curve only: along the curve, the DNA content is scaled to [0, 1] between its minimum and maximum. G1 is below 0.25, S is rising between 0.25 and 0.75, G2/M is at or above 0.75, and the descending part between 0.25 and 0.75 after the peak is the "M-to-G1 exit":

  | sector | phase interval |
  |---|---|
  | G1 | 0 to 0.253 and 0.881 to 1 |
  | S | 0.253 to 0.447 |
  | G2/M | 0.447 to 0.681 |
  | M-to-G1 exit | 0.681 to 0.881 |

  In NTC cells: G1 63%, S 7%, G2/M 25%, exit 5%;
- the NTC 95th and 99th percentile of the distance to the curve (3.52 at the 99th), used to flag **cells off the NTC cycle**.

The refit uses the same procedure and random seed as the NTC score and **reproduces the existing NTC phase for all 35,633 cells to machine precision** (maximum difference 1e-16 turns); this check is run every time the reference is built (`reference_summary.json`).

### Step 2 – two predictors for a new cell (`apply_cc_transfer.py`)

| score | uses | how |
|---|---|---|
| `cc_phase_stain` | p21, Rb, pRb nuclear intensities and DNA content | normalise with the **frozen NTC per-well constants**, then project onto the **frozen NTC curve** (no refit) |
| `cc_phase_phase` | phase-contrast features only | the NTC-trained ridge circular regression of Part 3, retrained on all labelled NTC cells with the label `cc_phase_dna` (the quantity the stain projection reproduces); tile-grouped cross-validation on NTC: mean error 0.103 turns |
| `cc_phase` | | `cc_phase_stain` where available, else `cc_phase_phase` (`cc_phase_source`) |

The PC layer is **not** used for transfer (the PCA models of the NTC PCs were not saved, and the PCs changed the phase by only 0.028 turns on average).

Each cell also gets its NTC-defined sector, `out_of_ntc_support` (further from the NTC cycle than 99% of NTC cells), and `phase_confidence`.

### What the new-cell file must contain
- `layers['raw']`: the **raw, unscaled** 4i features (same columns as the NTC 4i file). The pre-scaled `X` of the knockout file used another NTC file's constants and is not used;
- `obsm['X_phase_raw']` and `uns['phase_var_names']`: raw phase features;
- `obs['well']` (`A/1/0`-style), `obs['dna_content_log2']`, `obs['dna_reliable']`;
- the input's own derived columns (`cycle_stage`, `phase_predicted_*`) are carried as `ref_*` for **comparison only** and never used as inputs.

Cells in a well not seen in the NTC reference would be normalised with pooled NTC constants and flagged (`unseen_well`).

### Pickled model
The phase-contrast model is a scikit-learn pipeline saved with `joblib` (`phase_ridge_model.joblib`). It is not part of this repository because a pickle must only be loaded from a trusted source; notebook 18 re-creates it (about 10 minutes).

## Output columns

`reports/ops0077_phase_predicts_cc_phase/phase_predicted_cc_phase.csv` (NTC):

| column | meaning |
|---|---|
| `cc_phase` | the label |
| `pred_phase`, `pred_sin`, `pred_cos` | predicted phase of the best model (out-of-fold for labelled cells) |
| `pred_confidence` | length of the predicted (sin, cos) vector |
| `pred_is_out_of_fold` | true when the cell was predicted by a model that did not see its tile |
| `error_turns` | signed circular error |
| `pred_*_ridge`, `pred_*_hgb_pca` | per-model predictions |

`reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv` (knockout cells): `cc_phase_stain`, `cc_phase_phase`, `cc_phase`, `cc_phase_source`, `sector_stain`, `sector_phase`, `sector_final`, `cc_curve_dist`, `out_of_ntc_support`, `phase_confidence`, the normalised axes `z_*`, the carried metadata (`gene`, `sgRNA`, `well`, ...) and the `ref_*` comparison columns.
