# How the results are evaluated

There is **no orthogonal ground truth** for a cell's cell cycle phase in these data, and DNA content is itself an input of the score. Evaluation therefore rests on four kinds of evidence: *internal consistency* (does the score behave like the biology it is built from), *held-out prediction* (does it generalise to cells and wells the model never saw), *agreement of independent measurements* (different instruments, disjoint cells, a gating shipped in the input file), and *comparison with baselines and confounds* (does the model learn more than cell size or position).

Contents: [Part 2: the NTC score](#part-2-checking-the-ntc-score) · [Part 3: the phase-contrast prediction](#part-3-checking-the-phase-contrast-prediction) · [Part 4: the knockout cells](#part-4-evaluating-the-knockout-cells) · [Figure and table index](#figure-and-table-index) · [Known corrections](#known-corrections-and-pitfalls)

## Part 2: checking the NTC score

Results: `reports/ops0077_cell_cycle_score/`, explored in `notebooks/14_ops0077_cell_cycle_score_exploration.ipynb`.

| check | where | result |
|---|---|---|
| Do the markers and DNA content follow the expected curve around the loop? | `marker_profiles.png`, `phase_planes.png` | DNA peaks at phase ~0.55 to 0.6, p21 is lowest around S, pRb and phosphorylated Rb rise in late G1 |
| How much of each axis does the phase explain? (two-harmonic circular fit) | script log | DNA 94%, p21 54%, phosphorylated Rb 27% |
| Is the score a batch artefact? | `score_distributions.png`, script log | well effect η² = 0.0009 |
| Is it driven by nuclear size? | script log | circular R² with nuclear area 0.22; size alone is a weak explanation |
| What does DNA content add? | `umap_overview.png`, `score_distributions.png` | the markers-only loop recovers the DNA-anchored phase poorly (circular agreement 0.34): DNA content is essential |
| Do the PCs help? | `pc_screening.csv` | two PCs pass the screen; they change the phase by 0.028 turns on average (agreement 0.97) |

The open-ended markers-only analysis (`reports/ops0077_cell_cycle_score_only_markers/`) is evaluated the same way: score vs the mixture model's expected stage (Spearman 0.86), intensity-only vs PC-combined (0.95), no relation to nuclear area, negligible well effect.

## Part 3: checking the phase-contrast prediction

Results: `reports/ops0077_phase_predicts_cc_phase/`, explored in `notebooks/17_ops0077_phase_prediction_exploration.ipynb`.

### Metrics (all circular, computed on out-of-fold predictions)

| metric | meaning | random guess | perfect |
|---|---|---|---|
| `mae_turns` | mean absolute circular error, in fractions of a cycle | 0.25 | 0 |
| `within_0p10` | fraction of cells predicted within 0.10 of a cycle | 0.20 | 1 |
| `acc_4bins`, `acc_8bins` | accuracy after cutting the loop into 4 or 8 equal phase bins | 0.25, 0.125 | 1 |
| `circ_agreement` | length of the mean error vector (1 = all errors the same, 0 = unrelated) | 0 | 1 |

### Validation design
Grouped by tile (5-fold) and by well (leave one well out), so neighbouring cells never straddle train and test. Baselines: size only, registration distance and cell number only, random guess. A feature-group ablation and a comparison on the subset that has CellProfiler phase texture are also run.

### Results

| model | MAE (turns) | within 0.10 | accuracy, 8 bins |
|---|---|---|---|
| ridge, tile-grouped CV | **0.105** | 0.63 | 0.43 |
| ridge, leave-one-well-out | 0.113 | 0.60 | 0.40 |
| gradient boosting, tile-grouped | 0.124 | 0.57 | 0.39 |
| baseline: cell size only | 0.174 | 0.42 | 0.29 |
| baseline: registration distance + cell number | 0.188 | 0.41 | 0.26 |
| random guess | 0.25 | 0.20 | 0.125 |

- **Confidence is informative:** the error falls from 0.20 in the least confident decile to 0.065 in the most confident (`qc_error_drivers.png`).
- **Error is highest where cells are rare:** 0.08 to 0.09 in G1 and G2/M, about 0.23 around phase 0.75 to 0.875 (`error_by_true_phase.csv`, `qc_error_by_phase.png`); predictions are pulled towards the dense G1 and G2/M clusters.
- **Registration distance matters:** error rises from 0.097 (best-registered fifth) to 0.119 (worst) (`error_by_registration_distance.csv`); part of the "error" may be real progression between live imaging and fixation, or a wrong cell match.
- **No well effect** worth noting (`error_by_well.csv`).
- **Feature groups are redundant:** the 2-D phase and the focus-stack groups each score about the same alone (MAE 0.114 / 0.108) and dropping either costs little (`ablation_by_group.csv`).
- **CellProfiler phase texture adds nothing** on the 8,341 cells that have it (0.091 with, 0.090 without).
- **Biology is preserved:** the median DNA and marker profiles along the *predicted* phase look like those along the true phase (`qc_profiles_true_vs_predicted.png`).

## Part 4: evaluating the knockout cells

Results: `reports/ops0077_cc_transfer/` (`ko3_evaluation/`, `ko3_guide_analysis/`, `ko3_density_profiles/`), explored in notebooks 20 and 21.

### Design principles
- **Sectors fixed first.** G1, S, G2/M and exit are defined from the NTC curve before any knockout cell is looked at (see [PREDICTION.md](PREDICTION.md)); nothing is tuned to the result.
- **The expectation is stated, not tested as pass or fail.** From the experimental design: CCNB1 → G1, MCM3 → S, CENPC → G2/M. It is drawn on the figures; effect sizes come with bootstrap intervals. (Loss of cyclin B1 would classically be expected to arrest cells in G2/M, so a G2/M shift for CCNB1 would not by itself be a failure of the score.)
- **Two independent scores** are always compared: the stain score (markers + DNA, precise, 1,093 of 1,603 cells) and the phase-contrast score (all cells, less precise).

### Checks of the transfer itself (is the score meaningful for knockout cells?)

| check | file | result |
|---|---|---|
| Do knockout cells lie on the NTC cycle? DNA and marker profiles along the phase, knockout vs NTC | `fig5_profiles_knockout_vs_ntc.png`, density figures | the profiles coincide with NTC |
| How far are they from the NTC curve? | `fig6_support_and_distance_to_ntc_curve.png`, `coverage_by_gene.csv` | 2.2% (35 of 1,603) are beyond the NTC 99th percentile (CCNB1 3.3%, CENPC 3.3%, MCM3 0.8%): large, p21-high cells with DNA beyond 4N |
| Do the two scores agree? | `fig7_stain_vs_phase_contrast_scores.png`, `stain_vs_phase_agreement.csv` | mean error 0.09 to 0.10 turns, about as for NTC cells |
| Does it agree with the gating shipped in the input file (`cycle_stage`, DNA × pRb with an NTC threshold) and the file's own phase-contrast stage? These were not used as inputs | `fig8_comparison_with_file_reference_columns.png`, `crosstab_*.csv` | 94% / 86% of its G0/G1-early / G1-late cells fall in our G1, 100% of its G2/M cells in our G2/M |
| Is the DNA content on the same reference as the NTC? The lowest-DNA (2N) peak exists in every population, so it should sit at the same place | `fig9_dna_reference_check.png`, `dna_g1_peak_check.csv` | peak positions agree within about 0.03 log2 units per well (QC only, not used by any score) |

### Gene-level result

Sector fractions of the stain-scored cells (NTC: G1 0.63, S 0.07, G2/M 0.25, exit 0.05). Change against NTC, with 95% bootstrap intervals; the last column is a permutation test of a shift of the mean phase.

| gene | cells scored | G1 | S | G2/M | permutation p |
|---|---|---|---|---|---|
| CCNB1 | 279 | +0.02 [−0.03, +0.08] | −0.03 [−0.05, −0.01] | +0.02 [−0.03, +0.07] | 0.0002 |
| MCM3 | 515 | −0.00 [−0.05, +0.04] | −0.01 [−0.03, +0.02] | +0.03 [−0.01, +0.07] | 0.15 |
| CENPC | 299 | +0.06 [+0.01, +0.12] | −0.02 [−0.05, +0.01] | −0.03 [−0.08, +0.01] | 0.004 |

At the gene level the expected arrests are **not clearly visible**: no S enrichment for MCM3, no G2/M enrichment for CENPC (a small move towards G1), and CCNB1 mostly unchanged apart from fewer S cells. The phase-contrast score on all cells gives the same picture (`fig2_sector_fractions.png`, `ko3_density_profiles/`), so the result is not a selection effect of which cells have a stain score.

### Guide-level result (`guide_level_analysis.py`, notebook 21)

A **guide is the 20-nt protospacer** (`sgRNA`): 4 per gene. For each guide the change in sector fractions is computed against the NTC fractions **re-weighted to the wells the guide's cells sit in** (so a different well mix cannot create a difference), with bootstrap intervals; a mean-phase shift test against NTC with Benjamini–Hochberg correction across the 12 guides; and a heterogeneity test across the guides of a gene, with the null made by shuffling guide labels **within wells**.

- **Guides of one gene differ** (stain score, permutation p): CCNB1 0.016, MCM3 0.046, CENPC 0.027 (with all cells: 0.10, 0.006, 0.06).
- **Two patterns coexist.** More G2/M: CCNB1 `ATATTTGC…` (+15 points) and `GTCAGACC…` (+15), MCM3 `CAGTATGA…` (+12) and `ATATCCAC…` (+7). More G1: CCNB1 `TGAGGCCA…` (+12 G1; q = 0.002) and CENPC `ATCAACTT…` (+12.5; q = 0.004). Other guides sit close to NTC.
- **The guide effects are reproducible.** The same guide's phase composition is measured in two **disjoint sets of cells by different instruments** (cells with a stain score, scored by the stain; and the remaining cells, scored by phase contrast). Across the 12 guides the G1 fraction correlates at Spearman 0.82 (p = 0.001), G2/M at 0.50 (p = 0.095). The phase-contrast sectors have a different absolute calibration, so only the ranking is compared. The narrow S sector agrees poorly (0.11). Guide ranking agrees even better on the same cells (0.83 and 0.87).
- **A shift shared by all three genes.** The end of the loop (phase ≥ 0.92: p21-high, low DNA, just after division) holds 20 to 26% of the knockout cells against 13% in NTC (well-matched differences +0.07 to +0.13, intervals excluding 0), driven mostly by CCNB1 `TGAGGCCA…` and CENPC `ATCAACTT…` (about 40%). A shift common to different targets is not gene-specific. A general p53 / p21 response to Cas9 cutting, which non-targeting guides do not trigger, is one possible explanation; **it is not tested here** (a cutting control, for example intergenic or safe-targeting guides, would test it). Notebook 21 section 7 reproduces these numbers.
- **Other phenotypes.** Across guides the G2/M shift tracks the fraction of cells with DNA beyond 4N (Spearman 0.87) and the phosphorylated Rb level (0.66), but both are inputs of the score, so this is expected by construction and not independent confirmation; nuclear area (0.25) and p21 (about 0) show no clear relation.

### Power
A guide of n scored cells gives a 95% interval of about ±1.96·√(p(1−p)/n) on a fraction p; for p = 0.25 that is ±0.12 at n = 50 and ±0.07 at n = 150 (`power_by_guide.csv`). The 4 guides per gene have 36 to 185 stain-scored cells.

### Caveats
- One run per script; there is no replicate-level analysis beyond the per-guide and per-well tables.
- The knockout cells are compared with non-targeting controls: general effects of Cas9 cutting, guide efficiency, the time between transduction and fixation, and partial or mixed knockouts are not separated from gene-specific effects. The input file contains no expression or knockout-efficiency data.
- About one third of the knockout cells have no stain score (missing markers or no reliable DNA value) and only the less precise phase-contrast score.
- 12 guides are tested; three are significant after correction and two more are borderline (q ≈ 0.06).

## Figure and table index

**`reports/ops0077_cell_cycle_clock/`** (part 1): `embedding_by_marker.png`, `embedding_by_well.png` (2-D marker embedding); `per_cell_pca_comparison.png`, `per_cell_pca_pairwise_comparison.png` (per-channel and pairwise embeddings); `per_cell_pca_embeddings.csv`, `per_cell_pca_variance.csv`, `per_cell_feature_qc.csv`, `embedding_variance.csv`, `embedding_marker_correlations.csv`, `well_association_tests.csv`, `channel_intensity_correlations.csv`, `phase_label_feature_audit.csv`, `analysis_summary.json`.

**`reports/ops0077_cell_cycle_score/`** (part 2, markers + DNA): `cell_cycle_scores.csv` (all cells), `umap_overview.png` (UMAP coloured by phase and axes), `phase_planes.png` (pairs of axes coloured by phase), `marker_profiles.png` (profiles around the cycle and cell density), `score_distributions.png` (by well; fits compared), `pc_screening.csv`, `summary.json`. The heuristic marker proxy adds `per_cell_marker_stage_scores.csv`, `score_summary.json`, `marker_score_pca_visualization.png`.

**`reports/ops0077_cell_cycle_score_only_markers/`** (part 2, no DNA): the same file types with the suffix `_only_markers`.

**`reports/ops0077_phase_predicts_cc_phase/`** (part 3): `qc_pred_vs_true.png` (predicted vs true, error histograms), `qc_error_by_phase.png` (error along the cycle and confusion of 8 phase bins), `qc_error_drivers.png` (confidence, registration distance, well, label quality, cell number), `qc_model_comparison.png` (models vs baselines), `qc_profiles_true_vs_predicted.png`, `qc_feature_groups_and_top_features.png`; tables `phase_predicted_cc_phase.csv`, `summary.json`, `error_by_*.csv`, `ablation_by_group.csv`, `feature_importance.csv`, `feature_table.csv`.

**`reports/ops0077_cc_transfer/reference/`** (frozen NTC reference): `ntc_cc_reference.json`, `ntc_reference_scores.csv`, `ntc_profiles_by_phase.csv`, `phase_model_features.json`, `reference_summary.json` (reproduction check, model cross-validation).

**`reports/ops0077_cc_transfer/ko3_scores/`**: `cc_phase_new_cells.csv`, `summary.json`.

**`reports/ops0077_cc_transfer/ko3_evaluation/`**: `fig1` phase distributions vs NTC; `fig2` sector fractions with intervals; `fig3` change vs NTC heat map; `fig4` per guide and per well; `fig5` profiles vs NTC; `fig6` distance to the NTC curve; `fig7` stain vs phase-contrast score; `fig8` comparison with the file's own gating; `fig9` DNA reference check; tables `sector_shift_vs_ntc.csv`, `sector_fractions_with_ci.csv`, `fractions_per_guide.csv`, `fractions_per_well.csv`, `coverage_by_gene.csv`, `stain_vs_phase_agreement.csv`, `dna_g1_peak_check.csv`, `crosstab_*.csv`, `evaluation_summary.json`.

**`reports/ops0077_cc_transfer/ko3_guide_analysis/`**: `g1` per-guide change vs NTC (forest plots), `g2` density along the cycle per guide, `g3` reproducibility across disjoint cells, `g4` other phenotypes ranked by the G2/M shift, `g5` heterogeneity tests, `g6` well composition of each guide; tables `guide_sector_shifts.csv`, `heterogeneity_across_guides.csv`, `guide_effect_concordance.csv`, `phenotypes_by_guide.csv`, `phenotype_vs_g2m_shift_across_guides.csv`, `power_by_guide.csv`.

**`reports/ops0077_cc_transfer/ko3_density_profiles/`**: `density_profiles_by_gene_stain_score.png` and `density_profiles_by_gene_phase_contrast_score.png` (marker and DNA profiles with the cell density along the cycle, per gene, NTC density as a grey step), `sector_fractions_by_gene.csv`.

## Known corrections and pitfalls

- **Use `sgRNA`, not `barcode`, to identify a guide.** The `barcode` column of the knockout file is truncated to 9 or 10 nucleotides, so the same guide appears under two barcodes. An early version of the per-guide plots used the barcode and doubled the number of "guides"; the scripts and notebooks in this repository use the 20-nt `sgRNA`.
- **Missing sectors become `NaN` after saving to CSV** (they are empty strings in memory); the evaluation code drops both when computing fractions.
- **Never renormalise new cells on themselves** (see [PREDICTION.md](PREDICTION.md)).
