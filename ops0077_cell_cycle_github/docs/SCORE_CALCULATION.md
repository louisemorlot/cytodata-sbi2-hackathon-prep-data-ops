# How the scores are calculated (NTC cells)

All steps in this document use the **non-targeting-control (NTC) cells only**: 52,220 cells in three wells (`A/1/0`, `A/2/0`, `A/3/0`). The results are in `reports/`; the code is in `scripts/` and `notebooks/`.

Contents: [1. Per-cell inputs](#1-per-cell-inputs) · [2. PCA embeddings](#2-pca-embeddings-notebook-13) · [3. Markers + DNA: the circular score](#3-markers--dna-the-circular-cell-cycle-phase-main-score) · [4. Markers only: the open-ended score](#4-markers-only-the-open-ended-score-no-dna) · [5. Heuristic marker proxy](#5-the-earlier-heuristic-marker-proxy) · [6. Output columns](#6-output-columns) · [7. Limitations](#7-limitations)

## Biology behind the markers

| state | p21 | pRb (phosphorylated Rb) | DNA content |
|---|---|---|---|
| G0 / early G1 | high | low | 2N |
| late G1 | falling | rising (CDK4/6 and CDK2 activity) | 2N |
| S | lowest (p21 is degraded during replication) | high | rising from 2N to 4N |
| G2 / M | rising again | highest | 4N |

The most informative marker quantity is the **phosphorylated fraction of Rb** (pRb relative to total Rb): raw pRb mostly follows cell size and staining, the ratio follows CDK activity. p21 separates early G1 (p21 high, pRb low) from G2 (p21 moderate, pRb high). Marker intensities alone do not separate S from G2/M well; **DNA content** does.

## 1. Per-cell inputs

**Nuclear marker intensity.** The nuclear mean intensity of each marker (pRb, Rb, p21) is derived from the CellProfiler features of the cell and its cytoplasm:

- nuclear area (px) = cell area − cytoplasm area, where each area = integrated intensity / mean intensity;
- nuclear mean intensity = (cell integrated − cytoplasm integrated) / nuclear area;
- cells with a nuclear area of 50 px or less are not scored.

**DNA content** (`dna_rel_R234`, from the DAPI file): per cell, the mean over imaging rounds R2 to R4 of log2(integrated DAPI − non-nuclear background × area), minus that round's median. It is a *relative* log2 DNA content (a doubling is +1). Round 1 saturates and is not used. Only cells flagged `dna_reliable` keep a DNA value (36,014 of 42,940 cells in the DAPI table).

**Normalisation, per well** (the same recipe is applied to every axis):

1. background = 0.5th percentile of the well's nuclear mean for that marker, subtracted; the result is floored at the 1st percentile of the positive remainder so dim cells cannot collapse to `log(0)`;
2. natural log;
3. robust z-score: (x − median) / (1.4826 × MAD), per well, then clipped to ±4.

**The four axes** used by the scores:

| axis | definition |
|---|---|
| `phos` | log(pRb) − log(Rb), the phosphorylated fraction of Rb |
| `p21` | log(p21) |
| `pRb` | log(pRb) |
| `dna` | relative DNA content |

Cells with all three markers: 42,410. Cells that also have a reliable DNA value: **35,633**.

## 2. PCA embeddings (notebook 13)

`notebooks/13_ops0077_marker_pca_cell_cycle_clock.ipynb` builds two kinds of embeddings of the NTC cells from the p21, Rb and pRb features and checks what they capture. The results are in `reports/ops0077_cell_cycle_clock/`.

- **2-D marker embedding** (read from the `..._marker_pca.h5ad` file): PCA of the log1p, standardised mean intensities of the three markers. PC1 explains 49% and PC2 28% of the variance (`embedding_variance.csv`). PC1 correlates with all three markers (Spearman about −0.6 to −0.75), so it mostly tracks overall staining intensity; there is a clear well effect (`well_association_tests.csv`).
- **Per-cell feature-block embeddings** (`per_cell_pca_embeddings.csv.gz`, two PCs each):
  - one PCA per channel (`pRb_*`, `Rb_*`, `p21_*`) from that channel's OrganelleProfiler feature block (intensity, morphology, localisation, network; 194 / 221 / 194 features);
  - for each pair (pRb–Rb, pRb–p21, Rb–p21): the two channels' blocks together (`*_channels_*`), only the pair's CellProfiler colocalisation features (`*_coloc_*`, 16 features) and both together (`*_integrated_*`);
  - features with < 80% finite values (among cells with measurements for that block) or no variation are dropped (`per_cell_feature_qc.csv`); remaining gaps are median-imputed; features are median-centred, IQR-scaled and clipped; the PCA is fit on a random sample of 12,000 cells and applied to all cells;
  - explained variance is in `per_cell_pca_variance.csv`.
- The notebook also audits the data for phase labels and finds none (`phase_label_feature_audit.csv`): this analysis **predates** the DNA content table, and its conclusion that no supervised clock can be trained refers to that situation.

These PCs are an optional extra layer for the scores below (section 3); they are not cell cycle scores themselves.

## 3. Markers + DNA: the circular cell cycle phase (main score)

Script: `scripts/ops0077_cell_cycle_score.py`. Runner: `notebooks/15_run_ops0077_cell_cycle_score_markers_dna.ipynb`. Explore: `notebooks/14_ops0077_cell_cycle_score_exploration.ipynb`. Results: `reports/ops0077_cell_cycle_score/`.

### The idea
Each cell is a point in the 4-D space of the axes `[phos, p21, pRb, dna]`. Cells that cycle form a **closed loop** in that space. The score of a cell is *where along the loop* it lies, as a fraction of the loop length.

### The steps
1. **Weights.** The axes are weighted `[1, 1, 1, 1.5]`: DNA content gets more weight because it is the cleanest separator of S and G2/M.
2. **Closed principal curve.** A closed curve with 16 nodes is fitted through the 35,633 cells (on a random subsample of 20,000, then all cells are projected onto it):
   - initialise each cell's phase from the angle in the plane of the first two principal components of the weighted axes;
   - repeat 30 times: cut the cells into 16 equal-count bins along the current phase, set each node to the mean of its bin, smooth the nodes with a periodic (wrap-around) filter, then re-project every cell onto the closed polyline;
   - the phase of a cell is the arc length of its projection, divided by the loop length (so it lies in [0, 1)).
3. **Anchoring** (the loop has no natural start or direction):
   - phase 0 is the node with the **lowest DNA content** (the start of G1);
   - the forward direction is the one in which the **p21 minimum (S) is reached before the DNA maximum (G2/M)**.
4. **Distance to the curve** (`cc_curve_dist_dna`) is kept for each cell as a quality measure.

Check on the fitted loop (`reports/ops0077_cell_cycle_score/marker_profiles.png`): DNA content peaks at phase ~0.55 to 0.6, p21 is lowest around S, pRb and the phosphorylated fraction rise in late G1; the phase explains 94% of the DNA-content variance, 54% of p21 and 27% of the phosphorylated fraction (two-harmonic fit), and has essentially no well dependence (η² = 0.0009).

### The optional PC layer
Per-channel and pairwise PCs (section 2) can refine the ordering. Each PC is screened and kept only if it
- is available for more than half of the cells,
- has a circular association with the DNA-only phase (R² of a two-harmonic fit ≥ 0.03),
- has little dependence on the well (η² < 0.10),
- is not dominated by nuclear size (|Spearman ρ| with nuclear area < 0.6).

Two PCs pass (`prb_p21_coloc_PC1`, `prb_p21_coloc_PC2`; `reports/ops0077_cell_cycle_score/pc_screening.csv`). They are added to the loop fit with a combined weight of 1 (so they refine and do not dominate), starting from the DNA-only phase. The result agrees closely with the DNA-only phase (circular agreement 0.97, mean absolute difference 0.028 turns).

### The final score
`cc_phase` = the markers + DNA + PCs phase where available, otherwise the markers + DNA phase (`cc_phase_source` records which): **35,633 cells**. Cells without three markers and a reliable DNA value get **no** `cc_phase`.

### Markers only on the same loop (QC, not a score)
The same closed-curve fit on `[phos, p21, pRb]` alone (`cc_phase_markers`) is aligned to the DNA-anchored phase afterwards. It agrees poorly (circular agreement 0.34): the three markers alone do not recover the cycle. It is kept as a QC column only and is **not** used as a fallback. This quantifies what DNA content adds.

## 4. Markers only: the open-ended score (no DNA)

Script: `scripts/ops0077_cell_cycle_score_only_markers.py`. Results: `reports/ops0077_cell_cycle_score_only_markers/` (all file names carry the suffix `_only_markers`). It was written to show later what DNA content adds.

- Same axes without DNA (`[phos, p21, pRb]`), same normalisation, 42,456 cells.
- An **open** principal curve (12 nodes, 30 iterations, 4 smoothing passes), oriented so that pRb and the phosphorylated fraction increase along it; the score is the percentile rank of the position along the curve (0 to 1, ends do **not** join).
- **Soft states.** A 4-component Gaussian mixture is initialised from the biology table above (prior means for G0/early G1, late G1, S, G2/M in the z-space of the three axes) and gives a probability per state and an expected stage.
- **PC layer** with the same screening idea (|ρ| with the score ≥ 0.15, well η² < 0.10, |ρ| with area < 0.6): three PCs kept; 34,716 cells have the combined score.
- Checks: Spearman 0.86 between the score and the mixture's expected stage, 0.95 between the intensity-only and the PC-combined scores, about zero relation to nuclear area (ρ = 0.065), negligible well effect (η² = 0.003).
- Its limitation: without DNA the ends of the cycle (G2/M and G1) cannot be joined, so it is an ordering of cells, not a circular phase.

## 5. The earlier heuristic marker proxy

Script: `scripts/score_ops0077_cell_cycle_markers.py`; notebook `14_ops0077_marker_cell_cycle_score.ipynb`; results `reports/ops0077_cell_cycle_score/per_cell_marker_stage_scores.csv.gz`, `score_summary.json`, `marker_score_pca_visualization.png`.

An exploratory, hypothesis-based **G1-to-S progression proxy** from the OrganelleProfiler mean-intensity features of pRb, Rb and p21 (17,161 cells with all three). Per well it takes the percentile of `0.65 × (robust z(log1p pRb) − robust z(log1p Rb)) + 0.35 × robust z(−log1p p21)` and labels cells as G1-like, intermediate or G1/S-entry-like (and "ambiguous" when the two evidence terms disagree). Its own `score_summary.json` lists its limitations; it is not a validated phase classifier and is not circular.

## 6. Output columns

`reports/ops0077_cell_cycle_score/cell_cycle_scores.csv.gz` (one row per NTC cell):

| column | meaning |
|---|---|
| `cc_phase` | final circular phase in [0, 1) (see above) |
| `cc_phase_source` | `dna+pcs`, or `dna`; empty if the cell has no score |
| `cc_phase_dna` | markers + DNA phase |
| `cc_phase_dna_pcs` | markers + DNA + PCs phase |
| `cc_phase_markers` | markers-only phase, aligned to the DNA phase (QC only) |
| `cc_curve_dist_dna`, `cc_curve_dist_markers` | distance of the cell from the fitted loop |
| `z_phos`, `z_p21`, `z_pRb`, `z_dna` | the normalised axes |
| `cc_phase_angle_rad` | `cc_phase` as an angle in radians |
| `nuc_area_px`, `dapi_log_area`, `dna_rel_R234`, `dna_reliable`, `dense_field`, `nucleus_found_R234`, `dna_R1_any_saturated`, `x_cp`, `y_cp`, `well`, `cell_key` | QC and position covariates |

**Using the phase.** It is circular. Compare two phases with the circular difference `((a − b + 0.5) mod 1) − 0.5`, or use `sin(2π·phase)` and `cos(2π·phase)`; never subtract raw phases.

## 7. Limitations

- DNA content is an input of the score, so it cannot be used as independent validation of it.
- The anchoring (phase 0 at the lowest-DNA node, forward direction by the p21 minimum) is a convention; it is approximate, not an exact G1/S boundary.
- Cells are unevenly spread around the loop (many in G1 and G2/M, few in S), so S is the least well resolved part.
- The score needs all three markers and a reliable DNA value; about one third of the NTC cells have no score.
