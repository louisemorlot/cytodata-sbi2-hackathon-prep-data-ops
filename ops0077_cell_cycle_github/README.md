# OPS0077: a cell cycle phase score from 4i markers and DNA content, and its transfer to phase contrast and gene-knockout cells

This repository contains the code, notebooks and results for estimating where every single cell sits in the cell cycle in the Biohub OPS0077 optical pooled screen (experiment `ops0144_20260406`), using only features that were already extracted per cell.

The work has four parts:

| part | question | main outputs |
|---|---|---|
| **1. PCA embeddings of the NTC cells** | how do the p21 / Rb / pRb features (per channel and per channel pair) organise the non-targeting-control (NTC) cells? | per-cell PCs, variance, QC |
| **2. Cell cycle score in NTC cells** | can p21, Rb, pRb and DNA content be turned into a continuous, *circular* cell cycle phase per cell? | `cc_phase` per cell (markers only, and markers + DNA) |
| **3. Prediction from phase contrast (NTC)** | how well do label-free phase-contrast features predict that score? | circular regression, out-of-fold QC |
| **4. Transfer to gene-knockout cells** | what does the NTC-derived score say about cells with CCNB1, MCM3 or CENPC targeted? | scores for 1,603 knockout cells, per-gene and per-guide evaluation |

The cell cycle is circular (the end of mitosis is the start of G1), so the score is a **phase in [0, 1)** where 0 and 1 are the same point: 0 is the start of G1, it increases through S to G2/M, and wraps back to G1.

> **Status.** This is exploratory analysis. The score is derived from marker and DNA-content measurements and is checked against internal consistency and an independent gating in the input file. It is not validated against an orthogonal ground truth (for example live-cell cycle reporters), and the knockout results are descriptive.

## Documentation

| file | read it to understand |
|---|---|
| [docs/SCORE_CALCULATION.md](docs/SCORE_CALCULATION.md) | how the PCA embeddings and the cell cycle scores are calculated |
| [docs/PREDICTION.md](docs/PREDICTION.md) | how the phase-contrast model is built and how the NTC score is transferred to new cells |
| [docs/EVALUATION.md](docs/EVALUATION.md) | how every result is checked, what the QC figures show, and the main findings and caveats |
| [figures/README.md](figures/README.md) | every figure, with a caption of what it shows |
| [data/README.md](data/README.md) | which input files are needed (they are **not** in this repository) and where they go |

## Repository layout

```
scripts/                          analysis code (run from the repository root)
  ops0077_cell_cycle_score.py             part 2: circular score, markers + DNA (+ PCs)
  ops0077_cell_cycle_score_only_markers.py part 2: open-ended score, markers only (no DNA)
  score_ops0077_cell_cycle_markers.py     part 2: earlier heuristic G1-to-S marker proxy
  ops0077_phase_predict_cc_phase.py       part 3: phase contrast -> cc_phase, circular regression
  ops0077_cc_transfer/                    part 4: freeze the NTC score, apply it to new cells, evaluate
  collect_figures.py                      copies every figure from reports/ into figures/
notebooks/                        runner notebooks (re-run an analysis) and exploration notebooks
reports/                          results of every analysis (tables and figures), one folder per part
figures/                          all figures in one place, grouped by part, with captions (figures/README.md)
docs/                             the explanations (see above)
data/                             placeholder only; put the input .h5ad files here (never committed)
```

### What is where

| part | code | notebooks | results |
|---|---|---|---|
| 1. PCA embeddings | in the notebook | `13_ops0077_marker_pca_cell_cycle_clock` | `reports/ops0077_cell_cycle_clock/` |
| 2. Score: markers only, open-ended | `scripts/ops0077_cell_cycle_score_only_markers.py` | | `reports/ops0077_cell_cycle_score_only_markers/` |
| 2. Score: heuristic marker proxy | `scripts/score_ops0077_cell_cycle_markers.py` | `14_ops0077_marker_cell_cycle_score` | `reports/ops0077_cell_cycle_score/` (`per_cell_marker_stage_scores.csv.gz`, `score_summary.json`, `marker_score_pca_visualization.png`) |
| 2. Score: markers + DNA, circular | `scripts/ops0077_cell_cycle_score.py` | `15_run_...` (run), `14_ops0077_cell_cycle_score_exploration` (explore) | `reports/ops0077_cell_cycle_score/` (the other files) |
| 3. Phase-contrast prediction (NTC) | `scripts/ops0077_phase_predict_cc_phase.py` | `16_run_...` (run), `17_ops0077_phase_prediction_exploration` (explore) | `reports/ops0077_phase_predicts_cc_phase/` |
| 4. Knockout cells | `scripts/ops0077_cc_transfer/` | `18_...` (freeze reference), `19_...` (score new cells), `20_...` (evaluate), `21_...` (guides and density) | `reports/ops0077_cc_transfer/` |

### Figures

All 36 figures are in [figures/](figures/README.md), grouped by part and captioned. They are copies of the figures in `reports/` (the notebooks read them there); `python scripts/collect_figures.py .` refreshes the copies.

### Compressed result tables

The large per-cell tables are stored **gzip-compressed** (`*.csv.gz`) to keep the repository small (these tables went from 70 MB to 27 MB). The scripts read and write them directly and pandas does it transparently from the file extension, for example `pd.read_csv("reports/ops0077_cell_cycle_score/cell_cycle_scores.csv.gz", index_col=0)`. No extra step is needed. The tables are:

| table | content | read by (inputs of) |
|---|---|---|
| `reports/ops0077_cell_cycle_clock/per_cell_pca_embeddings.csv.gz` (4 MB) | per-cell PCs of the channel and channel-pair feature blocks | **scoring**: both score scripts (optional PC layer) |
| `reports/ops0077_cell_cycle_score/cell_cycle_scores.csv.gz` (6 MB) | the NTC cell cycle phase of every cell (markers + DNA) | **prediction**: `ops0077_phase_predict_cc_phase.py` (the label) and `build_ntc_reference.py` (reproduction check); notebooks 14, 17, 18 |
| `reports/ops0077_cc_transfer/reference/ntc_reference_scores.csv.gz` (3 MB) | NTC cells re-scored with the frozen reference | evaluation scripts and notebooks 20, 21 |
| `reports/ops0077_cc_transfer/ko3_scores/cc_phase_new_cells.csv.gz` (0.2 MB) | scores of the 1,603 knockout cells | evaluation scripts and notebooks 20, 21 |
| `reports/ops0077_phase_predicts_cc_phase/phase_predicted_cc_phase.csv.gz` (6 MB) | out-of-fold phase-contrast predictions (NTC) | notebook 17 |
| `reports/ops0077_cell_cycle_score_only_markers/cell_cycle_scores_only_markers.csv.gz` (5 MB) | open-ended markers-only score | results |
| `reports/ops0077_cell_cycle_score/per_cell_marker_stage_scores.csv.gz` (3 MB) | heuristic marker proxy | results |

The smaller tables (QC summaries, tests, per-guide tables) stay as plain CSV.

## How to run

1. Install the dependencies (Python 3.10 or newer):
   ```bash
   pip install -r requirements.txt
   ```
2. Put the input files in `data/altair/` as described in [data/README.md](data/README.md).
3. Run the parts in order, either with the runner notebooks or the scripts (always from the repository root):

   | step | notebook | or script |
   |---|---|---|
   | PCA embeddings | `notebooks/13_ops0077_marker_pca_cell_cycle_clock.ipynb` | n/a |
   | score, markers only | n/a | `python scripts/ops0077_cell_cycle_score_only_markers.py .` |
   | score, markers + DNA | `notebooks/15_run_ops0077_cell_cycle_score_markers_dna.ipynb` | `python scripts/ops0077_cell_cycle_score.py .` |
   | phase-contrast prediction | `notebooks/16_run_ops0077_phase_predict_cc_phase.ipynb` | `python scripts/ops0077_phase_predict_cc_phase.py .` |
   | freeze the NTC reference | `notebooks/18_build_ops0077_ntc_cc_phase_reference.ipynb` | `python scripts/ops0077_cc_transfer/build_ntc_reference.py .` |
   | score the knockout cells | `notebooks/19_run_ops0077_cc_phase_transfer_to_new_cells.ipynb` | `python scripts/ops0077_cc_transfer/apply_cc_transfer.py .` |
   | evaluate | `notebooks/20_ops0077_ko_cc_phase_evaluation.ipynb` | `python scripts/ops0077_cc_transfer/evaluate_cc_transfer.py .` |
   | guides and density plots | `notebooks/21_ops0077_ko_guide_level_and_density_profiles.ipynb` | `python scripts/ops0077_cc_transfer/guide_level_analysis.py .` and `python scripts/ops0077_cc_transfer/plot_density_profiles_by_gene.py .` |

   Notebooks 14 (both), 17 and 20 are for exploring existing results.

All results are already in `reports/`, so notebooks 17, 20 and 21 can be explored without the large input files. Notebook 19 needs the **pickled phase-contrast model**, which is not stored in this repository (a pickle must only be loaded from a trusted source): create it first with notebook 18 (about 10 minutes).

### Which notebooks need the large `.h5ad` inputs

| needs the `.h5ad` files | can run from `reports/` alone |
|---|---|
| 13, 14 (marker score), 15, 16, 18, 19, and the "link to the full feature matrix" section of 14 (exploration) | 20, 21 (21 re-runs its two scripts, which only read CSV files), 17 except its last section |

## The main numbers

- **NTC score.** 35,633 of 52,220 NTC cells have all four inputs (three markers and a reliable DNA measurement) and get a circular `cc_phase`. Phase 0 is anchored at the lowest-DNA part of the loop; the DNA content peaks at phase ~0.55 to 0.6 and p21 is lowest around S.
- **Phase contrast alone** predicts `cc_phase` with a mean circular error of **0.105 turns** on cells never seen in training (a random guess gives 0.25), and 0.113 on a held-out well.
- **Knockout cells (1,603)** lie on the NTC cycle. At the gene level no clear arrest is seen. At the guide level the four guides of each gene differ, and the guide ranking reproduces in disjoint sets of cells. See [docs/EVALUATION.md](docs/EVALUATION.md).

## Limitations to keep in mind

- No orthogonal ground truth for the phase; DNA content is itself an input of the score.
- Cells need all three markers and a reliable DNA measurement to get the full (stain) score; the others only get the less precise phase-contrast score in the knockout analysis.
- The knockout analysis compares targeting guides with non-targeting controls; a general effect of Cas9 cutting is not separated from gene-specific effects (no cutting control).
- Per-guide estimates rest on 36 to 185 scored cells each.

## License and data

Code is released under the MIT license ([LICENSE](LICENSE)). The input data come from the Biohub OPS0077 screen (public S3 bucket `ops-explorer-public`); please cite and follow the terms of the original data source.
