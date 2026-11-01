# SBI2 Hackathon: processed imaging data starter

This repository stages processed phenotypic profiles for the SBI2 hackathon. It includes guide-level phase and fluorescence profiles from Altair, processed morphology profiles from POSH, and both aggregate and guide-level profiles from Aconcagua. The data are intended for participant analysis on a shared AWS VM; the scripts download only the selected processed tables and metadata, not raw microscopy images.

The starter flow is modeled on the [Cytodata hackathon prep OPS repository](https://github.com/timtreis/cytodata-hackathon-prep-ops): preview what will be fetched, download selected data, run first-pass checks, then explore profile overlap before designing an analysis.

## Start here

1. Install the Python dependencies from the repository root: `python -m pip install -r requirements.txt`.
2. Open [`notebooks/01_download_data_started.ipynb`](notebooks/01_download_data_started.ipynb). Select the datasets needed for your project, inspect the size preview, and enable the download cell when ready. In a terminal, the same downloader supports `python scripts/download_data.py --datasets altair-phase altair-4i posh aconcagua --dry-run` followed by the same command without `--dry-run`.
3. In notebook 01, stage the supplied Aconcagua guide/gene Parquet bundle from `aconcagua_guide_gene_features/` after reviewing its checksum-verified preview. On a VM, transfer that folder first and set `ACONCAGUA_GUIDE_SOURCE_DIR` if needed.
4. Run [`notebooks/02_data_inventory_and_overlap.ipynb`](notebooks/09_preflight_dataset_inventory_and_overlap.ipynb) to inventory staged files, inspect schemas, compare Altair phase and 4i guide identifiers, compare Aconcagua guide IDs across cell-stage profiles, and measure candidate gene-symbol overlap. It also works with partial data.
5. For a compact administrative integrity/schema report, run `python scripts/run_quality_checks.py --data-dir data --allow-incomplete`. The QC report, remote download manifest, and local staging manifest help the VM operator confirm what was staged.

All data files are excluded from Git. Keep the code, notebooks, and descriptive metadata in the public repository; stage the downloaded data on the VM's mounted storage. For VM setup, see [`docs/AWS_VM.md`](docs/AWS_VM.md). For exact object inventory and source byte counts, see [`docs/DATASETS.md`](docs/DATASETS.md).

## Dataset overview

The default `core` download is approximately **5,235,632,383 bytes (5.24 GB / 4.88 GiB) transferred**. The separately supplied six-file Aconcagua guide/gene bundle adds **654,509,313 bytes (624.2 MiB / 0.61 GiB)** when staged locally. The core download includes both the per-guide phase profile and the much larger per-guide-per-screen phase profile. Extracting the 10 OPS0077 guide-by-screen H5AD profiles adds about **96.3 MB (91.9 MiB)** from the already-downloaded marker archive, for approximately **5,986,455,879 bytes (5.99 GB / 5.58 GiB) of total staged content** including the new local bundle. These are measured source sizes; the downloader checks live sizes where possible and records actual downloaded sizes and SHA-256 checksums in its manifest.

| Dataset | Screen and experimental context | Downloaded processed data | Approximate size |
|---|---|---|---:|
| **Altair phase** | A549 optical pooled CRISPR perturbation screen. Includes one phase profile per guide (median of screen medians) and one per guide and screen (median over cells). | `phase_guide.h5ad`, `phase_guide_x_screen.h5ad`, screen map, unexpressed-gene metadata | 48.0 MB + 3.13 GB (2.92 GiB) + metadata |
| **Altair Biohub_OPS0077 4i** | Altair's Biohub_OPS0077 screen maps to experiment `ops0144_20260406`. The four-dimensional imaging assay includes ten antibody marker readouts: β-catenin, c-Myc, γH2AX, NFκB, p21, p53, pRB, pS6, RB, and RPS6. | Marker archive plus ten extracted guide-by-screen H5AD profiles and marker/manifest metadata | 2.04 GB (1.90 GiB) transfer; extracted profiles 96.3 MB (91.9 MiB) |
| **POSH morphology-124** | A549 CRISPR imaging screen with processed cell-painting morphology profiles for a 124-gene subset. | Zenodo `morphology_124.pq` Parquet profile table | 4.69 MB (4.47 MiB) |
| **Aconcagua aggregate profiles** | HeLa optical pooled CRISPR screen. Public aggregate profiles are provided separately for interphase and mitotic populations; experiment metadata, feature definitions, and the perturbation library support interpretation. | Two aggregate H5AD files plus three metadata files | About 12.7 MB total |
| **Aconcagua guide/gene feature bundle** | Six brieflow-processed HeLa morphology profiles: guide-level (20,444 all-cell/interphase; 20,406 mitotic) and gene-level (5,299 per population), for all-cell, interphase, and mitotic groups. Channels include DAPI, tubulin, γH2AX, and phalloidin. | Six Parquet tables; guide rows use `cell_barcode_0` and `gene_symbol_0`; feature values are robust z-scores | 654.5 MB (624.2 MiB), copied from the supplied local bundle |

Altair is an A549 optical pooled CRISPR screen of 1,000 knockout targets, each with four guides, plus 211 non-targeting guides, across 88 separate screens. Phase and fluorescence measurements within a screen come from the same cells; cells are not shared across screens. Altair's marker archive is large because it contains marker-level profiles and sampled-cell files. The additional `phase_guide_x_screen.h5ad` is about 3.13 GB and has one observation per guide and screen, enabling screen-to-screen variation checks and avoiding aggregation across screens. Only the ten guide-by-screen profiles for the selected Biohub_OPS0077 experiment are extracted into `data/altair/ops0077_4i/`; no extra network transfer is needed for those profiles. Individual sampled-cell files are not part of this starter bundle. Aconcagua's raw `cell_data.parquet` (about 323 GB in the inspected public inventory) and raw/image data are excluded. POSH's optional MoA-300 and Druggable Genome-1640 profile tables can be added with `--preset full`; see the inventory for sizes.

## What can be compared?

- **Altair phase versus Altair 4i:** this is the clearest guide-level integration opportunity. The prepared phase table has 4,211 distinct guides; the ten OPS0077 profiles each have 4,208. The staged quality report found 4,208 shared guides, three phase-only guides, no 4i-only guides, and no target mismatches among shared guides. Compare phenotypes by shared guide (or target after checking guide consistency). Measurements are from different assay modalities and feature spaces, so do not concatenate their columns as if they were the same measurements. The shared guide index supports cross-modal comparisons and multimodal modeling.
- **Altair and POSH:** both use A549 cells, which makes shared-target phenotype comparisons a reasonable exploratory path. POSH's library is smaller and its Parquet row semantics and identifiers should be verified before joining. Use the preflight notebook to inspect candidate fields and quantify symbol overlap; compare at shared targets only after mapping names and aggregation levels.
- **Altair/POSH and Aconcagua:** the supplied Aconcagua Parquet bundle now provides guide-level as well as gene-level profiles. Cross-atlas conservation can be explored at shared target genes, but Aconcagua uses HeLa while Altair and POSH use A549, so agreement remains a cross-cell-line, cross-assay hypothesis. Its `cell_barcode_0` guide barcode must not be equated with Altair guide sequence identifiers without an explicit library mapping. The separate OPS Explorer aggregate H5ADs remain gene/perturbation-level and cannot be disaggregated into guides.
- **Feature overlap:** feature names and units are assay-specific. A matching feature label alone does not guarantee equivalent measurements. Normalize and compare within each dataset first, then compare target-level effects or representation similarity with assay-aware methods.

On the supplied Aconcagua guide bundle, the current preflight finds 20,444 shared guide barcodes between all-cell and interphase, and 20,406 shared with mitotic (38 absent from mitotic); no target-symbol conflicts were found among shared barcodes. Against the currently staged Altair OPS0077 profiles, the preflight finds 809 shared target symbols between Aconcagua and Altair. These counts are recomputed from the staged files. The exact overlaps depend on the particular staged files and identifier columns. `02_data_inventory_and_overlap.ipynb` reports measured overlaps from the files available in `SBI2_DATA_DIR` (or `data/` by default); it does not infer missing values or download files.

## File layout

```text
data/
  altair/
    phase_guide.h5ad
    phase_guide_x_screen.h5ad
    screen_map.csv
    unexpressed_genes.csv
    markers.csv
    manifest.csv
    markers_guide.zip
    ops0077_4i/                 # ten extracted guide-by-screen H5ADs
  posh/
    morphology_124.pq
  aconcagua_guide_gene_features/   # six locally supplied Parquet profiles
  aconcagua/
    interphase_aggregated_data.h5ad
    mitotic_aggregated_data.h5ad
    metadata/
      experimental_metadata.yaml
      feature_definitions.csv
      perturbation_library.csv
  download_manifest.json
```

Remote downloads are recorded in `data/download_manifest.json`; the staged Aconcagua Parquet tables have a separate `staging_manifest.json` with MD5 and SHA-256. The six supplied MD5 values are checked before and after copying. The source and contents inventory, optional POSH files, and known limitations are described in [`docs/DATASETS.md`](docs/DATASETS.md). Initial 4i guide-overlap checks are documented in [`docs/OPS0077_QC.md`](docs/OPS0077_QC.md).

## Sources and experimental references

- **Altair prepared profiles and schema:** [Cytodata hackathon prep OPS](https://github.com/timtreis/cytodata-hackathon-prep-ops); the release base is listed in [`docs/DATASETS.md`](docs/DATASETS.md). Altair's `screen_map.csv` links the human-readable screen label to the experiment identifier used by the marker profiles.
- **POSH processed profiles:** [Zenodo record 17470307](https://doi.org/10.5281/zenodo.17470307) and [insitro/cp-posh](https://github.com/insitro/cp-posh). The prepared Zenodo Parquet files are the source used by this downloader. Consult the record and repository for attribution and use conditions; the source materials state CC BY-NC-SA 4.0.
- **Aconcagua profiles and metadata:** [OPS Explorer collection](https://biohub.ai/ops-explorer?collection=3b32c693-a8b6-46d2-8f42-d4511a07f2f7), [Aconcagua screen portal](https://screens.wi.mit.edu/aconcagua/), and the public object prefix documented in [`docs/DATASETS.md`](docs/DATASETS.md). The downloader retrieves the available interphase and mitotic aggregate profiles; it does not synthesize an all-population matrix.

To stage the supplied bundle outside the notebook, use `python scripts/stage_aconcagua_guide_data.py --source-dir ../aconcagua_guide_gene_features --data-dir /mnt/hackathon-data --dry-run`, then rerun without `--dry-run` to copy it. The supplied folder describes these profiles as brieflow-processed (commit `2e41c7f`); it is not downloaded from a public URL by this project.

Before redistribution, check the source records for dataset-specific attribution and license requirements. For Aconcagua object details and the current download commands, consult [`docs/DATASETS.md`](docs/DATASETS.md); for AWS disk sizing and setup, use [`docs/AWS_VM.md`](docs/AWS_VM.md).
