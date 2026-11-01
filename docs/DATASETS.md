# Dataset inventory, contents, and sizes

Object sizes below are source object sizes observed on 2026-10-29 where available. The downloader queries current HTTP/S3 sizes at run time and records actual downloaded sizes and SHA-256 values in `data/download_manifest.json`. Some small Altair CSV endpoints do not return a `Content-Length` header; those sizes are determined from the completed download.

## Default `core` bundle

Altair uses A549 cells with inducible Cas9 and 88 screens. The library targets 1,000 genes with four guides each, plus 211 non-targeting guides. Phase and marker measurements within a screen come from the same cells; different screens use different cells. `phase_guide.h5ad` aggregates screen medians per guide, while `phase_guide_x_screen.h5ad` retains a separate profile per guide and screen.

| Dataset | File(s) | Contents / unit | Source size | Download destination |
|---|---|---|---:|---|
| Altair | `phase_guide.h5ad` | A549 phase-imaging profile; one profile per guide (median of per-screen medians) | 47,992,569 B (45.8 MiB) | `altair/` |
| Altair | `phase_guide_x_screen.h5ad` | A549 phase-imaging profile; one profile per guide and screen, median over cells | 3,132,434,691 B (2.92 GiB) | `altair/` |
| Altair | `markers_guide.zip` | Guide-level fluorescence marker profiles and guide-by-screen profiles for marker assays | 2,037,832,363 B (1.90 GiB) | `altair/` |
| Altair Biohub_OPS0077 / 4i | Ten `*_guide_x_screen.h5ad` files for b-catenin, c-myc, γH2AX, NFκB, p21, p53, pRB, pS6, RB, and RPS6 | Screen-specific guide-by-screen profiles for experiment `ops0144_20260406`; extracted from the marker ZIP | 96,314,183 B (91.85 MiB) uncompressed; no additional network bytes beyond the ZIP | `altair/ops0077_4i/` |
| Altair metadata | `screen_map.csv`, `unexpressed_genes.csv`, `markers.csv`, `manifest.csv` | Screen IDs, expression context, marker inventory, prepared-release file manifest | 24,052 B combined (about 23.5 KiB) | `altair/` |
| POSH | `morphology_124.pq` | Published processed 124-gene CRISPR KO morphology profile table; no imagery | 4,686,984 B (4.47 MiB) | `posh/` |
| Aconcagua | `interphase/aggregated_data.h5ad` | HeLa aggregate feature profiles, interphase population | 5,428,564 B (5.18 MiB) | `aconcagua/interphase_aggregated_data.h5ad` |
| Aconcagua | `mitotic/aggregated_data.h5ad` | HeLa aggregate feature profiles, mitotic population | 5,423,724 B (5.17 MiB) | `aconcagua/mitotic_aggregated_data.h5ad` |
| Aconcagua metadata | `experimental_metadata.yaml`, `feature_definitions.csv`, `perturbation_library.csv` | Experiment description, feature names/definitions, guide and target library | 1,809,436 B combined (1.73 MiB) | `aconcagua/metadata/` |

The Aconcagua perturbation library contains 20,445 records and 5,299 distinct gene symbols in the copy already inspected for this project. It includes both targeting and control guides.

## Supplied Aconcagua guide/gene feature bundle

This additional processed bundle is supplied in the sibling folder `aconcagua_guide_gene_features/`; no public download URL was provided. `scripts/stage_aconcagua_guide_data.py` verifies the included `MD5SUMS.txt` values and copies the six Parquet files into `<data-dir>/aconcagua_guide_gene_features/`. A staging manifest records file sizes, MD5, and SHA-256 checksums. The source includes README processing notes and checksums.

| File | Rows | Columns | Size | Profile level / population |
|---|---:|---:|---:|---|
| `aconcagua_guide_level_all_cells.parquet` | 20,444 | 1,638 | 173,524,236 B (165.5 MiB) | Guide / all cells |
| `aconcagua_guide_level_interphase_cells.parquet` | 20,444 | 1,638 | 173,608,922 B (165.5 MiB) | Guide / interphase |
| `aconcagua_guide_level_mitotic_cells.parquet` | 20,406 | 1,646 | 172,071,372 B (164.1 MiB) | Guide / mitotic |
| `aconcagua_gene_level_all_cells.parquet` | 5,299 | 1,637 | 45,393,979 B (43.3 MiB) | Gene / all cells |
| `aconcagua_gene_level_interphase_cells.parquet` | 5,299 | 1,637 | 45,395,758 B (43.3 MiB) | Gene / interphase |
| `aconcagua_gene_level_mitotic_cells.parquet` | 5,299 | 1,645 | 44,515,046 B (42.5 MiB) | Gene / mitotic |
| **Total** | — | — | **654,509,313 B (624.2 MiB)** | Six Parquet files |

Guide tables identify perturbations with `cell_barcode_0` and targets with `gene_symbol_0`; gene-level tables use `gene_symbol_0`. `cell_count`, `cell_stage`, and `cell_stage_confidence` are metadata, not morphology measurements. Remaining numeric columns are robust-z feature values from DAPI, tubulin, γH2AX, and phalloidin assays. The source README reports brieflow commit `2e41c7f`, 5,299 genes, 20,445 designed guides, and 250 non-targeting controls. Mitotic missing-value filtering leaves eight more features than all-cell/interphase profiles.

The supplied guide bundle supports guide-level comparisons among Aconcagua's own all-cell/interphase/mitotic populations and gene-level overlap with other datasets. In the supplied tables, all-cell and interphase each have 20,444 guide barcodes; mitotic has 20,406, all shared with the other two profiles, with no shared-barcode target-symbol conflicts. The current staged Altair OPS0077 profiles share 809 target symbols with the Aconcagua tables; the notebooks recompute this from the files. Aconcagua uses barcode identifiers; do not join them directly to Altair guide sequences without a validated mapping. Cross-atlas target-symbol comparisons are possible, with the A549 versus HeLa and assay differences kept in mind. These Parquet files complement the separate OPS Explorer aggregate H5ADs; they are not derived by expanding the aggregate H5ADs.

Preview and stage the files with:

```bash
python scripts/stage_aconcagua_guide_data.py --source-dir ../aconcagua_guide_gene_features --data-dir /mnt/hackathon-data --dry-run
python scripts/stage_aconcagua_guide_data.py --source-dir ../aconcagua_guide_gene_features --data-dir /mnt/hackathon-data
```

The extra bundle is not included in the HTTP `core` / `full` downloader presets because it is currently distributed as a local source folder. Keep its source folder available to the VM operator or transfer it separately before staging.

Known-size core network transfer is about **5,235,632,383 bytes (5.24 GB / 4.88 GiB)**, including `phase_guide_x_screen.h5ad`. Including the extracted Biohub_OPS0077 profiles, local storage is about **5,331,946,566 bytes (5.33 GB / 4.97 GiB)**. The 4i subset is already a member of the marker ZIP, so its extraction does not increase the download total. Exact sizes vary if source files are updated.

## Participant download selections

`notebooks/01_download_data_started.ipynb` is the participant-facing selector. The same choices are available at the command line:

```bash
python scripts/download_data.py --datasets altair-phase --dry-run
python scripts/download_data.py --datasets altair-4i --dry-run
python scripts/download_data.py --datasets posh aconcagua --dry-run
python scripts/download_data.py --datasets altair-phase altair-4i posh aconcagua --dry-run
```

Remove `--dry-run` after reviewing the preview to download. The `--datasets` names map to the four groups listed above. `--preset full` adds POSH's MoA-300 and Druggable Genome-1640 profile tables when POSH is selected.

## Optional files in `full`

| Dataset | File | Contents / unit | Source size |
|---|---|---|---:|
| POSH | `moa_300.pq` | Processed 300-gene MoA screen profile table | 31,421,848 B (29.97 MiB) |
| POSH | `druggable_genome_1640.pq` | Processed 1,640-gene screen profile table | 54,547,608 B (52.02 MiB) |

## Source locations

- **Altair prepared profiles and metadata:** [Cytodata hackathon prep OPS release](https://exampledata.scverse.org/sbi2-cytodata/ops/) and its schema/guide at [cytodata-hackathon-prep-ops](https://github.com/timtreis/cytodata-hackathon-prep-ops). The marker ZIP is optional in analyses that do not need fluorescence; it is included in `core` for the phase-versus-fluorescence workstream.
- **Biohub_OPS0077 mapping:** Altair's `screen_map.csv` maps `Biohub_OPS0077` to `ops0144_20260406` with `match=marker_identity`. Altair's `markers.csv` lists 10 marker assays for that same experiment (b-catenin, c-myc, γH2AX, NFκB, p21, p53, pRB, pS6, RB, RPS6), and the `*_guide_x_screen.h5ad` members in `markers_guide.zip` provide the small profile-level data used here. The inventory lists sampled cell files separately; those files are large (about 1.2–1.6 GB each) and are not part of this guide-profile subset.

For an OPS0077-only staging run, `python scripts/download_data.py --only-ops0077-4i --dry-run` previews the archive download needed to extract that subset. The individual H5AD members do not have standalone URLs in the Altair release.

See [the measured Biohub_OPS0077 checks](OPS0077_QC.md) for the successful download and guide-overlap results.
- **POSH:** [Zenodo record 17470307](https://doi.org/10.5281/zenodo.17470307) and [insitro/cp-posh](https://github.com/insitro/cp-posh). The Zenodo record provides the processed Parquet tables and is a more stable source for this repo than an S3 aggregate key that currently returned 404. Do not confuse these small processed profiles with the much larger POSH single-cell/tile imagery in the S3 bucket.
- **Aconcagua:** public bucket prefix `s3://ops-explorer-public/aconcagua/aconcagua/`; exact profile keys are under `visualizations/{interphase,mitotic}/aggregated_data.h5ad`; metadata is under `metadata/`. The [OPS Explorer collection](https://biohub.ai/ops-explorer?collection=3b32c693-a8b6-46d2-8f42-d4511a07f2f7) also describes the collection. The portal is [screens.wi.mit.edu/aconcagua](https://screens.wi.mit.edu/aconcagua/).

## Contents and analysis limits

- Altair `phase_guide.h5ad` is per guide; `phase_guide_x_screen.h5ad` is per guide and screen. Both carry guide/target annotations and phase feature matrices. Marker archive contents include profile-level tables; it is not an image archive.
- The Biohub_OPS0077 subset uses screen-specific guide-by-screen summaries. It is suitable for profile-level comparisons against Altair phase guides for the same perturbations. The ten corresponding `*_cells_n100.h5ad` files hold sampled cell-level measurements and are much larger; this repo does not download them by default.
- POSH Zenodo profiles are processed numeric profiles. Inspect the Parquet schema and library identifiers with `run_quality_checks.py` before treating their rows as guide-level or assuming their target column is a gene symbol.
- Aconcagua `aggregated_data.h5ad` is a gene/perturbation aggregate representation, not one row per guide. `perturbation_library.csv` maps target symbols to guide/barcode records, but this mapping cannot split an already gene-aggregated profile back into guide-specific profiles.
- The Aconcagua portal exposes an `all` population in its interactive interface. No `all/aggregated_data.h5ad` was present at the public OPS Explorer prefix when checked. This repo therefore downloads both available public split profiles (interphase and mitotic); it does not construct an `all` profile by averaging them.
- Aconcagua raw cell data (`cell_data.parquet`, observed at 323,471,309,291 B) and image/Zarr content are intentionally excluded from `core` and `full`.
- POSH is CC BY-NC-SA 4.0 according to its source materials. Confirm the applicable terms for the hackathon and any later data redistribution. Keep data files out of GitHub.

### Previously staged POSH preflight files

The parent workspace contains `well_normalized_profiles.pq` (about 1.5 GB) and `124_gene_library.csv`. Their original download URL/source manifest is not recorded in the files available here, so the new downloader does not silently claim to reproduce them. Direct requests to the previously documented POSH S3 profile keys returned 404 during preparation. This repository uses the smaller, source-identified POSH Parquet files in the Zenodo record. If the 1.5 GB local profile is required for a particular analysis, transfer it separately and record its original source and checksum in the VM's `download_manifest.json` before distribution.
