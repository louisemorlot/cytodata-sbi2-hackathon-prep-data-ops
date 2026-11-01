# Notebook guide

`01_download_data_started.ipynb` is the participant-facing starter. It lets participants select Altair phase profiles per guide and per guide-by-screen, Altair Biohub_OPS0077 4i, POSH, and Aconcagua; previews sizes; downloads when explicitly enabled; and runs first-pass QC. Install dependencies using `python -m pip install -r requirements.txt` from the repository root before opening it.

`00_download_preview_and_qc.ipynb` is the administrative preview/QC entry point. For unattended VM preparation, use the documented scripts instead of an interactive notebook.

The default download also extracts the 10 `*_guide_x_screen.h5ad` profiles associated with the Biohub_OPS0077 / 4i screen from Altair's marker-profile ZIP into `data/altair/ops0077_4i/`.

`02_data_inventory_and_overlap.ipynb` inventories the repository's staged files and checks OPS0077 guide overlap against both Altair phase summaries, plus candidate gene-level overlap across the available profile datasets. It is safe to run with partial downloads.

The Aconcagua guide/gene Parquet tables are supplied separately from the HTTP downloader. Notebook 01 can verify and stage them from `aconcagua_guide_gene_features/` after the folder has been transferred to the participant VM. Notebook 09 checks their schema, guide-barcode overlap across cell-stage profiles, and target-gene overlap with the other staged profiles.
