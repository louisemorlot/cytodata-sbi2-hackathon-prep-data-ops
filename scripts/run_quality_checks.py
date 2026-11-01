#!/usr/bin/env python3
"""Check staged file integrity, schemas, metadata coverage, and tentative gene overlap."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import zipfile

GENE_FIELDS = ("gene_symbol", "gene_symbol_0", "target_gene", "gene_name", "target", "gene")
GUIDE_FIELDS = ("barcode", "cell_barcode_0", "guide_id", "guide_sequence", "sgrna", "sgRNA", "protospacer_sequence")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def normalize(values) -> set[str]:
    out = set()
    for value in values:
        if value is None:
            continue
        text = str(value).strip().upper()
        if text and text not in {"NAN", "NONE", "NULL"}:
            out.add(text)
    return out


def best_gene_field(frame, known_genes: set[str] | None = None):
    candidates = [c for c in frame.columns if str(c).lower() in {x.lower() for x in GENE_FIELDS}]
    if known_genes:
        candidates += [c for c in frame.columns if str(c).lower() in {"aggregate_id", "gene_id"}]
        scored = []
        for col in dict.fromkeys(candidates):
            values = normalize(frame[col].dropna().astype(str).head(50000))
            scored.append((len(values & known_genes), col))
        if scored:
            score, col = max(scored)
            if score:
                return col, score
    return (candidates[0], None) if candidates else (None, None)


def inspect_h5ad(path: Path, known_genes: set[str] | None = None) -> dict:
    import anndata as ad
    data = ad.read_h5ad(path, backed="r")
    obs = data.obs
    gene_field, matched = best_gene_field(obs, known_genes)
    genes = normalize(obs[gene_field].dropna()) if gene_field else set()
    screen_cols = [c for c in obs.columns if "screen" in str(c).lower() or "experiment" in str(c).lower()]
    screen_values = {}
    for col in screen_cols:
        values = obs[col].dropna().astype(str).unique()[:50]
        screen_values[str(col)] = list(map(str, values))
    result = {"format": "H5AD", "observations": int(data.n_obs), "features": int(data.n_vars),
              "obs_columns": list(map(str, obs.columns)), "var_columns": list(map(str, data.var.columns)),
              "gene_field": str(gene_field) if gene_field is not None else None,
              "unique_gene_labels": len(genes) if gene_field else None,
              "gene_labels_matching_library": matched,
              "screen_metadata_values_sample": screen_values,
              "duplicate_observation_names": int(obs.index.duplicated().sum()),
              "duplicate_feature_names": int(data.var.index.duplicated().sum())}
    data.file.close()
    if genes:
        result["_gene_symbols"] = sorted(genes)
    return result


def inspect_parquet(path: Path, known_genes: set[str] | None = None) -> dict:
    import pyarrow.parquet as pq
    parquet = pq.ParquetFile(path)
    names = parquet.schema.names
    lookup = {name.lower(): name for name in names}
    result = {"format": "Parquet", "rows": parquet.metadata.num_rows,
              "columns": len(names), "row_groups": parquet.metadata.num_row_groups,
              "column_names": names}
    gene_col = next((lookup[x] for x in GENE_FIELDS if x in lookup), None)
    guide_col = next((lookup[x.lower()] for x in GUIDE_FIELDS if x.lower() in lookup), None)
    result["gene_field"] = gene_col
    result["guide_field"] = guide_col
    if gene_col or guide_col:
        columns = list(dict.fromkeys(c for c in (gene_col, guide_col, lookup.get("cell_count")) if c))
        frame = parquet.read(columns=columns).to_pandas()
        genes = normalize(frame[gene_col].dropna()) if gene_col else set()
        if gene_col:
            result["unique_gene_labels"] = len(genes)
            if known_genes:
                result["gene_labels_matching_library"] = len(genes & known_genes)
            if genes:
                result["_gene_symbols"] = sorted(genes)
        if guide_col:
            guide_values = frame[guide_col].dropna().astype(str)
            result["unique_guide_ids"] = int(guide_values.nunique())
            result["duplicate_guide_id_rows"] = int(guide_values.duplicated().sum())
            if gene_col:
                pair = frame[[guide_col, gene_col]].dropna().drop_duplicates()
                per_guide_targets = pair.groupby(guide_col)[gene_col].nunique()
                result["guides_with_multiple_target_symbols"] = int((per_guide_targets > 1).sum())
        if "cell_count" in frame:
            counts = frame["cell_count"].dropna()
            result["cell_count_min"] = int(counts.min()) if len(counts) else None
            result["cell_count_median"] = float(counts.median()) if len(counts) else None
            result["cell_count_max"] = int(counts.max()) if len(counts) else None
            result["nonpositive_cell_count_rows"] = int((counts <= 0).sum())
    metadata_names = {"cell_barcode_0", "gene_symbol_0", "cell_count", "cell_stage", "cell_stage_confidence"}
    result["numeric_feature_columns"] = sum(
        (field.name not in metadata_names) and
        (str(field.type).startswith(("int", "uint", "float", "double")))
        for field in parquet.schema_arrow
    )
    return result


def guide_target_map(path: Path) -> tuple[dict[str, str], list[str]]:
    import anndata as ad
    data = ad.read_h5ad(path, backed="r")
    obs = data.obs
    guide_field = next((c for c in ["guide", "sgRNA", "Metadata_sgRNA", "barcode", "guide_id"]
                        if c in obs.columns), None)
    target_field = next((c for c in ["target_gene", "gene", "gene_name", "Metadata_Gene"]
                         if c in obs.columns), None)
    if guide_field is None or target_field is None:
        data.file.close()
        raise ValueError(f"Missing guide or target field in {path.name}; columns={list(obs.columns)}")
    pairs = obs[[guide_field, target_field]].dropna().astype(str)
    targets_by_guide = pairs.groupby(guide_field)[target_field].agg(lambda values: sorted(set(values)))
    conflicting = sorted(guide for guide, targets in targets_by_guide.items() if len(targets) > 1)
    mapping = {guide: targets[0] for guide, targets in targets_by_guide.items()}
    data.file.close()
    return mapping, conflicting


def parquet_guide_target_map(path: Path) -> tuple[dict[str, str], list[str]]:
    import pyarrow.parquet as pq
    parquet = pq.ParquetFile(path)
    lookup = {name.lower(): name for name in parquet.schema.names}
    guide_field = lookup.get("cell_barcode_0")
    target_field = lookup.get("gene_symbol_0")
    if not guide_field or not target_field:
        raise ValueError(f"Expected cell_barcode_0 and gene_symbol_0 in {path.name}; columns={parquet.schema.names[:12]}")
    frame = parquet.read(columns=[guide_field, target_field]).to_pandas().dropna().astype(str)
    grouped = frame.drop_duplicates().groupby(guide_field)[target_field].agg(lambda values: sorted(set(values)))
    conflicting = sorted(guide for guide, targets in grouped.items() if len(targets) > 1)
    return {guide: targets[0] for guide, targets in grouped.items()}, conflicting


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default="data")
    parser.add_argument("--outdir", help="Defaults to <data-dir>/qc")
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="Report missing core files without returning failure (useful for checking a partial staging run).")
    parser.add_argument("--phase-reference-file", help="Optional single Altair phase H5AD to use for guide-overlap checks; by default checks both phase_guide.h5ad and phase_guide_x_screen.h5ad when present.")
    args = parser.parse_args()
    root = Path(args.data_dir).expanduser().resolve()
    out = Path(args.outdir).expanduser().resolve() if args.outdir else root / "qc"
    report = {"created_utc": datetime.now(timezone.utc).isoformat(), "data_dir": str(root), "files": [],
              "gene_overlap": {}, "limits": [
                  "Integrity and schema checks do not validate biological correctness or prove guide-to-target assignments.",
                  "Aconcagua is HeLa; Altair and POSH are A549. Cross-dataset agreement remains exploratory.",
                  "Aconcagua cell_barcode_0 values are not assumed to match guide sequence identifiers in other screens; cross-atlas comparisons use target symbols unless an explicit mapping is validated.",
                  "The OPS Explorer aggregate H5AD profiles are gene/perturbation-level and cannot be converted to guide-level profiles from the library map alone; the separately supplied Parquet bundle contains independently processed guide-level profiles."]}

    library = root / "aconcagua" / "metadata" / "perturbation_library.csv"
    known_genes = set()
    if library.exists():
        import pandas as pd
        lib = pd.read_csv(library, low_memory=False)
        gene_col = next((c for c in lib.columns if str(c).lower() == "gene_symbol"), None)
        if gene_col:
            known_genes = normalize(lib[gene_col].dropna())
            report["aconcagua_library"] = {"rows": len(lib), "unique_gene_symbols": len(known_genes),
                                             "columns": list(map(str, lib.columns))}

    files = sorted(p for p in root.rglob("*") if p.is_file() and ".part" not in p.name and "qc" not in p.parts)
    expected_core = [
        "altair/phase_guide.h5ad", "altair/phase_guide_x_screen.h5ad", "altair/markers_guide.zip", "altair/screen_map.csv",
        "altair/unexpressed_genes.csv", "altair/markers.csv", "altair/manifest.csv",
        "posh/morphology_124.pq", "aconcagua/interphase_aggregated_data.h5ad",
        "aconcagua/mitotic_aggregated_data.h5ad", "aconcagua/metadata/experimental_metadata.yaml",
        "aconcagua/metadata/feature_definitions.csv", "aconcagua/metadata/perturbation_library.csv",
    ]
    ops0077_profiles = [
        "b-catenin", "c-myc", "gh2ax", "nfkb", "p21", "p53", "prb", "ps6", "rb", "rsp6",
    ]
    expected_core.extend(f"altair/ops0077_4i/{marker}_guide_x_screen.h5ad" for marker in ops0077_profiles)
    report["missing_core_files"] = [relative for relative in expected_core if not (root / relative).is_file()]
    report["Biohub_OPS0077_4i"] = {"mapped_experiment": "ops0144_20260406",
                                    "expected_guide_by_screen_profiles": len(ops0077_profiles),
                                    "profiles_found": sum((root / "altair" / "ops0077_4i" / f"{marker}_guide_x_screen.h5ad").is_file()
                                                          for marker in ops0077_profiles)}
    gene_sets: dict[str, set[str]] = {}
    for path in files:
        rel = str(path.relative_to(root))
        item = {"path": rel, "exists": True, "size_bytes": path.stat().st_size, "sha256": digest(path)}
        try:
            if path.suffix == ".h5ad":
                item.update(inspect_h5ad(path, known_genes or None))
            elif path.suffix == ".pq" or path.suffix == ".parquet":
                item.update(inspect_parquet(path, known_genes or None))
            elif path.suffix.lower() == ".zip":
                with zipfile.ZipFile(path) as archive:
                    members = [m for m in archive.infolist() if not m.is_dir()]
                    item.update({"format": "ZIP", "members": len(members),
                                 "h5ad_members": sum(m.filename.lower().endswith(".h5ad") for m in members),
                                 "uncompressed_bytes": sum(m.file_size for m in members),
                                 "crc_test": "not run (use unzip -t for full archive check)"})
            elif path.suffix.lower() == ".csv":
                import pandas as pd
                head = pd.read_csv(path, nrows=5)
                item.update({"format": "CSV", "columns": list(map(str, head.columns))})
            elif path.suffix.lower() in {".yaml", ".yml"}:
                item["format"] = "YAML"
            if "_gene_symbols" in item:
                gene_sets[rel] = set(item.pop("_gene_symbols"))
            if rel.startswith("altair/ops0077_4i/") and path.suffix == ".h5ad":
                item["expected_experiment_id"] = "ops0144_20260406"
                item["biohub_dataset_id"] = "Biohub_OPS0077"
                screen_columns = [c for c in item.get("obs_columns", [])
                                  if "screen" in str(c).lower() or "experiment" in str(c).lower()]
                item["screen_metadata_columns"] = screen_columns
                screen_evidence = json.dumps(item.get("screen_metadata_values_sample", {})).lower()
                item["experiment_id_in_screen_metadata"] = ("ops0144_20260406" in screen_evidence or
                                                             "biohub_ops0077" in screen_evidence)
            if item.get("size_bytes", 0) == 0:
                item["error"] = "empty file"
        except Exception as exc:
            item["error"] = f"{type(exc).__name__}: {exc}"
        report["files"].append(item)

    labels = list(gene_sets)
    for i, left in enumerate(labels):
        for right in labels[i + 1:]:
            shared = sorted(gene_sets[left] & gene_sets[right])
            report["gene_overlap"][f"{left} vs {right}"] = {"shared_count": len(shared), "examples": shared[:50]}

    acon_guide_paths = sorted((root / "aconcagua_guide_gene_features").glob("aconcagua_guide_level_*.parquet"))
    if acon_guide_paths:
        population_maps = {}
        for path in acon_guide_paths:
            population_maps[path.name], conflicts = parquet_guide_target_map(path)
            report.setdefault("Aconcagua_guide_level_profiles", {})[path.name] = {
                "unique_guides": len(population_maps[path.name]),
                "guides_with_multiple_target_symbols": len(conflicts),
            }
        population_overlap = {}
        labels = list(population_maps)
        for i, left in enumerate(labels):
            for right in labels[i + 1:]:
                left_map, right_map = population_maps[left], population_maps[right]
                shared = set(left_map) & set(right_map)
                mismatches = sorted(g for g in shared if left_map[g].upper() != right_map[g].upper())
                population_overlap[f"{left} vs {right}"] = {
                    "shared_guide_ids": len(shared), "left_only": len(set(left_map) - set(right_map)),
                    "right_only": len(set(right_map) - set(left_map)),
                    "shared_guides_with_target_mismatch": len(mismatches),
                    "mismatch_examples": mismatches[:10],
                }
        report["Aconcagua_guide_population_overlap"] = population_overlap

    phase_paths = ([Path(args.phase_reference_file).expanduser().resolve()] if args.phase_reference_file else
                   [root / "altair" / "phase_guide.h5ad", root / "altair" / "phase_guide_x_screen.h5ad"])
    ops0077_paths = sorted((root / "altair" / "ops0077_4i").glob("*_guide_x_screen.h5ad"))
    if ops0077_paths:
        for phase_path in phase_paths:
            if not phase_path.is_file():
                continue
            phase_map, phase_conflicts = guide_target_map(phase_path)
            profile_overlaps = {}
            for path in ops0077_paths:
                target_map, target_conflicts = guide_target_map(path)
                shared = set(phase_map) & set(target_map)
                mismatches = sorted(guide for guide in shared if phase_map[guide].upper() != target_map[guide].upper())
                profile_overlaps[path.name] = {
                    "phase_unique_guides": len(phase_map), "4i_unique_guides": len(target_map),
                    "shared_guides": len(shared), "phase_only_guides": len(set(phase_map) - set(target_map)),
                    "4i_only_guides": len(set(target_map) - set(phase_map)),
                    "shared_guides_with_target_mismatch": len(mismatches),
                    "phase_guides_with_internal_target_conflicts": len(phase_conflicts),
                    "4i_guides_with_internal_target_conflicts": len(target_conflicts),
                    "target_mismatch_examples": mismatches[:10],
                }
            key = "guide_overlap_vs_" + phase_path.stem.replace("phase_", "altair_phase_")
            report["Biohub_OPS0077_4i"][key] = profile_overlaps
        report["Biohub_OPS0077_4i"]["phase_reference_files"] = [str(p) for p in phase_paths if p.is_file()]

    out.mkdir(parents=True, exist_ok=True)
    report_path = out / "quality_report.json"
    report_path.write_text(json.dumps(report, indent=2, default=str) + "\n")
    total = sum(x["size_bytes"] for x in report["files"])
    errors = sum("error" in x for x in report["files"])
    print(json.dumps({"report": str(report_path), "files_found": len(report["files"]),
                      "total_bytes": total, "inspection_errors": errors,
                      "missing_core_files": report["missing_core_files"],
                      "gene_overlap_pairs": len(report["gene_overlap"])}, indent=2))
    return 1 if errors or (report["missing_core_files"] and not args.allow_incomplete) else 0


if __name__ == "__main__":
    raise SystemExit(main())
