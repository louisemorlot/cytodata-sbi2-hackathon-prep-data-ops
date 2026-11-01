#!/usr/bin/env python3
"""Preview or download the selected processed SBI2 hackathon datasets."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import requests
from extract_altair_ops0077_4i import BIOHUB_DATASET_ID, EXPERIMENT_ID, extract_profiles

ALT_BASE = "https://exampledata.scverse.org/sbi2-cytodata/ops"
ACON_BASE = "https://ops-explorer-public.s3.amazonaws.com/aconcagua/aconcagua"
ZENODO_BASE = "https://zenodo.org/api/records/17470307/files"
USER_AGENT = "sbi2-hackathon-prep-data-ops/1.0 (+public dataset staging)"
ALT_SHA256 = {
    "phase_guide.h5ad": "8cffb86078d6ed9839821b9176c4ef744f1b9b35386b7710794e607188d7b985",
    "phase_guide_x_screen.h5ad": "3f1ae51ff5958b62c16c539160ae0cc5f51c91b5445659cfdf1dc6142f36e976",
    "screen_map.csv": "9dbd0486fbb685403ec26eff87dc12bb06e7ac773a2a0428ad3ff3482ad3d509",
    "unexpressed_genes.csv": "7ca185b0bc4b6f4c0f50f517846e9bd783ca0f61ad19fde04ac6ea4367b16acc",
    "markers_guide.zip": "8a7257138d26c2dfec7414c17ef9f3fee27f64e6b1f7178a957ea78c24311930",
}
ALT_SIZE_HINTS = {"phase_guide.h5ad": 47_992_569,
                  "phase_guide_x_screen.h5ad": 3_132_434_691,
                  "screen_map.csv": 8_346, "unexpressed_genes.csv": 1_585,
                  "markers.csv": 6_249, "manifest.csv": 7_872,
                  "markers_guide.zip": 2_037_832_363}
POSH_SIZE_HINTS = {"morphology_124.pq": 4_686_984, "moa_300.pq": 31_421_848,
                   "druggable_genome_1640.pq": 54_547_608}
ACON_METADATA_SIZE_HINTS = {"experimental_metadata.yaml": 1_967,
                             "feature_definitions.csv": 168_230,
                             "perturbation_library.csv": 1_639_239}
ACON_PROFILE_SIZE_HINTS = {"interphase": 5_428_564, "mitotic": 5_423_724}


def objects(full: bool) -> list[dict]:
    items = []
    altair_groups = {
        "phase_guide.h5ad": ["altair-phase"],
        "phase_guide_x_screen.h5ad": ["altair-phase"],
        "screen_map.csv": ["altair-phase", "altair-4i"],
        "unexpressed_genes.csv": ["altair-phase"],
        "markers.csv": ["altair-4i"],
        "manifest.csv": ["altair-4i"],
        "markers_guide.zip": ["altair-4i"],
    }
    for name, groups in altair_groups.items():
        items.append({"dataset": "Altair", "name": name, "url": f"{ALT_BASE}/{name}",
                      "path": f"altair/{name}", "sha256": ALT_SHA256.get(name),
                      "size_hint_bytes": ALT_SIZE_HINTS.get(name),
                      "groups": groups,
                      "description": "One phase profile per guide" if name == "phase_guide.h5ad" else
                      "One phase profile per guide and screen" if name == "phase_guide_x_screen.h5ad" else
                      "Guide-level and guide-by-screen fluorescence profiles" if name == "markers_guide.zip" else
                      "Altair release metadata"})
    for name in ["morphology_124.pq"] + (["moa_300.pq", "druggable_genome_1640.pq"] if full else []):
        items.append({"dataset": "CellPaint-POSH", "name": name,
                      "url": f"{ZENODO_BASE}/{quote(name)}?download=1", "path": f"posh/{name}",
                      "sha256": None, "groups": ["posh"],
                      "size_hint_bytes": POSH_SIZE_HINTS.get(name),
                      "description": "Published processed POSH Parquet profile"})
    for name in ["interphase", "mitotic"]:
        items.append({"dataset": "Aconcagua", "name": f"{name}_aggregated_data.h5ad",
                      "url": f"{ACON_BASE}/visualizations/{name}/aggregated_data.h5ad",
                      "path": f"aconcagua/{name}_aggregated_data.h5ad", "sha256": None,
                      "groups": ["aconcagua"],
                      "size_hint_bytes": ACON_PROFILE_SIZE_HINTS[name],
                      "description": f"Gene/perturbation aggregate H5AD; {name} population"})
    for name in ["experimental_metadata.yaml", "feature_definitions.csv", "perturbation_library.csv"]:
        items.append({"dataset": "Aconcagua metadata", "name": name,
                      "url": f"{ACON_BASE}/metadata/{name}", "path": f"aconcagua/metadata/{name}",
                      "sha256": None, "groups": ["aconcagua"],
                      "size_hint_bytes": ACON_METADATA_SIZE_HINTS.get(name),
                      "description": "Aconcagua experiment or feature/perturbation metadata"})
    return items


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def size_at(url: str, session: requests.Session) -> int | None:
    try:
        response = session.head(url, allow_redirects=True, timeout=(15, 45))
        if response.ok and response.headers.get("Content-Length"):
            return int(response.headers["Content-Length"])
        content_range = response.headers.get("Content-Range", "")
        if "/" in content_range and content_range.rsplit("/", 1)[1].isdigit():
            return int(content_range.rsplit("/", 1)[1])
    except requests.RequestException:
        pass
    try:
        # Some CDNs omit Content-Length on HEAD but provide a total size for a byte-range GET.
        with session.get(url, headers={"Range": "bytes=0-0"}, stream=True,
                         allow_redirects=True, timeout=(15, 45)) as response:
            content_range = response.headers.get("Content-Range", "")
            if "/" in content_range and content_range.rsplit("/", 1)[1].isdigit():
                return int(content_range.rsplit("/", 1)[1])
            content_length = response.headers.get("Content-Length")
            if response.status_code == 200 and content_length and content_length.isdigit():
                return int(content_length)
    except requests.RequestException:
        pass
    return None


def download(item: dict, target: Path, session: requests.Session, expected_size: int | None, force: bool) -> dict:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and not force:
        if expected_size is not None and target.stat().st_size != expected_size:
            raise IOError(f"Existing size mismatch for {target}: {target.stat().st_size} != {expected_size}; use --force to replace")
        digest = sha256(target)
        if item["sha256"] and digest != item["sha256"]:
            raise IOError(f"Existing checksum mismatch for {target}; use --force to replace")
        return {"path": str(target), "size_bytes": target.stat().st_size, "sha256": digest, "reused": True}

    partial = target.with_name(target.name + ".part")
    partial.unlink(missing_ok=True)
    last_error = None
    for attempt in range(3):
        try:
            count = 0
            with session.get(item["url"], stream=True, allow_redirects=True, timeout=(20, 300)) as response:
                response.raise_for_status()
                header_size = response.headers.get("Content-Length")
                reported_size = int(header_size) if header_size and header_size.isdigit() else expected_size
                with partial.open("wb") as stream:
                    for block in response.iter_content(chunk_size=2 * 1024 * 1024):
                        if block:
                            stream.write(block)
                            count += len(block)
            if reported_size is not None and reported_size != count:
                raise IOError(f"Incomplete download: received {count} bytes, source reported {reported_size}")
            digest = sha256(partial)
            if item["sha256"] and digest != item["sha256"]:
                raise IOError(f"SHA-256 mismatch: expected {item['sha256']}, received {digest}")
            if expected_size is not None and count != expected_size:
                raise IOError(f"Size mismatch: received {count}, expected {expected_size}")
            partial.replace(target)
            return {"path": str(target), "size_bytes": count, "sha256": digest, "reused": False}
        except (requests.RequestException, OSError) as exc:
            last_error = exc
            if attempt < 2:
                time.sleep(2 ** attempt)
    partial.unlink(missing_ok=True)
    raise IOError(f"Could not download {item['url']}: {last_error}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preset", choices=("core", "full"), default="core")
    parser.add_argument("--only-ops0077-4i", action="store_true",
                        help="Fetch only the Altair marker archive, then extract the Biohub_OPS0077 guide-by-screen profiles.")
    parser.add_argument("--datasets", nargs="+", choices=("altair-phase", "altair-4i", "posh", "aconcagua"),
                        help="Choose one or more participant dataset groups. Default: all four groups.")
    parser.add_argument("--data-dir", default="data", help="Data root (set to the mounted EBS path on AWS)")
    parser.add_argument("--dry-run", action="store_true", help="Show planned files and source sizes without writing")
    parser.add_argument("--yes", action="store_true", help="Skip the download confirmation prompt")
    parser.add_argument("--force", action="store_true", help="Replace existing files after reviewing the plan")
    args = parser.parse_args()
    out = Path(args.data_dir).expanduser().resolve()
    selected = objects(full=args.preset == "full")
    if args.only_ops0077_4i:
        if args.datasets:
            parser.error("Do not combine --only-ops0077-4i with --datasets; use --datasets altair-4i instead.")
        selected = [item for item in selected if item["path"] == "altair/markers_guide.zip"]
        selected.extend(item for item in objects(full=False)
                        if item["path"] in {"altair/screen_map.csv", "altair/markers.csv", "altair/manifest.csv"})
    elif args.datasets:
        selected = [item for item in selected if set(item["groups"]) & set(args.datasets)]
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    print("Planned processed data only (no raw images or cell_data.parquet):")
    total, sized = 0, 0
    for item in selected:
        live_size = size_at(item["url"], session)
        size = live_size if live_size is not None else item.get("size_hint_bytes")
        item["size_bytes"] = size
        if size is not None:
            total += size
            sized += 1
            size_text = f"{size / 1024**2:,.2f} MiB"
        else:
            size_text = "size unavailable from response headers"
        size_source = "" if live_size is not None or size is None else " (release inventory size)"
        print(f"{size_text:>31}{size_source}  {item['dataset']}: {item['path']}\n{'':33}{item['description']}\n{'':33}{item['url']}")
    print(f"Known-size total: {total / 1024**3:.3f} GiB ({sized}/{len(selected)} objects sized).")
    if args.dry_run:
        return 0
    if not args.yes and input("Proceed with these downloads? [y/N] ").strip().lower() not in {"y", "yes"}:
        print("Cancelled; no files downloaded.")
        return 0

    if args.only_ops0077_4i:
        selected_groups = ["altair-4i"]
    elif args.datasets:
        selected_groups = args.datasets
    else:
        selected_groups = ["altair-phase", "altair-4i", "posh", "aconcagua"]
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "preset": "ops0077-4i-only" if args.only_ops0077_4i else args.preset,
                "datasets_selected": selected_groups,
                "data_dir": str(out), "raw_images_downloaded": False, "aconcagua_cell_data_parquet_downloaded": False,
                "files": []}
    downloaded = {}
    for item in selected:
        target = out / item["path"]
        print(f"Fetching {item['dataset']}: {item['name']}")
        result = download(item, target, session, item["size_bytes"], args.force)
        downloaded[item["path"]] = target
        manifest["files"].append({"dataset": item["dataset"], "name": item["name"], "source": item["url"],
                                  "path": str(target.relative_to(out)), "size_bytes": result["size_bytes"],
                                  "sha256": result["sha256"], "expected_sha256_checked": bool(item["sha256"]),
                                  "reused": result["reused"]})
    marker_archive = downloaded.get("altair/markers_guide.zip")
    if marker_archive:
        ops0077_dir = out / "altair" / "ops0077_4i"
        extracted = extract_profiles(marker_archive, ops0077_dir, force=args.force)
        for item in extracted:
            target = Path(item["path"])
            manifest["files"].append({"dataset": f"Altair {BIOHUB_DATASET_ID} 4i", "name": item["name"],
                                      "source": f"{ALT_BASE}/markers_guide.zip#member={item['name']}",
                                      "path": str(target.relative_to(out)), "size_bytes": item["size_bytes"],
                                      "sha256": item["sha256"], "expected_sha256_checked": False,
                                      "zip_crc_checked": item["zip_crc_checked"], "reused": item["reused"]})
        (ops0077_dir / "extraction_manifest.json").write_text(json.dumps({
            "biohub_dataset_id": BIOHUB_DATASET_ID, "altair_experiment_id": EXPERIMENT_ID,
            "profiles": [{**item, "path": str(Path(item["path"]).relative_to(out))} for item in extracted],
            "profile_count": len(extracted), "total_bytes": sum(item["size_bytes"] for item in extracted)}, indent=2) + "\n")
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "download_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
