#!/usr/bin/env python3
"""Stage the locally supplied Aconcagua guide/gene feature Parquet files."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

FILES = (
    "aconcagua_guide_level_all_cells.parquet",
    "aconcagua_guide_level_interphase_cells.parquet",
    "aconcagua_guide_level_mitotic_cells.parquet",
    "aconcagua_gene_level_all_cells.parquet",
    "aconcagua_gene_level_interphase_cells.parquet",
    "aconcagua_gene_level_mitotic_cells.parquet",
)
AUXILIARY = ("README.md", "MD5SUMS.txt")


def digests(path: Path) -> tuple[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return md5.hexdigest(), sha256.hexdigest()


def expected_md5(source: Path) -> dict[str, str]:
    checksum_file = source / "MD5SUMS.txt"
    result = {}
    for line in checksum_file.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            result[fields[-1].lstrip("*")] = fields[0].lower()
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", required=True, help="Folder containing the six Parquet files, README.md, and MD5SUMS.txt")
    parser.add_argument("--data-dir", default="data", help="Staged data root; files go in aconcagua_guide_gene_features/")
    parser.add_argument("--dry-run", action="store_true", help="Verify and report files without copying them")
    parser.add_argument("--yes", action="store_true", help="Skip the copy confirmation prompt")
    parser.add_argument("--force", action="store_true", help="Replace existing destination files after verifying sources")
    args = parser.parse_args()
    source = Path(args.source_dir).expanduser().resolve()
    destination = Path(args.data_dir).expanduser().resolve() / "aconcagua_guide_gene_features"
    sums = expected_md5(source)
    items = []
    total = 0
    for name in FILES:
        path = source / name
        if not path.is_file():
            parser.error(f"Missing source file: {path}")
        md5, sha256 = digests(path)
        if name not in sums:
            parser.error(f"No expected MD5 for {name} in {source / 'MD5SUMS.txt'}")
        if md5 != sums[name]:
            parser.error(f"MD5 mismatch for {path}: expected {sums[name]}, got {md5}")
        size = path.stat().st_size
        total += size
        target = destination / name
        if target.exists() and not args.force:
            if target.stat().st_size != size or digests(target)[0] != md5:
                parser.error(f"Existing destination differs: {target}; use --force to replace")
            status = "already staged and verified"
        else:
            status = "will copy" if args.dry_run else "ready to copy"
        items.append({"name": name, "source": str(path), "destination": str(target),
                      "size_bytes": size, "md5": md5, "sha256": sha256, "status": status})
        print(f"{size / 1024**2:9.2f} MiB  {status:27} {name}  MD5 {md5}")
    print(f"Total: {total:,} bytes ({total / 1024**2:.2f} MiB, {total / 1024**3:.3f} GiB)")
    if args.dry_run:
        return 0
    if not args.yes and input(f"Copy verified files to {destination}? [y/N] ").strip().lower() not in {"y", "yes"}:
        print("Cancelled; no files copied.")
        return 0
    destination.mkdir(parents=True, exist_ok=True)
    for item in items:
        source_path = Path(item["source"])
        target = Path(item["destination"])
        if target.exists() and item["status"] == "already staged and verified":
            continue
        partial = target.with_name(target.name + ".part")
        shutil.copy2(source_path, partial)
        copied_md5, copied_sha256 = digests(partial)
        if copied_md5 != item["md5"] or copied_sha256 != item["sha256"]:
            partial.unlink(missing_ok=True)
            raise IOError(f"Copied file checksum mismatch: {target.name}")
        partial.replace(target)
    for name in AUXILIARY:
        src = source / name
        if src.is_file():
            shutil.copy2(src, destination / name)
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "source_dir": str(source), "dataset": "Aconcagua guide- and gene-level morphology profiles",
                "files": [{k: v for k, v in item.items() if k != "status"} for item in items]}
    (destination / "staging_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Staged and verified {len(items)} Parquet files; manifest: {destination / 'staging_manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
