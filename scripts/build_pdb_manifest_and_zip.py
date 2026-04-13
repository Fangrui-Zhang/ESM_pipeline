"""
build_pdb_manifest_and_zip.py

Scans all PDB files in esmfold_structures/, builds a CSV manifest with
columns: protein_id, bmrb_id, entity_id, pdb_filename, pdb_path, sequence
(sequence joined from new_manifest.csv if available), then zips the
esmfold_structures/ folder plus the manifest CSV into a single archive.
"""

import os
import re
import csv
import zipfile
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR        = Path(__file__).resolve().parent.parent          # esm_pipeline/
PDB_DIR         = BASE_DIR / "data" / "esmfold_structures"
MANIFEST_CSV    = BASE_DIR / "data" / "new_manifest.csv"
OUT_CSV         = BASE_DIR / "data" / "esmfold_pdb_manifest.csv"
OUT_ZIP         = BASE_DIR / "data" / "esmfold_structures_package.zip"


def load_sequence_lookup(manifest_path: Path) -> dict:
    """Return {protein_id: sequence} from new_manifest.csv."""
    lookup = {}
    if not manifest_path.exists():
        print(f"[WARN] Manifest not found: {manifest_path}")
        return lookup
    with open(manifest_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lookup[row["protein_id"].strip()] = row["sequence"].strip()
    print(f"[INFO] Loaded {len(lookup):,} sequences from {manifest_path.name}")
    return lookup


def collect_pdb_records(pdb_dir: Path, seq_lookup: dict) -> list[dict]:
    """Return a list of dicts, one per PDB file."""
    records = []
    pattern = re.compile(r"^bmr(\d+)_entity(\d+)\.pdb$", re.IGNORECASE)

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    print(f"[INFO] Found {len(pdb_files):,} PDB files in {pdb_dir}")

    for pdb_path in pdb_files:
        fname = pdb_path.name
        m = pattern.match(fname)
        if not m:
            print(f"[WARN] Skipping unexpected filename: {fname}")
            continue

        bmrb_id   = m.group(1)
        entity_id = m.group(2)
        protein_id = f"bmr{bmrb_id}_entity{entity_id}"
        sequence   = seq_lookup.get(protein_id, "")

        records.append({
            "protein_id":        protein_id,
            "bmrb_id":           bmrb_id,
            "entity_id":         entity_id,
            "pdb_filename":      fname,
            # relative to the project root (esm_pipeline/) – portable across machines
            "pdb_path_relative": f"data/esmfold_structures/{fname}",
            # path inside the zip archive – use this when loading from the zip
            "pdb_path_in_zip":   f"esmfold_structures/{fname}",
            "sequence":          sequence,
        })

    return records


def write_csv(records: list[dict], out_path: Path) -> None:
    fieldnames = ["protein_id", "bmrb_id", "entity_id",
                  "pdb_filename", "pdb_path_relative", "pdb_path_in_zip", "sequence"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    print(f"[INFO] Wrote manifest CSV → {out_path}  ({len(records):,} rows)")


def create_zip(pdb_dir: Path, manifest_csv: Path, out_zip: Path) -> None:
    """
    Zip all PDB files and the manifest CSV into a single archive.
    Structure inside the zip:
        esmfold_structures/<filename>.pdb
        esmfold_pdb_manifest.csv
    """
    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    total = len(pdb_files) + 1  # +1 for the manifest

    print(f"[INFO] Creating ZIP archive → {out_zip}")
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for i, pdb_path in enumerate(pdb_files, 1):
            arcname = f"esmfold_structures/{pdb_path.name}"
            zf.write(pdb_path, arcname=arcname)
            if i % 500 == 0 or i == len(pdb_files):
                print(f"  [{i}/{len(pdb_files)}] added PDB files …")

        # Add the manifest CSV at the top level of the zip
        zf.write(manifest_csv, arcname=manifest_csv.name)
        print(f"  [{total}/{total}] added manifest CSV")

    size_mb = out_zip.stat().st_size / (1024 * 1024)
    print(f"[INFO] ZIP created: {out_zip}  ({size_mb:.1f} MB)")


def main():
    seq_lookup = load_sequence_lookup(MANIFEST_CSV)
    records    = collect_pdb_records(PDB_DIR, seq_lookup)

    if not records:
        print("[ERROR] No PDB records found – check PDB_DIR path.")
        return

    write_csv(records, OUT_CSV)
    create_zip(PDB_DIR, OUT_CSV, OUT_ZIP)
    print("\n✅  Done!")
    print(f"   Manifest CSV : {OUT_CSV}")
    print(f"   ZIP archive  : {OUT_ZIP}")


if __name__ == "__main__":
    main()
