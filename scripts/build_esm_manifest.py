#!/usr/bin/env python3
"""
Build a protein manifest CSV for ESMFold fetching.

For each BMRB ID in an input CSV:
  1. Downloads the NMR-STAR .str file via makeshift
  2. Extracts all polypeptide entity sequences
  3. Writes one row per entity: protein_id, sequence, bmrb_id, entity_id

The output manifest is directly usable as input to fetch_esmfold.py.

Usage
-----
    python3 scripts/build_esm_manifest.py \
        --ids new_bmrb_ids.csv \
        --out data/new_manifest.csv \
        --str-cache data/str_cache \
        --resume

Merge with existing manifest:
    python3 scripts/build_esm_manifest.py \
        --ids new_bmrb_ids.csv \
        --out data/new_manifest.csv \
        --merge-into data/existing_manifest.csv
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

# makeshift lives in the s2_pipeline tree; add it to path
_SCRIPT_DIR = Path(__file__).resolve().parent
_S2_PIPELINE = _SCRIPT_DIR.parent.parent / "s2_pipeline"
if str(_S2_PIPELINE) not in sys.path:
    sys.path.insert(0, str(_S2_PIPELINE))

try:
    import makeshift as ms
except ImportError:
    print("ERROR: Could not import makeshift. Make sure s2_pipeline is at "
          f"{_S2_PIPELINE} or add it to PYTHONPATH.", file=sys.stderr)
    sys.exit(1)


def get_sequences_for_entry(bmrb_id: int, str_cache_dir: Path) -> list[dict]:
    """
    Download (or reuse cached) NMR-STAR file for bmrb_id and extract
    all polypeptide entity sequences.

    Returns a list of dicts with keys:
        protein_id, sequence, bmrb_id, entity_id
    """
    str_file = str_cache_dir / f"bmr{bmrb_id}_3.str"

    # Download if not cached
    if not str_file.exists():
        ms.fetch_nmrstar_file(bmrb_id, output_dir=str(str_cache_dir))

    if not str_file.exists():
        print(f"  ✗ [{bmrb_id}] Download failed — skipping")
        return []

    try:
        entry = ms.parse_nmr_star(str(str_file))
        seqs_df = ms.get_sequences(entry)
    except Exception as exc:
        print(f"  ✗ [{bmrb_id}] Parse error: {exc}")
        return []

    if seqs_df is None or len(seqs_df) == 0:
        print(f"  ✗ [{bmrb_id}] No sequences found")
        return []

    polypeptide = seqs_df[
        seqs_df["Polymer_type"].str.contains("polypeptide", case=False, na=False)
    ].copy()

    if len(polypeptide) == 0:
        print(f"  ✗ [{bmrb_id}] No polypeptide entity found")
        return []

    rows = []
    for _, row in polypeptide.iterrows():
        entity_id = int(row["ID"])
        seq = str(row["Polymer_seq_one_letter_code"]).strip().replace("\n", "").replace(" ", "").replace("\r", "")
        if not seq:
            continue
        protein_id = f"bmr{bmrb_id}_entity{entity_id}"
        rows.append({
            "protein_id": protein_id,
            "sequence": seq,
            "bmrb_id": bmrb_id,
            "entity_id": entity_id,
        })

    print(f"  ✓ [{bmrb_id}] {len(rows)} polypeptide entity/entities")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ids", required=True,
                        help="CSV file with a column named BMRB_entry or Entry_ID")
    parser.add_argument("--out", required=True,
                        help="Output manifest CSV path")
    parser.add_argument("--str-cache", default="data/str_cache",
                        help="Directory to cache downloaded .str files (default: data/str_cache)")
    parser.add_argument("--merge-into", default=None,
                        help="If provided, merge the new manifest into this existing manifest CSV "
                             "and write the merged result back to --out")
    parser.add_argument("--resume", action="store_true",
                        help="Skip BMRB IDs already present in --out")
    parser.add_argument("--pause", type=float, default=0.5,
                        help="Seconds to wait between BMRB downloads (default: 0.5)")
    args = parser.parse_args()

    # ── Load IDs ──────────────────────────────────────────────────────────────
    ids_df = pd.read_csv(args.ids)
    id_col = next((c for c in ids_df.columns if c in ("BMRB_entry", "Entry_ID")), None)
    if id_col is None:
        print(f"ERROR: CSV must have a column named 'BMRB_entry' or 'Entry_ID'. "
              f"Found: {list(ids_df.columns)}", file=sys.stderr)
        sys.exit(1)

    all_ids = ids_df[id_col].dropna().astype(int).tolist()
    print(f"Loaded {len(all_ids)} BMRB IDs from {args.ids}")

    # ── Resume: skip already-done IDs ─────────────────────────────────────────
    done_ids: set[int] = set()
    if args.resume and Path(args.out).exists():
        existing = pd.read_csv(args.out)
        done_ids = set(existing["bmrb_id"].dropna().astype(int).tolist())
        print(f"Resume: skipping {len(done_ids)} already-processed IDs")

    ids_to_process = [i for i in all_ids if i not in done_ids]
    print(f"Processing {len(ids_to_process)} IDs...\n")

    # ── Prepare cache dir ─────────────────────────────────────────────────────
    str_cache = Path(args.str_cache)
    str_cache.mkdir(parents=True, exist_ok=True)

    # ── Main loop ─────────────────────────────────────────────────────────────
    new_rows: list[dict] = []
    n_total = len(ids_to_process)
    for i, bmrb_id in enumerate(ids_to_process, 1):
        print(f"[{i}/{n_total}] BMRB {bmrb_id}", end="  ")
        rows = get_sequences_for_entry(bmrb_id, str_cache)
        new_rows.extend(rows)

        # Checkpoint: write partial results every 50 entries
        if i % 50 == 0:
            _save_checkpoint(new_rows, done_ids, args)
            print(f"  → Checkpoint saved ({len(new_rows)} rows so far)")

        if args.pause > 0:
            time.sleep(args.pause)

    # ── Final save ────────────────────────────────────────────────────────────
    _save_checkpoint(new_rows, done_ids, args)
    print(f"\nDone. {len(new_rows)} new manifest rows written to {args.out}")
    if args.merge_into:
        print(f"Merged with {args.merge_into}")


def _save_checkpoint(new_rows: list[dict], done_ids: set[int], args) -> None:
    """Save current rows (plus any previously existing rows if resuming/merging)."""
    new_df = pd.DataFrame(new_rows) if new_rows else pd.DataFrame(
        columns=["protein_id", "sequence", "bmrb_id", "entity_id"])

    out_path = Path(args.out)

    # Merge with an existing manifest if requested
    if args.merge_into and Path(args.merge_into).exists():
        base = pd.read_csv(args.merge_into)
        combined = pd.concat([base, new_df], ignore_index=True)
    elif args.resume and out_path.exists():
        # Append to what's already in --out
        existing = pd.read_csv(out_path)
        combined = pd.concat([existing, new_df], ignore_index=True)
    else:
        combined = new_df

    # Deduplicate by protein_id, keep latest
    if not combined.empty and "protein_id" in combined.columns:
        combined = combined.drop_duplicates(subset=["protein_id"], keep="last")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()
