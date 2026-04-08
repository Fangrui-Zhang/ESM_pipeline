#!/usr/bin/env python3
"""Fetch ESMFold structures for each protein sequence in a manifest."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
import requests

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def fetch_structure(sequence: str, timeout: int) -> str:
    response = requests.post(ESMFOLD_API_URL, data=sequence, timeout=timeout)
    response.raise_for_status()
    return response.text


def main(
    manifest_csv: str,
    out_dir: str,
    status_out: str,
    timeout: int,
    pause_seconds: float,
    skip_existing: bool,
) -> None:
    manifest = pd.read_csv(manifest_csv)
    required_columns = {"protein_id", "sequence"}
    missing = required_columns - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")

    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    status_rows = []
    for _, row in manifest.iterrows():
        protein_id = row["protein_id"]
        sequence = str(row["sequence"]).strip()
        output_path = output_dir / f"{protein_id}.pdb"

        if skip_existing and output_path.exists():
            status_rows.append(
                {
                    "protein_id": protein_id,
                    "bmrb_id": row.get("bmrb_id"),
                    "entity_id": row.get("entity_id"),
                    "pdb_path": str(output_path.resolve()),
                    "status": "skipped_existing",
                    "error": "",
                }
            )
            continue

        try:
            pdb_text = fetch_structure(sequence, timeout=timeout)
            output_path.write_text(pdb_text)
            status_rows.append(
                {
                    "protein_id": protein_id,
                    "bmrb_id": row.get("bmrb_id"),
                    "entity_id": row.get("entity_id"),
                    "pdb_path": str(output_path.resolve()),
                    "status": "ok",
                    "error": "",
                }
            )
        except Exception as exc:
            status_rows.append(
                {
                    "protein_id": protein_id,
                    "bmrb_id": row.get("bmrb_id"),
                    "entity_id": row.get("entity_id"),
                    "pdb_path": str(output_path.resolve()),
                    "status": "error",
                    "error": str(exc),
                }
            )

        if pause_seconds > 0:
            time.sleep(pause_seconds)

    status_df = pd.DataFrame.from_records(status_rows)
    status_path = Path(status_out)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_df.to_csv(status_path, index=False)
    ok_count = int((status_df["status"] == "ok").sum()) if not status_df.empty else 0
    print(f"Wrote {ok_count} ESMFold structures to {output_dir}")
    print(f"Wrote status table to {status_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="data/protein_manifest.csv",
        help="Protein manifest CSV with protein_id and sequence columns",
    )
    parser.add_argument(
        "--out-dir",
        default="data/esmfold_structures",
        help="Directory to store fetched PDB files",
    )
    parser.add_argument(
        "--status-out",
        default="data/esmfold_status.csv",
        help="CSV that records fetch status for each protein",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=300,
        help="HTTP timeout in seconds for each ESMFold request",
    )
    parser.add_argument(
        "--pause-seconds",
        type=float,
        default=0.0,
        help="Optional delay between API calls",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip proteins that already have a saved PDB file",
    )
    args = parser.parse_args()
    main(
        args.manifest,
        args.out_dir,
        args.status_out,
        args.timeout,
        args.pause_seconds,
        args.skip_existing,
    )
