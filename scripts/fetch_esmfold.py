#!/usr/bin/env python3
"""Fetch ESMFold structures for each protein sequence in a manifest."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
ALLOWED_SEQUENCE_CHARS = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")


def fetch_structure(sequence: str, timeout: int) -> str:
    response = requests.post(ESMFOLD_API_URL, data=sequence, timeout=timeout)
    response.raise_for_status()
    return response.text


def normalize_sequence(sequence: Any) -> str:
    return "".join(str(sequence).split()).upper()


def validate_sequence(sequence: str) -> str | None:
    if not sequence:
        return "empty sequence"
    invalid_chars = sorted({char for char in sequence if char not in ALLOWED_SEQUENCE_CHARS})
    if invalid_chars:
        return f"invalid sequence characters: {''.join(invalid_chars)}"
    return None


def clean_row_values(row: pd.Series) -> dict[str, str]:
    cleaned: dict[str, str] = {}
    for key, value in row.items():
        cleaned[key] = "" if pd.isna(value) else str(value)
    return cleaned


def resolve_output_path(
    row: pd.Series,
    manifest_path: Path,
    out_dir: str,
    pdb_path_template: str | None,
    path_root: str | None,
) -> Path:
    row_values = clean_row_values(row)
    protein_id = row_values.get("protein_id", "").strip()
    explicit_path = row_values.get("pdb_path", "").strip()

    if explicit_path:
        candidate = Path(explicit_path)
        base_dir = Path(path_root) if path_root else manifest_path.parent
    elif pdb_path_template:
        try:
            candidate = Path(pdb_path_template.format(**row_values))
        except KeyError as exc:
            raise ValueError(
                f"Missing column '{exc.args[0]}' required by --pdb-path-template"
            ) from exc
        base_dir = Path(path_root) if path_root else manifest_path.parent
    else:
        candidate = Path(out_dir) / f"{protein_id}.pdb"
        base_dir = None

    if candidate.is_absolute():
        return candidate
    if base_dir is not None:
        return base_dir / candidate
    return candidate


def write_status(status_rows: list[dict[str, Any]], status_out: str) -> None:
    status_df = pd.DataFrame.from_records(status_rows)
    status_path = Path(status_out)
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_df.to_csv(status_path, index=False)


def print_progress(status_rows: list[dict[str, Any]], completed: int, total: int) -> None:
    status_counts = pd.Series([row["status"] for row in status_rows]).value_counts()
    ok = int(status_counts.get("ok", 0))
    skipped = int(status_counts.get("skipped_existing", 0))
    invalid = int(status_counts.get("invalid_sequence", 0))
    errors = int(status_counts.get("error", 0))
    print(
        f"[{completed}/{total}] ok={ok} skipped_existing={skipped} "
        f"invalid_sequence={invalid} error={errors}"
    )


def main(
    manifest_csv: str,
    out_dir: str,
    status_out: str,
    timeout: int,
    pause_seconds: float,
    skip_existing: bool,
    pdb_path_template: str | None,
    path_root: str | None,
    checkpoint_every: int,
) -> None:
    manifest_path = Path(manifest_csv)
    manifest = pd.read_csv(manifest_path)
    required_columns = {"protein_id", "sequence"}
    missing = required_columns - set(manifest.columns)
    if missing:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing)}")

    output_dir = Path(out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    status_rows = []
    total_rows = len(manifest)
    last_progress_completed = 0
    for index, row in manifest.iterrows():
        row_values = clean_row_values(row)
        protein_id = row_values["protein_id"]
        sequence = normalize_sequence(row_values["sequence"])
        output_path = resolve_output_path(
            row=row,
            manifest_path=manifest_path,
            out_dir=out_dir,
            pdb_path_template=pdb_path_template,
            path_root=path_root,
        )

        invalid_reason = validate_sequence(sequence)
        if invalid_reason:
            status_rows.append(
                {
                    "protein_id": protein_id,
                    "bmrb_id": row.get("bmrb_id"),
                    "entity_id": row.get("entity_id"),
                    "pdb_path": str(output_path.resolve()),
                    "status": "invalid_sequence",
                    "error": invalid_reason,
                }
            )
            if checkpoint_every > 0 and (index + 1) % checkpoint_every == 0:
                write_status(status_rows, status_out)
                print_progress(status_rows, index + 1, total_rows)
                last_progress_completed = index + 1
            continue

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
            output_path.parent.mkdir(parents=True, exist_ok=True)
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
        if checkpoint_every > 0 and (index + 1) % checkpoint_every == 0:
            write_status(status_rows, status_out)
            print_progress(status_rows, index + 1, total_rows)
            last_progress_completed = index + 1

    write_status(status_rows, status_out)
    status_df = pd.DataFrame.from_records(status_rows)
    status_path = Path(status_out)
    ok_count = int((status_df["status"] == "ok").sum()) if not status_df.empty else 0
    if last_progress_completed != total_rows:
        print_progress(status_rows, total_rows, total_rows)
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
    parser.add_argument(
        "--pdb-path-template",
        default=None,
        help="Optional format string for output paths, e.g. "
        "'esmfold_structures/{protein_id}.pdb' or '{bmrb_id}/{protein_id}.pdb'",
    )
    parser.add_argument(
        "--path-root",
        default=None,
        help="Base directory for relative pdb_path values or --pdb-path-template outputs. "
        "Defaults to the manifest's parent directory.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="Write status CSV every N rows so long runs can resume safely (default: 25)",
    )
    args = parser.parse_args()
    main(
        args.manifest,
        args.out_dir,
        args.status_out,
        args.timeout,
        args.pause_seconds,
        args.skip_existing,
        args.pdb_path_template,
        args.path_root,
        args.checkpoint_every,
    )
