"""
retry_failed_esmfold.py

Analyses the combined ESMFold status logs, classifies every failed entry,
retries the ones that are worth retrying (timeouts / 5xx server errors),
and writes an updated status CSV + a summary report.

Error categories
----------------
  413 Too Large     → sequence exceeds ESMFold public API limit (~400 aa)
                      → SKIP permanently, log as 'too_large'
  invalid_sequence  → sequence contains non-amino-acid characters
                      → SKIP permanently, log as 'invalid_sequence'
  504 / timeout     → transient server / network issue  → RETRY
  500               → transient server error            → RETRY
  422               → malformed request                 → RETRY once
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import requests

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent
DATA_DIR    = BASE_DIR / "data"
PDB_DIR     = DATA_DIR / "esmfold_structures"
STATUS_FILES = [DATA_DIR / "esmfold_status.csv",
                DATA_DIR / "new_esmfold_status.csv"]

RETRY_STATUS_OUT   = DATA_DIR / "retry_esmfold_status.csv"
FAILED_REPORT_OUT  = DATA_DIR / "failed_entries_report.csv"

MANIFEST_CSV = DATA_DIR / "new_manifest.csv"

API_URL  = "https://api.esmatlas.com/foldSequence/v1/pdb/"
TIMEOUT  = 600   # seconds – generous for large-ish sequences
MAX_RETRIES = 3
RETRY_DELAY = 10  # seconds between retries


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_all_status() -> dict[str, dict]:
    """Merge both status files; later file wins on duplicates."""
    by_protein: dict[str, dict] = {}
    for sf in STATUS_FILES:
        if not sf.exists():
            continue
        with open(sf) as f:
            for row in csv.DictReader(f):
                by_protein[row["protein_id"]] = row
    return by_protein


def load_sequences() -> dict[str, str]:
    seq: dict[str, str] = {}
    with open(MANIFEST_CSV) as f:
        for row in csv.DictReader(f):
            seq[row["protein_id"]] = row["sequence"]
    return seq


def classify_error(error_msg: str) -> str:
    """Return 'retry', 'too_large', 'invalid_sequence', or 'skip_422'."""
    if not error_msg:
        return "retry"
    if "413" in error_msg:
        return "too_large"
    if "invalid sequence" in error_msg.lower():
        return "invalid_sequence"
    if "504" in error_msg or "timed out" in error_msg.lower() or "timeout" in error_msg.lower():
        return "retry"
    if "500" in error_msg:
        return "retry"
    if "422" in error_msg:
        return "retry"
    return "retry"


def fetch_structure(sequence: str) -> str:
    response = requests.post(API_URL, data=sequence, timeout=TIMEOUT)
    response.raise_for_status()
    return response.text


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    PDB_DIR.mkdir(parents=True, exist_ok=True)

    by_protein = load_all_status()
    sequences  = load_sequences()
    manifest_ids = {pid for pid, row in by_protein.items()
                    if row["status"] in ("success", "skipped_existing")}

    # Find all entries that are missing a PDB
    missing = {pid: row for pid, row in by_protein.items()
               if pid not in manifest_ids
               and (PDB_DIR / f"{pid}.pdb").exists() is False}

    # Classify
    to_retry     = {}
    permanent_fail = {}
    for pid, row in missing.items():
        cat = classify_error(row.get("error", ""))
        if cat == "retry":
            to_retry[pid] = row
        else:
            permanent_fail[pid] = {"protein_id": pid,
                                   "status": cat,
                                   "error": row.get("error", "")}

    print(f"\n{'='*60}")
    print(f"  Total missing PDBs          : {len(missing)}")
    print(f"  Permanently skipped (413)   : {sum(1 for r in permanent_fail.values() if r['status']=='too_large')}")
    print(f"  Permanently skipped (invalid): {sum(1 for r in permanent_fail.values() if r['status']=='invalid_sequence')}")
    print(f"  Will retry                  : {len(to_retry)}")
    print(f"{'='*60}\n")

    # ── Retry loop ────────────────────────────────────────────────────────────
    retry_results: list[dict] = []

    for i, (pid, row) in enumerate(to_retry.items(), 1):
        sequence = sequences.get(pid, "")
        if not sequence:
            print(f"[{i}/{len(to_retry)}] {pid} — no sequence in manifest, skipping")
            retry_results.append({"protein_id": pid, "bmrb_id": row["bmrb_id"],
                                   "entity_id": row["entity_id"],
                                   "pdb_path": "", "status": "no_sequence", "error": ""})
            continue

        pdb_path = PDB_DIR / f"{pid}.pdb"
        success  = False
        last_err = ""

        for attempt in range(1, MAX_RETRIES + 1):
            try:
                print(f"[{i}/{len(to_retry)}] {pid}  attempt {attempt}/{MAX_RETRIES} "
                      f"(seq len={len(sequence)}) …", end=" ", flush=True)
                pdb_text = fetch_structure(sequence)
                pdb_path.write_text(pdb_text)
                print("✅ success")
                retry_results.append({
                    "protein_id": pid,
                    "bmrb_id":   row["bmrb_id"],
                    "entity_id": row["entity_id"],
                    "pdb_path":  str(pdb_path),
                    "status":    "success",
                    "error":     "",
                })
                success = True
                break
            except Exception as e:
                last_err = str(e)
                print(f"❌ {last_err[:80]}")
                if attempt < MAX_RETRIES:
                    time.sleep(RETRY_DELAY)

        if not success:
            retry_results.append({
                "protein_id": pid,
                "bmrb_id":   row["bmrb_id"],
                "entity_id": row["entity_id"],
                "pdb_path":  "",
                "status":    "failed_retry",
                "error":     last_err,
            })

    # ── Write retry status CSV ────────────────────────────────────────────────
    fieldnames = ["protein_id", "bmrb_id", "entity_id", "pdb_path", "status", "error"]
    with open(RETRY_STATUS_OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(retry_results)

    # ── Write full failed report (permanent + still-failed retries) ───────────
    all_failed = list(permanent_fail.values()) + \
                 [r for r in retry_results if r["status"] not in ("success",)]

    report_fields = ["protein_id", "bmrb_id", "entity_id", "status", "error"]
    with open(FAILED_REPORT_OUT, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=report_fields, extrasaction="ignore")
        writer.writeheader()
        # Add bmrb/entity for permanent fails from original status
        for row in permanent_fail.values():
            orig = by_protein.get(row["protein_id"], {})
            row["bmrb_id"]   = orig.get("bmrb_id", "")
            row["entity_id"] = orig.get("entity_id", "")
            writer.writerow(row)
        writer.writerows([r for r in retry_results if r["status"] not in ("success",)])

    # ── Summary ───────────────────────────────────────────────────────────────
    n_success      = sum(1 for r in retry_results if r["status"] == "success")
    n_still_failed = sum(1 for r in retry_results if r["status"] == "failed_retry")

    print(f"\n{'='*60}")
    print(f"  Retry results:")
    print(f"    Newly succeeded  : {n_success}")
    print(f"    Still failed     : {n_still_failed}")
    print(f"    Permanent skips  : {len(permanent_fail)}")
    print(f"\n  Retry status CSV : {RETRY_STATUS_OUT}")
    print(f"  Full failure report: {FAILED_REPORT_OUT}")
    print(f"{'='*60}\n")

    if n_success > 0:
        print("  ⚠️  New PDB files were added — re-run build_pdb_manifest_and_zip.py")
        print("     to update the manifest and zip.")


if __name__ == "__main__":
    main()
