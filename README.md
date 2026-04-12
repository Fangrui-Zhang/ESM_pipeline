# esm_pipeline

Utilities for generating ESMFold protein structures from per-protein sequences.

This repo is intentionally limited in scope:
- `esm_pipeline` owns ESMFold structure generation only.
- Manifest building and final dataset collation live in a separate repo.

## Data model

Every protein is tracked with a stable `protein_id`:

`bmr<BMRB_ID>_entity<ENTITY_ID>`

Example:

`bmr6457_entity1`

This matters because one BMRB entry can contain more than one protein entity.

## Expected input

Provide a manifest CSV with at least:
- `protein_id`
- `sequence`

Optional columns such as `bmrb_id` and `entity_id` are preserved in the status output.

## Fetch ESMFold structures

```bash
python scripts/fetch_esmfold.py \
  --manifest /path/to/protein_manifest.csv \
  --out-dir data/esmfold_structures \
  --status-out data/esmfold_status.csv \
  --skip-existing
```

This writes one PDB per protein ID.

If the manifest includes a `pdb_path` column, that per-row path is used instead.
Relative `pdb_path` values are resolved relative to the manifest location by default.

You can also derive output paths without editing the manifest:

```bash
python scripts/fetch_esmfold.py \
  --manifest data/new_manifest.csv \
  --pdb-path-template 'esmfold_structures/{protein_id}.pdb' \
  --status-out data/new_esmfold_status.csv \
  --skip-existing
```

Long runs checkpoint the status CSV every 25 rows by default.

## Notes

The current script uses the public ESM Atlas fold API. If you later want local inference instead, we can swap the backend without changing the input manifest format.
