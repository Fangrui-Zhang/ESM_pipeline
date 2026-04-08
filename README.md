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

## Notes

The current script uses the public ESM Atlas fold API. If you later want local inference instead, we can swap the backend without changing the input manifest format.
