# Final Structure Bundle Report

- Source CSV: `<local_path_removed>/qc_final_dedup_by_sequence_and_structure.csv`
- Source structure root: `<local_path_removed>
- Generated at: `2026-03-15 22:46:45`
- Bundle root: `<local_path_removed>
- Bundled CSV with rewritten paths: `<local_path_removed>/qc_final_dedup_by_sequence_and_structure.csv`

## Directory Style

- The bundle preserves the original directory style under a new root, including:
  - `all_structures`
  - `all_sequences`
  - `Result_adnectin` ... `Result_phdfingerdomain`
  - `LocalPrediction_AllTargets_1595_pdb_only`

## Copy Summary

| category | copied_count |
| --- | --- |
| pdb_structures | 682 |
| pdb_fastas | 682 |
| uniprot_fastas | 316 |
| afdb_parent_proxy_pdbs | 404 |
| local_prediction_pdbs | 1327 |
| local_prediction_fastas | 2654 |
| subset_local_prediction_manifest_sequences | 1327 |
| missing_source_files | 0 |

## Notes

- Files were copied as a subset from the original structure repository, not regenerated from scratch.
- Relative layout under the new bundle root follows the original `<local_path_removed>/结构预测` organization style.
- Active path columns in the bundled CSV were rewritten from the old structure root to the new bundle root.
- A subset `all_local_prediction_targets.fasta` was generated under `all_structures/Local_Prediction_Input` for the current final main-version records.