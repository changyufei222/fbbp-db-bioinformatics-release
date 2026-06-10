# Clean FASTA Regeneration Report

- Source CSV: `<local_path_removed>/qc_final_dedup_by_sequence_and_structure.csv`
- Generated at: `2026-03-15 23:04:40`
- Per-protein FASTA directory: `<local_path_removed>
- Combined FASTA: `<local_path_removed>/all_protein_sequences.fasta`
- Total protein FASTA entries regenerated: `1996`
- Rows falling back from `Final_Tested_Sequence` to `original_Final_Tested_Sequence`: `56`

## Rule

- Default: use `Final_Tested_Sequence`.
- If `Final_Tested_Sequence` contains `X` and `original_Final_Tested_Sequence` is available, use the original clean sequence instead.