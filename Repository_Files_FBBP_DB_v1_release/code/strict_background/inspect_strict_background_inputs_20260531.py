from pathlib import Path

import pandas as pd


ROOT = Path(r"<local_path_removed>")
AIR_PDBS = ROOT / "柔性分析最终" / "ALL_conformations_extracted" / "ALL_conformations" / "ALL_conformations_PDBs.csv"
TBP_ENHANCED = ROOT / "柔性分析最终" / "tbp_vs_airs_enhanced_loop_level_v2.csv"
STRUCT = ROOT / "FBBP_DB_v1_frozen_release_20260531" / "data" / "schema_tables" / "structural_info.csv"


def main() -> None:
    air_cols = pd.read_csv(AIR_PDBS, nrows=0).columns.tolist()
    tbp_cols = pd.read_csv(TBP_ENHANCED, nrows=0).columns.tolist()
    struct_cols = pd.read_csv(STRUCT, nrows=0).columns.tolist()
    print("AIR_PDBS columns:", air_cols)
    print("TBP columns:", tbp_cols)
    print("STRUCT columns:", struct_cols)

    tbp = pd.read_csv(
        TBP_ENHANCED,
        usecols=[
            "loop_id",
            "structure_unique_sequence_id",
            "loop_sequence",
            "loop_length",
            "selected_structure_source",
            "itsflex_class",
            "tbp_binary_label",
            "analysis_cohort",
            "highest_priority_cohort",
        ],
    )
    print("\nTBP rows:", len(tbp))
    print("selected_structure_source:")
    print(tbp["selected_structure_source"].value_counts(dropna=False).to_string())
    print("\nanalysis_cohort:")
    print(tbp["analysis_cohort"].value_counts(dropna=False).to_string())
    print("\ntbp_binary_label:")
    print(tbp["tbp_binary_label"].value_counts(dropna=False).to_string())

    struct = pd.read_csv(STRUCT)
    print("\nSTRUCT method:")
    print(struct["method"].value_counts(dropna=False).head(20).to_string())
    print("\nSTRUCT resolution non-null:", int(struct["resolution"].notna().sum()))
    print(struct["resolution"].describe().to_string())

    air = pd.read_csv(AIR_PDBS, usecols=["pdb", "chain", "seq_loop", "loop_length", "structure_method", "resolution"])
    print("\nAIR rows:", len(air))
    print("AIR methods:")
    print(air["structure_method"].value_counts(dropna=False).head(20).to_string())
    print("\nAIR resolution non-null:", int(air["resolution"].notna().sum()))
    print(air["resolution"].describe().to_string())


if __name__ == "__main__":
    main()
