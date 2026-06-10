from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd

try:
    from scipy.stats import chi2_contingency
except Exception:  # pragma: no cover
    chi2_contingency = None


ROOT = Path(r"<local_path_removed>")
LOOP_RESULTS = ROOT / "loop_flexibility_results_long.csv"
WIDE_TABLE = Path(r"<local_path_removed>/main_with_all_results_wide.csv")
AIR_FLEX = ROOT / "ALL_conformations_extracted" / "ALL_conformations" / "ALL_conformations_flexibility.csv"


def read_csv_fallback(path: Path) -> pd.DataFrame:
    for encoding in (None, "utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding) if encoding else pd.read_csv(path)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("csv", b"", 0, 1, f"Unable to decode {path}")


def derive_loop_sequence(full_sequence: str, start: float, end: float) -> str | pd.NA:
    if pd.isna(full_sequence) or pd.isna(start) or pd.isna(end):
        return pd.NA
    sequence = str(full_sequence).strip()
    start_i = int(start)
    end_i = int(end)
    if start_i < 1 or end_i < start_i or end_i > len(sequence):
        return pd.NA
    return sequence[start_i - 1 : end_i]


def classify_tbp_binary(itsflex_class: str) -> str | pd.NA:
    if pd.isna(itsflex_class):
        return pd.NA
    rigid = {"high_confidence_rigid", "low_confidence_rigid"}
    flexible = {"low_confidence_flexible", "high_confidence_flexible"}
    if itsflex_class in rigid:
        return "rigid"
    if itsflex_class in flexible:
        return "flexible"
    return pd.NA


def build_outputs() -> dict[str, Path]:
    ts = datetime.now().isoformat(timespec="seconds")
    outputs = {
        "air_reference": ROOT / "tbp_vs_airs_air_reference_v1.csv",
        "loop_level": ROOT / "tbp_vs_airs_loop_level_comparison_v1.csv",
        "exact_matches": ROOT / "tbp_vs_airs_exact_sequence_matches_v1.csv",
        "length_summary": ROOT / "tbp_vs_airs_length_summary_v1.csv",
        "scaffold_summary": ROOT / "tbp_vs_airs_scaffold_summary_v1.csv",
        "dataset_summary": ROOT / "tbp_vs_airs_dataset_summary_v1.csv",
        "json_summary": ROOT / "tbp_vs_airs_summary_v1.json",
        "report_md": ROOT / "tbp_vs_airs_analysis_20260401.md",
    }
    outputs["timestamp"] = ts  # type: ignore[index]
    return outputs


def main() -> None:
    outputs = build_outputs()

    loop_df = read_csv_fallback(LOOP_RESULTS)
    wide_df = read_csv_fallback(WIDE_TABLE)
    air_df = pd.read_csv(AIR_FLEX)

    loop_df["protein_row_id"] = pd.to_numeric(loop_df["protein_row_id"], errors="coerce")
    wide_merge = wide_df[
        [
            "protein_row_id",
            "structure_unique_sequence_id",
            "Final_Tested_Sequence",
            "Scaffold_Category",
            "Targets_gene_name",
            "Target_PDB_ID",
        ]
    ].drop_duplicates()
    tbp = loop_df.merge(
        wide_merge,
        on=["protein_row_id", "structure_unique_sequence_id"],
        how="left",
    )
    tbp["loop_sequence"] = tbp.apply(
        lambda row: derive_loop_sequence(row.get("Final_Tested_Sequence"), row.get("loop_start"), row.get("loop_end")),
        axis=1,
    )
    tbp["loop_sequence"] = tbp["loop_sequence"].astype("string").str.upper()
    tbp["tbp_itsflex_binary_label"] = tbp["itsflex_class"].apply(classify_tbp_binary)
    tbp["tbp_consensus_binary_label"] = tbp["flexibility_consensus_label"].map(
        {"rigid": "rigid", "flexible": "flexible"}
    )
    tbp["tbp_has_dynamic_score"] = tbp["flexscore_ensemble_diversity"].notna()
    tbp["loop_sequence_length_check"] = tbp["loop_sequence"].str.len()
    tbp["loop_sequence_length_matches"] = tbp["loop_sequence_length_check"] == tbp["loop_length"]

    air = air_df.rename(columns={"Unnamed: 0": "air_record_id"}).copy()
    air["seq_loop"] = air["seq_loop"].astype("string").str.upper()
    air["air_group"] = "AIRs"
    air["is_pdb_set"] = air["set"].eq("PDB")
    air["is_immune_set"] = ~air["set"].eq("PDB")
    air_reference = air[
        [
            "air_record_id",
            "seq_loop",
            "flexibility_label",
            "loop_length",
            "set",
            "n_structures",
            "n_pdbs",
            "n_conformations_0.5A",
            "n_conformations_1.25A",
            "n_conformations_2A",
            "chains",
        ]
    ].copy()
    air_reference.to_csv(outputs["air_reference"], index=False)

    exact_matches = (
        tbp[["loop_id", "loop_sequence"]]
        .merge(
            air_reference,
            left_on="loop_sequence",
            right_on="seq_loop",
            how="inner",
        )
        .sort_values(["loop_id", "air_record_id"], kind="stable")
    )
    exact_matches.to_csv(outputs["exact_matches"], index=False)

    exact_summary = (
        exact_matches.groupby("loop_id", dropna=False)
        .agg(
            air_exact_seq_match_count=("air_record_id", "size"),
            air_exact_seq_match_flexible_count=("flexibility_label", lambda s: int((s == "flexible").sum())),
            air_exact_seq_match_rigid_count=("flexibility_label", lambda s: int((s == "rigid").sum())),
            air_exact_seq_match_sets=("set", lambda s: "|".join(sorted(pd.Series(s).dropna().astype(str).unique()))),
        )
        .reset_index()
    )

    air_length_summary = (
        air_reference.groupby("loop_length", dropna=False)
        .agg(
            air_total=("air_record_id", "size"),
            air_flexible_count=("flexibility_label", lambda s: int((s == "flexible").sum())),
            air_rigid_count=("flexibility_label", lambda s: int((s == "rigid").sum())),
        )
        .reset_index()
    )
    air_length_summary["air_flexible_fraction"] = (
        air_length_summary["air_flexible_count"] / air_length_summary["air_total"]
    )

    air_length_pdb = (
        air_reference[air_reference["set"] == "PDB"]
        .groupby("loop_length", dropna=False)
        .agg(
            air_pdb_total=("air_record_id", "size"),
            air_pdb_flexible_count=("flexibility_label", lambda s: int((s == "flexible").sum())),
        )
        .reset_index()
    )
    air_length_pdb["air_pdb_flexible_fraction"] = air_length_pdb["air_pdb_flexible_count"] / air_length_pdb["air_pdb_total"]

    air_length_immune = (
        air_reference[air_reference["set"] != "PDB"]
        .groupby("loop_length", dropna=False)
        .agg(
            air_immune_total=("air_record_id", "size"),
            air_immune_flexible_count=("flexibility_label", lambda s: int((s == "flexible").sum())),
        )
        .reset_index()
    )
    air_length_immune["air_immune_flexible_fraction"] = (
        air_length_immune["air_immune_flexible_count"] / air_length_immune["air_immune_total"]
    )

    tbp = tbp.merge(air_length_summary, on="loop_length", how="left")
    tbp = tbp.merge(air_length_pdb, on="loop_length", how="left")
    tbp = tbp.merge(air_length_immune, on="loop_length", how="left")
    tbp = tbp.merge(exact_summary, on="loop_id", how="left")
    tbp["air_exact_seq_match_count"] = tbp["air_exact_seq_match_count"].fillna(0).astype(int)
    tbp["air_exact_seq_match_flexible_count"] = tbp["air_exact_seq_match_flexible_count"].fillna(0).astype(int)
    tbp["air_exact_seq_match_rigid_count"] = tbp["air_exact_seq_match_rigid_count"].fillna(0).astype(int)
    tbp.to_csv(outputs["loop_level"], index=False)

    tbp_binary = tbp[tbp["tbp_itsflex_binary_label"].notna()].copy()
    length_summary = (
        tbp_binary.groupby("loop_length", dropna=False)
        .agg(
            tbp_total=("loop_id", "size"),
            tbp_flexible_count=("tbp_itsflex_binary_label", lambda s: int((s == "flexible").sum())),
            tbp_rigid_count=("tbp_itsflex_binary_label", lambda s: int((s == "rigid").sum())),
            tbp_primary_count=("highest_priority_cohort", lambda s: int((s == "primary_approved_v1").sum())),
            tbp_exact_seq_match_count=("air_exact_seq_match_count", lambda s: int((s > 0).sum())),
        )
        .reset_index()
        .merge(air_length_summary, on="loop_length", how="outer")
        .sort_values("loop_length", kind="stable")
    )
    length_summary["tbp_flexible_fraction"] = length_summary["tbp_flexible_count"] / length_summary["tbp_total"]
    length_summary["tbp_minus_air_flexible_fraction"] = (
        length_summary["tbp_flexible_fraction"] - length_summary["air_flexible_fraction"]
    )
    length_summary.to_csv(outputs["length_summary"], index=False)

    length_weights = air_length_summary.set_index("loop_length")["air_flexible_fraction"].to_dict()
    scaffold_summary = (
        tbp_binary.groupby("Scaffold_Category", dropna=False)
        .agg(
            tbp_total=("loop_id", "size"),
            tbp_flexible_count=("tbp_itsflex_binary_label", lambda s: int((s == "flexible").sum())),
            tbp_rigid_count=("tbp_itsflex_binary_label", lambda s: int((s == "rigid").sum())),
            mean_loop_length=("loop_length", "mean"),
            exact_seq_match_loops=("air_exact_seq_match_count", lambda s: int((s > 0).sum())),
        )
        .reset_index()
    )
    scaffold_summary["tbp_flexible_fraction"] = scaffold_summary["tbp_flexible_count"] / scaffold_summary["tbp_total"]

    scaffold_weighted = []
    for scaffold, group in tbp_binary.groupby("Scaffold_Category", dropna=False):
        weights = [length_weights.get(length) for length in group["loop_length"] if length in length_weights]
        scaffold_weighted.append(
            {
                "Scaffold_Category": scaffold,
                "air_length_weighted_flexible_fraction": (sum(weights) / len(weights)) if weights else pd.NA,
            }
        )
    scaffold_summary = scaffold_summary.merge(pd.DataFrame(scaffold_weighted), on="Scaffold_Category", how="left")
    scaffold_summary["tbp_minus_air_length_weighted_flexible_fraction"] = (
        scaffold_summary["tbp_flexible_fraction"] - scaffold_summary["air_length_weighted_flexible_fraction"]
    )
    scaffold_summary.to_csv(outputs["scaffold_summary"], index=False)

    dataset_summary = pd.DataFrame(
        [
            {
                "dataset": "TBP_all_loops",
                "n_loops": int(len(tbp)),
                "n_flexible": int((tbp["tbp_itsflex_binary_label"] == "flexible").sum()),
                "n_rigid": int((tbp["tbp_itsflex_binary_label"] == "rigid").sum()),
                "n_binary_comparable": int(tbp["tbp_itsflex_binary_label"].notna().sum()),
                "flexible_fraction": float((tbp["tbp_itsflex_binary_label"] == "flexible").mean(skipna=True)),
            },
            {
                "dataset": "TBP_primary_approved",
                "n_loops": int((tbp["highest_priority_cohort"] == "primary_approved_v1").sum()),
                "n_flexible": int(((tbp["highest_priority_cohort"] == "primary_approved_v1") & (tbp["tbp_itsflex_binary_label"] == "flexible")).sum()),
                "n_rigid": int(((tbp["highest_priority_cohort"] == "primary_approved_v1") & (tbp["tbp_itsflex_binary_label"] == "rigid")).sum()),
                "n_binary_comparable": int(((tbp["highest_priority_cohort"] == "primary_approved_v1") & tbp["tbp_itsflex_binary_label"].notna()).sum()),
                "flexible_fraction": float(
                    (
                        (tbp["highest_priority_cohort"] == "primary_approved_v1")
                        & (tbp["tbp_itsflex_binary_label"] == "flexible")
                    ).sum()
                    / max(1, ((tbp["highest_priority_cohort"] == "primary_approved_v1") & tbp["tbp_itsflex_binary_label"].notna()).sum())
                ),
            },
            {
                "dataset": "AIRs_all",
                "n_loops": int(len(air_reference)),
                "n_flexible": int((air_reference["flexibility_label"] == "flexible").sum()),
                "n_rigid": int((air_reference["flexibility_label"] == "rigid").sum()),
                "n_binary_comparable": int(len(air_reference)),
                "flexible_fraction": float((air_reference["flexibility_label"] == "flexible").mean()),
            },
            {
                "dataset": "AIRs_PDB",
                "n_loops": int((air_reference["set"] == "PDB").sum()),
                "n_flexible": int(((air_reference["set"] == "PDB") & (air_reference["flexibility_label"] == "flexible")).sum()),
                "n_rigid": int(((air_reference["set"] == "PDB") & (air_reference["flexibility_label"] == "rigid")).sum()),
                "n_binary_comparable": int((air_reference["set"] == "PDB").sum()),
                "flexible_fraction": float(
                    ((air_reference["set"] == "PDB") & (air_reference["flexibility_label"] == "flexible")).sum()
                    / max(1, (air_reference["set"] == "PDB").sum())
                ),
            },
        ]
    )
    dataset_summary.to_csv(outputs["dataset_summary"], index=False)

    chi_square_result = None
    if chi2_contingency is not None:
        contingency = pd.DataFrame(
            {
                "TBP": tbp_binary["tbp_itsflex_binary_label"].value_counts().reindex(["rigid", "flexible"], fill_value=0),
                "AIRs": air_reference["flexibility_label"].value_counts().reindex(["rigid", "flexible"], fill_value=0),
            }
        )
        chi2, p_value, dof, expected = chi2_contingency(contingency)
        chi_square_result = {
            "chi2": float(chi2),
            "p_value": float(p_value),
            "dof": int(dof),
            "contingency": contingency.to_dict(),
        }

    json_summary = {
        "timestamp": outputs["timestamp"],
        "tbp_total_loops": int(len(tbp)),
        "tbp_binary_comparable_loops": int(tbp["tbp_itsflex_binary_label"].notna().sum()),
        "tbp_exact_sequence_matches": int((tbp["air_exact_seq_match_count"] > 0).sum()),
        "air_total_loops": int(len(air_reference)),
        "air_set_counts": air_reference["set"].value_counts(dropna=False).to_dict(),
        "air_flexible_fraction": float((air_reference["flexibility_label"] == "flexible").mean()),
        "tbp_flexible_fraction_itsflex_binary": float((tbp_binary["tbp_itsflex_binary_label"] == "flexible").mean()),
        "chi_square_tbp_vs_airs_binary": chi_square_result,
    }
    with outputs["json_summary"].open("w", encoding="utf-8") as fh:
        json.dump(json_summary, fh, indent=2, ensure_ascii=False)

    report = f"""# TBP vs AIRs Comparison v1

生成时间：{outputs["timestamp"]}

## 数据源

- TBP loop 结果：`{LOOP_RESULTS}`
- TBP 主表：`{WIDE_TABLE}`
- AIR 参考集：`{AIR_FLEX}`

## 本版比较能回答什么

这是一版首轮对照，不是假设 TBP 与 AIRs 有逐条一一对应关系。
本版主要做三件事：

1. 清洗 AIR `ALL_conformations_flexibility.csv` 成可复用参考表。
2. 给每个 TBP loop 接上同长度 AIR 柔性基线。
3. 检查 TBP loop 是否与 AIR `seq_loop` 存在 exact sequence overlap。

## 关键现状

- TBP loop 总数：{len(tbp)}
- TBP 可做二元柔/刚比较的 loop 数：{int(tbp['tbp_itsflex_binary_label'].notna().sum())}
- AIR loop 总数：{len(air_reference)}
- TBP exact sequence match 数：{int((tbp['air_exact_seq_match_count'] > 0).sum())}

## 说明

- TBP 二元柔/刚标签来自 `itsflex_class`：
  - `high/low_confidence_flexible -> flexible`
  - `high/low_confidence_rigid -> rigid`
  - `ambiguous / NA -> excluded from binary comparison`
- AIR 侧直接使用官方 `flexibility_label`
- `tbp_vs_airs_length_summary_v1.csv` 是首选长度分层对照表
- `tbp_vs_airs_scaffold_summary_v1.csv` 提供了 TBP scaffold 与 AIR 长度加权基线的比较
- 当前 TBP 动态分数字段仍来自现有 `loop_flexibility_results_long.csv`；等 `full_v2_masked` 跑完并重新 merge 后，这组对照表建议再刷新一次
- `tbp_vs_airs_exact_sequence_matches_v1.csv` 中的 exact match 多数是短 loop 或短 motif 的序列重合，不能直接等同于同源或同功能位点映射
"""
    outputs["report_md"].write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
