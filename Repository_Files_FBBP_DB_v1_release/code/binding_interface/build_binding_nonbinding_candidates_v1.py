from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(r"<local_path_removed>")
LOOP_RESULTS = ROOT / "loop_flexibility_results_long.csv"
WIDE_TABLE = Path(r"<local_path_removed>/main_with_all_results_wide.csv")


def read_csv_fallback(path: Path) -> pd.DataFrame:
    for encoding in (None, "utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding) if encoding else pd.read_csv(path)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("csv", b"", 0, 1, f"Unable to decode {path}")


def parse_target_pdb_ids(value: object) -> list[str]:
    if pd.isna(value):
        return []
    parts = re.split(r"[;|,\s]+", str(value).strip())
    return [p.upper() for p in parts if p]


def infer_tbp_chain_id(loop_id: object) -> str | pd.NA:
    if pd.isna(loop_id):
        return pd.NA
    parts = str(loop_id).split(":")
    if len(parts) >= 2 and parts[1].strip():
        return parts[1].strip()
    return pd.NA


def infer_target_chain_candidate(tbp_chain_id: object, actual_chain_ids: object) -> str | pd.NA:
    if pd.isna(tbp_chain_id) or pd.isna(actual_chain_ids) or not str(actual_chain_ids).strip():
        return pd.NA
    chains = [c for c in str(actual_chain_ids).split("|") if c]
    targets = [c for c in chains if c != str(tbp_chain_id)]
    if not targets:
        return pd.NA
    return "|".join(targets)


def choose_structure_path(row: pd.Series) -> str | pd.NA:
    for col in ("local_pdb_exact_structure_paths", "pdb_path"):
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return pd.NA


def count_chains_in_pdb(path: str) -> tuple[int | pd.NA, str | pd.NA]:
    p = Path(path)
    if not p.exists():
        return pd.NA, pd.NA
    chain_ids: list[str] = []
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) > 21:
                chain_id = line[21].strip() or "_"
                chain_ids.append(chain_id)
    if not chain_ids:
        return 0, ""
    unique = sorted(set(chain_ids))
    return len(unique), "|".join(unique)


def extract_pdb_code(row: pd.Series, structure_path: str | pd.NA) -> str | pd.NA:
    if pd.notna(row.get("pdb_code")) and str(row.get("pdb_code")).strip():
        return str(row.get("pdb_code")).strip().upper()
    if pd.isna(structure_path):
        return pd.NA
    match = re.search(r"_([0-9A-Za-z]{4})_", Path(str(structure_path)).name)
    if match:
        return match.group(1).upper()
    return pd.NA


def classify_route(
    has_target_pdb: bool,
    actual_chain_count: int | pd.NA,
    structure_matches_target_pdb: bool,
) -> tuple[str, str, str]:
    if has_target_pdb and pd.notna(actual_chain_count) and int(actual_chain_count) > 1 and structure_matches_target_pdb:
        return (
            "direct_complex_from_local_exact_pdb",
            "direct",
            "可以直接按复合物局部距离定义 binding residues",
        )
    if has_target_pdb and pd.notna(actual_chain_count) and int(actual_chain_count) > 1:
        return (
            "contextual_complex_needs_chain_mapping",
            "contextual",
            "本地是多链结构，但需要确认哪条链是真正 target chain",
        )
    if has_target_pdb:
        return (
            "target_pdb_available_but_local_structure_single_chain",
            "inferred",
            "有 target PDB 线索，但本地结构看起来是单链或已裁剪，需要补复合物结构或 target chain",
        )
    return (
        "manual_binding_annotation_needed",
        "inferred",
        "当前没有足够结构证据，需要人工提供 binding residue 列表或复合物结构",
    )


def build_binding_candidate_table(loop_df: pd.DataFrame, wide_df: pd.DataFrame) -> pd.DataFrame:
    wide_cols = [
        "protein_row_id",
        "structure_unique_sequence_id",
        "Scaffold_Category",
        "Targets_gene_name",
        "Target_PDB_ID",
        "local_pdb_exact_structure_paths",
        "pdb_path",
        "pdb_code",
        "chain_id",
    ]
    wide_subset = wide_df.reindex(columns=wide_cols)
    merged = loop_df.merge(
        wide_subset.drop_duplicates(),
        on=["protein_row_id", "structure_unique_sequence_id"],
        how="left",
    ).copy()
    merged["selected_local_structure_path"] = merged.apply(choose_structure_path, axis=1)
    merged["target_pdb_id_list"] = merged["Target_PDB_ID"].apply(parse_target_pdb_ids)
    merged["target_pdb_id_count"] = merged["target_pdb_id_list"].apply(len)
    merged["has_target_pdb"] = merged["target_pdb_id_count"] > 0
    merged["tbp_chain_id"] = merged["loop_id"].apply(infer_tbp_chain_id)
    merged["structure_pdb_code"] = merged.apply(
        lambda row: extract_pdb_code(row, row["selected_local_structure_path"]),
        axis=1,
    )
    chain_info = merged["selected_local_structure_path"].apply(
        lambda path: count_chains_in_pdb(str(path)) if pd.notna(path) else (pd.NA, pd.NA)
    )
    merged["actual_chain_count"] = chain_info.apply(lambda x: x[0])
    merged["actual_chain_ids"] = chain_info.apply(lambda x: x[1])
    merged["target_chain_candidate"] = merged.apply(
        lambda row: infer_target_chain_candidate(row["tbp_chain_id"], row["actual_chain_ids"]),
        axis=1,
    )
    merged["structure_matches_target_pdb"] = merged.apply(
        lambda row: pd.notna(row["structure_pdb_code"]) and row["structure_pdb_code"] in row["target_pdb_id_list"],
        axis=1,
    )

    routes = merged.apply(
        lambda row: classify_route(
            has_target_pdb=bool(row["has_target_pdb"]),
            actual_chain_count=row["actual_chain_count"],
            structure_matches_target_pdb=bool(row["structure_matches_target_pdb"]),
        ),
        axis=1,
    )
    merged["binding_candidate_route"] = routes.apply(lambda x: x[0])
    merged["binding_annotation_confidence"] = routes.apply(lambda x: x[1])
    merged["action_needed"] = routes.apply(lambda x: x[2])

    ordered_cols = [
        "protein_row_id",
        "structure_unique_sequence_id",
        "loop_id",
        "highest_priority_cohort",
        "Scaffold_Category",
        "Targets_gene_name",
        "Target_PDB_ID",
        "target_pdb_id_count",
        "selected_local_structure_path",
        "structure_pdb_code",
        "chain_id",
        "tbp_chain_id",
        "actual_chain_count",
        "actual_chain_ids",
        "target_chain_candidate",
        "structure_matches_target_pdb",
        "binding_candidate_route",
        "binding_annotation_confidence",
        "action_needed",
    ]
    return merged[ordered_cols].sort_values(
        ["binding_annotation_confidence", "binding_candidate_route", "loop_id"],
        kind="stable",
    ).reset_index(drop=True)


def main() -> None:
    timestamp = datetime.now().isoformat(timespec="seconds")
    loop_df = read_csv_fallback(LOOP_RESULTS)
    wide_df = read_csv_fallback(WIDE_TABLE)

    candidates = build_binding_candidate_table(
        loop_df=loop_df[
            [
                "protein_row_id",
                "structure_unique_sequence_id",
                "loop_id",
                "highest_priority_cohort",
            ]
        ].copy(),
        wide_df=wide_df,
    )
    out_csv = ROOT / "binding_nonbinding_complex_candidates_v1.csv"
    out_csv.write_text(candidates.to_csv(index=False), encoding="utf-8")

    summary = (
        candidates.groupby(["binding_candidate_route", "binding_annotation_confidence"], dropna=False)
        .agg(n_loops=("loop_id", "size"))
        .reset_index()
        .sort_values("n_loops", ascending=False, kind="stable")
    )
    summary_csv = ROOT / "binding_nonbinding_candidate_summary_v1.csv"
    summary_csv.write_text(summary.to_csv(index=False), encoding="utf-8")

    prefilled_input = candidates.copy()
    prefilled_input["complex_structure_path"] = prefilled_input["selected_local_structure_path"]
    prefilled_input["apo_structure_path"] = pd.NA
    prefilled_input["target_chain_id"] = prefilled_input["target_chain_candidate"]
    prefilled_input["binding_residue_list"] = pd.NA
    prefilled_input["binding_annotation_source"] = prefilled_input["binding_candidate_route"]
    prefilled_input["binding_definition_method"] = pd.NA
    prefilled_input["notes"] = prefilled_input["action_needed"]
    input_cols = [
        "protein_row_id",
        "structure_unique_sequence_id",
        "loop_id",
        "complex_structure_path",
        "apo_structure_path",
        "tbp_chain_id",
        "target_chain_id",
        "binding_residue_list",
        "binding_annotation_source",
        "binding_definition_method",
        "notes",
    ]
    input_csv = ROOT / "loop_subregion_annotations_input_v1.csv"
    input_csv.write_text(prefilled_input[input_cols].to_csv(index=False), encoding="utf-8")

    report_md = ROOT / "binding_nonbinding_next_steps_20260402.md"
    report = f"""# Binding / Non-binding 子区分析下一步

生成时间：{timestamp}

本文件的目的不是直接给出 residue-level binding 结论，而是把当前库里哪些 loop 已经具备“可以继续往界面定义推进”的证据梳理清楚。当前生成的 [binding_nonbinding_complex_candidates_v1.csv](<local_path_removed>/binding_nonbinding_complex_candidates_v1.csv) 按 loop 列出了 target PDB 线索、本地结构路径、实际链数以及建议的后续路线。

如果某个 loop 被标成 `direct_complex_from_local_exact_pdb`，说明当前本地 exact PDB 同时满足两个条件：一是它本身是多链结构，二是它的 PDB 编号已经出现在该条目的 `Target_PDB_ID` 列里。这类条目最适合优先推进，因为理论上已经具备了按距离规则自动提取 binding residues 的必要前提。

如果条目被标成 `contextual_complex_needs_chain_mapping`，说明本地结构是多链的，但我们还不知道哪条链是 target chain。这类条目通常只差最后一步链映射，补上 `target_chain_id` 后就能进入自动距离判定。

如果条目被标成 `target_pdb_available_but_local_structure_single_chain`，说明宽表里虽然已经有 target PDB 线索，但当前本地结构看起来是单链、裁剪结构或非复合物结构。这类条目不是完全没法做，而是需要你后续补其中任意一种：复合物结构路径、target chain、或者人工 binding residue 列表。

如果条目被标成 `manual_binding_annotation_needed`，说明现有表里没有足够的结构证据来定义 interface。这类条目最直接的推进方式就是人工补一个 `binding_residue_list`。为了减少你手工整理的负担，我已经把可填写模板预生成到 [loop_subregion_annotations_input_v1.csv](<local_path_removed>/loop_subregion_annotations_input_v1.csv)，你只需要在有证据的条目上补 `target_chain_id` 或 `binding_residue_list` 即可。

建议实际执行顺序是：先优先处理 `direct` 和 `contextual` 条目，把能自动算 interface 的那部分先拿下来；然后再看是否有必要对剩余高价值 loop 人工补 hotspot 残基。这样可以尽快把 `binding vs non-binding` 子区分析从“全库都缺注释”推进到“至少一批高置信条目可做配对统计”。
"""
    report_md.write_text(report, encoding="utf-8")

    json_path = ROOT / "binding_nonbinding_candidate_summary_v1.json"
    payload = {
        "timestamp": timestamp,
        "n_loops": int(len(candidates)),
        "route_counts": candidates["binding_candidate_route"].value_counts(dropna=False).to_dict(),
        "confidence_counts": candidates["binding_annotation_confidence"].value_counts(dropna=False).to_dict(),
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


if __name__ == "__main__":
    main()
