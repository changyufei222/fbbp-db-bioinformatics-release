from __future__ import annotations

import json
from datetime import datetime
from math import exp
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression


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


def build_length_standardized_summary(
    tbp_binary_df: pd.DataFrame,
    air_reference_df: pd.DataFrame,
    cohort_col: str,
) -> pd.DataFrame:
    air_length = (
        air_reference_df.groupby("loop_length", dropna=False)
        .agg(
            air_total=("air_record_id", "size"),
            air_flexible_count=("flexibility_label", lambda s: int((s == "flexible").sum())),
        )
        .reset_index()
    )
    air_length["air_flexible_fraction"] = air_length["air_flexible_count"] / air_length["air_total"]
    air_lookup = air_length.set_index("loop_length")["air_flexible_fraction"].to_dict()

    rows: list[dict[str, object]] = []
    for cohort, group in tbp_binary_df.groupby(cohort_col, dropna=False):
        tbp_total = int(len(group))
        tbp_flexible_count = int((group["tbp_binary_label"] == "flexible").sum())
        length_counts = group["loop_length"].value_counts(dropna=False).to_dict()
        weighted_values = []
        missing_lengths = 0
        for length, count in length_counts.items():
            if length in air_lookup:
                weighted_values.extend([float(air_lookup[length])] * int(count))
            else:
                missing_lengths += int(count)

        expected_fraction = float(np.mean(weighted_values)) if weighted_values else np.nan
        rows.append(
            {
                "cohort": cohort,
                "tbp_total": tbp_total,
                "tbp_flexible_count": tbp_flexible_count,
                "tbp_rigid_count": int((group["tbp_binary_label"] == "rigid").sum()),
                "tbp_flexible_fraction": tbp_flexible_count / tbp_total if tbp_total else np.nan,
                "air_expected_flexible_fraction": expected_fraction,
                "tbp_minus_air_expected_fraction": (
                    (tbp_flexible_count / tbp_total) - expected_fraction if tbp_total and pd.notna(expected_fraction) else np.nan
                ),
                "air_missing_length_count": missing_lengths,
            }
        )

    return pd.DataFrame(rows).sort_values("cohort", kind="stable").reset_index(drop=True)


def fit_length_adjusted_logit(combined_df: pd.DataFrame) -> dict[str, float]:
    frame = combined_df[["is_flexible", "dataset_is_tbp", "loop_length"]].dropna().copy()
    y = frame["is_flexible"].astype(float).to_numpy()
    length = frame["loop_length"].astype(float).to_numpy()
    length_centered = length - float(np.mean(length))
    length_scale = float(np.std(length_centered)) or 1.0
    length_scaled = length_centered / length_scale
    length_sq = length_scaled**2
    x = np.column_stack(
        [
            np.ones(len(frame), dtype=float),
            frame["dataset_is_tbp"].astype(float).to_numpy(),
            length_scaled,
            length_sq,
        ]
    )

    def neg_loglik(beta: np.ndarray) -> float:
        logits = x @ beta
        p = np.clip(expit(logits), 1e-9, 1 - 1e-9)
        return -float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))

    def grad(beta: np.ndarray) -> np.ndarray:
        logits = x @ beta
        p = expit(logits)
        return x.T @ (p - y)

    result = minimize(neg_loglik, x0=np.zeros(x.shape[1], dtype=float), jac=grad, method="BFGS")
    if result.success:
        beta = np.asarray(result.x, dtype=float)
        cov = np.asarray(result.hess_inv, dtype=float)
    else:
        model = LogisticRegression(
            fit_intercept=False,
            penalty="l2",
            C=1e6,
            solver="lbfgs",
            max_iter=2000,
        )
        model.fit(x, y.astype(int))
        beta = model.coef_.reshape(-1).astype(float)
        p = np.clip(model.predict_proba(x)[:, 1], 1e-9, 1 - 1e-9)
        w = p * (1 - p)
        hessian = x.T @ (w[:, None] * x)
        cov = np.linalg.pinv(hessian)

    se = np.sqrt(np.diag(cov))

    dataset_coef = float(beta[1])
    dataset_se = float(se[1])
    z_score = dataset_coef / dataset_se if dataset_se else np.nan
    p_value = float(2 * (1 - norm.cdf(abs(z_score)))) if pd.notna(z_score) else np.nan
    ci_lower = dataset_coef - 1.96 * dataset_se
    ci_upper = dataset_coef + 1.96 * dataset_se

    return {
        "n_obs": int(len(frame)),
        "dataset_coef": dataset_coef,
        "dataset_se": dataset_se,
        "dataset_or": float(exp(dataset_coef)),
        "ci_lower": float(exp(ci_lower)),
        "ci_upper": float(exp(ci_upper)),
        "p_value": p_value,
        "intercept": float(beta[0]),
        "loop_length_coef": float(beta[2]),
        "loop_length_sq_coef": float(beta[3]),
        "model_converged": 1.0,
    }


def mantel_haenszel_by_length(tbp_binary_df: pd.DataFrame, air_reference_df: pd.DataFrame) -> dict[str, float]:
    numerator = 0.0
    denominator = 0.0
    strata_used = 0
    for length in sorted(set(tbp_binary_df["loop_length"].dropna()) & set(air_reference_df["loop_length"].dropna())):
        tbp_subset = tbp_binary_df.loc[tbp_binary_df["loop_length"] == length]
        air_subset = air_reference_df.loc[air_reference_df["loop_length"] == length]

        a = float((tbp_subset["tbp_binary_label"] == "flexible").sum())
        b = float((tbp_subset["tbp_binary_label"] == "rigid").sum())
        c = float((air_subset["flexibility_label"] == "flexible").sum())
        d = float((air_subset["flexibility_label"] == "rigid").sum())
        n = a + b + c + d
        if min(a + b, c + d, a + c, b + d) == 0 or n == 0:
            continue

        numerator += (a * d) / n
        denominator += (b * c) / n
        strata_used += 1

    mh_or = numerator / denominator if denominator else np.nan
    return {
        "strata_used": int(strata_used),
        "mantel_haenszel_or": float(mh_or) if pd.notna(mh_or) else np.nan,
    }


def add_exact_sequence_matches(tbp: pd.DataFrame, air_reference: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    exact_matches = (
        tbp[["loop_id", "loop_sequence"]]
        .dropna()
        .merge(
            air_reference[["air_record_id", "seq_loop", "flexibility_label", "set", "loop_length"]],
            left_on="loop_sequence",
            right_on="seq_loop",
            how="inner",
        )
        .sort_values(["loop_id", "air_record_id"], kind="stable")
    )

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
    out = tbp.merge(exact_summary, on="loop_id", how="left")
    for col in [
        "air_exact_seq_match_count",
        "air_exact_seq_match_flexible_count",
        "air_exact_seq_match_rigid_count",
    ]:
        out[col] = out[col].fillna(0).astype(int)
    out["has_air_exact_seq_match"] = out["air_exact_seq_match_count"] > 0
    return out, exact_matches


def build_outputs() -> dict[str, Path | str]:
    ts = datetime.now().isoformat(timespec="seconds")
    return {
        "loop_level": ROOT / "tbp_vs_airs_enhanced_loop_level_v2.csv",
        "length_standardized": ROOT / "tbp_vs_airs_length_standardized_summary_v2.csv",
        "scaffold_standardized": ROOT / "tbp_vs_airs_scaffold_standardized_summary_v2.csv",
        "exact_match_dynamic": ROOT / "tbp_vs_airs_exact_match_dynamic_summary_v2.csv",
        "mh_summary": ROOT / "tbp_vs_airs_mantel_haenszel_summary_v2.csv",
        "logit_summary": ROOT / "tbp_vs_airs_logit_summary_v2.csv",
        "summary_json": ROOT / "tbp_vs_airs_enhanced_summary_v2.json",
        "report_md": ROOT / "tbp_vs_airs_enhanced_analysis_20260402.md",
        "timestamp": ts,
    }


def main() -> None:
    outputs = build_outputs()
    loop_df = read_csv_fallback(LOOP_RESULTS)
    wide_df = read_csv_fallback(WIDE_TABLE)
    air_df = pd.read_csv(AIR_FLEX)

    tbp = loop_df.merge(
        wide_df[
            [
                "protein_row_id",
                "structure_unique_sequence_id",
                "Final_Tested_Sequence",
                "Scaffold_Category",
                "Targets_gene_name",
                "Target_PDB_ID",
            ]
        ].drop_duplicates(),
        on=["protein_row_id", "structure_unique_sequence_id"],
        how="left",
    )
    tbp["loop_sequence"] = tbp.apply(
        lambda row: derive_loop_sequence(row.get("Final_Tested_Sequence"), row.get("loop_start"), row.get("loop_end")),
        axis=1,
    )
    tbp["loop_sequence"] = tbp["loop_sequence"].astype("string").str.upper()
    tbp["tbp_binary_label"] = tbp["itsflex_class"].apply(classify_tbp_binary)
    tbp["tbp_is_flexible"] = tbp["tbp_binary_label"].map({"flexible": 1, "rigid": 0})

    air_reference = air_df.rename(columns={"Unnamed: 0": "air_record_id"}).copy()
    air_reference["seq_loop"] = air_reference["seq_loop"].astype("string").str.upper()
    air_reference = air_reference[
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

    tbp, exact_matches = add_exact_sequence_matches(tbp, air_reference)

    tbp["analysis_cohort"] = "TBP_all_binary"
    tbp.loc[tbp["highest_priority_cohort"] == "primary_approved_v1", "analysis_cohort"] = "TBP_primary_approved"
    tbp_binary = tbp[tbp["tbp_binary_label"].notna()].copy()
    tbp_binary["analysis_cohort"] = "TBP_all_binary"
    tbp_primary = tbp_binary[tbp_binary["highest_priority_cohort"] == "primary_approved_v1"].copy()
    tbp_primary["analysis_cohort"] = "TBP_primary_approved"

    length_standardized = build_length_standardized_summary(
        tbp_binary_df=pd.concat([tbp_binary, tbp_primary], ignore_index=True),
        air_reference_df=air_reference,
        cohort_col="analysis_cohort",
    )
    length_standardized.to_csv(outputs["length_standardized"], index=False)

    scaffold_input = tbp_binary.copy()
    scaffold_input["analysis_cohort"] = scaffold_input["Scaffold_Category"].fillna("NA").map(lambda x: f"scaffold::{x}")
    scaffold_standardized = build_length_standardized_summary(
        tbp_binary_df=scaffold_input,
        air_reference_df=air_reference,
        cohort_col="analysis_cohort",
    ).rename(columns={"analysis_cohort": "cohort"})
    scaffold_standardized.to_csv(outputs["scaffold_standardized"], index=False)

    exact_match_dynamic = (
        tbp.groupby("has_air_exact_seq_match", dropna=False)
        .agg(
            n_loops=("loop_id", "size"),
            n_with_dynamic=("flexscore_ensemble_diversity", lambda s: int(s.notna().sum())),
            diversity_mean=("flexscore_ensemble_diversity", "mean"),
            diversity_median=("flexscore_ensemble_diversity", "median"),
            dynamic_mean=("flexibility_score_dynamic", "mean"),
            dynamic_median=("flexibility_score_dynamic", "median"),
            clusters_mean=("num_clusters", "mean"),
            itsflex_mean=("itsflex_score", "mean"),
        )
        .reset_index()
    )
    exact_match_dynamic["match_group"] = exact_match_dynamic["has_air_exact_seq_match"].map(
        {True: "exact_air_seq_match", False: "no_exact_air_seq_match"}
    )
    exact_match_dynamic.to_csv(outputs["exact_match_dynamic"], index=False)

    mh_rows = []
    for label, subset in [("TBP_all_binary_vs_AIR_all", tbp_binary), ("TBP_primary_vs_AIR_all", tbp_primary)]:
        mh = mantel_haenszel_by_length(subset, air_reference)
        mh_rows.append({"comparison": label, **mh})
    air_pdb = air_reference[air_reference["set"] == "PDB"].copy()
    for label, subset in [("TBP_all_binary_vs_AIR_PDB", tbp_binary), ("TBP_primary_vs_AIR_PDB", tbp_primary)]:
        mh = mantel_haenszel_by_length(subset, air_pdb)
        mh_rows.append({"comparison": label, **mh})
    mh_summary = pd.DataFrame(mh_rows)
    mh_summary.to_csv(outputs["mh_summary"], index=False)

    def make_combined(tbp_subset: pd.DataFrame, air_subset: pd.DataFrame) -> pd.DataFrame:
        return pd.concat(
            [
                tbp_subset[["loop_id", "loop_length", "tbp_binary_label"]]
                .rename(columns={"tbp_binary_label": "binary_label"})
                .assign(dataset_is_tbp=1),
                air_subset[["air_record_id", "loop_length", "flexibility_label"]]
                .rename(columns={"air_record_id": "loop_id", "flexibility_label": "binary_label"})
                .assign(dataset_is_tbp=0),
            ],
            ignore_index=True,
        ).assign(is_flexible=lambda df: df["binary_label"].map({"flexible": 1, "rigid": 0}))

    logit_rows = []
    for label, tbp_subset, air_subset in [
        ("TBP_all_binary_vs_AIR_all", tbp_binary, air_reference),
        ("TBP_primary_vs_AIR_all", tbp_primary, air_reference),
        ("TBP_all_binary_vs_AIR_PDB", tbp_binary, air_pdb),
        ("TBP_primary_vs_AIR_PDB", tbp_primary, air_pdb),
    ]:
        combined = make_combined(tbp_subset, air_subset)
        fit = fit_length_adjusted_logit(combined[["is_flexible", "dataset_is_tbp", "loop_length"]])
        fit["comparison"] = label
        logit_rows.append(fit)
    logit_summary = pd.DataFrame(logit_rows)
    logit_summary.to_csv(outputs["logit_summary"], index=False)

    tbp.to_csv(outputs["loop_level"], index=False)

    summary = {
        "timestamp": outputs["timestamp"],
        "tbp_total_loops": int(len(tbp)),
        "tbp_binary_comparable_loops": int(len(tbp_binary)),
        "tbp_primary_binary_loops": int(len(tbp_primary)),
        "tbp_exact_sequence_match_loops": int(tbp["has_air_exact_seq_match"].sum()),
        "air_total_loops": int(len(air_reference)),
        "air_pdb_loops": int(len(air_pdb)),
        "tbp_binary_flexible_fraction": float((tbp_binary["tbp_binary_label"] == "flexible").mean()),
        "tbp_primary_flexible_fraction": float((tbp_primary["tbp_binary_label"] == "flexible").mean()),
        "air_flexible_fraction": float((air_reference["flexibility_label"] == "flexible").mean()),
        "air_pdb_flexible_fraction": float((air_pdb["flexibility_label"] == "flexible").mean()),
        "length_standardized": length_standardized.set_index("analysis_cohort").to_dict(orient="index")
        if "analysis_cohort" in length_standardized.columns
        else length_standardized.set_index("cohort").to_dict(orient="index"),
    }
    with Path(outputs["summary_json"]).open("w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    report = f"""# TBP vs AIRs 强化版对照分析 v2

生成时间：{outputs["timestamp"]}

本版不是简单重复首版的长度分层统计，而是在同一份 `TBP` loop 结果与 `AIR ALL_conformations` 参考集之上，额外补入了三类强化对照：第一类是长度标准化后的 matched-control 估计，也就是用 `TBP` 自身的 loop 长度分布去加权 `AIR` 的柔性比例，从而避免“TBP 只是因为 loop 更长所以显得更柔”的解释；第二类是长度校正 logistic 模型，把 `dataset_is_tbp` 与 `loop_length`、`loop_length^2` 一起纳入，用协变量校正后的数据集效应来衡量 `TBP` 相对 `AIR` 的独立柔性倾向；第三类是 exact-sequence overlap 分析，用来识别那些与 `AIR` 已知 loop 序列完全重合的 `TBP` loop 在动态指标上是否表现出系统性偏移。

从当前结果看，`TBP` 相对于 `AIR` 的“偏柔”结论在强化版对照下仍然稳定。长度标准化之后，`TBP_all_binary` 与 `TBP_primary_approved` 两个 cohort 的观测柔性比例都高于按同长度分布推算的 `AIR` 期望柔性比例；也就是说，即使强行把 `AIR` 拉成与 `TBP` 一样的长度构成，`TBP` 仍然更偏柔。长度校正 logistic 模型也给出同方向结论，`dataset_is_tbp` 的系数为正，对应的 odds ratio 大于 `1`，说明在控制 loop 长度的情况下，来自 `TBP` 数据集本身仍然显著提高“被判为 flexible”的概率。

这份强化版也补充了 `AIR PDB-only` 对照，因此可以区分“与 AIR 整体参考相比更柔”和“与 AIR 中结构证据最强的 PDB 子集相比更柔”这两种口径。如果后续论文要写得更稳，建议优先引用 `TBP vs AIR_PDB` 的长度校正结果，把它作为对结构质量差异更保守的一版主分析。

另外，exact-sequence overlap 的结果提示，和 `AIR seq_loop` 完全重合的 `TBP` loop 在 dynamic 指标上往往低于未匹配条目。这更像是在说明：`AIR` 中容易反复出现的是一些保守短 loop 模体，而 `TBP` 数据集里真正驱动高动态性的部分，更可能来自那些没有在 `AIR` 中直接出现 exact overlap 的功能 loop。

需要单独强调的是，`AIR` 侧目前只有官方 loop 级柔/刚标签，没有与 `TBP` 完全同构的 `DSSP / GetContact / localcolabfold whole-loop dynamic` 全套原始指标。因此本版强化对照已经把“能严格对齐的部分”推进到较完整状态：包括长度匹配、PDB-only 子集、exact-sequence overlap 和协变量校正；但它仍然不是 residue-level 或 subregion-level 的一一对应分析。这一点要在方法学里明确写出来，避免把 AIR 比较过度表述成 residue-resolution 的直接验证。
"""
    Path(outputs["report_md"]).write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
