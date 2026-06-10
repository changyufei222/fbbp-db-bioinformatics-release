from __future__ import annotations

import json
from math import erfc, exp, sqrt
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"<local_path_removed>")
OUT_DIR = ROOT / "论文写作" / "BMC修订_20260531"
TBP_ENHANCED = ROOT / "柔性分析最终" / "tbp_vs_airs_enhanced_loop_level_v2.csv"
STRUCTURAL_INFO = ROOT / "FBBP_DB_v1_frozen_release_20260531" / "data" / "schema_tables" / "structural_info.csv"
AIR_FLEX = ROOT / "柔性分析最终" / "ALL_conformations_extracted" / "ALL_conformations" / "ALL_conformations_flexibility.csv"
AIR_PDBS = ROOT / "柔性分析最终" / "ALL_conformations_extracted" / "ALL_conformations" / "ALL_conformations_PDBs.csv"

RESOLUTION_THRESHOLD = 3.0
IDENTITY_THRESHOLD = 0.90
N_BOOT = 1000
RANDOM_SEED = 20260531


def normalize_seq(value: object) -> str | None:
    if pd.isna(value):
        return None
    seq = str(value).strip().upper()
    return seq if seq else None


def sequence_identity(a: str, b: str) -> float:
    if len(a) != len(b) or not a:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


def cluster_by_length_identity(df: pd.DataFrame, seq_col: str, label_col: str, prefix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = df.copy()
    rows[seq_col] = rows[seq_col].map(normalize_seq)
    rows = rows.dropna(subset=[seq_col, "loop_length", label_col]).copy()
    rows["is_flexible"] = rows[label_col].map({"flexible": 1.0, "rigid": 0.0})
    rows = rows.dropna(subset=["is_flexible"]).sort_values(["loop_length", seq_col], kind="stable")

    assigned: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    next_id = 1
    for length, group in rows.groupby("loop_length", sort=True):
        representatives: list[dict[str, object]] = []
        for _, row in group.iterrows():
            seq = row[seq_col]
            match = None
            for rep in representatives:
                if sequence_identity(seq, rep["representative_sequence"]) >= IDENTITY_THRESHOLD:
                    match = rep
                    break
            if match is None:
                match = {
                    "cluster_id": f"{prefix}_CL{next_id:06d}",
                    "loop_length": int(length),
                    "representative_sequence": seq,
                    "members": [],
                }
                representatives.append(match)
                next_id += 1
            match["members"].append(row)
            assigned.append(
                {
                    "cluster_id": match["cluster_id"],
                    "source_sequence": seq,
                    "loop_length": int(length),
                    "source_label": row[label_col],
                    "source_is_flexible": float(row["is_flexible"]),
                }
            )

        for rep in representatives:
            member_df = pd.DataFrame(rep["members"])
            flexible_fraction = float(member_df["is_flexible"].mean())
            n_members = int(len(member_df))
            cluster_rows.append(
                {
                    "cluster_id": rep["cluster_id"],
                    "loop_length": int(length),
                    "representative_sequence": rep["representative_sequence"],
                    "n_members": n_members,
                    "cluster_flexible_fraction": flexible_fraction,
                    "cluster_binary_label": "flexible" if flexible_fraction >= 0.5 else "rigid",
                    "mixed_label_cluster": int(member_df["is_flexible"].nunique(dropna=True) > 1),
                }
            )

    return pd.DataFrame(cluster_rows), pd.DataFrame(assigned)


def load_tbp_strict() -> pd.DataFrame:
    tbp = pd.read_csv(
        TBP_ENHANCED,
        usecols=[
            "protein_row_id",
            "loop_id",
            "loop_sequence",
            "loop_length",
            "selected_structure_source",
            "tbp_binary_label",
            "analysis_cohort",
            "highest_priority_cohort",
        ],
    )
    tbp["domain_id"] = tbp["protein_row_id"].map(lambda x: f"TBP-{int(x):05d}" if pd.notna(x) else pd.NA)

    struct = pd.read_csv(STRUCTURAL_INFO)
    struct["resolution"] = pd.to_numeric(struct["resolution"], errors="coerce")
    struct = struct.rename(columns={"method": "structure_method", "resolution": "structure_resolution"})

    merged = tbp.merge(
        struct[["domain_id", "pdb_id", "chain_id", "structure_method", "structure_resolution"]],
        on="domain_id",
        how="left",
    )
    strict = merged.loc[
        (merged["selected_structure_source"] == "experimental_pdb")
        & (merged["structure_method"].astype(str).str.contains("X-RAY", case=False, na=False))
        & (merged["structure_resolution"].le(RESOLUTION_THRESHOLD))
        & (merged["tbp_binary_label"].isin(["flexible", "rigid"]))
    ].copy()
    strict["loop_sequence"] = strict["loop_sequence"].map(normalize_seq)
    strict = strict.dropna(subset=["loop_sequence", "loop_length"]).drop_duplicates("loop_id")
    return strict


def load_air_strict(tbp_sequences: set[str]) -> pd.DataFrame:
    air_flex = pd.read_csv(
        AIR_FLEX,
        usecols=["seq_loop", "loop_length", "flexibility_label", "set", "n_structures", "n_pdbs"],
    )
    air_flex = air_flex.loc[
        (air_flex["set"] == "PDB") & (air_flex["flexibility_label"].isin(["flexible", "rigid"]))
    ].copy()
    air_flex["seq_loop"] = air_flex["seq_loop"].map(normalize_seq)
    air_flex = air_flex.dropna(subset=["seq_loop", "loop_length"])

    air_pdbs = pd.read_csv(
        AIR_PDBS,
        usecols=["pdb", "chain", "seq_loop", "loop_length", "structure_method", "resolution"],
    )
    air_pdbs["resolution"] = pd.to_numeric(air_pdbs["resolution"], errors="coerce")
    air_pdbs["seq_loop"] = air_pdbs["seq_loop"].map(normalize_seq)
    air_pdbs = air_pdbs.loc[
        air_pdbs["structure_method"].astype(str).str.contains("x-ray", case=False, na=False)
        & air_pdbs["resolution"].le(RESOLUTION_THRESHOLD)
    ].dropna(subset=["seq_loop", "loop_length"])

    high_quality_sequences = air_pdbs[["seq_loop", "loop_length"]].drop_duplicates()
    strict = air_flex.merge(high_quality_sequences, on=["seq_loop", "loop_length"], how="inner")
    strict = strict.loc[~strict["seq_loop"].isin(tbp_sequences)].drop_duplicates(["seq_loop", "loop_length"])
    return strict


def length_matched_bootstrap(tbp_clusters: pd.DataFrame, air_clusters: pd.DataFrame) -> dict[str, float | int]:
    rng = np.random.default_rng(RANDOM_SEED)
    tbp = tbp_clusters.copy()
    air = air_clusters.copy()
    air_by_length = {
        int(length): group["cluster_flexible_fraction"].to_numpy(dtype=float)
        for length, group in air.groupby("loop_length")
    }

    usable_mask = tbp["loop_length"].map(lambda length: int(length) in air_by_length)
    usable_tbp = tbp.loc[usable_mask].copy()
    missing_tbp = int((~usable_mask).sum())
    observed = float(usable_tbp["cluster_flexible_fraction"].mean()) if len(usable_tbp) else np.nan

    boot_values = []
    lengths = usable_tbp["loop_length"].astype(int).to_numpy()
    for _ in range(N_BOOT):
        sampled = [rng.choice(air_by_length[int(length)]) for length in lengths]
        boot_values.append(float(np.mean(sampled)) if sampled else np.nan)
    boot = np.asarray(boot_values, dtype=float)
    boot = boot[~np.isnan(boot)]

    expected = float(np.mean(boot)) if len(boot) else np.nan
    p_one_sided = float((np.sum(boot >= observed) + 1) / (len(boot) + 1)) if len(boot) and pd.notna(observed) else np.nan
    return {
        "n_fbbp_clusters_length_matched": int(len(usable_tbp)),
        "n_fbbp_clusters_missing_air_length": missing_tbp,
        "observed_fbbp_flexible_fraction": observed,
        "strict_air_length_matched_expected_fraction": expected,
        "observed_minus_expected": float(observed - expected) if pd.notna(observed) and pd.notna(expected) else np.nan,
        "bootstrap_ci_lower": float(np.quantile(boot, 0.025)) if len(boot) else np.nan,
        "bootstrap_ci_upper": float(np.quantile(boot, 0.975)) if len(boot) else np.nan,
        "bootstrap_p_one_sided_air_ge_fbbp": p_one_sided,
        "n_bootstrap": int(len(boot)),
    }


def mantel_haenszel_or(tbp_clusters: pd.DataFrame, air_clusters: pd.DataFrame) -> dict[str, float | int]:
    numerator = 0.0
    denominator = 0.0
    strata = 0
    common_lengths = sorted(set(tbp_clusters["loop_length"].astype(int)) & set(air_clusters["loop_length"].astype(int)))
    for length in common_lengths:
        t = tbp_clusters.loc[tbp_clusters["loop_length"].astype(int) == length]
        a = float((t["cluster_binary_label"] == "flexible").sum())
        b = float((t["cluster_binary_label"] == "rigid").sum())
        r = air_clusters.loc[air_clusters["loop_length"].astype(int) == length]
        c = float((r["cluster_binary_label"] == "flexible").sum())
        d = float((r["cluster_binary_label"] == "rigid").sum())
        n = a + b + c + d
        if min(a + b, c + d, a + c, b + d) == 0 or n == 0:
            continue
        numerator += (a * d) / n
        denominator += (b * c) / n
        strata += 1
    return {
        "mh_strata_used": strata,
        "mantel_haenszel_or": float(numerator / denominator) if denominator else np.nan,
    }


def fit_length_adjusted_logit(tbp_clusters: pd.DataFrame, air_clusters: pd.DataFrame) -> dict[str, float | int]:
    tbp = tbp_clusters[["cluster_binary_label", "loop_length"]].copy()
    tbp["dataset_is_fbbp"] = 1.0
    air = air_clusters[["cluster_binary_label", "loop_length"]].copy()
    air["dataset_is_fbbp"] = 0.0
    combined = pd.concat([tbp, air], ignore_index=True)
    combined["is_flexible"] = combined["cluster_binary_label"].map({"flexible": 1.0, "rigid": 0.0})
    frame = combined.dropna(subset=["is_flexible", "loop_length", "dataset_is_fbbp"]).copy()
    y = frame["is_flexible"].to_numpy(dtype=float)
    length = frame["loop_length"].to_numpy(dtype=float)
    length_scaled = (length - length.mean()) / (length.std() or 1.0)
    x = np.column_stack([np.ones(len(frame)), frame["dataset_is_fbbp"].to_numpy(dtype=float), length_scaled, length_scaled**2])

    beta = np.zeros(x.shape[1], dtype=float)
    converged = False
    ridge = 1e-8
    for _ in range(100):
        eta = np.clip(x @ beta, -35, 35)
        p = 1.0 / (1.0 + np.exp(-eta))
        w = np.clip(p * (1.0 - p), 1e-8, None)
        z = eta + (y - p) / w
        xtw = x.T * w
        hessian = xtw @ x
        hessian.flat[:: hessian.shape[0] + 1] += ridge
        beta_new = np.linalg.pinv(hessian) @ (xtw @ z)
        if np.max(np.abs(beta_new - beta)) < 1e-7:
            beta = beta_new
            converged = True
            break
        beta = beta_new
    eta = np.clip(x @ beta, -35, 35)
    p = 1.0 / (1.0 + np.exp(-eta))
    w = np.clip(p * (1.0 - p), 1e-8, None)
    hessian = (x.T * w) @ x
    hessian.flat[:: hessian.shape[0] + 1] += ridge
    cov = np.linalg.pinv(hessian)
    se = np.sqrt(np.diag(cov))
    coef = float(beta[1])
    coef_se = float(se[1]) if len(se) > 1 else np.nan
    z_score = coef / coef_se if coef_se else np.nan
    p_value = float(erfc(abs(z_score) / sqrt(2))) if pd.notna(z_score) else np.nan
    return {
        "logit_n_obs": int(len(frame)),
        "length_adjusted_dataset_or": float(exp(coef)),
        "length_adjusted_or_ci_lower": float(exp(coef - 1.96 * coef_se)) if pd.notna(coef_se) else np.nan,
        "length_adjusted_or_ci_upper": float(exp(coef + 1.96 * coef_se)) if pd.notna(coef_se) else np.nan,
        "length_adjusted_p_value": p_value,
        "logit_converged": int(converged),
    }


def fmt(value: float, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}f}"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tbp_strict = load_tbp_strict()
    tbp_sequences = set(tbp_strict["loop_sequence"].dropna().astype(str))
    air_strict = load_air_strict(tbp_sequences)

    tbp_clusters, tbp_cluster_members = cluster_by_length_identity(
        tbp_strict.rename(columns={"loop_sequence": "seq"}),
        seq_col="seq",
        label_col="tbp_binary_label",
        prefix="FBBP",
    )
    air_clusters, air_cluster_members = cluster_by_length_identity(
        air_strict.rename(columns={"seq_loop": "seq"}),
        seq_col="seq",
        label_col="flexibility_label",
        prefix="AIR",
    )

    bootstrap = length_matched_bootstrap(tbp_clusters, air_clusters)
    mh = mantel_haenszel_or(tbp_clusters, air_clusters)
    logit = fit_length_adjusted_logit(tbp_clusters, air_clusters)

    summary = {
        "analysis": "strict_combined_background_pdb_resolution_dedup_length_matched",
        "resolution_threshold_angstrom": RESOLUTION_THRESHOLD,
        "sequence_identity_dedup_threshold": IDENTITY_THRESHOLD,
        "exact_fbbp_air_sequence_overlap_removed_from_air": 1,
        "n_fbbp_records_after_pdb_xray_resolution_filter": int(len(tbp_strict)),
        "n_fbbp_sequence_clusters_90id": int(len(tbp_clusters)),
        "n_air_records_after_pdb_xray_resolution_filter": int(len(air_strict)),
        "n_air_sequence_clusters_90id": int(len(air_clusters)),
        "fbbp_cluster_mixed_label_count": int(tbp_clusters["mixed_label_cluster"].sum()) if len(tbp_clusters) else 0,
        "air_cluster_mixed_label_count": int(air_clusters["mixed_label_cluster"].sum()) if len(air_clusters) else 0,
        **bootstrap,
        **mh,
        **logit,
    }

    summary_df = pd.DataFrame([summary])
    summary_path = OUT_DIR / "strict_combined_background_summary_20260531.csv"
    tbp_cluster_path = OUT_DIR / "strict_combined_background_fbbp_clusters_20260531.csv"
    air_cluster_path = OUT_DIR / "strict_combined_background_air_clusters_20260531.csv"
    member_path = OUT_DIR / "strict_combined_background_cluster_members_20260531.csv"
    report_path = OUT_DIR / "strict_combined_background_report_20260531.md"

    summary_df.to_csv(summary_path, index=False)
    tbp_clusters.to_csv(tbp_cluster_path, index=False)
    air_clusters.to_csv(air_cluster_path, index=False)
    pd.concat(
        [
            tbp_cluster_members.assign(dataset="FBBP"),
            air_cluster_members.assign(dataset="AIR"),
        ],
        ignore_index=True,
    ).to_csv(member_path, index=False)

    methods_text = (
        "Strict combined background analysis used only experimental PDB/X-ray structures with resolution "
        f"<= {RESOLUTION_THRESHOLD:.1f} A, collapsed loop sequences into greedy same-length sequence clusters at "
        f">= {IDENTITY_THRESHOLD:.0%} identity, removed exact FBBP-loop sequences from the AIR background, and "
        "then estimated the AIR expectation by matching the AIR clusters to the FBBP cluster length distribution."
    )
    result_text = (
        "After this combined restriction, the FBBP strict subset retained "
        f"{summary['n_fbbp_records_after_pdb_xray_resolution_filter']} records and "
        f"{summary['n_fbbp_sequence_clusters_90id']} sequence clusters; the AIR background retained "
        f"{summary['n_air_records_after_pdb_xray_resolution_filter']} records and "
        f"{summary['n_air_sequence_clusters_90id']} clusters. The FBBP cluster-level ITsFlexible-label flexible fraction "
        f"was {fmt(summary['observed_fbbp_flexible_fraction'])}, compared with a length-matched strict AIR expectation of "
        f"{fmt(summary['strict_air_length_matched_expected_fraction'])} "
        f"(bootstrap 95% interval {fmt(summary['bootstrap_ci_lower'])}-{fmt(summary['bootstrap_ci_upper'])}; "
        f"one-sided resampling P={fmt(summary['bootstrap_p_one_sided_air_ge_fbbp'], 4)}). The length-stratified "
        f"Mantel-Haenszel OR was {fmt(summary['mantel_haenszel_or'], 2)}, and the length-adjusted logistic OR was "
        f"{fmt(summary['length_adjusted_dataset_or'], 2)} "
        f"({fmt(summary['length_adjusted_or_ci_lower'], 2)}-{fmt(summary['length_adjusted_or_ci_upper'], 2)})."
    )

    report_path.write_text(
        "\n".join(
            [
                "# Strict combined background analysis",
                "",
                methods_text,
                "",
                result_text,
                "",
                "## Machine-readable summary",
                "",
                "```json",
                json.dumps(summary, ensure_ascii=False, indent=2),
                "```",
            ]
        ),
        encoding="utf-8",
    )
    print(result_text)
    print(f"Wrote: {summary_path}")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    main()
