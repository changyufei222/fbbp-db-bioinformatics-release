from __future__ import annotations

import copy
import math
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio.PDB import Model, PDBParser, Structure
from Bio.PDB.SASA import ShrakeRupley

try:
    from scipy.stats import mannwhitneyu
except Exception:  # pragma: no cover
    mannwhitneyu = None


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
ANNOTATION_INPUT = ROOT / "loop_subregion_annotations_input_v1.csv"
LOOP_RESULTS = ROOT / "loop_flexibility_results_long.csv"
DSSP_RESIDUE_LONG = PROJECT_ROOT / "loop分析最终" / "dssp_residue_long.csv"
WIDE_TABLE = PROJECT_ROOT / "总表" / "normalized" / "main_with_all_results_wide.csv"
DISTANCE_CUTOFF = 5.0
DELTA_ASA_SUPPORTED_THRESHOLD = 1.0
CORE_INTERFACE_THRESHOLD = 5.0
CORE_INTERFACE_STRICT_THRESHOLD = 10.0


def read_csv_fallback(path: Path) -> pd.DataFrame:
    for encoding in (None, "utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding, low_memory=False) if encoding else pd.read_csv(path, low_memory=False)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("csv", b"", 0, 1, f"Unable to decode {path}")


def parse_pdb_structure(path: Path | str) -> dict[str, dict[int, list[tuple[float, float, float]]]]:
    structure: dict[str, dict[int, list[tuple[float, float, float]]]] = {}
    p = Path(path)
    with p.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            chain_id = line[21].strip() or "_"
            try:
                resnum = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError:
                continue
            structure.setdefault(chain_id, {}).setdefault(resnum, []).append((x, y, z))
    return structure


def parse_pdb_model(path: Path | str):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("strict_full", str(path))
    return next(structure.get_models())


def subset_structure(model, chain_ids: Iterable[str], structure_id: str = "subset"):
    structure = Structure.Structure(structure_id)
    model_copy = Model.Model(0)
    structure.add(model_copy)
    for chain_id in chain_ids:
        if chain_id in model:
            model_copy.add(copy.deepcopy(model[chain_id]))
    return structure


def compute_chain_pair_delta_asa(path: Path | str, tbp_chain_id: str, target_chain_id: str) -> dict[int, dict[str, float]]:
    model = parse_pdb_model(path)
    mono = subset_structure(model, [tbp_chain_id], "mono_proxy")
    pair = subset_structure(model, [tbp_chain_id, target_chain_id], "pair_complex")

    sr = ShrakeRupley()
    sr.compute(mono, level="R")
    sr.compute(pair, level="R")

    mono_chain = mono[0][tbp_chain_id]
    pair_chain = pair[0][tbp_chain_id]
    pair_lookup = {rid: residue for rid, residue in pair_chain.child_dict.items() if rid[0] == " "}

    result: dict[int, dict[str, float]] = {}
    for rid, mono_residue in mono_chain.child_dict.items():
        if rid[0] != " " or rid not in pair_lookup:
            continue
        resnum = int(rid[1])
        mono_asa = float(getattr(mono_residue, "sasa", math.nan))
        pair_asa = float(getattr(pair_lookup[rid], "sasa", math.nan))
        result[resnum] = {
            "mono_proxy_asa": mono_asa,
            "complex_pair_asa": pair_asa,
            "delta_asa_proxy": mono_asa - pair_asa if not (math.isnan(mono_asa) or math.isnan(pair_asa)) else math.nan,
        }
    return result


def _min_distance(
    atoms_a: Iterable[tuple[float, float, float]],
    atoms_b: Iterable[tuple[float, float, float]],
) -> float:
    best = math.inf
    for ax, ay, az in atoms_a:
        for bx, by, bz in atoms_b:
            dist = math.dist((ax, ay, az), (bx, by, bz))
            if dist < best:
                best = dist
    return best


def _normalize_target_chain_ids(target_chain_ids: Iterable[str] | str | None) -> list[str]:
    if target_chain_ids is None:
        return []
    if isinstance(target_chain_ids, str):
        parts = [p.strip() for p in target_chain_ids.split("|")]
        return [p for p in parts if p]
    return [str(p).strip() for p in target_chain_ids if str(p).strip()]


def label_loop_residues_from_structure(
    loop_row: pd.Series,
    parsed_structure: dict[str, dict[int, list[tuple[float, float, float]]]],
    target_chain_ids: Iterable[str] | str | None,
    distance_cutoff: float = 5.0,
) -> list[dict[str, object]]:
    loop_chain_id = str(loop_row["loop_chain_id"]).strip()
    loop_chain = parsed_structure.get(loop_chain_id, {})
    loop_resnums = [r for r in sorted(loop_chain) if int(loop_row["loop_start"]) <= r <= int(loop_row["loop_end"])]
    candidates = _normalize_target_chain_ids(target_chain_ids)
    if not loop_resnums or not candidates:
        return []

    per_chain_scores: list[tuple[int, int, float, str]] = []
    residue_distances_by_chain: dict[str, dict[int, float]] = {}
    for chain_id in candidates:
        target_chain = parsed_structure.get(chain_id, {})
        if not target_chain:
            continue
        flat_target_atoms = [atom for atoms in target_chain.values() for atom in atoms]
        residue_distances: dict[int, float] = {}
        binding_residues = 0
        min_dist = math.inf
        for resnum in loop_resnums:
            dist = _min_distance(loop_chain[resnum], flat_target_atoms)
            residue_distances[resnum] = dist
            if dist <= distance_cutoff:
                binding_residues += 1
            min_dist = min(min_dist, dist)
        residue_distances_by_chain[chain_id] = residue_distances
        total_contacts = sum(1 for d in residue_distances.values() if d <= distance_cutoff)
        per_chain_scores.append((binding_residues, total_contacts, -min_dist, chain_id))

    if not per_chain_scores:
        return []

    selected_chain_id = sorted(per_chain_scores, reverse=True)[0][3]
    selected_distances = residue_distances_by_chain[selected_chain_id]
    rows: list[dict[str, object]] = []
    for resnum in loop_resnums:
        rows.append(
            {
                "residue_seq_num": resnum,
                "binding_label": 1 if selected_distances[resnum] <= distance_cutoff else 0,
                "selected_target_chain_id": selected_chain_id,
                "min_target_distance": selected_distances[resnum],
            }
        )
    return rows


def select_strict_annotation_rows(annotation_df: pd.DataFrame) -> pd.DataFrame:
    required = annotation_df.copy()
    for column in ["complex_structure_path", "tbp_chain_id", "target_chain_id"]:
        required = required.loc[required[column].notna() & required[column].astype(str).str.strip().ne("")]

    path_source = required["complex_structure_path"].astype(str)
    if "structure_path_hint" in required.columns:
        path_source = required["structure_path_hint"].fillna(required["complex_structure_path"]).astype(str)

    is_experimental = path_source.str.contains(r"PDB_Structures", case=False, regex=True)
    not_anomaly = ~path_source.str.contains(r"7BOF", case=False, regex=True)
    exists = required["complex_structure_path"].astype(str).map(lambda p: Path(p).exists())

    strict = required.loc[is_experimental & not_anomaly & exists].copy()
    return strict.reset_index(drop=True)


def _attach_scaffold(loop_df: pd.DataFrame, wide_df: pd.DataFrame | None) -> pd.DataFrame:
    if "Scaffold_Category" in loop_df.columns:
        return loop_df
    if wide_df is None:
        raise ValueError("wide_df is required when loop_df lacks Scaffold_Category")
    cols = ["protein_row_id", "structure_unique_sequence_id", "Scaffold_Category"]
    return loop_df.merge(wide_df[cols].drop_duplicates(), on=["protein_row_id", "structure_unique_sequence_id"], how="left")


def build_strict_residue_labels(
    loop_df: pd.DataFrame,
    annotation_df: pd.DataFrame,
    dssp_df: pd.DataFrame,
    wide_df: pd.DataFrame | None = None,
    distance_cutoff: float = 5.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    strict_annotations = select_strict_annotation_rows(annotation_df)
    loop_with_scaffold = _attach_scaffold(loop_df.copy(), wide_df)
    loop_map = loop_with_scaffold.set_index("loop_id")

    dssp_cols = ["protein_row_id", "structure_unique_sequence_id", "chain_id", "residue_seq_num", "aa", "ss", "rel_asa", "phi", "psi"]
    available_dssp_cols = [c for c in dssp_cols if c in dssp_df.columns]
    dssp_lookup = dssp_df[available_dssp_cols].copy()

    residue_rows: list[dict[str, object]] = []
    loop_summary_rows: list[dict[str, object]] = []
    parsed_cache: dict[str, dict[str, dict[int, list[tuple[float, float, float]]]]] = {}
    delta_asa_cache: dict[tuple[str, str, str], dict[int, dict[str, float]]] = {}

    for _, annotation_row in strict_annotations.iterrows():
        loop_id = annotation_row["loop_id"]
        if loop_id not in loop_map.index:
            continue
        loop_row = loop_map.loc[loop_id]
        if isinstance(loop_row, pd.DataFrame):
            loop_row = loop_row.iloc[0]

        complex_path = str(annotation_row["complex_structure_path"]).strip()
        if complex_path not in parsed_cache:
            parsed_cache[complex_path] = parse_pdb_structure(complex_path)

        loop_series = pd.Series(
            {
                "loop_id": loop_id,
                "loop_chain_id": annotation_row.get("tbp_chain_id", loop_row["loop_chain_id"]),
                "loop_start": loop_row["loop_start"],
                "loop_end": loop_row["loop_end"],
            }
        )
        residue_labels = label_loop_residues_from_structure(
            loop_row=loop_series,
            parsed_structure=parsed_cache[complex_path],
            target_chain_ids=annotation_row["target_chain_id"],
            distance_cutoff=distance_cutoff,
        )
        if not residue_labels:
            continue

        residue_df = pd.DataFrame(residue_labels)
        residue_df["protein_row_id"] = loop_row["protein_row_id"]
        residue_df["structure_unique_sequence_id"] = loop_row["structure_unique_sequence_id"]
        residue_df["loop_id"] = loop_id
        residue_df["Scaffold_Category"] = loop_row.get("Scaffold_Category", pd.NA)
        residue_df["tbp_chain_id"] = annotation_row["tbp_chain_id"]
        residue_df["target_chain_id"] = annotation_row["target_chain_id"]
        residue_df["complex_structure_path"] = complex_path
        residue_df["binding_annotation_source"] = annotation_row.get("binding_annotation_source", pd.NA)

        delta_key = (
            complex_path,
            str(annotation_row["tbp_chain_id"]).strip(),
            str(residue_df["selected_target_chain_id"].iloc[0]).strip(),
        )
        if delta_key not in delta_asa_cache:
            delta_asa_cache[delta_key] = compute_chain_pair_delta_asa(
                path=complex_path,
                tbp_chain_id=delta_key[1],
                target_chain_id=delta_key[2],
            )
        delta_lookup = delta_asa_cache[delta_key]
        residue_df["mono_proxy_asa"] = residue_df["residue_seq_num"].map(
            lambda n: delta_lookup.get(int(n), {}).get("mono_proxy_asa", math.nan)
        )
        residue_df["complex_pair_asa"] = residue_df["residue_seq_num"].map(
            lambda n: delta_lookup.get(int(n), {}).get("complex_pair_asa", math.nan)
        )
        residue_df["delta_asa_proxy"] = residue_df["residue_seq_num"].map(
            lambda n: delta_lookup.get(int(n), {}).get("delta_asa_proxy", math.nan)
        )
        residue_df["distance_defined_label"] = pd.to_numeric(residue_df["binding_label"], errors="coerce").fillna(0).astype(int)
        residue_df["delta_asa_supported_label"] = (
            (pd.to_numeric(residue_df["binding_label"], errors="coerce") == 1)
            & (pd.to_numeric(residue_df["delta_asa_proxy"], errors="coerce") > DELTA_ASA_SUPPORTED_THRESHOLD)
        ).astype(int)
        residue_df["core_interface_label"] = (
            (pd.to_numeric(residue_df["binding_label"], errors="coerce") == 1)
            & (pd.to_numeric(residue_df["delta_asa_proxy"], errors="coerce") >= CORE_INTERFACE_THRESHOLD)
        ).astype(int)
        residue_df["core_interface_strict_label"] = (
            (pd.to_numeric(residue_df["binding_label"], errors="coerce") == 1)
            & (pd.to_numeric(residue_df["delta_asa_proxy"], errors="coerce") >= CORE_INTERFACE_STRICT_THRESHOLD)
        ).astype(int)

        residue_df = residue_df.merge(
            dssp_lookup,
            left_on=["protein_row_id", "structure_unique_sequence_id", "tbp_chain_id", "residue_seq_num"],
            right_on=["protein_row_id", "structure_unique_sequence_id", "chain_id", "residue_seq_num"],
            how="left",
        )
        residue_df = residue_df.drop(columns=["chain_id"], errors="ignore")
        residue_rows.extend(residue_df.to_dict(orient="records"))

        binding_mask = residue_df["binding_label"] == 1
        loop_summary_rows.append(
            {
                "protein_row_id": loop_row["protein_row_id"],
                "structure_unique_sequence_id": loop_row["structure_unique_sequence_id"],
                "loop_id": loop_id,
                "Scaffold_Category": loop_row.get("Scaffold_Category", pd.NA),
                "complex_structure_path": complex_path,
                "tbp_chain_id": annotation_row["tbp_chain_id"],
                "target_chain_id": annotation_row["target_chain_id"],
                "selected_target_chain_id": residue_df["selected_target_chain_id"].iloc[0],
                "loop_start": loop_row["loop_start"],
                "loop_end": loop_row["loop_end"],
                "loop_length": loop_row.get("loop_length", len(residue_df)),
                "binding_residue_count": int(binding_mask.sum()),
                "nonbinding_residue_count": int((~binding_mask).sum()),
                "binding_fraction": float(binding_mask.mean()),
                "mean_binding_rel_asa": pd.to_numeric(residue_df.loc[binding_mask, "rel_asa"], errors="coerce").mean(),
                "mean_nonbinding_rel_asa": pd.to_numeric(residue_df.loc[~binding_mask, "rel_asa"], errors="coerce").mean(),
                "mean_binding_delta_asa_proxy": pd.to_numeric(residue_df.loc[binding_mask, "delta_asa_proxy"], errors="coerce").mean(),
                "mean_nonbinding_delta_asa_proxy": pd.to_numeric(residue_df.loc[~binding_mask, "delta_asa_proxy"], errors="coerce").mean(),
                "delta_asa_supported_residue_count": int(residue_df["delta_asa_supported_label"].sum()),
                "core_interface_residue_count": int(residue_df["core_interface_label"].sum()),
                "core_interface_strict_residue_count": int(residue_df["core_interface_strict_label"].sum()),
                "dssp_matched_residue_count": int(residue_df["rel_asa"].notna().sum()) if "rel_asa" in residue_df.columns else 0,
            }
        )

    return pd.DataFrame(residue_rows), pd.DataFrame(loop_summary_rows)


def _metric_stats(binding: pd.Series, nonbinding: pd.Series) -> dict[str, object]:
    result = {
        "n_binding": int(binding.notna().sum()),
        "n_nonbinding": int(nonbinding.notna().sum()),
        "binding_mean": pd.to_numeric(binding, errors="coerce").mean(),
        "nonbinding_mean": pd.to_numeric(nonbinding, errors="coerce").mean(),
        "binding_median": pd.to_numeric(binding, errors="coerce").median(),
        "nonbinding_median": pd.to_numeric(nonbinding, errors="coerce").median(),
        "delta_mean": pd.to_numeric(binding, errors="coerce").mean() - pd.to_numeric(nonbinding, errors="coerce").mean(),
        "delta_median": pd.to_numeric(binding, errors="coerce").median() - pd.to_numeric(nonbinding, errors="coerce").median(),
        "mannwhitney_p": pd.NA,
    }
    clean_binding = pd.to_numeric(binding, errors="coerce").dropna()
    clean_nonbinding = pd.to_numeric(nonbinding, errors="coerce").dropna()
    if mannwhitneyu is not None and len(clean_binding) > 0 and len(clean_nonbinding) > 0:
        result["mannwhitney_p"] = float(mannwhitneyu(clean_binding, clean_nonbinding, alternative="two-sided").pvalue)
    return result


def build_residue_level_stats(residue_df: pd.DataFrame) -> pd.DataFrame:
    if residue_df.empty:
        return pd.DataFrame()

    usable = residue_df.loc[pd.to_numeric(residue_df["binding_label"], errors="coerce").isin([0, 1])].copy()
    if usable.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    groups = {"overall_strict": usable}
    for scaffold, group_df in usable.groupby("Scaffold_Category", dropna=False):
        groups[str(scaffold)] = group_df

    metric_specs = [
        ("distance_defined_label", "rel_asa"),
        ("distance_defined_label", "delta_asa_proxy"),
        ("delta_asa_supported_label", "rel_asa"),
        ("core_interface_label", "rel_asa"),
        ("core_interface_strict_label", "rel_asa"),
    ]

    for group_name, group_df in groups.items():
        for label_column, metric_column in metric_specs:
            if label_column not in group_df.columns or metric_column not in group_df.columns:
                continue
            binding = group_df.loc[pd.to_numeric(group_df[label_column], errors="coerce") == 1, metric_column]
            nonbinding = group_df.loc[pd.to_numeric(group_df[label_column], errors="coerce") == 0, metric_column]
            if binding.dropna().empty or nonbinding.dropna().empty:
                continue
            row = {"group": group_name, "metric": metric_column, "label_layer": label_column}
            row.update(_metric_stats(binding, nonbinding))
            rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(
    residue_df: pd.DataFrame,
    loop_summary_df: pd.DataFrame,
    residue_stats_df: pd.DataFrame,
    output_prefix: str = "strict_v1",
) -> dict[str, Path]:
    residue_path = ROOT / f"binding_residue_labels_{output_prefix}.csv"
    summary_path = ROOT / f"binding_loop_residue_summary_{output_prefix}.csv"
    stats_path = ROOT / f"binding_vs_nonbinding_residue_stats_{output_prefix}.csv"
    report_path = ROOT / f"binding_residue_report_{output_prefix}.md"

    residue_df.to_csv(residue_path, index=False, encoding="utf-8")
    loop_summary_df.to_csv(summary_path, index=False, encoding="utf-8")
    residue_stats_df.to_csv(stats_path, index=False, encoding="utf-8")

    n_loops = int(loop_summary_df["loop_id"].nunique()) if not loop_summary_df.empty else 0
    n_structures = int(loop_summary_df["complex_structure_path"].nunique()) if not loop_summary_df.empty else 0
    n_binding_residues = int((pd.to_numeric(residue_df.get("binding_label"), errors="coerce") == 1).sum()) if not residue_df.empty else 0
    n_nonbinding_residues = int((pd.to_numeric(residue_df.get("binding_label"), errors="coerce") == 0).sum()) if not residue_df.empty else 0
    dssp_coverage = int(residue_df["rel_asa"].notna().sum()) if "rel_asa" in residue_df.columns and not residue_df.empty else 0
    delta_asa_supported_count = int(pd.to_numeric(residue_df.get("delta_asa_supported_label"), errors="coerce").fillna(0).sum()) if not residue_df.empty else 0
    core_interface_count = int(pd.to_numeric(residue_df.get("core_interface_label"), errors="coerce").fillna(0).sum()) if not residue_df.empty else 0
    core_interface_strict_count = int(pd.to_numeric(residue_df.get("core_interface_strict_label"), errors="coerce").fillna(0).sum()) if not residue_df.empty else 0

    lines = [
        "# Strict Residue-Level Binding Report",
        "",
        "This strict result only includes loops with explicit `target_chain_id`, existing experimental multichain structures, and successful residue-level contact extraction.",
        "",
        f"- Eligible loops with residue-level labels: `{n_loops}`",
        f"- Eligible structures: `{n_structures}`",
        f"- Binding residues: `{n_binding_residues}`",
        f"- Non-binding residues: `{n_nonbinding_residues}`",
        f"- Residues with DSSP rel_asa matched: `{dssp_coverage}`",
        f"- ΔASA-supported residues (`ΔASA_proxy > {DELTA_ASA_SUPPORTED_THRESHOLD:.1f} Å²`): `{delta_asa_supported_count}`",
        f"- Core-interface residues (`ΔASA_proxy >= {CORE_INTERFACE_THRESHOLD:.1f} Å²`): `{core_interface_count}`",
        f"- Core-interface strict residues (`ΔASA_proxy >= {CORE_INTERFACE_STRICT_THRESHOLD:.1f} Å²`): `{core_interface_strict_count}`",
        "",
        "## Output Files",
        f"- `{residue_path}`",
        f"- `{summary_path}`",
        f"- `{stats_path}`",
        "",
        "## Overall Strict Statistics",
    ]
    overall = residue_stats_df.loc[residue_stats_df["group"] == "overall_strict"]
    if overall.empty:
        lines.append("- No residue-level statistics available.")
    else:
        for _, row in overall.iterrows():
            lines.append(
                f"- `{row['label_layer']} / {row['metric']}`: binding_mean=`{row['binding_mean']}`, nonbinding_mean=`{row['nonbinding_mean']}`, delta_mean=`{row['delta_mean']}`, p=`{row['mannwhitney_p']}`"
            )
    report_path.write_text("\n".join(lines), encoding="utf-8")

    return {
        "residue_labels": residue_path,
        "loop_summary": summary_path,
        "residue_stats": stats_path,
        "report": report_path,
    }


def main() -> None:
    annotation_df = read_csv_fallback(ANNOTATION_INPUT)
    loop_df = read_csv_fallback(LOOP_RESULTS)
    dssp_df = read_csv_fallback(DSSP_RESIDUE_LONG)
    wide_df = read_csv_fallback(WIDE_TABLE)

    residue_df, loop_summary_df = build_strict_residue_labels(
        loop_df=loop_df,
        annotation_df=annotation_df,
        dssp_df=dssp_df,
        wide_df=wide_df,
    )
    residue_stats_df = build_residue_level_stats(residue_df)
    write_outputs(residue_df, loop_summary_df, residue_stats_df)


if __name__ == "__main__":
    main()
