from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Iterable

import pandas as pd

try:
    from scipy.stats import mannwhitneyu
except Exception:  # pragma: no cover - fallback if scipy is unavailable
    mannwhitneyu = None


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parent
PRIOR_TABLE = ROOT / "2026-04-17_scaffold_binding_prior_definitions.csv"
LOOP_RESULTS = ROOT / "loop_flexibility_results_long.csv"
CANDIDATE_TABLE = ROOT / "binding_nonbinding_complex_candidates_v1.csv"
ANNOTATION_INPUT = ROOT / "loop_subregion_annotations_input_v1.csv"
MANUAL_FINAL_REVIEW = ROOT / "binding_subregion_manual_review_primary_binding_all_v1_reviewed_round2_with_ibody_7BOF_remap.csv"
WIDE_TABLE = PROJECT_ROOT / "总表" / "normalized" / "main_with_all_results_wide.csv"
GETCONTACT_STRUCTURES = ROOT / "getcontact_stage_v1" / "structures"
PRIMARY_LOOPCENTRIC = {"adnectin", "affimer", "cyclotide", "ibody", "knottin", "kunitz"}
SECONDARY_HYBRID = {"obody", "centyrin"}

TEXT_COLUMNS_DEFAULT = (
    "Sources_title",
    "Domains_domain_name",
    "Domains_evidence_quote",
    "Proteins_description",
    "manual_notes",
    "Notes",
)

SCAFFOLD_PATTERN_RULES: dict[str, list[tuple[re.Pattern[str], set[int]]]] = {
    "adnectin": [
        (re.compile(r"\bbc/?fg\b|\bbc.*fg\b|\bfg.*bc\b"), {1, 3}),
        (re.compile(r"\bfg loop\b|\bfg-loop\b"), {3}),
        (re.compile(r"\bbc loop\b|\bbc-loop\b"), {1}),
        (re.compile(r"\bde loop\b|\bde-loop\b"), {2}),
    ],
    "affimer": [
        (re.compile(r"\bvr1\b|\bloop ?1\b"), {1}),
        (re.compile(r"\bvr2\b|\bloop ?2\b"), {2}),
    ],
    "knottin": [
        (re.compile(r"\bloop ?1\b"), {1}),
        (re.compile(r"\bloop ?3\b"), {3}),
        (re.compile(r"\bloop ?4\b"), {4}),
        (re.compile(r"\btrypsin loop\b|\breactive loop\b"), {1}),
    ],
    "cyclotide": [
        (re.compile(r"\bloop ?1\b"), {1}),
        (re.compile(r"\bloop ?3\b"), {3}),
        (re.compile(r"\bloop ?6\b"), {6}),
    ],
    "ibody": [
        (re.compile(r"\bcdr1\b"), {1}),
        (re.compile(r"\bcdr3\b"), {3}),
    ],
    "obody": [
        (re.compile(r"\bl4\b"), {4}),
    ],
    "kunitz": [
        (re.compile(r"\breactive loop\b|\binhibitory loop\b|\bp1/?p1'?\b"), {1}),
    ],
}


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
                "loop_id": loop_row["loop_id"],
                "chain_id": loop_chain_id,
                "resnum": resnum,
                "binding_label": 1 if selected_distances[resnum] <= distance_cutoff else 0,
                "selected_target_chain_id": selected_chain_id,
                "min_target_distance": selected_distances[resnum],
                "label_source": "direct_structure_contact",
            }
        )
    return rows


def aggregate_document_text(row: pd.Series, text_columns: Iterable[str] = TEXT_COLUMNS_DEFAULT) -> str:
    parts: list[str] = []
    for col in text_columns:
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            parts.append(str(value).strip())
    return " | ".join(parts).lower()


def build_empirical_loop_priors(high_evidence_df: pd.DataFrame) -> dict[str, dict[tuple[object, object], dict[str, float]]]:
    if high_evidence_df.empty:
        return {"slot_rates": {}, "flank_rates": {}}

    df = high_evidence_df.copy()
    df["scaffold_key"] = df["Scaffold_Category"].map(_normalize_scaffold_name)

    slot_rates: dict[tuple[object, object], dict[str, float]] = {}
    flank_rates: dict[tuple[object, object], dict[str, float]] = {}

    for (scaffold, slot), group in df.groupby(["scaffold_key", "loop_slot"], dropna=False):
        slot_rates[(scaffold, int(slot))] = {
            "binding_rate": float(pd.to_numeric(group["binding_loop_label"], errors="coerce").mean()),
            "n": int(len(group)),
        }

    for (scaffold, flank), group in df.groupby(["scaffold_key", "loop_flank_pattern"], dropna=False):
        flank_rates[(scaffold, str(flank))] = {
            "binding_rate": float(pd.to_numeric(group["binding_loop_label"], errors="coerce").mean()),
            "n": int(len(group)),
        }

    return {"slot_rates": slot_rates, "flank_rates": flank_rates}


def infer_prior_loop_label(
    loop_row: pd.Series,
    wide_row: pd.Series,
    prior_row: pd.Series,
    active_loop_count: int,
    empirical_priors: dict[str, dict[tuple[object, object], dict[str, float]]] | None = None,
) -> dict[str, object]:
    scaffold = str(wide_row.get("Scaffold_Category", prior_row.get("scaffold", ""))).strip().lower()
    text = aggregate_document_text(wide_row)
    slot = int(loop_row.get("loop_slot", 0) or 0)
    architecture = str(prior_row.get("binding_architecture_type", "")).strip().lower()

    score = 0.0
    source = "prior_default"
    evidence = "family_default"

    for pattern, slots in SCAFFOLD_PATTERN_RULES.get(scaffold, []):
        if pattern.search(text):
            score = 0.95 if slot in slots else 0.05
            source = "prior_document_override"
            evidence = pattern.pattern
            break

    if score == 0.0:
        empirical_components: list[tuple[float, int]] = []
        if empirical_priors:
            slot_prior = empirical_priors.get("slot_rates", {}).get((scaffold, slot))
            flank_prior = empirical_priors.get("flank_rates", {}).get((scaffold, str(loop_row.get("loop_flank_pattern", ""))))
            if slot_prior and slot_prior["n"] >= 2:
                empirical_components.append((float(slot_prior["binding_rate"]), int(slot_prior["n"])))
            if flank_prior and flank_prior["n"] >= 2:
                empirical_components.append((float(flank_prior["binding_rate"]), int(flank_prior["n"])))
        if empirical_components:
            total_weight = sum(weight for _, weight in empirical_components)
            empirical_score = sum(value * weight for value, weight in empirical_components) / total_weight
            if total_weight >= 4 and abs(empirical_score - 0.5) >= 0.2:
                confidence = "high" if total_weight >= 10 else "medium"
                return {
                    "binding_loop_score": round(float(empirical_score), 4),
                    "binding_loop_label": 1 if empirical_score >= 0.5 else 0,
                    "label_source": "prior_empirical_calibrated",
                    "evidence_type": "empirical_slot_flank_prior",
                    "confidence_tier": confidence,
                }

        if architecture in {"loop-centric", "reactive-loop-centric"}:
            if active_loop_count == 1:
                score = 0.8
            elif scaffold == "affimer":
                score = 0.7 if slot in {1, 2} else 0.2
            elif scaffold == "adnectin":
                score = 0.6 if slot in {1, 3} else 0.35 if slot == 2 else 0.1
            elif scaffold in {"knottin", "cyclotide"}:
                score = 0.55 if slot in {1, 3, 4, 6} else 0.2
            elif scaffold == "ibody":
                score = 0.65 if slot in {1, 3} else 0.2
            elif scaffold == "kunitz":
                score = 0.65 if slot == 1 else 0.2
            else:
                score = 0.55 if slot <= max(active_loop_count, 1) else 0.1
        elif architecture in {"surface-plus-loop-centric", "loop-surface-hybrid"}:
            score = 0.45 if active_loop_count == 1 else 0.35
        else:
            score = 0.2

    confidence = "high" if source == "prior_document_override" else "low" if score < 0.7 else "medium"
    return {
        "binding_loop_score": round(float(score), 4),
        "binding_loop_label": 1 if score >= 0.5 else 0,
        "label_source": source,
        "evidence_type": evidence,
        "confidence_tier": confidence,
    }


def compact_residue_labels_to_subregions(residue_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    if not residue_rows:
        return []

    rows = sorted(residue_rows, key=lambda x: int(x["resnum"]))
    subregions: list[dict[str, object]] = []
    current = {
        "loop_id": rows[0]["loop_id"],
        "binding_label": rows[0]["binding_label"],
        "subregion_start": rows[0]["resnum"],
        "subregion_end": rows[0]["resnum"],
        "label_source": rows[0].get("label_source"),
        "selected_target_chain_id": rows[0].get("selected_target_chain_id"),
    }

    for row in rows[1:]:
        contiguous = int(row["resnum"]) == int(current["subregion_end"]) + 1
        same_label = row["binding_label"] == current["binding_label"]
        if contiguous and same_label:
            current["subregion_end"] = row["resnum"]
            continue
        subregions.append(current.copy())
        current = {
            "loop_id": row["loop_id"],
            "binding_label": row["binding_label"],
            "subregion_start": row["resnum"],
            "subregion_end": row["resnum"],
            "label_source": row.get("label_source"),
            "selected_target_chain_id": row.get("selected_target_chain_id"),
        }
    subregions.append(current.copy())
    return subregions


def read_csv_fallback(path: Path) -> pd.DataFrame:
    for encoding in (None, "utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding) if encoding else pd.read_csv(path)
        except UnicodeDecodeError:
            continue
    raise UnicodeDecodeError("csv", b"", 0, 1, f"Unable to decode {path}")


def _normalize_scaffold_name(value: object) -> str:
    return str(value).strip().lower().replace("-", "").replace("_", "")


def assign_analysis_tier(scaffold: object) -> str:
    key = _normalize_scaffold_name(scaffold)
    if key in PRIMARY_LOOPCENTRIC:
        return "primary_loopcentric"
    if key in SECONDARY_HYBRID:
        return "secondary_hybrid"
    return "unsupported_nonloopcentric"


def _prior_lookup(prior_df: pd.DataFrame) -> dict[str, pd.Series]:
    records: dict[str, pd.Series] = {}
    for _, row in prior_df.iterrows():
        records[_normalize_scaffold_name(row.get("scaffold", ""))] = row
    return records


def _build_whole_loop_subregion(loop_row: pd.Series, binding_label: int, label_source: str) -> dict[str, object]:
    return {
        "loop_id": loop_row["loop_id"],
        "binding_label": int(binding_label),
        "subregion_start": int(loop_row["loop_start"]),
        "subregion_end": int(loop_row["loop_end"]),
        "label_source": label_source,
        "selected_target_chain_id": pd.NA,
    }


def _pick_review_value(review_row: pd.Series, *columns: str) -> object:
    for column in columns:
        if column in review_row.index and pd.notna(review_row[column]) and str(review_row[column]).strip():
            return review_row[column]
    return pd.NA


def _decision_to_binding_label(decision: object) -> object:
    mapping = {
        "confirm_binding_window": 1,
        "downgrade_to_nonbinding_window": 0,
        "needs_more_evidence": pd.NA,
    }
    return mapping.get(str(decision).strip(), pd.NA)


def apply_manual_review_overrides(
    loop_labels: pd.DataFrame,
    subregions: pd.DataFrame,
    review_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if review_df.empty:
        final_loop_labels = loop_labels.copy()
        final_subregions = subregions.copy()
        final_loop_labels["final_binding_loop_label"] = pd.array(
            pd.to_numeric(final_loop_labels["binding_loop_label"], errors="coerce"),
            dtype="Int64",
        )
        final_loop_labels["final_review_decision"] = "auto_retained"
        final_loop_labels["final_review_basis"] = final_loop_labels.get("label_source", pd.Series(index=final_loop_labels.index, dtype="object"))
        final_loop_labels["final_review_note"] = pd.NA
        final_loop_labels["final_label_source"] = final_loop_labels.get("label_source", pd.Series(index=final_loop_labels.index, dtype="object"))
        final_loop_labels["manual_review_applied"] = False
        final_subregions["final_binding_label"] = pd.array(
            pd.to_numeric(final_subregions["binding_label"], errors="coerce"),
            dtype="Int64",
        )
        final_subregions["final_label_source"] = final_subregions.get("label_source", pd.Series(index=final_subregions.index, dtype="object"))
        return final_loop_labels, final_subregions, final_loop_labels.iloc[0:0].copy()

    review_unique = review_df.drop_duplicates(subset=["loop_id"]).copy()
    review_cols = ["loop_id", *[c for c in review_unique.columns if c != "loop_id" and c not in loop_labels.columns]]
    final_loop_labels = loop_labels.merge(review_unique[review_cols], on="loop_id", how="left")
    final_subregions = subregions.copy()

    final_loop_labels["final_binding_loop_label"] = pd.array(
        pd.to_numeric(final_loop_labels["binding_loop_label"], errors="coerce"),
        dtype="Int64",
    )
    final_loop_labels["final_review_decision"] = "auto_retained"
    final_loop_labels["final_review_basis"] = final_loop_labels["label_source"]
    final_loop_labels["final_review_note"] = pd.NA
    final_loop_labels["final_label_source"] = final_loop_labels["label_source"]
    final_loop_labels["manual_review_applied"] = False

    final_subregions["final_binding_label"] = pd.array(
        pd.to_numeric(final_subregions["binding_label"], errors="coerce"),
        dtype="Int64",
    )
    final_subregions["final_label_source"] = final_subregions["label_source"]

    for idx, review_row in review_unique.iterrows():
        loop_id = review_row["loop_id"]
        decision = _pick_review_value(review_row, "round2_review_decision", "manual_review_decision")
        basis = _pick_review_value(review_row, "round2_review_basis", "manual_review_basis")
        note = _pick_review_value(review_row, "round2_review_note", "manual_review_note")
        final_label = _decision_to_binding_label(decision)

        row_mask = final_loop_labels["loop_id"] == loop_id
        if not row_mask.any():
            continue

        final_loop_labels.loc[row_mask, "final_review_decision"] = decision
        final_loop_labels.loc[row_mask, "final_review_basis"] = basis
        final_loop_labels.loc[row_mask, "final_review_note"] = note
        final_loop_labels.loc[row_mask, "final_binding_loop_label"] = final_label
        final_loop_labels.loc[row_mask, "final_label_source"] = (
            f"manual_review::{basis}" if pd.notna(basis) and str(basis).strip() else "manual_review"
        )
        final_loop_labels.loc[row_mask, "manual_review_applied"] = True

        sub_mask = final_subregions["loop_id"] == loop_id
        if sub_mask.any():
            final_subregions.loc[sub_mask, "final_binding_label"] = final_label
            final_subregions.loc[sub_mask, "final_label_source"] = (
                f"manual_review::{basis}" if pd.notna(basis) and str(basis).strip() else "manual_review"
            )

    unresolved = final_loop_labels.loc[final_loop_labels["final_review_decision"] == "needs_more_evidence"].copy()
    return final_loop_labels, final_subregions, unresolved


def build_binding_subregion_outputs(
    loop_df: pd.DataFrame,
    candidate_df: pd.DataFrame,
    wide_df: pd.DataFrame,
    prior_df: pd.DataFrame,
    annotation_df: pd.DataFrame | None = None,
    distance_cutoff: float = 5.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prior_map = _prior_lookup(prior_df)
    candidate_map = candidate_df.set_index("loop_id")
    wide_map = wide_df.set_index(["protein_row_id", "structure_unique_sequence_id"])
    annotation_map = annotation_df.set_index("loop_id") if annotation_df is not None and not annotation_df.empty else None

    loop_label_rows: list[dict[str, object]] = []
    subregion_rows: list[dict[str, object]] = []
    manual_rows: list[dict[str, object]] = []
    parsed_cache: dict[str, dict[str, dict[int, list[tuple[float, float, float]]]]] = {}

    active_counts = (
        loop_df.groupby(["protein_row_id", "structure_unique_sequence_id"], dropna=False)["loop_id"]
        .size()
        .to_dict()
    )
    loop_lookup = {row["loop_id"]: row for _, row in loop_df.iterrows()}

    for _, loop_row in loop_df.iterrows():
        loop_id = loop_row["loop_id"]
        candidate_row = candidate_map.loc[loop_id] if loop_id in candidate_map.index else pd.Series(dtype=object)
        annotation_row = (
            annotation_map.loc[loop_id]
            if annotation_map is not None and loop_id in annotation_map.index
            else pd.Series(dtype=object)
        )
        wide_key = (loop_row["protein_row_id"], loop_row["structure_unique_sequence_id"])
        wide_row = wide_map.loc[wide_key] if wide_key in wide_map.index else pd.Series(dtype=object)
        if isinstance(wide_row, pd.DataFrame):
            wide_row = wide_row.iloc[0]

        scaffold_key = _normalize_scaffold_name(wide_row.get("Scaffold_Category", ""))
        prior_row = prior_map.get(scaffold_key, pd.Series(dtype=object))
        route = str(candidate_row.get("binding_candidate_route", "")).strip()
        structure_path = str(candidate_row.get("selected_local_structure_path", "")).strip()
        if (not structure_path or not Path(structure_path).exists()) and pd.notna(loop_row.get("structure_job_id")):
            fallback = GETCONTACT_STRUCTURES / f"{loop_row['structure_job_id']}.pdb"
            if fallback.exists():
                structure_path = str(fallback)

        residue_labels: list[dict[str, object]] = []
        loop_label = 0
        binding_score = 0.0
        label_source = "manual_needed"
        confidence_tier = "low"
        evidence_type = "unresolved"

        target_chain_ids = candidate_row.get("target_chain_candidate")
        if pd.notna(annotation_row.get("target_chain_id")) and str(annotation_row.get("target_chain_id")).strip():
            target_chain_ids = annotation_row.get("target_chain_id")

        if structure_path and pd.notna(annotation_row.get("complex_structure_path")) and str(annotation_row.get("complex_structure_path")).strip():
            ann_path = str(annotation_row.get("complex_structure_path")).strip()
            if Path(ann_path).exists():
                structure_path = ann_path

        if structure_path and (
            route in {"direct_complex_from_local_exact_pdb", "contextual_complex_needs_chain_mapping"}
            or (pd.notna(annotation_row.get("target_chain_id")) and str(annotation_row.get("target_chain_id")).strip())
        ):
            if structure_path not in parsed_cache:
                parsed_cache[structure_path] = parse_pdb_structure(structure_path)
            residue_labels = label_loop_residues_from_structure(
                loop_row=loop_row,
                parsed_structure=parsed_cache[structure_path],
                target_chain_ids=target_chain_ids,
                distance_cutoff=distance_cutoff,
            )

        if residue_labels:
            subregions = compact_residue_labels_to_subregions(residue_labels)
            subregion_rows.extend(subregions)
            loop_label = 1 if any(int(r["binding_label"]) == 1 for r in residue_labels) else 0
            binding_score = round(sum(int(r["binding_label"]) for r in residue_labels) / max(len(residue_labels), 1), 4)
            label_source = (
                "annotation_target_chain"
                if pd.notna(annotation_row.get("target_chain_id")) and str(annotation_row.get("target_chain_id")).strip()
                else "direct_structure_contact"
            )
            confidence_tier = "high" if route == "direct_complex_from_local_exact_pdb" else "medium"
            evidence_type = "structure_contact"
        else:
            prior_result = infer_prior_loop_label(
                loop_row=loop_row,
                wide_row=wide_row,
                prior_row=prior_row,
                active_loop_count=active_counts.get(wide_key, 1),
            )
            loop_label = int(prior_result["binding_loop_label"])
            binding_score = float(prior_result["binding_loop_score"])
            label_source = str(prior_result["label_source"])
            confidence_tier = str(prior_result["confidence_tier"])
            evidence_type = str(prior_result["evidence_type"])
            subregion_rows.append(_build_whole_loop_subregion(loop_row, loop_label, label_source))

            if label_source == "prior_default" and route == "manual_binding_annotation_needed":
                manual_rows.append(
                    {
                        "loop_id": loop_id,
                        "protein_row_id": loop_row["protein_row_id"],
                        "structure_unique_sequence_id": loop_row["structure_unique_sequence_id"],
                        "Scaffold_Category": wide_row.get("Scaffold_Category", pd.NA),
                        "binding_candidate_route": route,
                        "reason": "No direct complex evidence; prior-only labeling should be reviewed if used in high-confidence analyses.",
                    }
                )

        loop_label_rows.append(
            {
                "protein_row_id": loop_row["protein_row_id"],
                "structure_unique_sequence_id": loop_row["structure_unique_sequence_id"],
                "loop_id": loop_id,
                "loop_slot": loop_row.get("loop_slot", pd.NA),
                "Scaffold_Category": wide_row.get("Scaffold_Category", pd.NA),
                "analysis_tier": assign_analysis_tier(wide_row.get("Scaffold_Category", pd.NA)),
                "binding_loop_label": loop_label,
                "binding_loop_score": binding_score,
                "label_source": label_source,
                "confidence_tier": confidence_tier,
                "evidence_type": evidence_type,
                "binding_residue_count": sum(int(r["binding_label"]) for r in residue_labels) if residue_labels else (int(loop_row["loop_length"]) if loop_label else 0),
                "nonbinding_residue_count": (len(residue_labels) - sum(int(r["binding_label"]) for r in residue_labels)) if residue_labels else (0 if loop_label else int(loop_row["loop_length"])),
                "subregion_mode": "residue_compacted" if residue_labels else "whole_loop_prior",
                "selected_target_chain_id": residue_labels[0]["selected_target_chain_id"] if residue_labels else pd.NA,
            }
        )

    loop_labels_df = pd.DataFrame(loop_label_rows)
    subregions_df = pd.DataFrame(subregion_rows)

    high_evidence = loop_labels_df.loc[
        loop_labels_df["label_source"].isin(["annotation_target_chain", "direct_structure_contact"]),
        ["loop_id", "Scaffold_Category", "binding_loop_label"],
    ].merge(
        loop_df.reindex(columns=["loop_id", "loop_slot", "loop_flank_pattern"]),
        on="loop_id",
        how="left",
    )
    empirical_priors = build_empirical_loop_priors(high_evidence)

    for idx, row in loop_labels_df.loc[loop_labels_df["label_source"] == "prior_default"].iterrows():
        loop_id = row["loop_id"]
        loop_row = loop_lookup[loop_id]
        wide_key = (loop_row["protein_row_id"], loop_row["structure_unique_sequence_id"])
        wide_row = wide_map.loc[wide_key] if wide_key in wide_map.index else pd.Series(dtype=object)
        if isinstance(wide_row, pd.DataFrame):
            wide_row = wide_row.iloc[0]
        scaffold_key = _normalize_scaffold_name(wide_row.get("Scaffold_Category", ""))
        prior_row = prior_map.get(scaffold_key, pd.Series(dtype=object))
        recalibrated = infer_prior_loop_label(
            loop_row=loop_row,
            wide_row=wide_row,
            prior_row=prior_row,
            active_loop_count=active_counts.get(wide_key, 1),
            empirical_priors=empirical_priors,
        )
        if recalibrated["label_source"] == "prior_default":
            continue
        loop_labels_df.at[idx, "binding_loop_label"] = recalibrated["binding_loop_label"]
        loop_labels_df.at[idx, "binding_loop_score"] = recalibrated["binding_loop_score"]
        loop_labels_df.at[idx, "label_source"] = recalibrated["label_source"]
        loop_labels_df.at[idx, "confidence_tier"] = recalibrated["confidence_tier"]
        loop_labels_df.at[idx, "evidence_type"] = recalibrated["evidence_type"]
        loop_labels_df.at[idx, "binding_residue_count"] = int(loop_row["loop_length"]) if int(recalibrated["binding_loop_label"]) == 1 else 0
        loop_labels_df.at[idx, "nonbinding_residue_count"] = 0 if int(recalibrated["binding_loop_label"]) == 1 else int(loop_row["loop_length"])

        sub_idx = subregions_df.index[subregions_df["loop_id"] == loop_id].tolist()
        for sidx in sub_idx:
            if subregions_df.at[sidx, "label_source"] == "prior_default":
                subregions_df.at[sidx, "binding_label"] = recalibrated["binding_loop_label"]
                subregions_df.at[sidx, "label_source"] = recalibrated["label_source"]

    manual_rows = []
    for _, row in loop_labels_df.iterrows():
        loop_id = row["loop_id"]
        candidate_row = candidate_map.loc[loop_id] if loop_id in candidate_map.index else pd.Series(dtype=object)
        route = str(candidate_row.get("binding_candidate_route", "")).strip()
        if row["label_source"] == "prior_default" and route == "manual_binding_annotation_needed":
            manual_rows.append(
                {
                    "loop_id": loop_id,
                    "protein_row_id": row["protein_row_id"],
                    "structure_unique_sequence_id": row["structure_unique_sequence_id"],
                    "Scaffold_Category": row["Scaffold_Category"],
                    "binding_candidate_route": route,
                    "reason": "No direct complex evidence; prior-only labeling should be reviewed if used in high-confidence analyses.",
                }
            )

    return (
        loop_labels_df,
        subregions_df,
        pd.DataFrame(manual_rows),
    )


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


def build_primary_main_stats(loop_labels: pd.DataFrame, flex_df: pd.DataFrame, label_column: str = "binding_loop_label") -> pd.DataFrame:
    metrics = [
        "avg_rsa",
        "loop_constraint_score_total_norm",
        "itsflex_score",
        "flexscore_ensemble_diversity",
    ]
    merged = loop_labels.merge(flex_df[["loop_id", *metrics]], on="loop_id", how="left")
    primary = merged.loc[merged["analysis_tier"] == "primary_loopcentric"].copy()
    if primary.empty:
        return pd.DataFrame()
    primary = primary.loc[pd.to_numeric(primary[label_column], errors="coerce").isin([0, 1])].copy()
    if primary.empty:
        return pd.DataFrame()

    rows: list[dict[str, object]] = []
    groups = {"overall_primary": primary}
    for scaffold, group_df in primary.groupby("Scaffold_Category", dropna=False):
        groups[str(scaffold)] = group_df

    for group_name, group_df in groups.items():
        for metric in metrics:
            binding = group_df.loc[pd.to_numeric(group_df[label_column], errors="coerce") == 1, metric]
            nonbinding = group_df.loc[pd.to_numeric(group_df[label_column], errors="coerce") == 0, metric]
            row = {"group": group_name, "metric": metric}
            row.update(_metric_stats(binding, nonbinding))
            rows.append(row)
    return pd.DataFrame(rows)


def build_report_payload(
    loop_labels: pd.DataFrame,
    subregions: pd.DataFrame,
    manual_queue: pd.DataFrame,
    label_column: str = "binding_loop_label",
    label_source_column: str = "label_source",
    confidence_column: str = "confidence_tier",
) -> dict[str, object]:
    numeric_labels = pd.to_numeric(loop_labels.get(label_column), errors="coerce")
    return {
        "n_loops": int(len(loop_labels)),
        "n_binding_loops": int((numeric_labels == 1).sum()) if not loop_labels.empty else 0,
        "n_nonbinding_loops": int((numeric_labels == 0).sum()) if not loop_labels.empty else 0,
        "n_subregions": int(len(subregions)),
        "n_manual_queue": int(len(manual_queue)),
        "label_source_counts": loop_labels[label_source_column].value_counts(dropna=False).to_dict() if not loop_labels.empty else {},
        "confidence_tier_counts": loop_labels[confidence_column].value_counts(dropna=False).to_dict() if not loop_labels.empty else {},
        "scaffold_binding_counts": (
            loop_labels.assign(_binding_value=numeric_labels)
            .groupby(["Scaffold_Category", "_binding_value"], dropna=False)["loop_id"]
            .size()
            .reset_index(name="n")
            .rename(columns={"_binding_value": label_column})
            .to_dict(orient="records")
            if not loop_labels.empty
            else []
        ),
    }


def write_outputs(
    loop_labels: pd.DataFrame,
    subregions: pd.DataFrame,
    manual_queue: pd.DataFrame,
    primary_stats: pd.DataFrame | None = None,
    output_prefix: str = "v1",
    label_column: str = "binding_loop_label",
    label_source_column: str = "label_source",
    confidence_column: str = "confidence_tier",
) -> dict[str, Path]:
    loop_labels_path = ROOT / f"loop_binding_labels_{output_prefix}.csv"
    subregions_path = ROOT / f"loop_binding_subregions_long_{output_prefix}.csv"
    summary_path = ROOT / f"binding_subregion_summary_{output_prefix}.csv"
    manual_queue_path = ROOT / f"binding_subregion_manual_queue_{output_prefix}.csv"
    primary_stats_path = ROOT / f"binding_vs_nonbinding_primary_stats_{output_prefix}.csv"
    primary_report_path = ROOT / f"binding_vs_nonbinding_primary_report_{output_prefix}.md"
    report_path = ROOT / f"binding_subregion_report_{output_prefix}.md"
    json_path = ROOT / f"binding_subregion_summary_{output_prefix}.json"

    loop_labels.to_csv(loop_labels_path, index=False, encoding="utf-8")
    subregions.to_csv(subregions_path, index=False, encoding="utf-8")
    manual_queue.to_csv(manual_queue_path, index=False, encoding="utf-8")
    if primary_stats is not None:
        primary_stats.to_csv(primary_stats_path, index=False, encoding="utf-8")
        primary_lines = [
            "# Binding Vs Non-binding Primary Statistics",
            "",
            "This file summarizes the first-pass comparison between binding-labeled and non-binding-labeled loop windows within `primary_loopcentric` scaffolds.",
            "",
        ]
        overall = primary_stats.loc[primary_stats["group"] == "overall_primary"].copy()
        if not overall.empty:
            primary_lines.append("## Overall Primary")
            for _, row in overall.iterrows():
                primary_lines.append(
                    f"- `{row['metric']}`: binding_mean=`{row['binding_mean']}`, nonbinding_mean=`{row['nonbinding_mean']}`, delta_mean=`{row['delta_mean']}`, p=`{row['mannwhitney_p']}`"
                )
            primary_lines.append("")
        primary_lines.append("## Per Scaffold")
        for scaffold in sorted(x for x in primary_stats["group"].unique() if x != "overall_primary"):
            primary_lines.append(f"### {scaffold}")
            sub = primary_stats.loc[primary_stats["group"] == scaffold]
            for _, row in sub.iterrows():
                primary_lines.append(
                    f"- `{row['metric']}`: binding_mean=`{row['binding_mean']}`, nonbinding_mean=`{row['nonbinding_mean']}`, delta_mean=`{row['delta_mean']}`, p=`{row['mannwhitney_p']}`"
                )
            primary_lines.append("")
        primary_report_path.write_text("\n".join(primary_lines), encoding="utf-8")

    summary = (
        loop_labels.groupby([label_source_column, confidence_column, label_column], dropna=False)["loop_id"]
        .size()
        .reset_index(name="n_loops")
        .sort_values("n_loops", ascending=False, kind="stable")
    )
    summary.to_csv(summary_path, index=False, encoding="utf-8")

    payload = build_report_payload(
        loop_labels,
        subregions,
        manual_queue,
        label_column=label_column,
        label_source_column=label_source_column,
        confidence_column=confidence_column,
    )
    pd.Series(payload).to_json(json_path, force_ascii=False, indent=2)

    report_lines = [
        "# Binding Subregion Analysis Report",
        "",
        f"- Total loops: `{payload['n_loops']}`",
        f"- Binding loops: `{payload['n_binding_loops']}`",
        f"- Non-binding loops: `{payload['n_nonbinding_loops']}`",
        f"- Total subregions: `{payload['n_subregions']}`",
        f"- Manual queue rows: `{payload['n_manual_queue']}`",
        "",
        "## Label Sources",
    ]
    for key, value in payload["label_source_counts"].items():
        report_lines.append(f"- `{key}`: `{value}`")
    report_lines.extend(["", "## Confidence Tiers"])
    for key, value in payload["confidence_tier_counts"].items():
        report_lines.append(f"- `{key}`: `{value}`")
    report_lines.extend(
        [
            "",
            "## Output Files",
            f"- `{loop_labels_path}`",
            f"- `{subregions_path}`",
            f"- `{summary_path}`",
            f"- `{manual_queue_path}`",
            f"- `{primary_stats_path}`",
            f"- `{primary_report_path}`",
            f"- `{json_path}`",
        ]
    )
    report_path.write_text("\n".join(report_lines), encoding="utf-8")

    return {
        "loop_labels": loop_labels_path,
        "subregions": subregions_path,
        "summary": summary_path,
        "manual_queue": manual_queue_path,
        "primary_stats": primary_stats_path,
        "primary_report": primary_report_path,
        "report": report_path,
        "json": json_path,
    }


def main() -> None:
    loop_df = read_csv_fallback(LOOP_RESULTS)
    candidate_df = read_csv_fallback(CANDIDATE_TABLE)
    annotation_df = read_csv_fallback(ANNOTATION_INPUT)
    wide_df = read_csv_fallback(WIDE_TABLE)
    prior_df = read_csv_fallback(PRIOR_TABLE)

    loop_labels, subregions, manual_queue = build_binding_subregion_outputs(
        loop_df=loop_df,
        candidate_df=candidate_df,
        wide_df=wide_df,
        prior_df=prior_df,
        annotation_df=annotation_df,
    )
    primary_stats = build_primary_main_stats(loop_labels, loop_df)
    write_outputs(loop_labels, subregions, manual_queue, primary_stats=primary_stats)

    if MANUAL_FINAL_REVIEW.exists():
        review_df = read_csv_fallback(MANUAL_FINAL_REVIEW)
        final_loop_labels, final_subregions, unresolved = apply_manual_review_overrides(
            loop_labels=loop_labels,
            subregions=subregions,
            review_df=review_df,
        )
        final_primary_stats = build_primary_main_stats(
            final_loop_labels,
            loop_df,
            label_column="final_binding_loop_label",
        )
        write_outputs(
            final_loop_labels,
            final_subregions,
            unresolved,
            primary_stats=final_primary_stats,
            output_prefix="manual_final_v1",
            label_column="final_binding_loop_label",
            label_source_column="final_label_source",
            confidence_column="confidence_tier",
        )


if __name__ == "__main__":
    main()
