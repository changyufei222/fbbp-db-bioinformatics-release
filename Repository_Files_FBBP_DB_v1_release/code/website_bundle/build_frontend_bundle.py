import csv
import json
import os
import re
from collections import Counter, defaultdict


DB = r"<local_path_removed>"
OUT = os.path.join(os.path.dirname(__file__), "data")


def read_csv(name):
    path = os.path.join(DB, name)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def num(v):
    s = str(v or "").strip()
    if not s:
        return None
    try:
        if re.fullmatch(r"-?\d+", s):
            return int(s)
        return float(s)
    except Exception:
        return None


def truthy(v):
    s = str(v or "").strip().lower()
    return s in {"true", "1", "yes", "y"}


def norm(v):
    return str(v or "").strip()


def key_of(row):
    return (norm(row.get("protein_row_id") or row.get("protein_id")), norm(row.get("structure_unique_sequence_id")))


def make_link(kind, ident):
    ident = norm(ident)
    if not ident:
        return ""
    if kind == "uniprot":
        return f"https://www.uniprot.org/uniprotkb/{ident}/entry"
    if kind == "pdb":
        return f"https://www.rcsb.org/structure/{ident}"
    if kind == "alphafold":
        return ident if ident.startswith("http") else f"https://alphafold.ebi.ac.uk/entry/{ident}"
    if kind == "chembl":
        return ident if ident.startswith("http") else f"https://www.ebi.ac.uk/chembl/target_report_card/{ident}/"
    if kind == "drugbank":
        return ident if ident.startswith("http") else f"https://go.drugbank.com/unearth/q?searcher=targets&query={ident}"
    if kind == "kegg":
        return ident if ident.startswith("http") else f"https://www.kegg.jp/dbget-bin/www_bget?{ident}"
    return ident


def top_counts(rows, field, limit=10):
    counter = Counter()
    for row in rows:
        value = norm(row.get(field))
        if value:
            counter[value] += 1
    return [{"label": k, "count": v} for k, v in counter.most_common(limit)]


proteins = read_csv("proteins.csv")
sources = read_csv("sources.csv")
protein_identifiers = read_csv("protein_identifiers.csv")
domains = read_csv("domains.csv")
targets_conceptual = read_csv("targets_conceptual.csv")
target_species_variants = read_csv("target_species_variants.csv")
interactions = read_csv("interactions.csv")
affinity_data = read_csv("affinity_data.csv")
digestive_assays = read_csv("digestive_assays.csv")
functional_annotations = read_csv("functional_annotations.csv")
structural_source_info = read_csv("structural_source_info.csv")
active_structure_context = read_csv("active_structure_context.csv")

bcell = read_csv("bcell_epitope_results.csv")
mhci = read_csv("mhci_results.csv")
mhcii = read_csv("mhcii_results.csv")
immun = read_csv("immunogenicity_summary.csv")
ecoli = read_csv("ecoli_expression_results.csv")
oral = read_csv("oral_results.csv")
sol = read_csv("solubility_results.csv")
protparam = read_csv("protparam_results.csv")
ted = read_csv("ted_results.csv")
protrek = read_csv("protrek_results.csv")
foldseek = read_csv("foldseek_results.csv")
plmsearch = read_csv("plmsearch_results.csv")
protein_flex = read_csv("protein_flexibility_summary.csv")
loop_ann = read_csv("loop_annotations.csv")
loop_flex = read_csv("loop_flexibility_public_summary.csv")


table_indexes = {}
for rows in [bcell, mhci, mhcii, immun, ecoli, oral, sol, protparam, ted, protrek, foldseek, plmsearch, protein_flex]:
    table_indexes[id(rows)] = {key_of(r): r for r in rows}


sources_by_id = {norm(r["source_id"]): r for r in sources}
proteins_by_id = {norm(r["protein_id"]): r for r in proteins}
identifiers_by_protein = defaultdict(list)
for r in protein_identifiers:
    identifiers_by_protein[norm(r["protein_id"])].append(r)

func_by_protein = defaultdict(list)
for r in functional_annotations:
    func_by_protein[norm(r["protein_id"])].append(r)

domains_by_protein = defaultdict(list)
for r in domains:
    domains_by_protein[norm(r["protein_id"])].append(r)

struct_by_domain = defaultdict(list)
for r in structural_source_info:
    struct_by_domain[norm(r["domain_id"])].append(r)

inter_by_key = defaultdict(list)
for r in interactions:
    inter_by_key[(norm(r["protein_row_id"]), norm(r["structure_unique_sequence_id"]))].append(r)

aff_by_inter = defaultdict(list)
for r in affinity_data:
    aff_by_inter[norm(r["interaction_id"])].append(r)

digest_by_domain = defaultdict(list)
for r in digestive_assays:
    digest_by_domain[norm(r["domain_id"])].append(r)

variant_by_id = {norm(r["target_variant_id"]): r for r in target_species_variants}
concept_by_id = {norm(r["target_concept_id"]): r for r in targets_conceptual}


rows = []
for active in active_structure_context:
    key = key_of(active)
    pid = norm(active["protein_row_id"])
    sid = norm(active["structure_unique_sequence_id"])
    protein = proteins_by_id.get(pid, {})

    rec = {}
    rec["protein_row_id"] = pid
    rec["structure_unique_sequence_id"] = sid
    rec["_id"] = sid
    rec["File_Name"] = norm(protein.get("canonical_name")) or sid
    rec["Scaffold_Category"] = norm(protein.get("scaffold_category"))
    rec["Final_Tested_Sequence"] = norm(protein.get("final_tested_sequence"))
    rec["Original_Sequence"] = norm(protein.get("source_sequence"))
    rec["Sequence_Length"] = norm(protein.get("sequence_length"))
    rec["Proteins_canonical_name"] = norm(protein.get("canonical_name"))
    rec["Proteins_description"] = norm(protein.get("description"))
    rec["Proteins_organism"] = norm(protein.get("organism"))
    rec["Proteins_sequence"] = norm(protein.get("source_sequence"))
    rec["Has_Tag"] = "False"
    rec["Is_Fusion"] = "False"
    rec["Fusion_Type"] = ""

    protein_domains = domains_by_protein.get(pid, [])
    domain0 = protein_domains[0] if protein_domains else {}
    source0 = sources_by_id.get(norm(domain0.get("source_id")), {}) if domain0 else {}

    rec["Domains_domain_name"] = norm(domain0.get("domain_name"))
    rec["Domains_is_engineered"] = norm(domain0.get("is_engineered"))
    rec["Domains_scaffold_type"] = norm(domain0.get("scaffold_type"))
    rec["Domains_sequence"] = norm(domain0.get("sequence"))

    rec["Sources_title"] = norm(source0.get("title"))
    rec["Sources_identifier"] = norm(source0.get("identifier"))
    rec["Sources_publication_date"] = norm(source0.get("publication_date"))
    rec["Sources_authors"] = norm(source0.get("authors"))
    rec["Sources_journal_or_office"] = norm(source0.get("journal_or_office"))
    rec["Sources_source_type"] = norm(source0.get("source_type"))

    ints = inter_by_key.get(key, [])
    inter0 = ints[0] if ints else {}
    variant0 = variant_by_id.get(norm(inter0.get("target_variant_id")), {}) if inter0 else {}
    concept0 = concept_by_id.get(norm(variant0.get("target_concept_id")), {}) if variant0 else {}

    target_gene = norm(concept0.get("gene_name")) or norm(variant0.get("gene_name_species"))
    rec["Targets_gene_name"] = target_gene
    rec["Targets_species_name"] = norm(variant0.get("species_name"))
    rec["Targets_uniprot_id"] = norm(variant0.get("uniprot_id"))
    rec["Interactions_interaction_class"] = norm(inter0.get("interaction_class"))
    rec["Interactions_is_inhibitory"] = norm(inter0.get("is_inhibitory"))
    rec["Derived_Target_Label"] = target_gene or "未标注"

    affinities = []
    for interaction in ints:
        affinities.extend(aff_by_inter.get(norm(interaction.get("interaction_id")), []))
    affinity_values = [num(x.get("value")) for x in affinities if num(x.get("value")) is not None]
    min_kd = min(affinity_values) if affinity_values else None
    rec["Affinity_Data_value"] = "" if min_kd is None else str(min_kd)
    rec["Affinity_Data_value_type"] = norm(affinities[0].get("determination_method")) if affinities else ""
    rec["Derived_Kd_nM"] = "" if min_kd is None else min_kd
    if min_kd is None:
        kd_range = "未提供"
    elif min_kd < 1:
        kd_range = "<1 nM"
    elif min_kd < 10:
        kd_range = "1-10 nM"
    elif min_kd < 100:
        kd_range = "10-100 nM"
    elif min_kd < 1000:
        kd_range = "100-1000 nM"
    else:
        kd_range = ">1000 nM"
    rec["Derived_Kd_Range"] = kd_range

    digestive = digest_by_domain.get(norm(domain0.get("domain_id")), []) if domain0 else []
    digestive0 = digestive[0] if digestive else {}
    rec["Digestive_Assays_enzyme_name"] = norm(digestive0.get("enzyme_name"))
    rec["Digestive_Assays_result_value"] = norm(digestive0.get("result_value"))
    rec["Digestive_Assays_data_type"] = "Experimental" if digestive0 else ""

    identifiers = identifiers_by_protein.get(pid, [])
    seq_uniprot = ""
    seq_pdb = ""
    for item in identifiers:
        id_type = norm(item.get("id_type")).lower()
        id_value = norm(item.get("id_value"))
        if not seq_uniprot and "uniprot" in id_type:
            seq_uniprot = id_value
        if not seq_pdb and "pdb" in id_type:
            seq_pdb = id_value
    rec["Sequence_UniProt_ID"] = seq_uniprot
    rec["Sequence_PDB_ID"] = seq_pdb
    rec["CrossLink_Sequence_UniProt_URLs"] = make_link("uniprot", seq_uniprot)
    rec["CrossLink_Sequence_PDB_URLs"] = make_link("pdb", seq_pdb)

    ann_map = defaultdict(list)
    for annotation in func_by_protein.get(pid, []):
        ann_map[norm(annotation.get("annotation_type"))].append(norm(annotation.get("annotation_value")))
    rec["UniProt_GO_Terms"] = "; ".join(ann_map.get("GO_Terms", []))
    rec["UniProt_Keywords"] = "; ".join(ann_map.get("Keywords", []))
    rec["UniProt_Function"] = "; ".join(ann_map.get("Function", []))
    rec["UniProt_Pathway"] = "; ".join(ann_map.get("Pathway", []))
    rec["UniProt_Protein_Families"] = "; ".join(ann_map.get("Protein_Families", []))
    rec["UniProt_Similarity"] = "; ".join(ann_map.get("Similarity", []))
    rec["UniProt_Pfam"] = "; ".join(ann_map.get("Pfam", []))
    rec["UniProt_SMART"] = "; ".join(ann_map.get("SMART", []))
    rec["UniProt_Gene3D"] = "; ".join(ann_map.get("Gene3D", []))
    rec["UniProt_Protein_Name"] = norm(protein.get("canonical_name"))
    rec["UniProt_Gene_Names"] = target_gene
    rec["UniProt_Organism"] = norm(protein.get("organism"))
    rec["UniProt_Organism_ID"] = ""
    rec["UniProt_Annotation_Score"] = "reviewed" if seq_uniprot else ""
    rec["UniProt_Reviewed"] = "reviewed" if seq_uniprot else "unreviewed"
    rec["UniProt_Entry_Name"] = seq_uniprot
    rec["UniProt_AlphaFoldDB"] = ""
    rec["CrossLink_Sequence_AlphaFoldDB_URLs"] = ""
    rec["CrossLink_Target_AlphaFoldDB_URLs"] = ""
    rec["UniProt_ChEMBL"] = ""
    rec["CrossLink_ChEMBL_URLs"] = ""
    rec["UniProt_DrugBank"] = ""
    rec["CrossLink_DrugBank_URLs"] = ""
    rec["UniProt_KEGG"] = ""
    rec["CrossLink_KEGG_URLs"] = ""
    rec["CrossLink_Target_UniProt_URLs"] = make_link("uniprot", rec["Targets_uniprot_id"])
    rec["Target_PDB_ID"] = ""
    rec["CrossLink_Target_PDB_URLs"] = ""

    source_structs = struct_by_domain.get(norm(domain0.get("domain_id")), []) if domain0 else []
    source_struct0 = source_structs[0] if source_structs else {}
    rec["Structural_Info_method"] = norm(source_struct0.get("method"))
    rec["Structural_Info_resolution"] = norm(source_struct0.get("resolution"))
    rec["PDB_Experimental_Method"] = norm(source_struct0.get("method")) or norm(active.get("structure_source"))
    rec["PDB_Resolution_A"] = norm(source_struct0.get("resolution"))
    rec["PDB_Technique"] = norm(source_struct0.get("source_structure_type"))
    rec["PDB_Binding_Site_Details"] = norm(source_struct0.get("source_structure_note"))

    rec["Structure_Best_Available_Source"] = norm(active.get("structure_best_available_source"))
    rec["Structure_Source"] = norm(active.get("structure_source"))
    rec["Structure_Type"] = norm(active.get("structure_type"))
    rec["Structure_Locator"] = norm(active.get("structure_locator"))
    rec["Structure_Chain_ID"] = norm(active.get("chain_id"))
    rec["Structure_Resnum_Start"] = norm(active.get("resnum_start"))
    rec["Structure_Resnum_End"] = norm(active.get("resnum_end"))
    rec["Structure_Chain_Length"] = norm(active.get("chain_length"))
    rec["Structure_Final_Chain_Assignment"] = norm(active.get("final_chain_assignment"))
    rec["Structure_Complex_Chain_Annotation_Completed"] = norm(active.get("complex_chain_annotation_completed"))
    rec["Structure_PDB_Path"] = norm(active.get("pdb_path"))
    rec["Structure_PDB_Code"] = norm(active.get("pdb_code"))
    rec["Structure_Match_Type"] = norm(active.get("structure_match_type"))
    rec["Structure_Match_Detail"] = norm(active.get("structure_match_detail"))
    rec["Structure_Match_Identity"] = norm(active.get("structure_match_identity"))
    rec["Structure_Match_Coverage"] = norm(active.get("structure_match_coverage"))
    rec["Structure_Match_Chain_Coverage"] = norm(active.get("structure_match_chain_coverage"))

    rec["Structural_Info_backbone_desc"] = norm(active.get("structure_type")) or "not used in final database"
    rec["Structural_Info_flexibility_desc"] = "moved to loop/protein flexibility tables"
    rec["Structural_Info_confidence_score"] = norm(active.get("structure_match_identity"))

    for table in [bcell, mhci, mhcii, immun, ecoli, oral, sol, protparam, ted, protrek, foldseek, plmsearch, protein_flex]:
        hit = table_indexes[id(table)].get(key, {})
        rec.update(hit)

    if truthy(norm(domain0.get("is_engineered"))):
        rec_type = "工程化靶向蛋白"
    elif "fusion" in norm(protein.get("description")).lower():
        rec_type = "融合构建体"
    else:
        rec_type = "天然肽样骨架"
    rec["Derived_Record_Type"] = rec_type

    if norm(rec.get("Targets_gene_name")):
        if norm(rec.get("Affinity_Data_value")):
            scenario = "治疗开发"
        else:
            scenario = "靶向发现"
    elif norm(rec.get("Digestive_Assays_result_value")):
        scenario = "口服递送"
    else:
        scenario = "基础研究"
    rec["Derived_Application_Scenario"] = scenario

    oral_rating = norm(rec.get("Oral_Rating"))
    if oral_rating == "High":
        oral_potential = "已有稳定性数据"
    elif oral_rating == "Medium":
        oral_potential = "潜在口服"
    elif oral_rating == "Low":
        oral_potential = "口服受限"
    else:
        oral_potential = "待评估"
    rec["Derived_Oral_Potential"] = oral_potential

    if norm(rec.get("flexibility_summary_version")):
        workflow = "loop与柔性分析完成"
    elif norm(rec.get("overall_final_judgement")):
        workflow = "免疫与柔性扩展"
    elif norm(rec.get("Structure_Best_Available_Source")):
        workflow = "结构注释"
    elif norm(rec.get("Targets_gene_name")) or norm(rec.get("Affinity_Data_value")):
        workflow = "靶点与亲和力整理"
    else:
        workflow = "数据采集"
    rec["Derived_Workflow_Stage"] = workflow

    overall_conf = num(rec.get("overall_confidence_score_0_100"))
    star = "1"
    if overall_conf is not None:
        if overall_conf >= 90:
            star = "5"
        elif overall_conf >= 80:
            star = "4"
        elif overall_conf >= 65:
            star = "3"
        elif overall_conf >= 50:
            star = "2"
    rec["Derived_Star_Rating"] = star

    seq_len = num(rec.get("Sequence_Length"))
    rec["Length_Flag"] = "short" if seq_len and seq_len < 30 else "normal"

    display_name = norm(rec.get("Proteins_canonical_name")) or sid
    rec["_display_name"] = display_name
    search_parts = [
        pid,
        sid,
        display_name,
        rec.get("Scaffold_Category"),
        rec.get("Targets_gene_name"),
        rec.get("Targets_species_name"),
        rec.get("Sequence_UniProt_ID"),
        rec.get("Sequence_PDB_ID"),
        rec.get("Target_PDB_ID"),
        rec.get("Final_Tested_Sequence"),
        rec.get("Sources_identifier"),
        rec.get("Derived_Application_Scenario"),
        rec.get("Derived_Oral_Potential"),
        rec.get("overall_final_judgement"),
    ]
    rec["_search"] = " ".join(norm(x).lower() for x in search_parts if norm(x))
    rec["_metrics"] = {
        "sequenceLength": num(rec.get("Sequence_Length")),
        "kdNm": min_kd,
        "overallConfidence": num(rec.get("overall_confidence_score_0_100")),
        "bp3Confidence": num(rec.get("BP3_Confidence_Score_0_100")),
        "mhciConfidence": num(rec.get("mhci_confidence_score_0_100")),
        "mhciiConfidence": num(rec.get("mhcii_confidence_score_0_100")),
        "starRating": num(star),
    }
    rec["_flags"] = {
        "hasStructure": bool(norm(rec.get("Sequence_PDB_ID")) or norm(rec.get("PDB_Experimental_Method")) or norm(active.get("structure_best_available_source"))),
        "hasSequenceUniProt": bool(seq_uniprot),
        "hasAlphaFold": bool(norm(rec.get("UniProt_AlphaFoldDB")) or "prediction" in norm(active.get("structure_best_available_source")).lower()),
        "hasImmunogenicity": bool(norm(rec.get("overall_final_judgement"))),
        "reviewedState": norm(rec.get("UniProt_Reviewed")),
    }

    rows.append(rec)


rows.sort(key=lambda r: int(r["protein_row_id"]) if str(r["protein_row_id"]).isdigit() else 10**9)

loop_rows = []
for row in loop_ann:
    item = dict(row)
    item["_id"] = f"{norm(row['structure_unique_sequence_id'])}::{norm(row['loop_id'])}"
    loop_rows.append(item)

flex_rows = []
for row in loop_flex:
    item = dict(row)
    item["_id"] = f"{norm(row['structure_unique_sequence_id'])}::{norm(row['loop_id'])}"
    flex_rows.append(item)


schema_groups = [
    {
        "id": "overview",
        "title": "核心条目",
        "description": "蛋白主记录、骨架、序列与主状态字段。",
        "fields": [
            {"name": "protein_row_id", "label": "Protein Row ID", "kind": "number", "present": True, "planned": False, "description": "蛋白级主键之一"},
            {"name": "structure_unique_sequence_id", "label": "Structure Unique Sequence ID", "kind": "text", "present": True, "planned": False, "description": "蛋白级稳定业务键"},
            {"name": "_display_name", "label": "Entry", "kind": "text", "present": True, "planned": False, "description": "前端显示名"},
            {"name": "Scaffold_Category", "label": "Scaffold", "kind": "text", "present": True, "planned": False, "description": "骨架类别"},
            {"name": "Final_Tested_Sequence", "label": "Final Tested Sequence", "kind": "sequence", "present": True, "planned": False, "description": "最终测试序列"},
            {"name": "Sequence_Length", "label": "Sequence Length", "kind": "number", "present": True, "planned": False, "description": "序列长度"},
            {"name": "Derived_Record_Type", "label": "Record Type", "kind": "text", "present": True, "planned": False, "description": "前端派生记录类型"},
            {"name": "Derived_Star_Rating", "label": "Star Rating", "kind": "text", "present": True, "planned": False, "description": "前端派生星级"},
            {"name": "overall_final_judgement", "label": "Overall Final Judgement", "kind": "text", "present": True, "planned": False, "description": "综合免疫原性判断"},
            {"name": "overall_confidence_score_0_100", "label": "Overall Confidence", "kind": "number", "present": True, "planned": False, "description": "综合置信度得分"},
        ],
    },
    {
        "id": "source",
        "title": "来源与文献",
        "description": "文献、专利、来源机构与来源标识。",
        "fields": [
            {"name": "Sources_title", "label": "Source Title", "kind": "text", "present": True, "planned": False, "description": "来源标题"},
            {"name": "Sources_identifier", "label": "Source Identifier", "kind": "text", "present": True, "planned": False, "description": "来源标识符"},
            {"name": "Sources_publication_date", "label": "Publication Date", "kind": "text", "present": True, "planned": False, "description": "来源日期"},
            {"name": "Sources_authors", "label": "Authors", "kind": "text", "present": True, "planned": False, "description": "作者"},
            {"name": "Sources_journal_or_office", "label": "Journal / Office", "kind": "text", "present": True, "planned": False, "description": "期刊或专利机构"},
            {"name": "Sources_source_type", "label": "Source Type", "kind": "text", "present": True, "planned": False, "description": "文献 / 专利 / UniProt / PDB"},
        ],
    },
    {
        "id": "target_affinity",
        "title": "靶点与亲和力",
        "description": "靶点、物种、结合关系与 Kd 摘要。",
        "fields": [
            {"name": "Targets_gene_name", "label": "Target Gene", "kind": "text", "present": True, "planned": False, "description": "标准化靶点基因"},
            {"name": "Targets_species_name", "label": "Target Species", "kind": "text", "present": True, "planned": False, "description": "靶点物种"},
            {"name": "Targets_uniprot_id", "label": "Target UniProt", "kind": "text", "present": True, "planned": False, "description": "靶点 UniProt"},
            {"name": "Interactions_interaction_class", "label": "Interaction Class", "kind": "text", "present": True, "planned": False, "description": "相互作用类别"},
            {"name": "Interactions_is_inhibitory", "label": "Is Inhibitory", "kind": "text", "present": True, "planned": False, "description": "是否抑制型"},
            {"name": "Affinity_Data_value", "label": "Affinity Value", "kind": "number", "present": True, "planned": False, "description": "亲和力数值摘要"},
            {"name": "Affinity_Data_value_type", "label": "Affinity Type", "kind": "text", "present": True, "planned": False, "description": "亲和力方法 / 单位"},
            {"name": "Derived_Kd_nM", "label": "Derived Kd nM", "kind": "number", "present": True, "planned": False, "description": "前端派生 Kd 数值"},
            {"name": "Derived_Kd_Range", "label": "Derived Kd Range", "kind": "text", "present": True, "planned": False, "description": "前端派生 Kd 分层"},
        ],
    },
    {
        "id": "structure",
        "title": "结构层",
        "description": "区分来源结构事实与当前活跃结构上下文。",
        "fields": [
            {"name": "Structural_Info_method", "label": "Source Structure Method", "kind": "text", "present": True, "planned": False, "description": "来源结构实验方法"},
            {"name": "Structural_Info_resolution", "label": "Source Structure Resolution", "kind": "text", "present": True, "planned": False, "description": "来源结构分辨率"},
            {"name": "PDB_Experimental_Method", "label": "PDB Experimental Method", "kind": "text", "present": True, "planned": False, "description": "结构来源或采用方式"},
            {"name": "PDB_Resolution_A", "label": "PDB Resolution A", "kind": "text", "present": True, "planned": False, "description": "来源分辨率"},
            {"name": "PDB_Binding_Site_Details", "label": "Source Structure Note", "kind": "text", "present": True, "planned": False, "description": "来源结构说明"},
            {"name": "Structure_Best_Available_Source", "label": "Best Available Source", "kind": "text", "present": True, "planned": False, "description": "当前数据库实际采用的结构来源"},
            {"name": "Structure_Source", "label": "Active Structure Source", "kind": "text", "present": True, "planned": False, "description": "活跃结构来源"},
            {"name": "Structure_Type", "label": "Active Structure Type", "kind": "text", "present": True, "planned": False, "description": "活跃结构类型"},
            {"name": "Structure_Locator", "label": "Structure Locator", "kind": "text", "present": True, "planned": False, "description": "结构定位信息"},
            {"name": "Structure_Chain_ID", "label": "Chain ID", "kind": "text", "present": True, "planned": False, "description": "当前采用链 ID"},
            {"name": "Structure_Resnum_Start", "label": "Resnum Start", "kind": "number", "present": True, "planned": False, "description": "起始残基号"},
            {"name": "Structure_Resnum_End", "label": "Resnum End", "kind": "number", "present": True, "planned": False, "description": "终止残基号"},
            {"name": "Structure_Final_Chain_Assignment", "label": "Final Chain Assignment", "kind": "text", "present": True, "planned": False, "description": "最终链归属"},
            {"name": "Structure_Complex_Chain_Annotation_Completed", "label": "Complex Chain Annotation Completed", "kind": "text", "present": True, "planned": False, "description": "复合物链标注是否完成"},
        ],
    },
    {
        "id": "immunogenicity",
        "title": "免疫原性与表达/口服",
        "description": "BP3、MHC、E. coli、口服与溶解性结果。",
        "fields": [
            {"name": "BP3_BCell_Level", "label": "B-cell Level", "kind": "text", "present": True, "planned": False, "description": "B 细胞表位分层"},
            {"name": "BP3_Confidence_Score_0_100", "label": "B-cell Confidence", "kind": "number", "present": True, "planned": False, "description": "B 细胞置信度"},
            {"name": "mhci_binding_level", "label": "MHC-I Binding Level", "kind": "text", "present": True, "planned": False, "description": "MHC-I 分层"},
            {"name": "mhci_confidence_score_0_100", "label": "MHC-I Confidence", "kind": "number", "present": True, "planned": False, "description": "MHC-I 置信度"},
            {"name": "mhcii_binding_level", "label": "MHC-II Binding Level", "kind": "text", "present": True, "planned": False, "description": "MHC-II 分层"},
            {"name": "mhcii_best_combined_score", "label": "MHC-II Combined Score", "kind": "number", "present": True, "planned": False, "description": "MHC-II 组合分数"},
            {"name": "overall_final_judgement", "label": "Overall Final Judgement", "kind": "text", "present": True, "planned": False, "description": "综合免疫原性结论"},
            {"name": "EcoliExpr_screen_track", "label": "E. coli Track", "kind": "text", "present": True, "planned": False, "description": "E. coli 表达轨道"},
            {"name": "EcoliExpr_recommended_ecoli_route", "label": "Recommended Route", "kind": "text", "present": True, "planned": False, "description": "推荐表达路线"},
            {"name": "Oral_Rating", "label": "Oral Rating", "kind": "text", "present": True, "planned": False, "description": "口服潜力评分"},
            {"name": "QCBundleServerOfficial_final_assessment", "label": "Solubility Final Assessment", "kind": "text", "present": True, "planned": False, "description": "溶解性/开发性总结"},
        ],
    },
    {
        "id": "annotations",
        "title": "注释与理化",
        "description": "UniProt 注释、ProtParam、蛋白宇宙结果。",
        "fields": [
            {"name": "Sequence_UniProt_ID", "label": "Sequence UniProt ID", "kind": "text", "present": True, "planned": False, "description": "序列 UniProt"},
            {"name": "UniProt_Protein_Name", "label": "UniProt Protein Name", "kind": "text", "present": True, "planned": False, "description": "蛋白名"},
            {"name": "UniProt_Gene_Names", "label": "UniProt Gene Names", "kind": "text", "present": True, "planned": False, "description": "基因名"},
            {"name": "UniProt_Reviewed", "label": "UniProt Reviewed", "kind": "text", "present": True, "planned": False, "description": "UniProt 审核状态"},
            {"name": "UniProt_GO_Terms", "label": "GO Terms", "kind": "text", "present": True, "planned": False, "description": "GO 注释"},
            {"name": "UniProt_Keywords", "label": "Keywords", "kind": "text", "present": True, "planned": False, "description": "关键词"},
            {"name": "UniProt_Function", "label": "Function", "kind": "text", "present": True, "planned": False, "description": "功能说明"},
            {"name": "ProtParam_MW_Da", "label": "ProtParam MW", "kind": "number", "present": True, "planned": False, "description": "分子量"},
            {"name": "ProtParam_pI", "label": "ProtParam pI", "kind": "number", "present": True, "planned": False, "description": "等电点"},
            {"name": "TED_best_cath_superfamily_name", "label": "TED Superfamily", "kind": "text", "present": True, "planned": False, "description": "TED 最佳超家族"},
            {"name": "protrek_confidence", "label": "ProTrek Confidence", "kind": "text", "present": True, "planned": False, "description": "ProTrek 置信度"},
            {"name": "foldseek_confidence", "label": "Foldseek Confidence", "kind": "text", "present": True, "planned": False, "description": "Foldseek 置信度"},
            {"name": "PLMSearch_top1_target", "label": "PLMSearch Top1", "kind": "text", "present": True, "planned": False, "description": "PLMSearch top1 命中"},
        ],
    },
    {
        "id": "crosslinks",
        "title": "外部资源链接",
        "description": "外部数据库链接与跳转入口。",
        "fields": [
            {"name": "CrossLink_Sequence_UniProt_URLs", "label": "Sequence UniProt URL", "kind": "link", "present": True, "planned": False, "description": "序列 UniProt 外链"},
            {"name": "CrossLink_Target_UniProt_URLs", "label": "Target UniProt URL", "kind": "link", "present": True, "planned": False, "description": "靶点 UniProt 外链"},
            {"name": "CrossLink_Sequence_PDB_URLs", "label": "Sequence PDB URL", "kind": "link", "present": True, "planned": False, "description": "序列 PDB 外链"},
            {"name": "CrossLink_Target_PDB_URLs", "label": "Target PDB URL", "kind": "link", "present": True, "planned": False, "description": "靶点 PDB 外链"},
            {"name": "CrossLink_Sequence_AlphaFoldDB_URLs", "label": "Sequence AlphaFoldDB URL", "kind": "link", "present": True, "planned": False, "description": "序列 AlphaFoldDB 外链"},
            {"name": "CrossLink_Target_AlphaFoldDB_URLs", "label": "Target AlphaFoldDB URL", "kind": "link", "present": True, "planned": False, "description": "靶点 AlphaFoldDB 外链"},
            {"name": "CrossLink_ChEMBL_URLs", "label": "ChEMBL URL", "kind": "link", "present": True, "planned": False, "description": "ChEMBL 外链"},
            {"name": "CrossLink_DrugBank_URLs", "label": "DrugBank URL", "kind": "link", "present": True, "planned": False, "description": "DrugBank 外链"},
            {"name": "CrossLink_KEGG_URLs", "label": "KEGG URL", "kind": "link", "present": True, "planned": False, "description": "KEGG 外链"},
        ],
    },
    {
        "id": "frontend_derivatives",
        "title": "前端派生字段",
        "description": "仅服务前端筛选、搜索、排序和展示的派生字段。",
        "fields": [
            {"name": "Derived_Target_Label", "label": "Derived Target Label", "kind": "text", "present": True, "planned": False, "description": "前端派生靶点标签"},
            {"name": "Derived_Kd_nM", "label": "Derived Kd nM", "kind": "number", "present": True, "planned": False, "description": "前端派生 Kd 数值"},
            {"name": "Derived_Kd_Range", "label": "Derived Kd Range", "kind": "text", "present": True, "planned": False, "description": "前端派生 Kd 区间"},
            {"name": "Derived_Application_Scenario", "label": "Derived Application Scenario", "kind": "text", "present": True, "planned": False, "description": "前端派生应用场景"},
            {"name": "Derived_Oral_Potential", "label": "Derived Oral Potential", "kind": "text", "present": True, "planned": False, "description": "前端派生口服标签"},
            {"name": "Derived_Workflow_Stage", "label": "Derived Workflow Stage", "kind": "text", "present": True, "planned": False, "description": "前端派生流程阶段"},
            {"name": "Derived_Star_Rating", "label": "Derived Star Rating", "kind": "text", "present": True, "planned": False, "description": "前端派生星级"},
        ],
    },
    {
        "id": "planned_flexibility",
        "title": "loop / 柔性子数据",
        "description": "正式 loop 注释与 loop 柔性结果已经独立接入，不再是待接入预留位。",
        "fields": [
            {"name": "flexibility_summary_version", "label": "Protein Flexibility Version", "kind": "text", "present": True, "planned": False, "description": "蛋白级柔性摘要版本"},
            {"name": "flexibility_dynamic_score_mean", "label": "Flexibility Dynamic Score Mean", "kind": "number", "present": True, "planned": False, "description": "蛋白级动态柔性平均分"},
            {"name": "flexibility_consensus_flexible_count", "label": "Flexible Loop Count", "kind": "number", "present": True, "planned": False, "description": "蛋白级 flexible loop 数量"},
            {"name": "loop_annotations_link", "label": "Loop Annotation Link", "kind": "text", "present": True, "planned": False, "description": "loop 级正式注释表已在独立子数据结构中接入"},
            {"name": "loop_flexibility_public_summary_link", "label": "Loop Flexibility Link", "kind": "text", "present": True, "planned": False, "description": "loop 级柔性轻量结论表已在独立子数据结构中接入"},
        ],
    },
]


row_count = len(rows)
with_structure = sum(1 for r in rows if r["_flags"]["hasStructure"])
with_uniprot = sum(1 for r in rows if r["_flags"]["hasSequenceUniProt"])
with_alphafold = sum(1 for r in rows if r["_flags"]["hasAlphaFold"])
with_immun = sum(1 for r in rows if r["_flags"]["hasImmunogenicity"])


def facet(field):
    return sorted({norm(r.get(field)) for r in rows if norm(r.get(field))})


stats = {
    "rowCount": row_count,
    "fieldCount": len(rows[0]) if rows else 0,
    "scaffoldCount": len(facet("Scaffold_Category")),
    "recordTypeCount": len(facet("Derived_Record_Type")),
    "scenarioCount": len(facet("Derived_Application_Scenario")),
    "workflowStageCount": len(facet("Derived_Workflow_Stage")),
    "topScaffolds": top_counts(rows, "Scaffold_Category", 10),
    "overallDistribution": top_counts(rows, "overall_final_judgement", 10),
    "withStructureCount": with_structure,
    "withUniProtCount": with_uniprot,
    "withAlphaFoldCount": with_alphafold,
    "withImmunogenicityCount": with_immun,
    "plannedFlexFieldCount": 2,
}


facet_values = {
    "scaffold": facet("Scaffold_Category"),
    "overall": facet("overall_final_judgement"),
    "bp3": facet("BP3_BCell_Level"),
    "mhci": facet("mhci_binding_level"),
    "mhcii": facet("mhcii_binding_level"),
    "reviewed": facet("UniProt_Reviewed"),
    "recordType": facet("Derived_Record_Type"),
    "kdRange": facet("Derived_Kd_Range"),
    "scenario": facet("Derived_Application_Scenario"),
    "oral": facet("Derived_Oral_Potential"),
    "workflow": facet("Derived_Workflow_Stage"),
    "stars": facet("Derived_Star_Rating"),
}


table_columns = [
    {"field": "protein_row_id", "label": "ID"},
    {"field": "_display_name", "label": "Entry"},
    {"field": "Scaffold_Category", "label": "Scaffold"},
    {"field": "Derived_Target_Label", "label": "Target"},
    {"field": "Derived_Kd_nM", "label": "Kd nM"},
    {"field": "overall_final_judgement", "label": "Overall"},
    {"field": "Derived_Star_Rating", "label": "Stars"},
]


meta = {
    "title": "FBBP 柔性主链结合蛋白数据库",
    "subtitle": "基于当前正式数据库表的蛋白级研究工作台，区分事实层、结构层、预测层、loop/柔性层与 provenance 层。",
    "generatedAt": "2026-05-14T20:00:00+08:00",
    "sourceCsv": r"<local_path_removed>",
    "schemaDoc": r"<local_path_removed>/database_schema_description.md",
    "ruleConfigVersion": "2026-05-14-final-db",
    "designReferences": [
        {"name": "UniProt", "url": "https://www.uniprot.org/uniprotkb"},
        {"name": "RCSB PDB", "url": "https://www.rcsb.org/"},
        {"name": "AlphaFold DB", "url": "https://alphafold.ebi.ac.uk/"},
    ],
    "workflowStages": [
        {"title": "1. 数据采集", "description": "多源文献、专利、UniProt、PDB 等来源事实层整合。"},
        {"title": "2. 靶点与亲和力整理", "description": "规范化靶点概念、物种变体、交互关系与亲和力信息。"},
        {"title": "3. 结构上下文确定", "description": "区分来源结构事实与数据库实际采用的活跃结构上下文。"},
        {"title": "4. 预测与扩展分析", "description": "整合 B-cell、MHC、E. coli、口服、溶解性、蛋白宇宙和柔性摘要。"},
    ],
    "featureModules": [
        {"title": "多维筛选", "description": "Scaffold、靶点、Kd、结构、免疫原性、口服潜力等联合筛选。"},
        {"title": "结构分层", "description": "明确区分来源结构事实与当前活跃结构上下文。"},
        {"title": "loop / 柔性分离", "description": "蛋白级主视图与 loop / 柔性子数据结构分开维护。"},
        {"title": "外部资源跳转", "description": "聚合 UniProt、PDB、AlphaFoldDB、ChEMBL、DrugBank、KEGG 入口。"},
    ],
    "ragDemoPrompts": [
        {"query": "找到人源高亲和力 EGFR 结合蛋白", "intent": "靶向发现", "suggestedFilters": {"target": "EGFR", "stars": "4", "sortBy": "kd"}},
        {"query": "筛选具有消化稳定性证据的 cyclotide", "intent": "口服递送", "suggestedFilters": {"scaffold": "cyclotide", "oral": "已有稳定性数据"}},
        {"query": "找适合后续问答训练的高质量记录", "intent": "数据集下载", "suggestedFilters": {"stars": "5", "overall": "High"}},
    ],
}

schema_tables = [
    {"name": "proteins", "layer": "A.实验/事实层", "grain": "protein", "description": "蛋白核心实体表"},
    {"name": "sources", "layer": "A.实验/事实层", "grain": "source", "description": "文献/专利/数据库来源表"},
    {"name": "protein_identifiers", "layer": "A.实验/事实层", "grain": "identifier", "description": "UniProt / PDB / Entry 等标识表"},
    {"name": "domains", "layer": "A.实验/事实层", "grain": "domain", "description": "domain / construct / scaffold 实体表"},
    {"name": "targets_conceptual", "layer": "A.实验/事实层", "grain": "target concept", "description": "靶点概念层"},
    {"name": "target_species_variants", "layer": "A.实验/事实层", "grain": "target variant", "description": "靶点物种变体层"},
    {"name": "interactions", "layer": "A.实验/事实层", "grain": "interaction", "description": "相互作用关系表"},
    {"name": "affinity_data", "layer": "A.实验/事实层", "grain": "affinity", "description": "亲和力数值表"},
    {"name": "digestive_assays", "layer": "A.实验/事实层", "grain": "assay", "description": "来源层消化实验事实表"},
    {"name": "functional_annotations", "layer": "A.实验/事实层", "grain": "annotation", "description": "功能注释表"},
    {"name": "structural_source_info", "layer": "B.Structure层", "grain": "source structure", "description": "来源结构事实"},
    {"name": "active_structure_context", "layer": "B.Structure层", "grain": "protein", "description": "当前活跃结构上下文"},
    {"name": "bcell_epitope_results", "layer": "C.预测/结果层", "grain": "protein", "description": "B-cell 结果"},
    {"name": "mhci_results", "layer": "C.预测/结果层", "grain": "protein", "description": "MHC-I 结果"},
    {"name": "mhcii_results", "layer": "C.预测/结果层", "grain": "protein", "description": "MHC-II 结果"},
    {"name": "immunogenicity_summary", "layer": "C.预测/结果层", "grain": "protein", "description": "综合免疫原性"},
    {"name": "ecoli_expression_results", "layer": "C.预测/结果层", "grain": "protein", "description": "E. coli 表达预测"},
    {"name": "oral_results", "layer": "C.预测/结果层", "grain": "protein", "description": "口服预测"},
    {"name": "solubility_results", "layer": "C.预测/结果层", "grain": "protein", "description": "溶解性/开发性结果"},
    {"name": "protparam_results", "layer": "C.预测/结果层", "grain": "protein", "description": "ProtParam 理化性质"},
    {"name": "ted_results", "layer": "C.预测/结果层", "grain": "protein", "description": "TED 蛋白宇宙结果"},
    {"name": "protrek_results", "layer": "C.预测/结果层", "grain": "protein", "description": "ProTrek 结果"},
    {"name": "foldseek_results", "layer": "C.预测/结果层", "grain": "protein", "description": "Foldseek 结果"},
    {"name": "plmsearch_results", "layer": "C.预测/结果层", "grain": "protein", "description": "PLMSearch 结果"},
    {"name": "loop_annotations", "layer": "D.loop/柔性层", "grain": "loop", "description": "正式 loop 注释表"},
    {"name": "loop_flexibility_public_summary", "layer": "D.loop/柔性层", "grain": "loop", "description": "loop 柔性轻量结论表"},
    {"name": "protein_flexibility_summary", "layer": "D.loop/柔性层", "grain": "protein", "description": "蛋白级柔性摘要"},
]

scaffold_descriptions = {
    "cyclotide": "小型环肽骨架，常与消化稳定性、口服潜力和特殊 loop 解释相关。",
    "adnectin": "工程化 FN3 类骨架，常见于靶向开发和蛋白工程场景。",
    "knottin": "富二硫键紧致骨架，结构稳定但表达/工程化路线常需单独评估。",
    "kunitz": "经典抑制剂样骨架，常具明确结构与功能线索。",
    "obody": "工程化小蛋白骨架，适合做结合与筛选展示。",
}


bundle = {
    "meta": meta,
    "stats": stats,
    "fieldOrder": list(rows[0].keys()) if rows else [],
    "tableColumns": table_columns,
    "facetValues": facet_values,
    "schemaGroups": schema_groups,
    "schemaTables": schema_tables,
    "scaffoldDescriptions": scaffold_descriptions,
    "plannedFields": [field for group in schema_groups for field in group["fields"] if field.get("planned")],
    "rows": rows,
}

loop_bundle = {"meta": {"title": "FBBP loop bundle", "generatedAt": meta["generatedAt"]}, "rows": loop_rows}
flex_bundle = {"meta": {"title": "FBBP flexibility bundle", "generatedAt": meta["generatedAt"]}, "rows": flex_rows}

os.makedirs(OUT, exist_ok=True)

with open(os.path.join(OUT, "database-data.js"), "w", encoding="utf-8", newline="\n") as f:
    f.write("window.__FLEX_DB__ = ")
    json.dump(bundle, f, ensure_ascii=False, separators=(",", ":"))
    f.write(";\n")

with open(os.path.join(OUT, "loop-data.js"), "w", encoding="utf-8", newline="\n") as f:
    f.write("window.__FBBP_LOOP_DATA__ = ")
    json.dump(loop_bundle, f, ensure_ascii=False, separators=(",", ":"))
    f.write(";\n")

with open(os.path.join(OUT, "flexibility-data.js"), "w", encoding="utf-8", newline="\n") as f:
    f.write("window.__FBBP_FLEX_DATA__ = ")
    json.dump(flex_bundle, f, ensure_ascii=False, separators=(",", ":"))
    f.write(";\n")

print("wrote database-data.js, loop-data.js, flexibility-data.js")
print("rows", len(rows), "loops", len(loop_rows), "flex", len(flex_rows))
