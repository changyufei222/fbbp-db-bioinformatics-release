# 数据库 Schema 表结构清单

更新时间：2026-05-14

## 1. 设计原则

- 实验/事实层 与 预测/结果层严格分开
- Structure 层单独拆开，区分来源事实结构和当前活跃结构上下文
- loop/柔性层保持 long table，不直接硬拼进蛋白级母表
- 审计/Provenance 层单独保留，避免污染主业务表

## 2. 表结构清单

### A.实验/事实层

| table_name | primary_key | foreign_keys | 是否从当前总表直接拆 | 备注 |
| --- | --- | --- | --- | --- |
| proteins | protein_id | none | 否，需从主总表拆分并去重 | 蛋白核心实体表，1 protein 1 row |
| sources | source_id | none | 否，需从主总表拆分并去重 | 来源文献/专利/数据库表 |
| protein_identifiers | identifier_id | protein_id -> proteins.protein_id | 否，需从主总表拆分并去重 | UniProt/PDB/Entry 等标识表 |
| domains | domain_id | protein_id -> proteins.protein_id; source_id -> sources.source_id | 否，需从主总表拆分并去重 | domain / construct / scaffold 实体表 |
| targets_conceptual | target_concept_id | none | 否，需从主总表拆分并去重 | 靶点概念层 |
| target_species_variants | target_variant_id | target_concept_id -> targets_conceptual.target_concept_id | 否，需从主总表拆分并去重 | 靶点物种变体层 |
| interactions | interaction_id | domain_id -> domains.domain_id; target_variant_id -> target_species_variants.target_variant_id; source_id -> sources.source_id | 否，需从主总表拆分并去重 | 结合对象与作用关系表 |
| affinity_data | affinity_id | interaction_id -> interactions.interaction_id | 否，需从主总表拆分并去重 | 亲和力数值表 |
| digestive_assays | assay_id | domain_id -> domains.domain_id | 否，需从主总表拆分并去重 | 来源层消化实验事实表 |
| functional_annotations | annotation_id | protein_id -> proteins.protein_id | 否，需从主总表拆分并去重 | GO/Keyword/Function 等注释表 |

### B.Structure层

| table_name | primary_key | foreign_keys | 是否从当前总表直接拆 | 备注 |
| --- | --- | --- | --- | --- |
| structural_source_info | structure_id | domain_id -> domains.domain_id | 否，需从主总表拆分并去重 | 来源层结构事实：pdb_id/method/resolution |
| active_structure_context | protein_row_id + structure_unique_sequence_id | none | 是，直接来自预测前主表 | 数据库实际采用的活跃结构上下文 |

### C.预测/结果层

| table_name | primary_key | foreign_keys | 是否从当前总表直接拆 | 备注 |
| --- | --- | --- | --- | --- |
| bcell_epitope_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | B 细胞表位结果 |
| mhci_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | MHC-I 结果 |
| mhcii_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | MHC-II 结果 |
| immunogenicity_summary | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | 综合免疫原性汇总 |
| ecoli_expression_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | E. coli 表达预测结果 |
| oral_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | 口服预测结果 |
| solubility_results | protein_row_id + structure_unique_sequence_id | none | 是，已在主总表中完成映射 | QCBundle / solubility 结果 |
| protparam_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | ProtParam 理化性质结果 |
| ted_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | TED 蛋白宇宙结果 |
| protrek_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | ProTrek 结果 |
| foldseek_results | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | Foldseek 结果 |
| plmsearch_results | protein_row_id + structure_unique_sequence_id | none | 是，已在主总表中完成映射 | PLMSearch 结果 |

### D.loop/柔性层

| table_name | primary_key | foreign_keys | 是否从当前总表直接拆 | 备注 |
| --- | --- | --- | --- | --- |
| loop_annotations | protein_row_id + structure_unique_sequence_id + loop_id | none | 是，独立 long table 直接使用 | 数据库正式 loop 表 |
| loop_flexibility_public_summary | protein_row_id + structure_unique_sequence_id + loop_id | none | 是，独立 long table 直接使用 | loop 级柔性轻量结论表 |
| protein_flexibility_summary | protein_row_id + structure_unique_sequence_id | none | 是，独立结果表直接使用 | 蛋白级柔性摘要 |

### E.审计/Provenance层

| table_name | primary_key | foreign_keys | 是否从当前总表直接拆 | 备注 |
| --- | --- | --- | --- | --- |
| trusted_source_qc | protein_row_id + structure_unique_sequence_id 或 audit_id | none | 部分已有，建议后续单独整理 | 可信来源抽样 QC 结果 |
| manual_review_audit | protein_row_id + structure_unique_sequence_id 或 audit_id | none | 部分已有，建议后续单独整理 | 人工审核痕迹表 |
| structure_refinement_audit | protein_row_id + structure_unique_sequence_id | none | 否，需从 QC 真实最终结构细化层单独整理 | 结构细化前后对照审计 |
| pipeline_versions | pipeline_name + source_table | none | 否，建议后续手工整理 | 流程版本表 |
| result_metadata | key | none | 否，建议后续手工整理 | 结果元数据表 |
