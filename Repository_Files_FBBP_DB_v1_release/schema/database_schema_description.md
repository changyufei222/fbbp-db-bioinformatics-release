# 数据库 Schema 说明

更新时间：2026-05-14

## 1. 文档定位

本文件是当前 `<local_path_removed>

它整合并替代了此前的阶段性拆表说明，后续以本文件为准。

## 2. 当前设计总览

当前数据库采用 5 层设计：

1. 实验 / 事实层
2. Structure 层
3. 预测 / 结果层
4. loop / 柔性层
5. 审计 / Provenance 层

设计原则：

- 实验/事实层 与 预测/结果层严格分开
- Structure 层单独拆开，区分来源事实结构和当前活跃结构上下文
- loop/柔性层保持 long table，不直接硬拼进蛋白级母表
- 审计/Provenance 层单独保留，避免污染主业务表

## 3. 当前正式表概况

### A.实验/事实层

- `sources.csv`：`1214` 行，`8` 列
- `proteins.csv`：`1996` 行，`9` 列
- `protein_identifiers.csv`：`1270` 行，`5` 列
- `domains.csv`：`1996` 行，`9` 列
- `targets_conceptual.csv`：`196` 行，`3` 列
- `target_species_variants.csv`：`213` 行，`5` 列
- `interactions.csv`：`1996` 行，`8` 列
- `affinity_data.csv`：`1996` 行，`4` 列
- `digestive_assays.csv`：`17` 行，`4` 列
- `functional_annotations.csv`：`1753` 行，`4` 列

### B.Structure层

- `structural_source_info.csv`：`1996` 行，`9` 列
- `active_structure_context.csv`：`1996` 行，`20` 列

### C.预测/结果层

- `bcell_epitope_results.csv`：`1996` 行，`12` 列
- `mhci_results.csv`：`1996` 行，`22` 列
- `mhcii_results.csv`：`1996` 行，`17` 列
- `immunogenicity_summary.csv`：`1996` 行，`5` 列
- `ecoli_expression_results.csv`：`1996` 行，`24` 列
- `oral_results.csv`：`1996` 行，`21` 列
- `solubility_results.csv`：`1996` 行，`32` 列
- `protparam_results.csv`：`1996` 行，`9` 列
- `ted_results.csv`：`1996` 行，`41` 列
- `protrek_results.csv`：`1996` 行，`16` 列
- `foldseek_results.csv`：`1996` 行，`17` 列
- `plmsearch_results.csv`：`1996` 行，`27` 列

### D.loop/柔性层

- `loop_annotations.csv`：`3383` 行，`16` 列
- `loop_flexibility_public_summary.csv`：`3383` 行，`27` 列
- `protein_flexibility_summary.csv`：`1996` 行，`22` 列

### E.审计/Provenance层

- `trusted_source_qc.csv`：`1996` 行，`10` 列
- `manual_review_audit.csv`：`1996` 行，`4` 列
- `structure_refinement_audit.csv`：`97` 行，`16` 列
- `pipeline_versions.csv`：`5` 行，`4` 列
- `result_metadata.csv`：`7` 行，`2` 列

## 4. 关键结构细化结果

本轮最终刷新后，以下优化已写入正式表与 DDL：

- `sources.source_type` 已细化为 `Literature / Patent / UniProt / PDB`
- `targets_conceptual / target_species_variants` 已按规范化后的 gene + uniprot 体系重建
- `interactions` 已按新的 target variant 映射重建
- `digestive_assays` 保持长表模式，保留来源层消化实验语义
- `structural_source_info` 已补充：
  - `source_structure_identifier`
  - `source_structure_type`
  - `source_structure_note`

## 5. Schema 清单文件说明

当前保留两份清单：

### `database_schema_table_inventory.csv`

作用：

- 机器可读清单
- 适合程序读取、自动核对、后续生成脚本

字段包括：

- 表名
- 层级
- 主键
- 主要外键
- 是否从当前总表直接拆
- 备注

### `database_schema_table_inventory.md`

作用：

- 人类可读清单
- 适合直接阅读、汇报和沟通

所以这两份不是冲突版本，而是：

- 一份给程序
- 一份给人看

## 6. ER / SQL 正式资产

### Core ER

- `protein_database_er_core.sql`
- `protein_database_er_core.mmd`

### Full Schema

- `protein_database_er_full.sql`
- `protein_database_er_full_sqlite.sql`
- `protein_database_er_full_mysql.sql`
- `protein_database_er_full.mmd`

说明：

- `protein_database_er_core.sql`：Core 层的**通用 ER 解析友好版**，适合 AI 检查、ER 工具解析、画图，不作为最终导库文件
- `protein_database_er_full.sql`：Full Schema 的**通用 ER 解析友好版**，适合 AI 检查、ER 工具解析、画图，不作为最终导库文件
- `protein_database_er_full_sqlite.sql`：SQLite 真正建库推荐版
- `protein_database_er_full_mysql.sql`：MySQL 真正建库推荐版
- `protein_database_er_full.mmd`：Full Schema 关系图

## 6.1 SQL 文件分工

为避免把“画图/检查用 SQL”和“真正导库用 SQL”混用，当前正式分工如下：

### 用于 AI / ER 工具 / 画图

- `protein_database_er_core.sql`
- `protein_database_er_full.sql`

特点：

- 更偏中性 DDL
- 更适合 ER 解析器、AI 检查器、画图工具读取
- 不强调 SQLite / MySQL 方言细节

### 用于真正建 SQLite 库

- `protein_database_er_full_sqlite.sql`

特点：

- 已按 SQLite 方言整理
- 已实际通过 SQLite 解析验证

### 用于真正建 MySQL 库

- `protein_database_er_full_mysql.sql`

特点：

- 已按 MySQL 方言整理
- 使用反引号、VARCHAR、INT、DOUBLE、TINYINT(1) 等更适合 MySQL 的定义

## 7. 当前主键 / 外键约定

- 蛋白级结果层统一使用：`protein_row_id + structure_unique_sequence_id`
- loop 级结果层统一使用：`protein_row_id + structure_unique_sequence_id + loop_id`
- `active_structure_context` 是蛋白级结果层的核心锚点表
- `loop_annotations` 是 loop 级结果层的核心锚点表