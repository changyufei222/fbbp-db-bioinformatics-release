# FBBP Supplementary Data / repository 最终包

本包用于把编号 Supplementary Tables 之外的机器可读文件统一打包。

目录结构：
- `Supplementary_Data_1_Field_Level_Data_Dictionary`：完整字段字典。
- `Supplementary_Data_2_LLM_Prompts_and_Curation_QC`：完整 LLM prompts 和整理/QC 指标。
- `Supplementary_Data_3_RAG_Agent_Benchmark`：固定 120 题 RAG/agent benchmark 题集及评测表。
- `Supplementary_Data_4_Figure_Source_Data_Code_Manifest`：去本机路径后的 figure source-data/code manifest。
- `Repository_Files_FBBP_DB_v1_release`：数据库 release tables、schema、分析结果、代码脚本、manifest 和 checksums。

本包中新生成或清理的 manifest 已去除本地电脑绝对路径。最终校验以根目录 `checksums_sha256.txt` 为准。
