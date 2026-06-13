# FBBP Bioinformatics 数据与复现发布包

[English](./README.md) | **中文**

这是 FBBP 项目面向论文投稿、同行评审和结果复核的公开发布仓库。仓库把补充数据、字段字典、LLM 整理质控、RAG/Agent 基准、图源数据、分析代码和数据库模式分开管理，避免把“论文补充材料”和“软件工程仓库”混为一体。

## 建议阅读顺序

1. 查看 [supplementary_data_repository_inventory.csv](./supplementary_data_repository_inventory.csv)，了解全部发布对象。
2. 阅读 [Repository_Files_FBBP_DB_v1_release/README_SUPPLEMENTARY_REPOSITORY_FILES.md](./Repository_Files_FBBP_DB_v1_release/README_SUPPLEMENTARY_REPOSITORY_FILES.md)，理解主发布目录。
3. 使用 [DATA_DICTIONARY.tsv](./Repository_Files_FBBP_DB_v1_release/DATA_DICTIONARY.tsv) 解释字段。
4. 运行 [alidate_release.py](./Repository_Files_FBBP_DB_v1_release/scripts/validate_release.py) 做发布完整性检查。
5. 使用根目录 [checksums_sha256.txt](./checksums_sha256.txt) 核对文件。

## 内容分层

| 区域 | 用途 |
|---|---|
| Supplementary_Data_1_* | 字段级数据字典 |
| Supplementary_Data_2_* | LLM 提示词与整理质控 |
| Supplementary_Data_3_* | 固定 120 问题的 RAG/Agent 基准 |
| Supplementary_Data_4_* | 图源数据与代码清单 |
| Repository_Files_FBBP_DB_v1_release/ | 表、模式、分析结果、脚本和版本元数据 |

## 复现与边界

- 本仓库是投稿与复核入口，不等同于在线数据库服务。
- 本机绝对路径和敏感运行配置已从公开文件中移除。
- 结果解释应同时参考数据字典、分析状态说明和具体方法文件。
- 尚未分配 DOI；当前引用方式见 [CITATION.cff](./CITATION.cff)。

## 相关仓库

- [FBBP Research Workbench](https://github.com/changyufei222/fbbp-research-workbench)
- [FBBP MCP RAG Server](https://github.com/changyufei222/fbbp-mcp-rag-server)
- [LLM RAG Knowledge Base](https://github.com/changyufei222/llm-rag-knowledge-base)
- [LLM Eval Benchmark](https://github.com/changyufei222/llm-eval-benchmark)

## 工程作品集证据

[`benchmarks/vllm_qwen25_0p5b_20260613`](./benchmarks/vllm_qwen25_0p5b_20260613/) 目录提供一个脱敏后的 vLLM OpenAI-compatible 推理 benchmark，用作私有化模型后端接入能力的工程作品集证据。该目录与 FBBP Bioinformatics 补充数据 DOI 包保持分离：它是 single-GPU smoke/formal benchmark for private OpenAI-compatible inference backend integration; not a high-concurrency serving or kernel optimization benchmark.

详细界面与公开边界说明见 [INTERFACE_GUIDE_CN.md](./INTERFACE_GUIDE_CN.md)。
