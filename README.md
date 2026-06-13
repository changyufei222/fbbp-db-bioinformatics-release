[**English**](./README.md) | [中文](./README_CN.md)

# FBBP Bioinformatics Data and Reproducibility Release

[![Repository quality](https://github.com/changyufei222/fbbp-db-bioinformatics-release/actions/workflows/repository-quality.yml/badge.svg)](https://github.com/changyufei222/fbbp-db-bioinformatics-release/actions/workflows/repository-quality.yml)
[![Release](https://img.shields.io/github/v/release/changyufei222/fbbp-db-bioinformatics-release?display_name=tag)](https://github.com/changyufei222/fbbp-db-bioinformatics-release/releases)
[![Citation](https://img.shields.io/badge/citation-CITATION.cff-2f6f9f.svg)](./CITATION.cff)

Publication-facing supplementary data, repository files, source-data manifests, schemas, and reproducibility assets for the FBBP project.

**Status:** Submission-oriented public release | **Public release:** 2026-06-10

| Start here | Resource |
|---|---|
| Primary documentation | [Release inventory](./supplementary_data_repository_inventory.csv) |
| Reproducibility / implementation | [Repository files](./Repository_Files_FBBP_DB_v1_release/) |
| Verified outcomes | [Checksums](./checksums_sha256.txt) |

![PCSK9 case-study evidence chain](./Repository_Files_FBBP_DB_v1_release/figure_source_data/pcsk9_case_figure/Figure_5_PCSK9_1459D05_evidence_chain.png)

---
Package structure:
- `Supplementary_Data_1_Field_Level_Data_Dictionary`: field-level data dictionary.
- `Supplementary_Data_2_LLM_Prompts_and_Curation_QC`: full LLM prompt archive and sanitized LLM curation-QC metrics.
- `Supplementary_Data_3_RAG_Agent_Benchmark`: fixed 120-question RAG/agent benchmark set and evaluation tables.
- `Supplementary_Data_4_Figure_Source_Data_Code_Manifest`: sanitized figure source-data/code manifest.
- `Repository_Files_FBBP_DB_v1_release`: release tables, schema, analysis outputs, code/scripts and release metadata.

Absolute local filesystem paths were removed from generated public-facing files. Sanitized build provenance is provided without local source paths, and package-level checksums are provided in `checksums_sha256.txt`.

## Engineering Portfolio Evidence

The [`benchmarks/vllm_qwen25_0p5b_20260613`](./benchmarks/vllm_qwen25_0p5b_20260613/) directory provides a sanitized vLLM OpenAI-compatible inference benchmark used as engineering portfolio evidence for private model-backend integration. This directory is intentionally separate from the FBBP Bioinformatics supplementary-data DOI package: it is a single-GPU smoke/formal benchmark for private OpenAI-compatible inference backend integration, not a high-concurrency serving or kernel optimization benchmark.

