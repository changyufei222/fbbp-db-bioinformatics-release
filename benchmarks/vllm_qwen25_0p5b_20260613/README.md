# vLLM Qwen2.5-0.5B Benchmark Evidence

This directory provides a sanitized engineering benchmark record for validating private LLM inference backend integration with the FBBP Agent/RAG workflow.

Scope boundary:

> single-GPU smoke/formal benchmark for private OpenAI-compatible inference backend integration; not a high-concurrency serving or kernel optimization benchmark.

## What This Shows

- A vLLM OpenAI-compatible inference server was deployed on a single NVIDIA GeForce RTX 4090 node.
- The local service was exercised through OpenAI-compatible chat-completion requests.
- Qwen2.5-0.5B-Instruct completed a warm-up-excluded formal benchmark with 50 measured requests.
- The benchmark completed with 50 successful requests and 0 errors.

## What This Does Not Claim

- This is not a multi-GPU, tensor-parallel, or production-scale serving benchmark.
- This is not an FP8, quantization, CUDA kernel, or vLLM engine optimization result.
- This does not claim high-concurrency throughput or production SLO readiness.

## Files

| File | Description |
|---|---|
| `vllm_benchmark_report.md` | Human-readable benchmark report and interpretation boundary. |
| `vllm_benchmark_summary.csv` | Sanitized one-row benchmark summary table. |

The raw `vllm_benchmark_details.jsonl` file is intentionally not published here because raw per-request records may include request text, internal execution paths, host metadata, or other environment-specific details. A redacted sample can be added later after reviewing the original remote file.

## Result Summary

| Metric | Value |
|---|---:|
| Measured requests | 50 |
| OK requests | 50 |
| Errors | 0 |
| p50 latency | 329.93 ms |
| p95 latency | 332.88 ms |
| Total throughput | 536.37 tokens/s |

