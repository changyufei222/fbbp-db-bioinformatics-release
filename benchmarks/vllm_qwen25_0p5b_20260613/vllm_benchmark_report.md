# vLLM OpenAI-Compatible Inference Benchmark

## Purpose

This benchmark validates that the FBBP Agent/RAG workflow can connect to a private OpenAI-compatible LLM inference backend, instead of relying only on external hosted model APIs.

Scope boundary:

> single-GPU smoke/formal benchmark for private OpenAI-compatible inference backend integration; not a high-concurrency serving or kernel optimization benchmark.

## Environment

| Item | Value |
|---|---|
| Inference server | vLLM OpenAI-compatible server |
| Model | qwen2.5-0.5b-instruct |
| GPU | NVIDIA GeForce RTX 4090 |
| Benchmark type | Warm-up-excluded formal request benchmark |
| Warm-up requests excluded | 3 |
| Measured requests | 50 |

## Formal Results

| Metric | Value |
|---|---:|
| OK requests | 50 |
| Errors | 0 |
| avg_latency_ms | 262.46 |
| p50_latency_ms | 329.93 |
| p95_latency_ms | 332.88 |
| min_latency_ms | 91.72 |
| max_latency_ms | 333.50 |
| completion_tokens | 4990 |
| total_tokens | 7040 |
| measured_wall_s | 13.13 |
| completion_tok_per_s_wall | 380.18 |
| total_tok_per_s_wall | 536.37 |

## Interpretation

The result confirms that the local/private vLLM service path is functional for an OpenAI-compatible inference backend:

- the service completed all 50 measured requests successfully;
- no benchmark requests failed;
- p50 and p95 latency were both approximately 330 ms in this single-GPU small-model setting;
- total wall-clock throughput reached 536.37 tokens/s.

This supports the engineering claim that the Agent/RAG workflow can be connected to a private model-serving backend. It should not be interpreted as evidence of production-scale high-concurrency serving, quantization optimization, tensor parallelism, or inference-kernel optimization.

## Public-Release Sanitization

Only the aggregate report and summary table are published in this repository. The raw per-request details file is not included because it has not been reviewed for possible request text, internal paths, host metadata, or private execution context.

