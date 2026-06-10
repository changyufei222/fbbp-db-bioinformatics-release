# Round C2

模型配置：

- `Extractor = gemini-3.1-pro-preview`
- `Reviewer = deepseek-v3.2`
- `Arbiter = gemini-3.1-pro-preview`

设计目的：

- 保留 Gemini 的召回优势
- 用 DeepSeek Reviewer 形成独立第二意见
- 让 Arbiter 回到 Gemini，看最终裁决是否会更完整

这是 Round C 的对照配置。
