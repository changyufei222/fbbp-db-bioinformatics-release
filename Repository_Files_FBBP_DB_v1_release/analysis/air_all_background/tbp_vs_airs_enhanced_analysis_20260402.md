# TBP vs AIRs 强化版对照分析 v2

生成时间：2026-04-02T21:38:07

本版不是简单重复首版的长度分层统计，而是在同一份 `TBP` loop 结果与 `AIR ALL_conformations` 参考集之上，额外补入了三类强化对照：第一类是长度标准化后的 matched-control 估计，也就是用 `TBP` 自身的 loop 长度分布去加权 `AIR` 的柔性比例，从而避免“TBP 只是因为 loop 更长所以显得更柔”的解释；第二类是长度校正 logistic 模型，把 `dataset_is_tbp` 与 `loop_length`、`loop_length^2` 一起纳入，用协变量校正后的数据集效应来衡量 `TBP` 相对 `AIR` 的独立柔性倾向；第三类是 exact-sequence overlap 分析，用来识别那些与 `AIR` 已知 loop 序列完全重合的 `TBP` loop 在动态指标上是否表现出系统性偏移。

从当前结果看，`TBP` 相对于 `AIR` 的“偏柔”结论在强化版对照下仍然稳定。长度标准化之后，`TBP_all_binary` 与 `TBP_primary_approved` 两个 cohort 的观测柔性比例都高于按同长度分布推算的 `AIR` 期望柔性比例；也就是说，即使强行把 `AIR` 拉成与 `TBP` 一样的长度构成，`TBP` 仍然更偏柔。长度校正 logistic 模型也给出同方向结论，`dataset_is_tbp` 的系数为正，对应的 odds ratio 大于 `1`，说明在控制 loop 长度的情况下，来自 `TBP` 数据集本身仍然显著提高“被判为 flexible”的概率。

这份强化版也补充了 `AIR PDB-only` 对照，因此可以区分“与 AIR 整体参考相比更柔”和“与 AIR 中结构证据最强的 PDB 子集相比更柔”这两种口径。如果后续论文要写得更稳，建议优先引用 `TBP vs AIR_PDB` 的长度校正结果，把它作为对结构质量差异更保守的一版主分析。

另外，exact-sequence overlap 的结果提示，和 `AIR seq_loop` 完全重合的 `TBP` loop 在 dynamic 指标上往往低于未匹配条目。这更像是在说明：`AIR` 中容易反复出现的是一些保守短 loop 模体，而 `TBP` 数据集里真正驱动高动态性的部分，更可能来自那些没有在 `AIR` 中直接出现 exact overlap 的功能 loop。

需要单独强调的是，`AIR` 侧目前只有官方 loop 级柔/刚标签，没有与 `TBP` 完全同构的 `DSSP / GetContact / localcolabfold whole-loop dynamic` 全套原始指标。因此本版强化对照已经把“能严格对齐的部分”推进到较完整状态：包括长度匹配、PDB-only 子集、exact-sequence overlap 和协变量校正；但它仍然不是 residue-level 或 subregion-level 的一一对应分析。这一点要在方法学里明确写出来，避免把 AIR 比较过度表述成 residue-resolution 的直接验证。
