# Binding / Non-binding 完整分析汇总

## 1. 分析目标

本轮分析的目标不是再做一版探索性自动标注，而是基于已经确认的人工终版结果，系统回答两个问题：

1. `binding` 与 `non-binding` loop window 在静态暴露、静态约束和高通量柔性指标上到底有什么差异。
2. 哪些 scaffold 家族在更严格的 residue-level 证据下，真正支持“功能界面更埋藏”这一结论。

因此，本报告明确分成两层证据：

- `loop-window level`：基于人工终版 `binding / non-binding` loop 标签，对整库主分析集做窗口级比较。
- `strict residue-level`：只保留已经有 `target_chain_id` 的实验多链结构，在 loop 内部按原子接触重新标记 `binding residue / non-binding residue`，用来判断“界面更埋藏”是否成立。

## 2. 使用的数据与文件

### 2.1 loop-window 主分析

- 终版标签表：[loop_binding_labels_manual_final_v1.csv](<local_path_removed>/loop_binding_labels_manual_final_v1.csv)
- 主统计表：[binding_vs_nonbinding_primary_stats_manual_final_v1.csv](<local_path_removed>/binding_vs_nonbinding_primary_stats_manual_final_v1.csv)
- 主结果报告：[binding_vs_nonbinding_primary_report_manual_final_v1.md](<local_path_removed>/binding_vs_nonbinding_primary_report_manual_final_v1.md)
- 人工复核来源：[binding_subregion_manual_review_primary_binding_all_v1_reviewed_round2_with_ibody_7BOF_remap.csv](<local_path_removed>/binding_subregion_manual_review_primary_binding_all_v1_reviewed_round2_with_ibody_7BOF_remap.csv)

### 2.2 strict residue-level 验证

- residue 标签表：[binding_residue_labels_strict_v1.csv](<local_path_removed>/binding_residue_labels_strict_v1.csv)
- loop 汇总表：[binding_loop_residue_summary_strict_v1.csv](<local_path_removed>/binding_loop_residue_summary_strict_v1.csv)
- residue 统计表：[binding_vs_nonbinding_residue_stats_strict_v1.csv](<local_path_removed>/binding_vs_nonbinding_residue_stats_strict_v1.csv)
- residue 报告：[binding_residue_report_strict_v1.md](<local_path_removed>/binding_residue_report_strict_v1.md)
- 严格版构建脚本：[build_binding_residue_strict_v1.py](<local_path_removed>/build_binding_residue_strict_v1.py)

### 2.3 结构与注释输入

- strict 入组表：[loop_subregion_annotations_input_v1.csv](<local_path_removed>/loop_subregion_annotations_input_v1.csv)
- loop 主结果表：[loop_flexibility_results_long.csv](<local_path_removed>/loop_flexibility_results_long.csv)
- residue DSSP 表：[dssp_residue_long.csv](<local_path_removed>/dssp_residue_long.csv)

## 3. 分析方案

### 3.1 loop-window level：人工终版主分析

这层分析的核心思想是：先把 loop 窗口定义为 `binding window` 或 `non-binding window`，再比较不同窗口在多个柔性/约束指标上的分布差异。

人工终版标签的生成规则是：

1. 先用自动标签生成初版 `binding_loop_label`。
2. 再把人工复核结果并回去，优先级为：
   - `round2_review_decision`
   - `manual_review_decision`
   - 原始自动标签
3. `confirm_binding_window` 记为 `binding`
4. `downgrade_to_nonbinding_window` 记为 `non-binding`
5. `needs_more_evidence` 记为 `unresolved`，不纳入二元比较

窗口级主分析使用的是 `primary_loopcentric` 主分析集，并在每个指标上对 `binding` 和 `non-binding` 做 complete-case 比较。

比较的 4 个指标为：

- `avg_rsa`：窗口平均相对溶剂暴露
- `loop_constraint_score_total_norm`：GetContact 静态约束分
- `itsflex_score`：ITsFlexible 柔性预测分
- `flexscore_ensemble_diversity`：localcolabfold 动态系综多样性

### 3.2 strict residue-level：高可信界面验证

这层分析不是沿用 window 标签，而是在 loop 内部重新做 residue 级别的接触判定。

strict 版入组标准全部满足才保留：

1. `target_chain_id` 非空
2. `tbp_chain_id` 非空
3. `complex_structure_path` 非空且本地结构文件存在
4. 结构路径属于实验结构目录 `PDB_Structures`
5. 排除已知异常映射 `7BOF`

满足条件的 loop 共 `235` 条，对应 `100` 个实验结构。

residue-level 标注规则：

1. 对每个 loop，读取 `tbp_chain_id` 对应链，并取 `loop_start ~ loop_end` 范围内全部残基。
2. 若 `target_chain_id` 有多个候选链，则逐个计算：
   - loop 内接触残基数
   - 接触总数
   - 最小原子距离
3. 选择“接触残基数最多，其次接触总数最多，其次最小距离最小”的 target 链作为 `selected_target_chain_id`。
4. 定义 residue-level 接触规则：
   - loop 内某残基任一原子到 target 链任一原子最小距离 `<= 5.0 A`，记为 `binding residue`
   - 否则记为 `non-binding residue`
5. 再把每个 residue 与 DSSP residue 表合并，读取 `rel_asa`

这一层只比较一个最直接的界面埋藏指标：

- `rel_asa`

## 4. 统计口径

所有 `binding vs non-binding` 比较都使用双侧 `Mann-Whitney U` 检验。

报告的核心数字包括：

- `n_binding`
- `n_nonbinding`
- `binding_mean`
- `nonbinding_mean`
- `binding_median`
- `nonbinding_median`
- `delta_mean = binding_mean - nonbinding_mean`
- `delta_median = binding_median - nonbinding_median`
- `mannwhitney_p`

本报告对“真正支持界面更埋藏”的判定口径是：

- 在 strict residue-level 中，`distance-defined residues` 的 `rel_asa` 低于同一批 loop 中的非界面残基
- `distance-defined residues` 的 `ΔASA_proxy` 显著高于非界面残基
- 在 `ΔASA-supported` 和 `core-interface` 层里，这个低暴露趋势仍然保留
- 且样本量不算过小
- 且统计上能站住脚

因此，最终判断主要以 strict residue-level 为准，而不是只看 loop-window 的 `avg_rsa`。

## 5. 数据规模与终版标签概况

人工终版 loop 标签总量如下：

- 总 loop 数：`3383`
- 最终 `binding loop`：`1044`
- 最终 `non-binding loop`：`2278`
- 最终 `unresolved`：`61`
- 实际应用了人工复核的条目：`518`

`61` 条 unresolved 的 scaffold 分布为：

- `knottin = 37`
- `kunitz = 18`
- `Ibody = 5`
- `adnectin = 1`

按 scaffold 统计的终版 loop 数如下：

| Scaffold | Total loops | Binding loops | Non-binding loops | Unresolved |
| --- | ---: | ---: | ---: | ---: |
| adnectin | 1020 | 587 | 432 | 1 |
| obody | 702 | 50 | 652 | 0 |
| cyclotide | 479 | 5 | 474 | 0 |
| knottin | 370 | 178 | 155 | 37 |
| kunitz | 276 | 127 | 131 | 18 |
| EVH1domain | 185 | 5 | 180 | 0 |
| centyrin | 162 | 2 | 160 | 0 |
| affimer | 121 | 60 | 61 | 0 |
| Ibody | 68 | 30 | 33 | 5 |

## 6. loop-window level 总体结果

在人工终版 `primary_loopcentric` 主分析集中，整体比较结果如下：

| Metric | n_binding | n_nonbinding | Binding mean | Non-binding mean | Delta mean | P value |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| avg_rsa | 987 | 1286 | 0.5090 | 0.5564 | -0.0474 | 1.68e-21 |
| loop_constraint_score_total_norm | 987 | 1286 | 3.2016 | 3.7358 | -0.5342 | 3.68e-02 |
| itsflex_score | 966 | 1268 | 0.1595 | 0.1483 | 0.0112 | 1.09e-02 |
| flexscore_ensemble_diversity | 937 | 1235 | 2.0624 | 2.4598 | -0.3974 | 4.66e-02 |

这组结果说明：

1. `binding window` 在整体上更不暴露，`avg_rsa` 明显更低。
2. `binding window` 的静态约束分更低，而不是更高。
3. `binding window` 的 `ITsFlexible` 略高，但 `ensemble diversity` 更低。

因此，主结论不是“binding 区域统一更柔”，而是：

> binding 区域更像是“受功能约束的局部可塑界面”，而不是最大尺度乱动的自由 loop。

## 7. loop-window level scaffold 分层结果

窗口级别中，`avg_rsa` 的 scaffold 分层结果最值得看：

| Scaffold | n_binding | n_nonbinding | Binding mean avg_rsa | Non-binding mean avg_rsa | Delta mean | P value | 解释 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| adnectin | 587 | 432 | 0.5113 | 0.5358 | -0.0245 | 1.80e-02 | 支持 binding 窗口更不暴露 |
| affimer | 60 | 61 | 0.5758 | 0.5704 | 0.0054 | 9.52e-01 | 不支持 |
| cyclotide | 5 | 474 | 0.5275 | 0.5769 | -0.0493 | 6.27e-01 | 方向支持，但 binding 窗口太少 |
| Ibody | 30 | 33 | 0.5503 | 0.5684 | -0.0181 | 3.82e-01 | 方向支持，但不显著 |
| knottin | 178 | 155 | 0.5955 | 0.5618 | 0.0337 | 1.84e-05 | 与“更埋藏”相反 |
| kunitz | 127 | 131 | 0.3348 | 0.5340 | -0.1992 | 1.02e-43 | 极强支持 binding 窗口更不暴露 |

窗口级别最清楚的家族是：

- `kunitz`：最强的低暴露信号
- `adnectin`：中等但稳定的低暴露信号

窗口级别最需要谨慎的家族是：

- `knottin`：方向相反，binding 窗口整体更暴露
- `cyclotide`：binding window 数量只有 `5`，窗口级结论非常不稳定

## 8. strict residue-level 总体结果

strict residue-level 的入组规模如下：

- eligible loops：`235`
- eligible structures：`100`
- loop 内 residue 总数：`1832`
- `distance-defined residues`：`316`
- `non-binding residue`：`1516`
- `DSSP rel_asa` 匹配：`1832 / 1832`
- `ΔASA-supported residues (> 1 Å²)`：`294`
- `core-interface residues (>= 5 Å²)`：`220`
- `core-interface strict residues (>= 10 Å²)`：`199`

整体 strict residue-level 结果为：

| Layer | Metric | n_binding residues | n_nonbinding residues | Binding mean | Non-binding mean | Delta mean | P value |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| distance-defined | rel_asa | 316 | 1516 | 0.3475 | 0.4427 | -0.0952 | 1.46e-07 |
| distance-defined | ΔASA_proxy | 316 | 1516 | 27.56 | 0.06 | 27.50 | ~0 |
| ΔASA-supported | rel_asa | 294 | 1538 | 0.3482 | 0.4412 | -0.0930 | 7.88e-07 |
| core-interface | rel_asa | 220 | 1612 | 0.3534 | 0.4362 | -0.0828 | 1.49e-04 |
| core-interface strict | rel_asa | 199 | 1633 | 0.3483 | 0.4358 | -0.0875 | 1.24e-04 |

这意味着在高可信实验复合物子集中，真正与 target 发生接触的 loop 残基不仅平均 `rel_asa` 明显更低，而且其 `complex-derived ΔASA/BSA proxy` 远高于同一批 loop 中的非界面残基。也就是说，这些 residue 并不只是“碰到了 partner”，而是伴随了明确的界面埋藏面积变化。

换句话说：

> 当分析从 window-level 收缩到 strict residue-level 之后，“界面更埋藏”不仅没有消失，反而变得更清楚。

## 9. 哪些 scaffold 真正支持“界面更埋藏”

### 9.1 strict residue-level 覆盖规模

| Scaffold | Strict loops | Structures | Binding residues | Non-binding residues |
| --- | ---: | ---: | ---: | ---: |
| obody | 119 | 40 | 127 | 782 |
| adnectin | 46 | 26 | 59 | 266 |
| cyclotide | 20 | 15 | 67 | 129 |
| knottin | 20 | 7 | 31 | 153 |
| EVH1domain | 17 | 6 | 13 | 91 |
| affimer | 5 | 2 | 2 | 32 |
| kunitz | 5 | 3 | 10 | 48 |
| centyrin | 3 | 3 | 7 | 15 |
| Ibody | 0 | 0 | 0 | 0 |

### 9.2 strict residue-level 支持度判断

| Scaffold | Binding mean rel_asa | Non-binding mean rel_asa | Delta mean | P value | 判断 |
| --- | ---: | ---: | ---: | ---: | --- |
| adnectin | 0.3169 | 0.4787 | -0.1618 | 9.38e-05 | 强支持 |
| cyclotide | 0.3569 | 0.4712 | -0.1143 | 1.16e-02 | 强支持 |
| obody | 0.3551 | 0.4363 | -0.0812 | 4.10e-03 | 强支持 |
| kunitz | 0.3043 | 0.4319 | -0.1275 | 2.20e-01 | 方向支持，但样本不足 |
| EVH1domain | 0.4031 | 0.4784 | -0.0753 | 3.53e-01 | 方向支持，但样本不足 |
| knottin | 0.3171 | 0.3768 | -0.0598 | 3.54e-01 | 方向支持，但样本不足 |
| centyrin | 0.3954 | 0.4114 | -0.0160 | 8.91e-01 | 方向极弱，暂不主张 |
| affimer | 0.6173 | 0.4293 | 0.1880 | 4.00e-01 | 不支持 |
| Ibody | NA | NA | NA | NA | 无 strict 证据 |

### 9.3 最终结论

如果问题是“哪些 scaffold 真正支持界面更埋藏”，当前最稳的答案是：

- `adnectin`
- `cyclotide`
- `obody`

这三个家族在 strict residue-level 子集中都表现为：

- `distance-defined residue rel_asa < non-interface residue rel_asa`
- 进入 `ΔASA-supported` 之后方向仍然保留
- 进入 `core-interface` 层后仍然保持更低暴露
- 效应量明确
- 统计上站得住

第二层答案是“方向支持，但 strict 样本还不够硬”的家族：

- `kunitz`
- `EVH1domain`
- `knottin`

其中：

- `kunitz` 在 loop-window level 已经非常强地支持更埋藏，而 strict residue-level 三层方向也一致，只是当前 strict 样本只有 `10 vs 48` 个 residue，证据仍偏薄。
- `knottin` 最有意思：window-level 上 binding 窗口更暴露，但 strict residue-level 内部接触残基方向却是更埋藏。更合理的解释不是“knottin 不埋藏”，而是“它的功能 loop 整体较暴露，但真正接触靶标的那部分残基并不一定更暴露”，同时也不能排除家族内部异质性和标签噪音。

当前不支持直接写成“界面更埋藏”的家族：

- `affimer`
- `Ibody`
- `centyrin`

## 10. 综合解释

把 loop-window 和 strict residue-level 两层放在一起看，最合理的总解释是：

1. 全库层面，`binding` 区域并不是简单“更柔、更乱、更暴露”。
2. 在窗口级别，binding 区域整体更不暴露，但不同柔性指标的方向并不一致。
3. 在真正有 `target_chain_id` 的 strict residue-level 子集里，真实接触残基的 `rel_asa` 明显更低，同时 `ΔASA_proxy` 显著更高，说明“距离定义的界面接触残基确实伴随界面埋藏”这个命题在高证据子集中是成立的。
4. 这种埋藏并不是所有 scaffold 完全一致的统一规律，而是 scaffold-dependent：
   - `adnectin / cyclotide / obody` 证据最完整
   - `kunitz` 很可能也成立，但 strict 样本偏少
   - `knottin` 呈现“窗口暴露但接触子位点未必暴露”的复杂模式

因此，不建议把结论写成“所有 binding loop 都更埋藏、更柔”。更准确的写法是：

> 功能性 binding 区域表现为 scaffold-dependent constrained plasticity；在高可信复合物子集中，真实界面残基总体上更埋藏，但不同 scaffold 的窗口尺度暴露和柔性模式并不相同。

## 11. 风险点与限制

本轮结果已经比纯 window-level 强很多，但仍然有 4 个限制必须明说：

1. strict residue-level 只覆盖 `235 / 3383` 个 loop，覆盖率有限。
2. strict 子集并不均衡，`obody`、`adnectin` 占比很高。
3. `Ibody` 当前没有 strict residue-level 证据，不能对其界面埋藏做强结论。
4. `knottin` 和 `kunitz` 的 strict 结果方向上有支持，但样本量仍偏小，不适合写成 definitive family rule。

## 12. 论文式结果段落

### 12.1 结果主段落

在人工复核后的 loop-window 主分析集中，共获得 `3383` 条功能 loop，其中 `1044` 条被标记为 binding windows，`2278` 条被标记为 non-binding windows，另有 `61` 条因证据不足保留为 unresolved。基于 `primary_loopcentric` scaffold 的窗口级比较显示，binding windows 的平均相对溶剂暴露显著低于 non-binding windows（`avg_rsa`: `0.5090` vs `0.5564`, Mann-Whitney `p = 1.68e-21`），同时其静态约束分也更低（`3.2016` vs `3.7358`, `p = 3.68e-02`）。相比之下，binding windows 的 `ITsFlexible` 分数略高（`0.1595` vs `0.1483`, `p = 1.09e-02`），但 `localcolabfold` 系综多样性反而更低（`2.0624` vs `2.4598`, `p = 4.66e-02`）。这些结果表明，功能性 binding 区域并不对应于最自由、最暴露的 loop，而更可能对应于受功能约束的局部可塑界面。

### 12.2 strict residue-level 验证段落

为了进一步验证“界面更埋藏”是否在 residue-level 成立，我们构建了一个严格的复合物子集，只保留具有显式 `target_chain_id`、真实实验多链结构且本地结构文件存在的 loop，并排除了异常映射结构。该 strict 子集共包含 `235` 条 loop、`100` 个实验结构和 `1832` 个 loop 内残基，其中 `316` 个残基被定义为 `distance-defined` 界面残基，`1516` 个残基被定义为同一批 loop 中的非界面残基。在此基础上，我们进一步计算了 `complex-derived ΔASA/BSA proxy`，并定义出 `294` 个 `ΔASA-supported` 残基和 `220` 个 `core-interface` 残基。基于 DSSP `rel_asa` 的比较显示，distance-defined 残基的平均相对暴露显著低于 non-interface residues（`0.3475` vs `0.4427`, `delta = -0.0952`, `p = 1.46e-07`），同时其 `ΔASA_proxy` 极高（`27.56 Å²` vs `0.06 Å²`）。而且，这种低暴露趋势在 `ΔASA-supported` 和 `core-interface` 层级中仍然保留，说明当分析尺度从 window-level 收缩到真实接触残基时，功能界面的埋藏特征不仅更加清晰，而且获得了面积埋藏证据的直接支撑。

### 12.3 scaffold 分层段落

scaffold 分层分析进一步表明，这一埋藏特征具有明显的家族依赖性。在 strict residue-level 子集中，`adnectin`（distance-defined `0.3169` vs `0.4787`, `p = 9.38e-05`；core-interface `0.3399` vs `0.4669`, `p = 6.18e-03`）、`cyclotide`（distance-defined `0.3569` vs `0.4712`, `p = 1.16e-02`）和 `obody`（distance-defined `0.3551` vs `0.4363`, `p = 4.10e-03`；core-interface `0.3437` vs `0.4358`, `p = 2.86e-03`）均明确支持“界面残基更埋藏”的模式。`kunitz`、`EVH1domain` 和 `knottin` 虽然在 strict residue-level 中也表现出更低的 `rel_asa`，并在 `ΔASA-supported/core-interface` 层保持相同方向，但由于样本量仍有限，当前更适合表述为方向支持而非最终定论。值得注意的是，`knottin` 在 loop-window level 中反而表现出 binding windows 更高的平均暴露，但在 strict residue-level 内部却呈现相反方向，提示该家族可能存在“整体 loop 暴露但真正接触子位点相对埋藏”的复杂界面模式。

### 12.4 结论段落

综合窗口级与 residue-level 证据，我们认为本研究不支持将功能性 binding loop 简单定义为“全局更柔或更暴露”的区域。相反，功能界面更符合一种 scaffold-dependent constrained plasticity 模式：即在不同 scaffold 中，binding 区域可以表现出不同的窗口尺度暴露和动态特征，但在高可信复合物子集中，真实接触残基总体上呈现更低的溶剂暴露，并伴随显著的 `complex-derived ΔASA/BSA proxy`，说明功能界面在分子尺度上更倾向于形成受约束的局部埋藏识别面。

