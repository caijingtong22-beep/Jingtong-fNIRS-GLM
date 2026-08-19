# Q8 Procedure 完成记录

更新日期：2026-08-19  
分析目录：`results/q8_hihe_vs_lile_procedure/cohort_specification`

## 研究问题与本次范围

Q8 的 Procedure 阶段比较高情绪强度与低情绪强度：**HIHE - LILE**。Procedure regressor 使用完整的 `LILEProcedure` 与 `HIHEProcedure` 时间段；按照既定决策，DOM/AOE 时间段没有从 Procedure 中剔除。DOM 与 AOE 将以独立 GLM 处理，而不是作为本模型中删去的时段。

本记录覆盖目前已完成的 Q8 Procedure 一阶 GLM、通道 contrast、n=30 组水平汇总、long-45 ROI 汇总和图像输出。

## 已完成流程

1. **Specification（n=30）**
   - 每位受试者、每种血红蛋白建立 canonical-HRF GLM。
   - 主要任务设计矩阵：`[LILEProcedure, HIHEProcedure, constant]`。
   - 受试者 001 使用已有 V2 输入；其他受试者采用既定输入策略。

2. **Beta estimation（n=30）**
   - 对 HbO 与 HbR 分别估计每个通道的 beta。

3. **Contrast（n=30）**
   - 对比向量：`[-1, +1, 0]`，即 **HIHEProcedure - LILEProcedure**。
   - 导出每位受试者、HbO/HbR、各通道的 beta、SE、t、p、BH-FDR q 等结果。

4. **组水平通道汇总（n=30）**
   - 每通道进行 one-sample t-test，并在每种血红蛋白的 61 个通道内作 BH-FDR 校正。
   - 结果：HbO 与 HbR 均没有通道达到 BH-FDR 显著。

5. **long-45 ROI 汇总与作图**
   - 只使用 45 条 long channels，保留全部 QC 标记；不重新按 QC 标记剔除通道。
   - broad ROI × laterality：以受试者的 ROI 平均 beta 为组水平检验单位。
   - cortical-system 版本：将 lateral/medial 分区等权平均为一个 ROI，并标记为 Prefrontal 或 Motor；该版本仅作辅助汇总，主图使用 broad ROI × laterality。
   - 结果：broad ROI × laterality 与 cortical-system ROI 均无 BH-FDR 显著结果。
   - 最小未校正 ROI 结果为 HbO Medial Frontopolar PFC：合并 ROI 后 `p = .039`、`q = .309`，因此不能报告为校正后的显著效应。

6. **Subject 001 敏感性检查（不作为主推断）**
   - 加入同一血红蛋白的 short-channel 均值、线性至三次漂移项，并使用 channel-wise feasible AR(1) GLS。
   - 该模型用于检查 beta 对 short-channel nuisance 与简单自相关处理的敏感性；并未重跑 n=30。
   - AR(1) 后仍存在残余自相关，故此结果不作为最终显著性推断；主结论仍基于上述全样本 NIRS-SPM 流程。

## 当前结论

在当前的 n=30 Q8 Procedure 分析中，没有足以通过多重比较校正的 HIHE - LILE 通道或 ROI 效应。因此，此阶段的合适表述是：**未发现经 BH-FDR 校正后稳定的 Procedure 阶段 HIHE 与 LILE 差异。**

HbO 与 HbR 在不少通道呈相同符号的现象仍存在；已确认 `structured_data` 的两个 plane 分别对应 OXY 和 DEOXY。简单的标签互换不能解释“同向”本身，仅会改变哪一项被命名为 HbO/HbR。该问题应作为数据质量/生理噪声与预白化限制讨论，不应据此作典型神经血管耦合的方向性解释。

## 主要输出

### 组水平 CSV

- `q8_HIHE_minus_LILE_procedure_group_channel_summary_n30.csv`
- `q8_HIHE_minus_LILE_procedure_group_overview_n30.csv`
- `q8_HIHE_minus_LILE_procedure_group_channel_long45_allQC_n30.csv`
- `q8_HIHE_minus_LILE_procedure_group_ROI_long45_allQC_n30.csv`
- `q8_HIHE_minus_LILE_procedure_group_ROI_cortical_system_long45_allQC_n30.csv`

### 主图（long-45，broad ROI × laterality）

- `figures_long45_allQC/q8_procedure_layer1_individual_roi_heatmaps_long45_allQC.png`
- `figures_long45_allQC/q8_procedure_layer2_group_roi_effects_long45_allQC.png`
- `figures_long45_allQC/q8_procedure_layer3_long_channel_effects_long45_allQC.png`

### 可复现 MATLAB 脚本

- `code/matlab/run_q8_procedure_cohort_specification.m`
- `code/matlab/run_q8_procedure_cohort_estimation.m`
- `code/matlab/run_q8_procedure_cohort_compute_contrast.m`
- `code/matlab/summarize_q8_procedure_group_contrasts.m`
- `code/matlab/summarize_q8_procedure_roi_long45_allqc.m`
- `code/matlab/summarize_q8_procedure_roi_cortical_system_long45_allqc.m`
- `code/matlab/plot_q8_procedure_design_and_fit.m`
- `code/matlab/plot_q8_procedure_roi_long45_allqc.m`
- `code/matlab/run_q8_procedure_subject001_shortmean_ar1_sensitivity.m`

## Git 版本控制范围

版本库包含本次 Q8 的 MATLAB 代码、CSV 结果和 PNG 图像。逐受试者 `.mat` 中间文件为约 805 MB 的可再生成产物，且由 `.gitignore` 排除，因此不会提交到 GitHub。
