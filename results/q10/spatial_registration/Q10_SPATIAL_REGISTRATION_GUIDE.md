# Q10 空间配准：完整操作与审计说明

## 1. 目的

本流程将 Q10 的 GLM channel-level 结果重新建立为可由 NIRS-SPM 呈现的空间图，同时保留原始实验分配表不变。它解决的是**通道的空间定位和可视化**，不是重新估计 GLM beta，也不是用 atlas 重新定义统计 ROI。

本次输出适用于 Q10：`HILoadHE - HILoadLE`，即固定高内在负荷时，高外在负荷相对低外在负荷的差异。

## 2. 原始资料与各自角色

|资料|用途|是否修改|
|---|---|---|
|`GLM/channel_allocation.csv`|原始 channel/ROI allocation；保留为历史原件|否|
|`GLM/amended_snirf/Subj001_Sess101_run001.snirf`|SNIRF 内实际的 source-detector measurement list 和探头三维坐标|否|
|`GLM/first_level/group_outputs/NIRS_SPM/channel_registration/subject001_probe_geometry.mat`|从 Subject 001 SNIRF 提取的 probe geometry|否|
|`q10_validated_channel_pair_mapping.csv`|本次经核验的、可用于空间配准的派生 mapping|新建|
|`q10_validated_channel_mni_coordinates.csv`|每个长通道 source/detector 中点和用于 NIRS-SPM 的投影坐标|新建|

因此，`channel_allocation.csv` 仍然是原始记录；新的 mapping 是一个有明确版本和审计轨迹的派生文件。

## 3. 为什么使用 VC source-detector 编号

allocation 表中的原始 source/detector 命名与 SNIRF 记录的 detector index 并非总是同一编号体系。例如原始 `Det17` 在虚拟通道（VC）编号中对应 SNIRF `Det01`。空间定位必须匹配 SNIRF measurement list，故以 allocation 表的 `VCSamplingSite.name` 作为正常情况下的 source-detector pair 来源。

每个 pair 都在 Subject 001 SNIRF 的 61 个唯一 source-detector pairs 中验证存在；只有被确认存在的 pair 才会被送入 NIRS-SPM registration。

## 4. Channel 40 的核验与修正

原始 allocation 的 ch40 为 `Src11-Det22`。这在 Subject 001 SNIRF 中不存在；与此同时，其他 allocation entries 显示原始 detector 到 VC/SNIRF detector 的规则为 `Det17 -> Det01`、...、`Det22 -> Det06`。

对于 source 11，SNIRF 中实际存在的 relevant pair 包括 `Src11-Det06`，而不存在 `Src11-Det22`。因此，派生 mapping 中仅对 channel 40 使用：

```text
allocation 中的原始文字：Src11-Det22
用于空间配准的经核验 pair：Src11-Det06
审计规则：audited_correction_raw_Det22_to_VC_Det06
```

这不是改写原始 allocation，也不改变 ch40 的现有 ROI label（Midline Medial PFC）；它仅修正了将该通道连接到 SNIRF 探头几何的编号。

## 5. 配准具体如何完成

1. 读取 `channel_allocation.csv`，识别 45 个 long channels 和 16 个 short channels。
2. 对每个 long channel，从 `VCSamplingSite.name` 取得 source/detector pair；ch40 应用上面的单一审计修正。
3. 读取 SNIRF 的 measurement list，提取所有唯一 `[sourceIndex, detectorIndex]`，并逐一验证 45 个 long-channel pairs 均存在且不重复。
4. 从 SNIRF probe geometry 读取 source 和 detector 的 3D 坐标；每个通道的位置取其 source 和 detector 的中点：

   `channel midpoint = (source position + detector position) / 2`

5. 使用 NIRS-SPM 的 `projection_CS` 将每个探头中点投影为 NIRS-SPM 的 cortical/MNI rendering coordinates；`render_MNI_coordinates` 再用于创建 registration 显示。
6. 保存为 NIRS-SPM 可读取的 `preproc_info`（61 个 measurement entries，其中 45 个 long channels 标记为 rendered）。short channels 本次不做空间结果图，但仍保留在 registration 结构中，避免伪造或删除原始测量。

本轮核验结果为：SNIRF 共 61 个唯一 pairs；45/45 long-channel pairs 可验证存在；ch40 的经核验 pair 为 `Src11-Det06`。

## 6. 重建 registration 的命令

在 MATLAB 中（已将 NIRS-SPM 加入 path）运行：

```matlab
cd('C:\aliceproject\GLM\first_level\group_outputs\NIRS_SPM');
build_q10_validated_spatial_registration;
```

主脚本：`build_q10_validated_spatial_registration.m`。

它会重建 `channel_registration_q10_validated` 中的 CSV、MAT 和审计报告，并且会在任意 long pair 缺失、重复或 ch40 未按预期修正时停止报错。

## 7. 用新 registration 绘制 Q10

Q10 的 first-level model 已估计，regressors 是：

```text
[HILoadLE, HILoadHE, constant]
contrast = [-1, +1, 0]
```

也就是 `HILoadHE - HILoadLE`。以 Subject 001 作为技术验证可运行：

```matlab
cd('C:\aliceproject\GLM\first_level\felipe_small_models');
run_q10_subject001_validated_spatial_maps;
```

该脚本使用新的 `q10_validated_channel_registration.mat`，并输出 HbO/HbR 各一张 dorsal 和 frontal view。输出目录为：

`q10_hiloadhe_vs_hiloadle/spatial_maps_subject001_validated/`

## 8. 如何解读这些空间图

- 这四张图是 Subject 001 的**个体第一层** Q10 contrast 的插值空间图，设置为 uncorrected `p = .05` 的可视化配置。
- 它们用于核验 registration、展示通道在 NIRS-SPM 框架内的空间分布；不能替代 30 名被试的正式群组推断。
- NIRS-SPM 的图是将离散通道统计量投影/插值到皮层表面，因此连续色块不是高分辨率体素测量，也不表示每个表面点都被直接采样。
- HbO 和 HbR 应分别呈现、并结合效应方向、置信区间与多重比较校正后的群组结果解释；单幅个体图不应作为显著性结论。

## 9. atlas 在此流程中的位置

这里的空间配准不依赖 Jess 当时使用的 atlas。NIRSite/FOLD 导出的 MNI-space montage 和 SNIRF probe geometry 已足以让 NIRS-SPM 建图。

atlas 只在需要把已验证的通道 MNI 坐标进一步赋予更细的解剖标签（例如 BA、AAL 或 Harvard-Oxford labels）时使用。它不会改变 Q10 的 beta、contrast 或 NIRS-SPM registration；也不应在等 Jess 确认 atlas 时阻塞 Q10/Q5-Q8 的 GLM 和空间图工作。

## 10. 可复现性与后续 Q5-Q8

Q5-Q8 可以复用**同一份** `q10_validated_channel_registration.mat` 和经核验 mapping。需要改变的只有对应模型的 `SPM.mat` 与 contrast 向量；不应重新手工编辑 channel allocation，也不应另行猜测 source/detector 编号。

每次画新研究问题的图前，应检查：

1. 模型的 regressor 名称和 contrast 与研究问题匹配；
2. 被调用的 registration 路径为 `channel_registration_q10_validated`；
3. 输出日志记录 subject、chromophore、view、contrast 和 threshold；
4. 群组统计与个体示意图分开保存和表述。
