# Extended Analyses

SIPEAC 扩展分析代码，与原有 Figure 脚本分开存放，不按投稿轮次命名。

在 RStudio 中将工作目录设为 SIPEAC 根目录（包含 Code 的那一层），新开 R 会话后打开对应脚本，按 `####` 分节从上到下运行。不要一次 source 全目录。

主脚本是实际分析代码，没有菜单、`main()` 或 step 包装。`R/` 只保留算法、统计和绘图中重复用到的函数。

| 文件 | 分析 |
|---|---|
| 01_Spatial_pairing_validation.R | 重复空间事件、随机配对与模型验证 |
| 02_RCTD_pairing.R | 单细胞参考、RCTD、singlet配对与比较 |
| 03_Bin1_Bin50_pairing.R | 固定GT的Bin1/Bin50配对比较 |
| 04_Clonotype_normalization.R | 匹配深度与重复rarefaction |
| 05_Repair_comparison.R | 从已有预测重建SIPEAC/repair比较，不是重跑repair |
| 06_Sequencing_QC.R | 原始输入预处理、覆盖率和基因饱和度 |
| 07_Spatial_BCR_GEX.R | BCR与GEX空间一致性 |
| 08_New_antibody_selection.R | 新增候选抗体筛选与合成序列提取 |
| 09_Bin50_pairing_examples.R | Bin50有/无AI的配对例图 |
| 10_Isotype_expansion.R | 跨组织亚型与克隆扩增 |
| 11_TCGA_antigen.R | TCGA TLS及肿瘤/癌旁比较 |
| 12_Antibody_TLS.R | 抗体TLS内外的bin比例 |
| 13_IFI30_distance.R | 抗体周围IFI30表达 |
| 14_Antigen_neighborhood.R | 抗原与抗体邻域 |
| 15_Antibody_phylogeny.R | 抗体家族与空间谱系树 |
| 16_ssDNA_prepare.R | h5ad转换与坐标导出 |
| 17_ssDNA_pairing.R | SIPEAC在ssDNA/scRNA参考下的比较 |
| 18_ssDNA_reference.R | ssDNA直接配对对照scRNA |
| 19_ssDNA_comparison.R | 汇总ssDNA与各分辨率结果 |

输入路径都在脚本开头。01–09使用 `data/raw/` 和 `data/metadata/`；10–19使用 `data/hzy/`。仓库只含代码，不含数据。

输出独立写到 `data/processed/Extended_Analyses/`、`results/figures/Extended_Analyses/`、`results/tables/Extended_Analyses/` 和 `results/stats_logs/Extended_Analyses/`，不覆盖旧Figure目录。

16–19是ssDNA分支。16的数据准备完成后，17按需要设置 `bin_size`（默认110，与原代码最后生效值相同）；如需三分辨率汇总，应分别运行20、50、110。18的参考默认bin20，保持原设置。19读取17/18的表。其他编号用于分类，不表示所有文件必须连续运行。

## 运行前注意

- 02的参考注释/UMAP输入需准备；原流程定义了SoupX函数但没有调用，本次没有补做。RCTD和06的原始测序预处理耗时较长，只运行所需分节。
- 10–19保留原稿的分析口径（10中NL/CLD扩增图样本已对调更正），尚未确认的科学问题见 [checks/notes.md](checks/notes.md)。
- 空间Python图需要 numpy、pandas、matplotlib、seaborn、opencv-python；在 `R/plot.R` 开头设置可用的 `PYTHON_BIN`。默认值为 `python3`。15另需PHYLIP的dnapars，16需已有hotspot conda环境。
- 已做语法检查和小型模拟数据测试，没有用真实数据重跑，结果复现情况未验证。12写xlsx需要 `writexl`，03的患者层面检验需要 `coin`。
