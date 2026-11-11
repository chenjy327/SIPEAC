# SIPEAC

## Code

- `Code/Figure1.R`、`Code/Figure E1.R`、`Code/Figure E2.R`、`Code/Figure2,E3.R`、`Code/Figure3.R`：在 RStudio 中把工作目录设为本仓库根目录后分节运行。脚本开头 source 的 `Code/figure_helpers.R` 提供共用函数和配色，空间区域外扩函数 `st_expand_n()` 在 `Code/Extended_Analyses/R/expand_n_layer.R`。
- 空间作图函数调用 `Code/Extended_Analyses/R/plot_cluster.py`，默认使用 `python3`，可用环境变量 `SIPEAC_PYTHON` 指定其他解释器。
- `Code/run.changeo.sh` 调用同目录下的 `add_info.R` 和 `rename_column.R`。
- 数据读写路径仍是原分析环境中的路径，运行前需按本地数据位置修改。
