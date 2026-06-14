---
name: r-lint
description: Check R scripts in PSS-Net for project coding-style compliance — base R only (no pipes / dplyr / tidyr), snake_case naming, script-header docstring, and the project plotting stack (ggplot2 > patchwork > igraph+grid). Use when the user asks to lint, review style, or check an R script before committing.
---

# r-lint — PSS-Net R 代码风格检查

对指定的 R 脚本（默认 `sim_script/` 与 `analysis_script/` 下全部 `.R`，或用户指定文件）做静态风格检查，**只查风格、不查算法正确性**（后者用 /code-review）。

## 检查项（按严重度）

### 🔴 必须修正（违反硬性规范）
1. **禁止管道符**：源码中不得出现 `%>%`、`|>`。
2. **禁止 tidyverse 数据操作**：不得 `library(dplyr)` / `library(tidyr)`，也不得调用其函数（`mutate` `filter` `select` `group_by` `summarise` `%in%`除外 `pivot_*` `gather` `spread` 等）。改用 base R（`subset` `transform` `[` `aggregate` `apply` 族 `reshape`）。
3. **命名 snake_case**：函数名、变量名用小写下划线；禁止 camelCase / dot.case（`T`/`F` 也应写 `TRUE`/`FALSE`）。
4. **脚本头清空环境**：每个 R 脚本第一行可执行代码必须是 `rm(list = ls())`，自动清空环境变量。
5. **脚本头注释**：紧随其后须有用途 + 输入 + 输出说明（见现有 `pss_net_v0.R` 的 banner 注释风格）。

### 🟡 建议修正（绘图栈优先级）
6. **绘图优先 ggplot2**：统计图一律 `ggplot2`，不用 base `plot()`/`hist()`/`barplot()` 出最终图；Python 端只输出数据交 R 画。
7. **拼图优先 patchwork**：多子图组合用 `patchwork`（`p1 + p2`、`/`、`plot_layout`），不用 `gridExtra::grid.arrange` / `cowplot`。
8. **网络图优先 igraph + grid**：网络可视化用 `igraph`；若需与其它图拼合，用 `grid`/`gridGraphics`（把 igraph base 绘图转 grob 再拼），不混用 patchwork 直接拼 base 网络图。
9. **网络配色规范**：促进边 `tomato3`、抑制边 `steelblue3`、FP 误报 `orange2`(虚线)、FN 遗漏 `grey60`(点虚线)。出现网络图但未用该配色 → 提示。

### 🟢 提示
10. 图形输出存为 PDF/PNG，文件名与脚本对应。
11. 模拟/推断脚本应在 `sim_script/`，分析脚本在 `analysis_script/`，结果写入 `results/`。

## 执行方式
- 用 Grep 扫描上述模式（如 `%>%|\|>`、`library\((dplyr|tidyr)\)`、`grid\.arrange|cowplot`、`\bplot\(|\bhist\(`、camelCase 标识符）；并检查每个脚本首行可执行代码是否为 `rm(list = ls())`，缺失则报 🔴。
- 逐文件汇报：按 🔴/🟡/🟢 分组，每条给 `file:line` + 问题 + 建议改法。
- 末尾给一句总评（是否可提交）。除非用户要求 `--fix`，否则只报告不改代码。
