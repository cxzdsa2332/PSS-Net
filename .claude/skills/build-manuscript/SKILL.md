---
name: build-manuscript
description: Compile the PSS-Net LaTeX manuscript to manuscript/main.pdf using tectonic (the available compiler in this environment), running enough passes to resolve citations and cross-references. Use when the user asks to build, compile, recompile, or 编译 the manuscript / paper / PDF, and automatically after any agent modifies manuscript sources or the manuscript bibliography.
---

# build-manuscript — 编译 PSS-Net 手稿 PDF

把 `manuscript/` 下的 LaTeX 源编译为 `manuscript/main.pdf`。

## 自动触发规则

- 修改 `manuscript/main.tex`、`manuscript/preamble.tex`、
  `manuscript/sections/*.tex`、手稿专用 figure/table 源，或手稿使用的
  `ref/references.bib` 后，必须调用本技能重新编译 `manuscript/main.pdf`。
- 如果编译器或网络/缓存环境导致无法编译，汇报实际命令、错误摘要和阻断原因；
  不要静默跳过。

## 环境说明

- 本地**没有** `latexmk` / `pdflatex` / `xelatex`，Makefile 的 `latexmk` 路径不可用。
- 可用编译器是 **`tectonic`**（`/opt/homebrew/bin/tectonic`），它自动处理多遍编译、
  BibTeX 与交叉引用，并自动下载缺失宏包。
- 参考文献来自根目录 `../ref/references.bib`（main.tex 内 `\bibliography{../ref/references}`）。

## 执行步骤

1. 进入手稿目录并编译（一条命令即可，tectonic 内部已跑足多遍）：

   ```bash
   cd /Users/angdong/Documents/R_proj/PSS-Net/manuscript && tectonic main.tex 2>&1 | tail -30
   ```

2. **检查输出**：
   - 成功标志：`Writing \`main.pdf\``。
   - tectonic 常见提示 `errors were issued by BibTeX, but were ignored` 一般来自
     `.bib` 中个别条目，不阻断 PDF 生成；若用户关心引用，用
     `tectonic main.tex --print --keep-logs` 看 BibTeX 细节。
   - 若某个 `\input` 的 section 报真正的 LaTeX 错误，定位到 `sections/*.tex`
     对应行修复后重新编译。

3. **汇报**：说明 PDF 是否生成、文件大小，以及任何未解决的引用/编译警告。
   不要把 `main.pdf` 或 LaTeX 中间件加入 git（它们被 `.gitignore` 忽略）。

## 注意

- 不要修改正文内容来"绕过"编译错误；先弄清错误来源再改。
- 改动正文文字时遵循 `manuscript-style` skill 的学术写作规范。
