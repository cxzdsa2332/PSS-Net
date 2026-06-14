# manuscript/

PSS-Net 论文 LaTeX 工程。

## 编译

```bash
cd manuscript
make          # latexmk -pdf -bibtex main.tex  → main.pdf
make clean    # 清理中间文件
make purge    # 连同 main.pdf 一起删除
```

需要本地有 TeX 发行版（TeX Live / MacTeX）和 `latexmk`。

## 结构

```
manuscript/
├── main.tex              # 主文件，\input 各部分
├── preamble.tex          # 宏包、记号宏（与 CLAUDE.md 记号保持一致）
├── Makefile
├── sections/
│   ├── 00_abstract.tex
│   ├── 01_introduction.tex
│   ├── 02_method.tex     # 对应 CLAUDE.md §1–§7
│   ├── 03_simulation.tex # 模拟实验（results/mcc_comparison.csv）
│   ├── 04_realdata.tex   # 实际数据（data/datasets.md）
│   ├── 05_discussion.tex
│   └── A1_appendix.tex
└── figures/              # 图（建议从 R 脚本输出 PDF 到此）
```

参考文献复用根目录 `../ref/references.bib`（`\bibliography{../ref/references}`）。
正文中以 TODO 注释标出各部分待撰写内容。
