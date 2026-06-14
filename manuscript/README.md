# manuscript/

PSS-Net LaTeX manuscript for perturbed steady-state inference in sparse nonlinear
dynamical systems.  Ecological and microbiome systems are treated as important
applications and simulation benchmarks, not as the only project scope.

## 编译

```bash
cd manuscript
make          # latexmk -pdf -bibtex main.tex  → main.pdf
make clean    # 清理中间文件
make purge    # 连同 main.pdf 一起删除
```

`make` requires a local TeX distribution and `latexmk`.  In the current local
environment, `tectonic` is available and can be used directly:

```bash
cd manuscript
tectonic main.tex
```

## 结构

```
manuscript/
├── main.tex              # 主文件，\input 各部分
├── preamble.tex          # 宏包、记号宏（与 CLAUDE.md 记号保持一致）
├── Makefile
├── sections/
│   ├── 00_abstract.tex
│   ├── 01_introduction.tex
│   ├── 02_method.tex     # General method for complex systems
│   ├── 03_simulation.tex # Simulation benchmarks, including gLV
│   ├── 04_realdata.tex   # 实际数据（data/datasets.md）
│   ├── 05_discussion.tex
│   └── A1_appendix.tex
└── figures/              # Manuscript-specific static figures if needed
```

参考文献复用根目录 `../ref/references.bib`（`\bibliography{../ref/references}`）。
Generated analysis figures are written to `../results/figure/` by
`analysis_script/` and are not tracked by Git while `results/` is ignored.
正文中以 TODO 注释标出各部分待撰写内容。
