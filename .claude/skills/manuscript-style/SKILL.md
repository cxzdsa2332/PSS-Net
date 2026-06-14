---
name: manuscript-style
description: Check and apply academic writing conventions when drafting or editing PSS-Net manuscript files (manuscript/**.tex, and prose in methods/*.md / README) — gene & species names in italics, define abbreviations on first use, consistent notation with preamble.tex, number/unit formatting, figure/table caption self-containment, and natbib \citep/\citet usage. Use whenever writing or revising paper text.
---

# manuscript-style — PSS-Net 学术写作规范检查

对 `manuscript/**.tex`（及 `methods/*.md`、`README` 中的论文性文字）做**写作/排版约定**检查与修订。
**只查写作规范，不查方法正确性**（后者另议）。逐条给 `file:line` + 问题 + 建议改法；
除非用户要求 `--fix`，否则只报告不改。

## 检查项（按严重度）

### 🔴 必须修正

1. **物种学名 / 基因名斜体**
   - 物种拉丁双名（属+种，如 *Bacteroides fragilis*、*Escherichia coli*）用斜体；
     属名首字母大写、种名小写。LaTeX 用 `\textit{...}`；Markdown 用 `*...*`。
   - 基因名斜体（如 *recA*）；其编码的**蛋白名正体**（RecA）。
   - 首次出现给全称后，后文可用 *E. coli* 形式（属名缩写仍斜体）。
2. **缩写首次出现定义**
   - 每个缩写在**正文首次出现**写"全称 (缩写)"，此后统一用缩写。
     例：`generalized Lotka--Volterra (gLV)`、`perturbed steady state (PSS)`、
     `Matthews correlation coefficient (MCC)`、`adaptive double sparse iterative
     hard thresholding (ADSIHT)`。摘要与正文各自独立定义一次。
   - 不在标题/小节标题里引入缩写定义。
3. **记号一致性**
   - 数学符号与 `manuscript/preamble.tex` 的宏、`methods/sindy_ss_method.md` 的记号
     保持一致（如 $\Psi_c$、$\theta_{ji}$、$f_{ji}$、$s, s_0$）；同一量勿用两种写法。

### 🟡 建议修正

4. **数字与单位**：句首不用阿拉伯数字；数字与单位间留空格（`5 mM`、`10^4`）；
   小数点前补 0（`0.91` 而非 `.91`）；统一有效数字位数。
5. **图表标题自洽**：每个 figure/table 的 caption 不依赖正文即可读懂（写明系统、$N$、
   重复数、误差棒含义）；正文须用 `\Cref{}`/`\ref{}` 引用每个图表。
6. **引用用法**：作者作句子成分用 `\citet{}`（"... \citet{x} showed"）；
   括号引用用 `\citep{}`。引用 key 必须在 `ref/references.bib` 中存在。
7. **连字符/破折号**：人名连用 en-dash（`Lotka--Volterra`、`Bickel--Ritov--Tsybakov`）；
   范围用 en-dash（`§3--5`、`10--20`）。

### 🟢 提示

8. 术语全文统一（如统一 "node/variable" 或 "species"，勿混用）。
9. 时态：方法/结果多用一般现在时或过去时，全文一致。
10. 美式或英式拼写择一统一。

## 执行方式

- 用 Grep 扫描常见违规：
  - 学名未斜体：搜物种名/属名裸文本（如 `Bacteroides`、`coli` 未在 `\textit`/`*` 内）。
  - 缩写未定义：搜 `gLV|PSS|MCC|ADSIHT|RIP|RE` 首次出现处是否有全称。
  - `\cite{`（应为 `\citep`/`\citet`）；`(^|\s)\.[0-9]`（缺前导 0）。
- 逐文件按 🔴/🟡/🟢 汇报，末尾给一句总评（是否可投/可合入）。
