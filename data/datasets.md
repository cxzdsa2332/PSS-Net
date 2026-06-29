# PSS-Net 候选真实数据集

本文档汇总 PSS-Net 真实数据分析与方法比较的候选数据，分两类：

- **微生物组群落数据**：多扰动条件下系统达到稳态后测量物种丰度，结构上完全契合
  PSS 框架（每个组合培养至稳态对应 $x^{*(k)}$，接种/药物量作为扰动 $u$）。
- **基因调控 / MRA 比较数据**：用于 PSS-Net 与 Modular Response Analysis（MRA）的公平
  比较，要求有明确干预靶点、测量时点可视为稳态或近稳态、覆盖多个模块，并有网络真值
  或可审计的外部参考。三套分别承担**合成金标准、小型真实机制案例、中型真实系统扰动**。

---

## 统一分析原则

1. 经典 MRA 使用官方 [`aiMeRA`](https://github.com/bioinfo-ircm/aiMeRA/) 0.99.0；包论文为
   Jimenez-Dominguez et al. (2021)，DOI: `10.1038/s41598-021-86544-0`。
2. `aiMeRA` 输出的是对角归一化局部响应矩阵，不应当作绝对 Jacobian 或非线性边函数。
3. KO/siRNA 数据通常只给出干预靶点，不给出等价的连续加性输入强度。PSS-Net 不能未经说明
   地把所有干预写成 `u=-1`；真实数据应使用 target-only 输入并把干预幅度作为 nuisance scale，
   或采用与 MRA 相同的响应归一化。
4. 有完整网络真值时报告 MCC/AUPRC/符号准确率；真实数据没有完整真值时，主指标改为留一干预
   预测、重复稳定性和 pathway/STRING enrichment，不能把数据库未收录的边一律记作假阳性。

---

# 一、微生物组群落数据集（按 PSS 适配度排序）

## ⭐⭐⭐ 首选（结构完全契合 PSS）

### 1. Clark et al. 2021 — 26株合成人肠道菌群
- **论文**: Clark RL et al. "Design of synthetic human gut microbiome assembly and butyrate production." *Nature Communications* 12:3254 (2021)
- **DOI**: https://doi.org/10.1038/s41467-021-22938-y
- **数据（Zenodo）**: https://zenodo.org/records/4642238
- **物种数**: 26株人肠道厌氧菌（有代表性的4个门）
- **条件数**: 1850种不同菌种组合（单菌、双菌、3–17菌）
- **扰动类型**: 物种存在/缺失（接种组合）→ 终点稳态丰度
- **数据类型**: 16S rRNA 扩增子测序（Illumina），批次培养终点
- **PSS适配**: 每个菌种组合独立培养至稳态，直接对应 $x^{*(k)}$；物种接种量可作为连续扰动 $u$
- **注意**: 16S 为组成型数据，需 CLR 或参考菌归一化处理

### 2. Venturelli et al. 2018 — 12株合成人肠道菌群
- **论文**: Venturelli OS et al. "Deciphering microbial interactions in synthetic human gut microbiome communities." *Molecular Systems Biology* 14:e8157 (2018)
- **DOI**: https://doi.org/10.15252/msb.20178157
- **数据（论文补充）**: https://www.embopress.org/doi/abs/10.15252/msb.20178157 → Supplementary Dataset EV1–EV4
- **代码**: https://github.com/VenturelliLab
- **物种数**: 12株（*Bacteroides*, *Clostridium*, *Lactobacillus* 等代表性菌株）
- **条件数**: ~200+（单菌66对双菌+高阶组合）
- **扰动类型**: 物种添加/去除；OD600校正的绝对丰度（非纯组成型）
- **数据类型**: 16S rRNA + OD600 总量校正
- **PSS适配**: 最接近理想 PSS 实验；OD600校正使丰度近似绝对值，减少组成型偏差

---

## ⭐⭐ 次选（稳态可识别，需额外预处理）

### 3. Maier et al. 2024 — 32株菌 × 药物扰动
- **论文**: Maier L et al. "Emergence of community behaviors in the gut microbiota upon drug treatment." *Cell* 187:1–17 (2024)
- **DOI**: https://doi.org/10.1016/j.cell.2024.09.014
- **配套多组学**: https://pmc.ncbi.nlm.nih.gov/articles/PMC10495815/
- **物种数**: 32株人肠道代表菌
- **条件数**: 30种药物（抗生素+非抗生素）× 多浓度梯度
- **扰动类型**: 药物浓度（连续标量） → 终点群落组成
- **数据类型**: 16S + 宏基因组 + 代谢组，4天培养终点
- **PSS适配**: 药物浓度直接映射为 $u_j$，扰动设计最接近框架；需确认4天是否达稳态

### 4. Stein et al. 2013 — 小鼠肠道 + 抗生素/C. diff 扰动
- **论文**: Stein RR et al. "Ecological modeling from time-series inference: insight into dynamics and stability of intestinal microbiota." *PLoS Comput Biol* 9:e1003388 (2013)
- **DOI**: https://doi.org/10.1371/journal.pcbi.1003388
- **数据（SRA）**: https://www.ncbi.nlm.nih.gov/sra/SRA026269
- **物种数**: 11个OTU（属级，小鼠盲肠）
- **条件数**: 未干预、clindamycin处理、*C. difficile* 感染（~3个稳态平台）
- **数据类型**: 16S qPCR（**绝对丰度**，非组成型）
- **PSS适配**: gLV 推断的标准基准数据；绝对丰度优势显著；条件数少（N≈3），推断欠定

### 5. Bucci et al. 2016 (MDSINE) — 无菌小鼠肠道 + 饮食切换
- **论文**: Bucci V et al. "MDSINE: Microbial Dynamical Systems INference Engine for microbiome time-series analyses." *Genome Biology* 17:121 (2016)
- **DOI**: https://doi.org/10.1186/s13059-016-0980-6
- **数据+代码**: https://github.com/gerberlab/MDSINE2
- **物种数**: 13–16株（VE-202益生菌组合 或 GnotoComplex群落）
- **条件数**: 高纤/低纤饮食2个稳态 × 多只小鼠
- **数据类型**: 16S rRNA（Illumina）+ 总量qPCR校正
- **PSS适配**: 饮食成分量作为 $u$；MDSINE2 GitHub 含完整处理流程，便于复现

### 微生物组数据获取优先级

| 优先级 | 数据集 | 理由 |
|--------|--------|------|
| 1 | Clark 2021 (Zenodo) | 直接下载，N=1850最大，26物种适中 |
| 2 | Venturelli 2018 (Supp) | 12物种可解释性强，OD600校正接近绝对丰度 |
| 3 | MDSINE2 (GitHub) | 代码齐全，有现成预处理流程 |
| 4 | Stein 2013 (SRA) | gLV基准，便于与已发表结果横向比较 |
| 5 | Maier 2024 (PMC) | 药物扰动设计最优，但多组学处理复杂 |

---

# 二、基因调控 / MRA 比较数据集

## 1. DREAM4 In Silico Network Challenge

### 数据与来源

- 网络：DREAM4 challenge 4，network 1，`p=10` 与 `p=100` 两个规模；
- 数据：逐基因 knockout 后的稳态响应；Santra et al. (2013) 明确使用其 single-knockout 子集
  比较 Bayesian MRA、stochastic MRA 等方法；
- 生成工具：GeneNetWeaver；
- 数据/生成入口：[GeneNetWeaver DREAM challenge resources](https://gnw.sourceforge.net/dreamchallenge.html)；
- 主要引用：Santra et al. (2013), DOI `10.1186/1752-0509-7-57`；Schaffter et al.
  (2011), DOI `10.1093/bioinformatics/btr373`。

### 评价

| 维度 | 评价 |
|---|---|
| MRA/PSS 适配 | **高**：单靶点扰动、稳态响应、10/100节点均被 MRA 文献直接使用。 |
| 真值 | **完整**：有有向网络 gold standard，适合 MCC、AUPRC、符号与排名评价。 |
| 生物真实性 | **低—中**：是 GeneNetWeaver 合成系统，不能代替真实实验验证。 |
| 重复/噪声 | 原 single-KO 数据没有生物或技术重复；不适合直接估计实验方差。 |
| 推荐角色 | 首要定量 benchmark；先以10节点调试数据管线，再冻结100节点结果。 |

### 建议分析

- 只用 single-knockout steady-state 子集作为 MRA/PSS 的第一轮公平比较；
- time-series、multifactorial、knockdown 和 double-KO 数据可作为额外信息实验，不与只使用
  single-KO 的方法混在同一主比较中；
- 报告 edge AUPRC/MCC、符号准确率和 held-out knockout response prediction。

## 2. ERBB–G1/S：HCC1954 乳腺癌细胞

### 数据与来源

- 原始实验：Sahin et al. (2009)，DOI `10.1186/1752-0509-3-1`；
- 系统：trastuzumab-resistant HCC1954 细胞；
- 干预：15个 ERBB/G1-S 相关基因分别进行 RNAi knockdown；
- 测量：EGF 刺激12小时后，测量 ERBB1、ERBB2、p21、p27、CDK2、CDK4、Cyclin-D1
  表达以及 ERK、AKT、pRB 磷酸化；
- MRA 子集：10个测量模块中9个有直接对应 siRNA，pRB 未被直接扰动；
- MRA 重分析与补充数据：Santra et al. (2013)，DOI `10.1186/1752-0509-7-57`。
- 开放全文与补充文件：[原始实验 PMC2652436](https://pmc.ncbi.nlm.nih.gov/articles/PMC2652436/)、
  [MRA 重分析 PMC3726398](https://pmc.ncbi.nlm.nih.gov/articles/PMC3726398/)。

### 评价

| 维度 | 评价 |
|---|---|
| MRA/PSS 适配 | **高**：小型、定向干预、12小时响应、直接用于 MRA，并有生物和技术重复。 |
| 真值 | **中等**：有文献整理的 ERBB–G1/S reference pathway，但不是完整实验真值。 |
| 主要限制 | siRNA off-target、干预幅度未知、一个测量模块没有直接干预；隐藏节点会表现为有效间接边。 |
| 推荐角色 | 最适合做小型真实机制案例和管线单元测试。 |

### 建议分析

- 保留重复层级，不只分析均值；用分层 bootstrap 给出边和预测误差区间；
- PSS-Net 与 `aiMeRA` 使用相同的9个直接靶向干预；
- 主结果使用留一干预预测、边符号、reference pathway 支持率；对未测量中介造成的有效边单独标注。

## 3. K61：HAP1 kinase knockout screen

### 数据与来源

- 原始研究：Gapp et al. (2016)，DOI `10.15252/msb.20166890`；
- 原始测序数据：[ENA ERP012914](https://www.ebi.ac.uk/ena/browser/view/ERP012914)；
- 系统：人 HAP1 单倍体细胞，55个酪氨酸激酶与6个非激酶逐一 full knockout；
- 环境：无刺激以及 FGF1、ACTA、BMP2、IFNβ、IFNγ、WNT3A、ionomycin、resveratrol、
  rotenone、deferoxamine，共11种条件；刺激6小时后按稳态或近稳态处理；
- MRA 应用：Mekedem et al. (2022)，DOI `10.1371/journal.pcbi.1009312`。该研究把数据限制
  到61个被扰动基因，对每种环境形成一个61×61响应矩阵。

### 评价

| 维度 | 评价 |
|---|---|
| MRA/PSS 适配 | **很高**：系统性逐基因 KO、WT 对照、近稳态，共11个环境。 |
| 规模 | **合适**：61节点足以检验稀疏高维方法，又远小于 L1000，便于完整复现。 |
| 真值 | **有限**：没有完整有向真值；STRING/pathway 只能作为不完全外部参考。 |
| 主要限制 | full KO 不等于已知连续加性 `u`；浅层 RNA-seq、环境差异和间接转录效应均会影响解释。 |
| 推荐角色 | 主要真实数据；检验网络是否随刺激环境发生可重复变化。 |

### 建议分析

1. 先冻结无刺激条件的 PSS-Net/aiMeRA 对比；
2. 再将11个环境作为独立网络分析，避免把刺激条件直接当作重复；
3. 比较跨环境稳定边与刺激特异边，并对 JAK/FGFR 等原论文重点模块进行定向复核；
4. 主指标采用留一 KO 预测、重复稳定性、STRING/Reactome enrichment；数据库匹配只作为
   支持证据，不作为完整 gold standard MCC。

## 暂不列入前三的数据

- **L1000/CMap**：已被大型 MRA 使用，但约1,000维、shRNA off-target 明显、表达为 landmark/
  推断特征且缺完整真值。适合作为后续 scalability 补充，不适合作为第一批真实数据。
- **Sachs phospho-flow**：是常用因果发现 benchmark，但干预覆盖与稳态假设不如上述三套数据
  清楚，优先级低于 ERBB 和 K61。

### MRA 比较实施顺序

1. DREAM4 `p=10`：验证数据方向、边方向与 MRA 归一化；
2. ERBB/HCC1954：验证真实重复、target-only 干预和隐藏中介处理；
3. DREAM4 `p=100` 与 K61：冻结定量 benchmark 和主要真实数据结果。
