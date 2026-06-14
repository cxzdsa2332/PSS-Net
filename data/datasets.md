# PSS-Net 候选真实数据集

按 PSS 框架适配度排序。核心要求：多扰动条件下系统达到稳态后测量物种丰度。

---

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

---

## 数据获取优先级

| 优先级 | 数据集 | 理由 |
|--------|--------|------|
| 1 | Clark 2021 (Zenodo) | 直接下载，N=1850最大，26物种适中 |
| 2 | Venturelli 2018 (Supp) | 12物种可解释性强，OD600校正接近绝对丰度 |
| 3 | MDSINE2 (GitHub) | 代码齐全，有现成预处理流程 |
| 4 | Stein 2013 (SRA) | gLV基准，便于与已发表结果横向比较 |
| 5 | Maier 2024 (PMC) | 药物扰动设计最优，但多组学处理复杂 |
