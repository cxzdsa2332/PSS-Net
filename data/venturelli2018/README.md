# Venturelli et al. 2018 — 12株合成人肠道菌群（下载数据）

**来源**: Venturelli OS et al. "Deciphering microbial interactions in synthetic
human gut microbiome communities." *Mol Syst Biol* 14:e8157 (2018).
DOI: https://doi.org/10.15252/msb.20178157

**下载日期**: 2026-06-26，自 EMBO Press / Springer 静态补充材料
(`static-content.springer.com/.../44320_2018_BFMSB178157_MOESM*`)。

## 文件清单（原始 MOESM 编号 → 本地名）

| 本地文件 | 原始 | 内容 |
|---|---|---|
| `Dataset_EV1.xlsx` | MOESM2 | 成对群落 PW1 的物种**相对丰度时间序列**（Fig 2A），sheet `PW1` + `Description` |
| `Dataset_EV2.xlsx` | MOESM3 | 成对群落 PW2 的相对丰度时间序列（Fig 2B） |
| `Dataset_EV3.xlsx` | MOESM4 | 多物种群落（含 11 株 leave-one-out 与全 12 株 "NONE"）相对丰度时间序列（Fig 3A） |
| `Dataset_EV4.xlsx` | MOESM5 | 单菌代谢物 log2 fold change（Fig 5A），97 代谢物 × 12 物种 |
| `Dataset_EV5.zip` | MOESM6 | `twelveSpeciesT1–T4.xml`（4 个时间点的群落组成 XML） |
| `Appendix.pdf` | MOESM1 | Appendix（含实验/建模细节，OD600 与 gLV 参数可能在此） |

## 12 物种缩写（EV3 中的排列顺序）

BH, CA, BU, PC, BO, BV, BT, EL, FP, CH, DP, ER

## 适配 PSS 的重要注意事项

1. **这些 EV 表是“图源数据”，不是完整组合实验矩阵。** 它们覆盖成对（Fig2）和单缺失/
   全群落（Fig3）等被选群落的时间序列，并非全部 ~200+ 组合条件。
2. **数据是相对丰度（compositional），不是 OD600 校正的绝对丰度。** 选 Venturelli 当主案例
   的主要理由之一（近绝对丰度、规避 closure loss）在这些 EV 表里**并不直接满足**——
   OD600 总量信息需到 `Appendix.pdf` 或作者 GitHub/原始数据核对，否则真实数据分析会落在
   组成型一档（见手稿 Fig3e 的 closure 退化）。
3. **是时间序列，不是单一稳态。** 用于 PSS 时，取每条件最后一个时间点近似稳态 $x^{*}$，
   并核对是否已收敛到平台。
4. VenturelliLab GitHub 无 2018 专用仓库（仅 `Clark_et_al_2021` 等）；如需完整组合
   终点 + OD600 数据，需进一步从 Appendix/作者处确认。

## 待办

- [ ] 解析 EV1–EV3，提取每条件末时间点作为稳态丰度，统一物种顺序。
- [ ] 在 Appendix.pdf 中定位 OD600 总量 / gLV 推断参数（用于绝对丰度换算与结果对照）。
- [ ] 确认是否需要从其他来源补全完整组合条件矩阵。
