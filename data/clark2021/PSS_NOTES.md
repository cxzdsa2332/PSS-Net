# Clark et al. 2021 — PSS-Net 使用说明

**来源**: clone 自 https://github.com/VenturelliLab/Clark_et_al_2021 （2026-06-27，已删除
嵌套 `.git`）。论文：Clark RL et al. *Nat Commun* 12:3254 (2021),
DOI: https://doi.org/10.1038/s41467-021-22938-y。

**为何改用本数据**: Venturelli 2018 的 EMBO 可下载补充材料只含相对丰度图源子集 + 拟合 gLV
模型（SBML），没有可用的 OD600 绝对丰度表（详见 `../venturelli2018/README.md`）。Clark 数据
**同时给出 OD 总量与每物种丰度**，可直接构造近绝对丰度，且为完整组合设计（631 组合 × 重复）。

## 已精简

原 clone 434 MB，多为 MCMC 训练/仿真输出（`RLC* Model Training/`、`RLC* Simulations/`），
PSS-Net 用不到，已删除；需要时可从 GitHub 重新 clone。保留 3.2 MB 核心：

| 文件 | 内容 |
|---|---|
| `commonfiles/2020_02_28_MasterDF.csv` | **主数据表**，1850 行 |
| `commonfiles/metadata_2019_06_17.py` | 物种代码、系统发育、code→拉丁名 `namedict` |
| `2020_12_16_README.pdf` / `.docx` | 官方数据说明 |
| `Python3CondaSpec.txt` / `JuliaSpec.txt` | 原分析环境 |

## 主数据表结构（`2020_02_28_MasterDF.csv`，1850 行）

- `Treatment`：群落组合 ID（631 个唯一组合）；`Rep`：重复 0–5。
- `OD`：群落 600 nm 吸光度 = **总生物量代理**（min −0.013, max 3.98, mean 1.85）。
- `Acetate, Butyrate, Lactate, Succinate`：代谢物浓度。
- 26 个物种丰度列：`ER FP AC HB CC RI DP BH CA PC DL CG BF EL CH BO BT BU BV BC BY PJ DF BL BP BA`
  （`metadata` 中 `numspecies=25`，主表多一列；以表头为准）。
- `* Fraction` 列：对应的相对丰度。

## PSS 映射

- **稳态 $x^{*(k)}$**：每个 `Treatment`（可对重复取均值）的终点群落组成。
- **绝对丰度**（规避 closure，对应手稿 Fig3e 的 absolute 档）：`OD` × 每物种相对丰度。
- **扰动 $u$**：物种存在/缺失的接种组合（target-only 设计输入）。
- **物种拉丁名**：用 `metadata_2019_06_17.py` 的 `namedict`（如 `BT` =
  *Bacteroides thetaiotaomicron*），手稿中首次出现给全名并斜体。

## 待办

- [ ] 写 `sim_script/` 或 `analysis_script/` 入口：读取 MasterDF、按 Treatment 聚合到稳态、
      构造绝对丰度与设计矩阵。
- [ ] 与已发表 gLV 交互（论文 Fig）做 sign/support 对照。
