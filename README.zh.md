# mrwin: 层次复合终点的因果赢统计

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## 本软件包的功能

`mrwin` 是一个 R 软件包，用于回答流行病学中的一个特定问题：

**由遗传预测的暴露是否会使人们在按严重程度进行临床排序的结局中，
整体结局轮廓变得更好或更差？**

在许多疾病中，患者并非仅经历单一事件。例如，心肾疾病患者可能死亡、因心力衰竭住院或出现肾功能下降。这些事件不可互换：死亡比住院更严重，住院又比实验室指标异常更严重。标准的分析方法将所有事件同等对待，或逐一单独分析，从而丢失了对患者和临床医生最为重要的临床顺序。

`mrwin` 整合了两个统计框架来解决这一问题：

1. **赢统计 (Win statistics)** —— 成对比较，尊重临床严重程度层次（死亡优先于住院，住院优先于生物标志物下降）。
2. **孟德尔随机化 (Mendelian randomization)** —— 遗传工具变量，即使在混杂因素未被测量的情况下，也能支持因果推断。

最终结果是一个**剂量标准化的因果赢比 (DS-CWR)**：一个单一数值，概括了由遗传预测的暴露是否使优先排序的结局轮廓变得更好或更差。

---

## 孟德尔随机化：简介

### 问题：未测量的混杂

在观察性流行病学中，我们经常想知道某种暴露（例如 LDL 胆固醇）是否会导致某种结局（例如心力衰竭）。困难在于，LDL 水平高的人群与 LDL 水平低的人群在许多其他方面也存在差异——饮食、运动、社会经济地位、其他用药情况。这些都是混杂因素。即便调整了已测量的混杂因素，未测量的混杂因素可能仍然存在。

### 解决方案：将遗传作为自然实验

在受孕时，每个人从父母双方各自随机继承一个等位基因。这被称为**孟德尔分离**。如果某个遗传变异影响某种暴露（例如，LDLR 基因中的某个变异会升高 LDL 胆固醇），那么携带该变异的人群，平均而言，一生中都暴露于较高的 LDL——不是因为他们的饮食或生活方式，而是因为他们的基因型。

这种随机分配构成了**孟德尔随机化 (MR)** 的基础。通过将遗传变异作为**工具变量**，MR 可以估计暴露对结局的因果效应，从而绕过未测量的混杂。

### MR 的关键假设

要使遗传工具有效，必须满足三个条件：

1. **相关性**：遗传变异必须与暴露相关。
2. **独立性**：遗传变异不得与混杂因素相关（在设计良好的研究中，这由孟德尔分离所保证）。
3. **排他性限制**：遗传变异必须仅通过暴露影响结局，而不能通过其他途径。

### MR 的数据类型

| 数据类型 | 含义 | 示例 |
|---|---|---|
| 个体水平数据 | 每人一行，包含基因型、暴露和结局 | 一个包含 10,000 人的队列，具有 SNP 数据、LDL 测量值和医院记录 |
| 汇总水平数据 | 来自 GWAS 的每个 SNP 效应估计值 | "SNP rs12345 对 LDL 的效应为 0.05（SE 0.01）" |
| 多基因风险评分 (PRS) | 多个 SNP 的加权和 | LDL 胆固醇的 PRS，由 100 个 SNP 计算得出 |

`mrwin` 使用**个体水平数据**和**PRS 工具**。PRS 作为一个连续的工具变量，捕捉了朝向该暴露的累积遗传倾向。

---

## 赢统计：目标估计量

### 什么是赢统计？

赢统计将研究中的每一对个体进行比较。对于每一对（A, B），比较遵循临床层次：

1. 首先，比较最重要的结局（例如，死亡）。如果 A 死亡而 B 未死亡，则 B 获胜。
2. 如果在第一个结局上打成平手，则进入下一个结局（例如，住院）。
3. 沿层次继续向下，直到某人获胜，或所有结局均打成平手。

**赢比 (win ratio)** 是暴露较高组的获胜次数除以暴露较低组的获胜次数。

### 为什么赢统计很重要

| 方法 | 同等对待事件？ | 尊重临床顺序？ | 处理竞争风险？ |
|---|---|---|---|
| 至首次事件时间分析 | 是 | 否 | 否 |
| 复合终点分析 | 是 | 否 | 部分 |
| 赢比 / 赢统计 | 否 | 是 | 是 |

赢统计由 Pocock 等人 (2012) 提出，在 2022 年至 2024 年间已被 36 项以上的随机临床试验采用。

### 因果赢比 (cCWR)

标准的赢统计需要随机化或强可忽略性假设。在观察性流行病学中，这一假设很少能够成立。

`mrwin` 将**连续因果赢比 (cCWR)** 定义为沿遗传工具分位数的梯度估计量：

- 按 PRS 将人群划分为有序层（例如，十分位数）。
- 在每一对相邻层之间计算赢比。
- 使用工具标准化梯度和 GLS 荟萃分析汇聚各层结果。

最终结果是一个**边际的、人群水平的因果效应**，不要求强可忽略性。

---

## 核心公式：通俗语言解释

### 1. 层次比较核函数

对于两个个体 *i* 和 *j*，按从高到低的优先级比较结局：

- 如果 *j* 在优先级 *k* 处发生了事件，而 *i* 未发生（或 *i* 存活时间更长），则 *i* 获胜 (+1)。
- 如果 *i* 在优先级 *k* 处发生了事件，而 *j* 未发生（或 *j* 存活时间更长），则 *i* 失败 (-1)。
- 如果打成平手，则进入下一优先级。

**简而言之**："这个人的结局轮廓是否比另一个人更好，遵循了临床严重程度的顺序？"

### 2. 层特异赢比

在某一层内（例如，第 6 个 PRS 十分位数与第 5 个相比）：

```
theta_d = （较高层的获胜次数）/（失败次数）
log_theta_d = log(theta_d)
```

**简而言之**："在遗传倾向略高的人群中，他们是否比失败更频繁地获胜？"

### 3. 工具标准化梯度 (ISG)

```
Delta_X_d = mean(层 d 的暴露) - mean(层 d-1 的暴露)
ISG_d = log_theta_d / Delta_X_d
```

**简而言之**："由遗传学预测的，每单位暴露增加对应的对数赢比变化有多大？"

### 4. GLS 汇聚

使用带有收缩协方差估计量的广义最小二乘法，将各层对的 ISG 值进行汇聚：

```
delta_GLS = 所有层汇聚后的 ISG
DS-CWR = exp(delta_GLS)
```

**简而言之**："暴露对优先排序结局轮廓的整体因果效应是什么？"

- **DS-CWR > 1**：暴露与更好的结局轮廓相关。
- **DS-CWR < 1**：暴露与更差的结局轮廓相关。
- **DS-CWR = 1**：无因果效应证据。

### 5. 乘数自助法

不确定性通过以下方式估计：
1. 扰动 GWAS 权重（外部不确定性）。
2. 抽取随机乘数权重（内部不确定性）。
3. 每次自助迭代重新计算完整流程。

由此得出的置信区间考虑了两种不确定性来源。

### 6. 逐步下降多效性诊断 (SDPD)

SDPD 检验遗传工具是否通过暴露以外的途径影响最高优先级结局（直接多效性）。这通过对每个 SNP 的暴露效应和结局效应拟合 MR-Egger 回归来实现。

**简而言之**："遗传工具是否有效，还是会通过其他途径影响结局？"

---

## 临床示例

### 示例 1：心血管疾病

**研究问题**：由遗传预测的 LDL 胆固醇是否会使整体心肾疾病轨迹恶化？

**终点层次**（优先级从高到低）：
1. 死亡
2. 心力衰竭住院
3. 肾功能下降

**解读**：如果 DS-CWR = 0.82（95% CI: 0.70 至 0.96），这意味着，每单位遗传预测 LDL 升高，优先排序的心肾结局轮廓恶化 18%。CI 不跨越 1，表明具有统计学显著的有害效应。

**临床意义**：较高的 LDL 不仅仅导致心脏病发作——它会将整个疾病轨迹推向更早的死亡、更多的住院和更快的肾脏衰退，完全遵循该临床顺序。

### 示例 2：痴呆症

**研究问题**：由遗传预测的睡眠障碍是否会将人们推向更早、更严重的神经认知结局？

**终点层次**：
1. 死亡
2. 痴呆诊断
3. 入住疗养院

**解读**：如果 DS-CWR = 0.91（95% CI: 0.78 至 1.06），点估计提示存在有害效应，但 CI 跨越 1。因此，没有强有力的因果效应证据。

**为什么这很重要**：在睡眠对痴呆的标准 MR 分析中，在发展为痴呆之前死亡的人会从分析中丢失。这是一种幸存者偏差：一种有害的睡眠变异可能显示为保护性，因为其携带者从未活到被诊断的年龄。`mrwin` 通过将死亡作为最高优先级事件来避免这一问题，因此差异生存被纳入估计量，而非被条件化排除。

### 示例 3：癌症流行病学

**研究问题**：由遗传预测的 BMI 是否会使优先排序的癌症结局恶化？

**终点层次**：
1. 癌症死亡
2. 进展或转移
3. 复发

**解读**：如果 DS-CWR = 0.75（95% CI: 0.60 至 0.94），这意味着，由遗传预测的较高 BMI 会在层次各个水平上推动癌症轨迹朝向更差的结局。

**临床意义**：BMI 不仅仅增加癌症风险——它会从复发乃至死亡，恶化整个病程。

---

## 安装

```r
# 从 GitHub 安装（开发版本）
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## 快速入门

### 模拟数据（无需真实队列）

```r
library(mrwin)

# 模拟心肾疾病数据集
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# 拟合模型
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# 查看结果
print(fit)
summary(fit)
plot(fit)
```

### 阅读输出结果

```
mrwin fit
  N: 500
  SNPs: 20
  Priorities: 3
  Strata: 5
  Bootstrap: 100 (100 valid)
  Backend: dense
  Adjustment: none
  delta_GLS: -0.2000
  DS-CWR: 0.8187
  95% CI: 0.7000 to 0.9600
  Q: 1.2345 ( df = 2 , p = 0.5394 )
```

| 字段 | 含义 | 需要关注的内容 |
|---|---|---|
| `delta_GLS` | 对数尺度的汇聚效应 | 负值 = 有害，正值 = 有益 |
| `DS-CWR` | 指数化的效应（赢比） | < 1 = 有害，> 1 = 有益，1 = 无效应 |
| `95% CI` | 自助法置信区间 | 是否跨越 1？ |
| `Q` | 异质性统计量 | p 值较小 = 效应在各层之间变化 |
| `Warnings` | 结构化警告 | 解读前务必检查 |

### 解读结果

如果 DS-CWR = 0.82（95% CI: 0.70 至 0.96）：

- 暴露与优先排序结局轮廓恶化 18% 相关。
- CI 不跨越 1，因此效应具有统计学显著性。
- 检查 Q 值对应的 p 值：如果显著，效应可能在暴露范围内有所不同。
- 检查警告：如果出现 `sdpd_rejected`，多效性可能会使估计无效。

---

## 真实队列模板

```r
snp_cols <- paste0("rs", 1:20)

fit <- mrwin(
  data = dat,
  endpoint = mrwin_endpoint(
    time = c("t_death", "t_hf", "t_renal"),
    status = c("d_death", "d_hf", "d_renal"),
    priority = c("death", "heart_failure", "renal_decline")
  ),
  genotype = snp_cols,
  exposure = "ldl_cholesterol",
  gwas = mrwin_gwas(beta = beta_hat, se = se_beta, snp = snp_cols),
  covariates = c("age", "sex", "PC1", "PC2"),
  controls = mrwin_controls(
    n_strata = 10,
    bootstrap = 500,
    seed = 20260510,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both"
  )
)

summary(fit)
plot(fit, type = "forest")
mrwin_report(fit, file = "analysis_report.md", format = "markdown")
tidy(fit)
```

---

## 输出格式

### 打印

```r
print(fit)
```

显示主要结果：样本量、SNP 数量、层数、自助法有效性、DS-CWR、CI、Q 统计量和警告。

### 摘要

```r
summary(fit)
```

生成格式化的报告，包含：
- 主要 DS-CWR 估计值、CI、p 值，以及可选的多效性约束 CI
- 相邻层梯度（log-theta、Delta-X、ISG、CWR）
- 异质性 Q 统计量
- Fieller 敏感性 CI
- SDPD 诊断（启用时）
- 结构化警告

### Tidy 输出（用于表格和进一步分析）

```r
tidy(fit)
#   term estimate    delta se_delta statistic p_value  ci_low ci_high ...
# DS-CWR    0.82 -0.2000   0.075    -2.667  0.0077  0.7000  0.9600 ...
```

### 报告

```r
mrwin_report(fit, format = "markdown")
mrwin_report(fit, file = "report.md", format = "markdown")
```

### 图表

```r
plot(fit, type = "isg")       # 相邻 ISG 梯度（默认）
plot(fit, type = "forest")    # 相邻 CWR 森林图
plot(fit, type = "bootstrap") # log-theta 的自助法分布
```

---

## 协变量调整

当混杂因素已被测量（年龄、性别、主成分），可使用有序 IPTW 调整：

```r
fit_adj <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  covariates = cbind(u_proxy = scale(dat$U)),
  controls = mrwin_controls(
    n_strata = 5,
    bootstrap = 100,
    seed = 2,
    adjustment = "ordinal_iptw",
    sdpd_scale = "both",
    pleiotropy_bias_radius = 0.05
  )
)
```

IPTW 调整通过在 PRS 层之间平衡协变量来提高精度。SDPD 用于检验直接多效性。

---

## 优势与局限性

### 优势

- **在未测量混杂下的因果推断**：使用遗传工具，而不仅仅是观察到的暴露。
- **尊重临床层次**：死亡优先于住院，住院优先于生物标志物下降。
- **自然处理竞争风险**：差异生存被纳入估计量，而非被条件化排除。
- **单一汇总度量**：一个 DS-CWR 取代多个分项分析。
- **诊断框架**：SDPD 检验排他性限制；Q 统计量检验异质性；结构化警告标记问题。
- **可复现性**：自助法由种子控制；所有输出在给定种子的情况下具有确定性。

### 局限性

- **需要个体水平数据**：本软件包目前需要每人一行的数据，包含基因型、暴露和结局。不支持仅汇总数据的 MR。
- **不可折叠性**：赢比不可折叠。边际 DS-CWR 与条件赢比不同。这是估计量的属性，而非软件错误。
- **线性 PRS 工具**：本软件包假设线性多基因风险评分。不模拟非线性基因-基因交互作用。
- **对角 GWAS 协方差**：当前的自助法使用逐个 SNP 的标准误，而非完整的 LD 协方差矩阵。当 SNP 处于连锁不平衡时，这可能低估不确定性。
- **多效性脆弱性**：cCWR 对层次污染型多效性敏感。死亡率水平的多效性即使小至 gamma = 0.05 也可能使覆盖率降至 12%。
- **样本量要求**：生物样本库级样本（N > 100,000）是可靠推断的严格统计先决条件。
- **计算可扩展性（开发中）**：当前的 win/loss 成对扫描在每次自助法迭代中关于 N 是二次复杂度，因此生物样本库规模是统计上的先决条件，但尚非计算上的默认。一个次二次（接近 N log N）后端、一个消除任意分层数的连续梯度估计量，以及一个解析方差均在第二阶段路线图（`inst/spec/acceleration-roadmap.md`）中规定，且均以可选后端形式加入，不改变现有结果。

---

## 假设

`mrwin` 框架依赖于以下假设：

1. **相关性**：PRS 必须与暴露相关。弱工具会产生不稳定估计（检查 `weak_instrument` 警告）。

2. **独立性**：PRS 不得与混杂因素相关。在同质人群中，孟德尔分离满足这一条件。人群分层可能违反这一假设；请使用主成分作为协变量。

3. **排他性限制**：PRS 必须仅通过暴露影响结局。SDPD 针对最高优先级终点检验此假设。如果被拒绝，估计值可能无效。

4. **单调性**：工具必须对所有个体以相同方向改变暴露。违反此假设可能导致估计偏差。

5. **无生存效应修饰**：因果效应不应在存活者与未存活者之间存在差异。这不可检验。

6. **正确的优先级排序**：临床层次必须在分析前指定。本软件包不会从数据中推断层次。

---

## 警告代码

| 代码 | 含义 | 应对措施 |
|---|---|---|
| `weak_instrument` | 表型偏移接近零或 Fieller CI 无界 | 将估计值视为不稳定 |
| `sdpd_rejected` | MR-Egger 截距提示多效性 | 添加敏感性讨论 |
| `sdpd_underpowered` | SNP 数量不足无法进行 SDPD | 报告功效限制 |
| `positivity_failure` | 某层中 IPTW ESS 太小 | 检查协变量重叠 |
| `bridged_strata` | 由于正向性失败而跳过的层 | 报告桥接方案 |
| `discordant_components` | 各组成部分 MR 方向冲突 | 不要过度宣称整体效应 |

---

## 模拟引擎

运行场景网格以验证方法：

```r
scenarios <- mrwin_scenarios(
  mrwin_config(n_outcome = 200, m_snps = 10, seed = 10),
  scenarios = c("A_null", "B_valid_IV", "C_pleiotropy", "D_hierarchy_discordant")
)

grid <- mrwin_run_simulation_grid(
  scenarios = scenarios, n_iter = 2,
  controls = mrwin_controls(n_strata = 4, bootstrap = 20, seed = 11)
)

mrwin_simulation_summary(grid)
```

| 场景 | 含义 | 预期行为 |
|---|---|---|
| A_null | 无因果效应 | DS-CWR 接近 1，拒绝率约 5% |
| B_valid_IV | 有效的因果效应 | DS-CWR 在预期方向 |
| C_pleiotropy | 直接多效性 | SDPD 拒绝率增加 |
| D_hierarchy_discordant | 各组成部分效应方向相反 | 发出警告 |

---

## 引用本软件包

如果您在研究中使用了 `mrwin`，请引用：

> Nguyen Thien Minh, N. Ahmad Aziz. "Causal Win Statistics: Integrating
> Instrumental Variable Estimation with Hierarchical Composite Endpoints."
> *arXiv preprint*, 2026.

BibTeX：

```bibtex
@article{NguyenAziz2026,
  title={Causal Win Statistics: Integrating Instrumental Variable
         Estimation with Hierarchical Composite Endpoints},
  author={Nguyen Thien Minh and N. Ahmad Aziz},
  journal={arXiv preprint},
  year={2026}
}
```

**作者：**
- Nguyen Thien Minh，胡志明市医药大学，越南 (minhnt@ump.edu.vn)
- N. Ahmad Aziz（通讯作者），德国神经退行性疾病中心 (DZNE) 和波恩大学，德国 (Ahmad.Aziz@dzne.de)

---

## 相关参考文献

**赢统计：**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**孟德尔随机化：**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**多基因风险评分：**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**多状态模型：**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## 软件包结构

```
mrwin/
  R/
    api.R              # 高层 mrwin() 工作流
    kernel.R           # 层次比较核函数
    estimate.R         # DS-CWR 估计量，GLS 汇聚
    bootstrap.R        # 乘数自助法推断
    adjustment.R       # 有序 IPTW，ESS，桥接
    sdpd.R             # 逐步下降多效性诊断
    simulate.R         # 模拟数据生成过程
    simulation_engine.R # 场景网格运行器
    validation.R       # 输入验证
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # 性能基准测试
    backend_sparse.R   # 稀疏后端
    config.R           # 模拟配置
    strata.R           # PRS 层分配
  tests/
    testthat/          # 273 个单元测试
  vignettes/           # 3 个小品文
  inst/spec/           # 实现规范
```

---

## 许可证

MIT
