# mrwin: 階層的複合エンドポイントのための因果的Win統計

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## このパッケージが行うこと

`mrwin`は、疫学における特定の問いに答えるためのRパッケージです。

**遺伝的に予測される曝露は、アウトカムが重症度によって臨床的にランク付けされている場合、人々をより良いまたはより悪い全体的なアウトカムプロファイルへとシフトさせるのか？**

多くの疾患では、患者は単一のイベントを経験するわけではありません。心腎疾患を持つ人は、死亡するか、心不全で入院するか、腎機能低下を発症する可能性があります。これらのイベントは互換性がありません：死亡は入院よりも悪く、入院は検査異常よりも悪いのです。標準的な分析では、すべてのイベントを同等に扱うか、一度に1つずつ分析するため、患者と臨床医にとって最も重要な臨床的順序が失われます。

`mrwin`は、これに対処するために2つの統計的フレームワークを組み合わせています。

1. **Win統計** — 臨床的重症度の階層を尊重するペアワイズ比較（死亡は入院に勝り、入院はバイオマーカー低下に勝る）。
2. **メンデルランダム化** — 未測定の交絡因子が存在する場合でも因果的主張を支持する遺伝的操作変数。

結果は**用量標準化因果的Win比（DS-CWR）**です。これは、遺伝的に予測される曝露が優先順位付けされたアウトカムプロファイルをより良くするか悪くするかを要約する単一の数値です。

---

## メンデルランダム化：簡単な紹介

### 問題：未測定の交絡

観察疫学では、曝露（例：LDLコレステロール）がアウトカム（例：心不全）を引き起こすかどうかを知りたいことがよくあります。難しいのは、高いLDLを持つ人々は低いLDLを持つ人々と多くの他の点で異なるということです — 食事、運動、社会経済的地位、他の薬剤。これらは交絡因子です。測定された交絡因子を調整した後でも、未測定の交絡因子が残る可能性があります。

### 解決策：自然実験としての遺伝学

受胎時において、各人は各遺伝子座で両親から1つずつの対立遺伝子をランダムに受け継ぎます。これは**メンデル分離**と呼ばれます。もし遺伝的変異が曝露に影響を与えるならば（例：LDLR遺伝子の変異がLDLコレステロールを上昇させる）、その変異を持つ人々は平均して生涯を通じて高いLDLに曝露されます — 食事や生活習慣のためではなく、遺伝子型のためです。

このランダムな割り当てが**メンデルランダム化（MR）**の基礎です。遺伝的変異を**操作変数**として使用することにより、MRは未測定の交絡を回避して、曝露がアウトカムに与える因果効果を推定します。

### MRの主要な仮定

遺伝的操作変数が有効であるためには、3つの条件が満たされなければなりません。

1. **関連性**：遺伝的変異は曝露と関連していなければならない。
2. **独立性**：遺伝的変異は交絡因子と関連していてはならない（これは、適切に設計された研究におけるメンデル分離から導かれる）。
3. **除外制約**：遺伝的変異は曝露を通じてのみアウトカムに影響を与え、他の経路を通じては影響を与えてはならない。

### MRにおけるデータの種類

| データタイプ | 内容 | 例 |
|---|---|---|
| 個人レベル | 遺伝子型、曝露、アウトカムを持つ1人1行のデータ | SNPデータ、LDL測定値、入院記録を持つ10,000人のコホート |
| 要約レベル | GWASからのSNPごとの効果推定値 | 「SNP rs12345はLDLに対して効果0.05（SE 0.01）を持つ」 |
| ポリジェニックリスクスコア（PRS） | 多数のSNPの加重和 | 100個のSNPから計算されたLDLコレステロールのPRS |

`mrwin`は**PRS操作変数**を用いた**個人レベルデータ**を使用します。PRSは、曝露への累積的な遺伝的素因を捉える連続的な操作変数として機能します。

---

## Win統計：ターゲット推定対象

### Win統計とは何か？

Win統計は、研究内のすべての個人のペアを比較します。各ペア（A, B）について、比較は臨床的階層に従います。

1. まず、最も重要なアウトカム（例：死亡）で比較します。Aが死亡しBが死亡しなかった場合、Bが勝ちます。
2. 最初のアウトカムで同点の場合、次へ進みます（例：入院）。
3. 一人が勝つか、すべてのアウトカムが同点になるまで階層を下り続けます。

**Win比**は、高曝露群の勝利数を低曝露群の勝利数で割ったものです。

### Win統計が重要な理由

| 方法 | イベントを同等に扱うか？ | 臨床的順序を尊重するか？ | 競合リスクを処理するか？ |
|---|---|---|---|
| 初回イベントまでの時間 | はい | いいえ | いいえ |
| 複合エンドポイント分析 | はい | いいえ | 部分的に |
| Win比 / Win統計 | いいえ | はい | はい |

Win統計はPocockら（2012）によって導入され、2022年から2024年の間に36以上のランダム化臨床試験で採用されています。

### 因果的Win比（cCWR）

標準的なWin統計はランダム化または強い無視可能性を必要とします。観察疫学では、この仮定が防御可能であることは稀です。

`mrwin`は、**連続的因果的Win比（cCWR）**を遺伝的操作変数の分位点にわたる勾配推定対象として定義します。

- PRS（例：十分位数）によって集団を順序付けられた層に分割します。
- 隣接する各層ペア内でWin比を計算します。
- 操作変数標準化勾配とGLSメタ分析を用いて層全体を統合します。

結果は、強い無視可能性を必要としない**周辺的で集団レベルの因果効果**です。

---

## 中核公式：平易な言葉での説明

### 1. 階層的比較カーネル

2人の個人 *i* と *j* について、最も高い優先度から最も低い優先度へとアウトカムを比較します。

- 優先度 *k* で *j* がイベントを経験したが *i* は経験しなかった場合（または *i* がより長く生存した場合）：*i* の勝ち（+1）。
- 優先度 *k* で *i* がイベントを経験したが *j* は経験しなかった場合（または *j* がより長く生存した場合）：*i* の負け（-1）。
- 同点の場合：次の優先度へ進む。

**要するに**：「この人は他方の人よりも、臨床的重症度の順序を尊重した上で、より良いアウトカムプロファイルを持っていたか？」

### 2. 層別Win比

層内（例：第6PRS十分位数 vs. 第5十分位数）では：

```
theta_d = （高層の勝利数） / （敗北数）
log_theta_d = log(theta_d)
```

**要するに**：「曝露に対する遺伝的素因がわずかに高い人々の間で、彼らは負けるよりも勝つ傾向があったか？」

### 3. 操作変数標準化勾配（ISG）

```
Delta_X_d = 層dにおける曝露の平均 - 層d-1における曝露の平均
ISG_d = log_theta_d / Delta_X_d
```

**要するに**：「遺伝学によって予測される曝露の1単位増加あたり、対数Win比はどれだけ変化するか？」

### 4. GLSプーリング

すべての層ペアにわたるISG値は、収縮共分散推定量を用いた一般化最小二乗法によって統合されます。

```
delta_GLS = 全層にわたるプールされたISG
DS-CWR = exp(delta_GLS)
```

**要するに**：「優先順位付けされたアウトカムプロファイルに対する曝露の全体的な因果効果は何か？」

- **DS-CWR > 1**：曝露はより良いアウトカムプロファイルと関連している。
- **DS-CWR < 1**：曝露はより悪いアウトカムプロファイルと関連している。
- **DS-CWR = 1**：因果効果の証拠はない。

### 5. マルチプライヤーブートストラップ

不確実性は以下によって推定されます。
1. GWASの重みを摂動させる（外部不確実性）。
2. ランダムなマルチプライヤー重みを抽出する（内部不確実性）。
3. 各ブートストラップ反復についてパイプライン全体を再計算する。

これにより、両方の不確実性の源を考慮した信頼区間が得られます。

### 6. ステップダウン多面発現診断（SDPD）

SDPDは、遺伝的操作変数が曝露以外の経路（直接的多面発現）を通じて最も優先度の高いアウトカムに影響を与えるかどうかを検定します。これは、SNPごとの曝露およびアウトカム効果に対してMR-Egger回帰を適合させることによって行われます。

**要するに**：「遺伝的操作変数は有効か、それとも他の経路を通じてアウトカムに影響を与えているか？」

---

## 臨床例

### 例1：心血管疾患

**研究上の問い**：遺伝的に予測されるLDLコレステロールは、全体的な心腎の軌跡を悪化させるか？

**エンドポイント階層**（優先度の高い順）：
1. 死亡
2. 心不全入院
3. 腎機能低下

**解釈**：DS-CWR = 0.82（95％CI：0.70〜0.96）の場合、これは遺伝的に予測されるLDLの1単位増加あたり、優先順位付けされた心腎アウトカムプロファイルが18％悪化することを意味します。CIは1をまたがず、統計的に有意な有害効果を示唆しています。

**臨床的意味**：高いLDLは単に心臓発作を引き起こすだけではありません — それは、その臨床的順序において、より早期の死亡、より多くの入院、より速い腎機能低下へと疾患の軌跡全体をシフトさせます。

### 例2：認知症

**研究上の問い**：遺伝的に予測される睡眠障害は、人々をより早期の重篤な神経認知アウトカムへとシフトさせるか？

**エンドポイント階層**：
1. 死亡
2. 認知症診断
3. 介護施設入所

**解釈**：DS-CWR = 0.91（95％CI：0.78〜1.06）の場合、点推定値は有害性を示唆しますが、CIは1をまたぎます。因果効果の強力な証拠はありません。

**これが重要な理由**：睡眠の認知症に対する標準的なMRでは、認知症を発症する前に死亡した人々は分析から失われます。これは生存者バイアスの一形態です：有害な睡眠変異は、その保有者が診断を受けるほど長く生きられないため、保護的に見える可能性があります。`mrwin`は死亡を最も優先度の高いイベントとして扱うことでこれを回避するため、差別的生存は条件付けされて除去されるのではなく、推定対象に吸収されます。

### 例3：がん疫学

**研究上の問い**：遺伝的に予測されるBMIは、優先順位付けされたがんアウトカムを悪化させるか？

**エンドポイント階層**：
1. がん死
2. 進行または転移
3. 再発

**解釈**：DS-CWR = 0.75（95％CI：0.60〜0.94）の場合、これは遺伝的に予測されるBMIが高いほど、階層のすべてのレベルでがんの軌跡をより悪いアウトカムへとシフトさせることを意味します。

**臨床的意味**：BMIは単にがんリスクを増加させるだけではありません — 再発から死亡に至るまで、疾患経過全体を悪化させます。

---

## インストール

```r
# GitHubから（開発版）
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## クイックスタート

### シミュレーションデータ（実際のコホートは不要）

```r
library(mrwin)

# 心腎データセットのシミュレーション
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# モデルの適合
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# 結果の確認
print(fit)
summary(fit)
plot(fit)
```

### 出力の読み方

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

| フィールド | 意味 | 注目すべき点 |
|---|---|---|
| `delta_GLS` | 対数スケールでのプールされた効果 | 負 = 有害、正 = 有益 |
| `DS-CWR` | 指数化された効果（Win比） | < 1 = 有害、> 1 = 有益、1 = 効果なし |
| `95% CI` | ブートストラップ信頼区間 | 1をまたぐか？ |
| `Q` | 異質性統計量 | 小さいp値 = 効果が層間で異なる |
| `Warnings` | 構造化された注意事項 | 解釈前に必ず確認する |

### 結果の解釈

DS-CWR = 0.82（95％CI：0.70〜0.96）の場合：

- 曝露は18％悪い優先順位付けされたアウトカムプロファイルと関連している。
- CIは1をまたがないため、効果は統計的に有意である。
- Qのp値を確認する：有意な場合、効果は曝露範囲にわたって異なる可能性がある。
- 警告を確認する：`sdpd_rejected`の場合、多面発現が推定値を無効にする可能性がある。

---

## 実際のコホートテンプレート

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

## 出力形式

### Print

```r
print(fit)
```

主な結果を表示します：サンプルサイズ、SNP、層、ブートストラップの妥当性、DS-CWR、CI、Q統計量、および警告。

### Summary

```r
summary(fit)
```

以下を含む整形されたレポートを生成します。
- 主要なDS-CWR推定値、CI、p値、およびオプションの多面発現境界CI
- 隣接層勾配（log-theta、Delta-X、ISG、CWR）
- 異質性Q統計量
- Fieller感度CI
- SDPD診断（有効な場合）
- 構造化された警告

### Tidy（表およびさらなる分析用）

```r
tidy(fit)
#   term estimate    delta se_delta statistic p_value  ci_low ci_high ...
# DS-CWR    0.82 -0.2000   0.075    -2.667  0.0077  0.7000  0.9600 ...
```

### Report

```r
mrwin_report(fit, format = "markdown")
mrwin_report(fit, file = "report.md", format = "markdown")
```

### プロット

```r
plot(fit, type = "isg")       # 隣接ISG勾配（デフォルト）
plot(fit, type = "forest")    # 隣接CWRのフォレストプロット
plot(fit, type = "bootstrap") # log-thetaのブートストラップ分布
```

---

## 共変量調整

交絡因子が測定されている場合（年齢、性別、主成分）、順序IPTW調整を使用します。

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

IPTW調整は、PRS層間で共変量のバランスを取ることで精度を向上させます。SDPDは直接的多面発現を検定します。

---

## 長所と限界

### 長所

- **未測定交絡下での因果推論**：観察された曝露ではなく、遺伝的操作変数を使用する。
- **臨床的階層を尊重する**：死亡は入院よりも優先され、入院はバイオマーカー低下よりも優先される。
- **競合リスクを自然に処理する**：差別的生存は推定対象に吸収され、条件付けされて除去されない。
- **単一の要約指標**：1つのDS-CWRが複数のコンポーネント別分析を置き換える。
- **診断フレームワーク**：SDPDは除外制約を検定し、Q統計量は異質性を検定し、構造化された警告は問題を指摘する。
- **再現可能**：ブートストラップはシード制御され、すべての出力はシードが与えられれば決定論的である。

### 限界

- **個人レベルデータが必要**：このパッケージは現在、遺伝子型、曝露、アウトカムを含む参加者ごとに1行のデータを必要とします。要約データのみのMRはサポートされていません。
- **非崩壊性**：Win比は崩壊可能ではありません。周辺DS-CWRは条件付きWin比とは異なります。これは推定対象の特性であり、バグではありません。
- **線形PRS操作変数**：このパッケージは線形ポリジェニックリスクスコアを仮定します。非線形の遺伝子間相互作用はモデル化されません。
- **対角GWAS共分散**：現在のブートストラップはSNPごとの標準誤差を使用し、完全なLD共分散行列は使用しません。これは、SNPが連鎖不平衡にある場合に不確実性を過小評価する可能性があります。
- **多面発現脆弱性**：cCWRは階層汚染的多面発現に敏感です。gamma = 0.05程度の死亡率レベルの多面発現でも、カバレッジが12％に低下する可能性があります。
- **サンプルサイズ要件**：バイオバンク規模のサンプル（N > 100,000）は、信頼できる推論のための厳格な統計的前提条件です。
- **計算スケーラビリティ（開発中）**：現在のwin/lossペア走査はブートストラップの各反復でNに関して2次であるため、バイオバンク規模は統計的な前提条件ですが、計算上のデフォルトではまだありません。準2次（ほぼ N log N）バックエンド、任意の層数を排除する連続勾配推定量、および解析的分散が第II期ロードマップ（`inst/spec/acceleration-roadmap.md`）で規定されており、いずれも既存の結果を変えないオプトインのバックエンドとして追加されます。

---

## 仮定

`mrwin`フレームワークは以下の仮定に依存しています。

1. **関連性**：PRSは曝露と関連していなければならない。弱い操作変数は不安定な推定値を生み出します（`weak_instrument`警告を確認してください）。

2. **独立性**：PRSは交絡因子と関連していてはならない。これは均質な集団におけるメンデル分離によって満たされます。集団層別化はこれを侵害する可能性があります。主成分を共変量として使用してください。

3. **除外制約**：PRSは曝露を通じてのみアウトカムに影響を与えなければならない。SDPDは最も優先度の高いエンドポイントについてこれを検定します。棄却された場合、推定値は無効である可能性があります。

4. **単調性**：操作変数はすべての個人に対して曝露を同じ方向にシフトさせなければならない。違反は推定値をバイアスさせる可能性があります。

5. **生存による効果修飾がないこと**：因果効果は生存する人と生存しない人との間で異なるべきではない。これは検定不可能です。

6. **正しい優先順位付け**：臨床的階層は分析前に指定されなければならない。パッケージはデータからそれを推論しません。

---

## 警告コード

| コード | 意味 | 対処法 |
|---|---|---|
| `weak_instrument` | 表現型シフトがゼロに近い、またはFieller CIが非有界 | 推定値を不安定として扱う |
| `sdpd_rejected` | MR-Egger切片が多面発現を示唆 | 感度に関する議論を追加する |
| `sdpd_underpowered` | SDPDに対してSNPが少なすぎる | 検出力の限界を報告する |
| `positivity_failure` | 層内でIPTW ESSが小さすぎる | 共変量の重なりを検査する |
| `bridged_strata` | 正値性失敗により層がスキップされた | ブリッジされたスキーマを報告する |
| `discordant_components` | コンポーネントMRの方向が矛盾している | 全体的効果を過大主張しない |

---

## シミュレーションエンジン

シナリオグリッドを実行してメソッドを検証します。

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

| シナリオ | 意味 | 期待される動作 |
|---|---|---|
| A_null | 因果効果なし | DS-CWRが1付近、棄却率が5％付近 |
| B_valid_IV | 有効な因果効果 | DS-CWRが期待される方向 |
| C_pleiotropy | 直接的多面発現 | SDPD棄却が増加 |
| D_hierarchy_discordant | 相反するコンポーネント効果 | 警告が発出される |

---

## このパッケージの引用

研究で`mrwin`を使用する場合は、以下を引用してください。

> Nguyen Thien Minh, N. Ahmad Aziz. "Causal Win Statistics: Integrating
> Instrumental Variable Estimation with Hierarchical Composite Endpoints."
> *arXiv preprint*, 2026.

BibTeX:

```bibtex
@article{NguyenAziz2026,
  title={Causal Win Statistics: Integrating Instrumental Variable
         Estimation with Hierarchical Composite Endpoints},
  author={Nguyen Thien Minh and N. Ahmad Aziz},
  journal={arXiv preprint},
  year={2026}
}
```

**著者：**
- Nguyen Thien Minh, ホーチミン市医科薬科大学、ベトナム (minhnt@ump.edu.vn)
- N. Ahmad Aziz（責任著者）, ドイツ神経変性疾患センター（DZNE）およびボン大学、ドイツ (Ahmad.Aziz@dzne.de)

---

## 関連参考文献

**Win統計：**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**メンデルランダム化：**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**ポリジェニックリスクスコア：**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**マルチステートモデル：**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## パッケージ構造

```
mrwin/
  R/
    api.R              # 高レベルmrwin()ワークフロー
    kernel.R           # 階層的比較カーネル
    estimate.R         # DS-CWR推定量、GLSプーリング
    bootstrap.R        # マルチプライヤーブートストラップ推論
    adjustment.R       # 順序IPTW、ESS、ブリッジング
    sdpd.R             # ステップダウン多面発現診断
    simulate.R         # シミュレーションデータ生成プロセス
    simulation_engine.R # シナリオグリッドランナー
    validation.R       # 入力検証
    methods.R          # print、summary、plot、tidy、report
    benchmark.R        # パフォーマンスベンチマーキング
    backend_sparse.R   # スパースバックエンド
    config.R           # シミュレーション設定
    strata.R           # PRS層割り当て
  tests/
    testthat/          # 273のユニットテスト
  vignettes/           # 3つのビネット
  inst/spec/           # 実装仕様
```

---

## ライセンス

MIT
