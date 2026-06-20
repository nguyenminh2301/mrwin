# mrwin: 계층적 복합 평가변수를 위한 인과 승리 통계량

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## 이 패키지가 하는 일

`mrwin`은 역학에서 특정 질문에 답하는 R 패키지입니다:

**유전적으로 예측된 노출이, 임상적으로 중증도에 따라 순위가 매겨진 결과를 기준으로
사람들을 더 나은 전반적 결과 프로파일로 이동시키는가, 아니면 더 나쁜 결과 프로파일로
이동시키는가?**

많은 질병에서 환자들은 단일 사건을 경험하지 않습니다. 심신질환(cardiorenal disease)을
가진 사람은 사망하거나, 심부전으로 입원하거나, 신장 기능 저하를 겪을 수 있습니다.
이러한 사건들은 서로 대체할 수 없습니다: 사망은 입원보다 더 나쁘고, 입원은 검사실
이상 소견보다 더 나쁩니다. 표준 분석은 모든 사건을 동등하게 취급하거나 한 번에
하나씩 분석하여, 환자와 임상의에게 가장 중요한 임상적 순서를 상실합니다.

`mrwin`은 두 가지 통계적 프레임워크를 결합하여 이 문제를 해결합니다:

1. **승리 통계량(Win statistics)** -- 임상적 중증도 계층(사망 > 입원 > 생체지표
   저하)을 존중하는 쌍별 비교.
2. **멘델식 무작위화(Mendelian randomization)** -- 미측정 교란요인이 있어도
   인과적 주장을 뒷받침하는 유전적 도구 변수.

그 결과는 **용량-표준화 인과 승리비(DS-CWR)** 입니다: 유전적으로 예측된 노출이
우선순위화된 결과 프로파일을 더 좋게 만드는지 나쁘게 만드는지를 요약하는 단일
숫자입니다.

---

## 멘델식 무작위화: 간략한 소개

### 문제: 미측정 교란

관찰 역학에서 우리는 종종 어떤 노출(예: LDL 콜레스테롤)이 어떤 결과(예: 심부전)를
유발하는지 알고 싶어 합니다. 어려운 점은 높은 LDL을 가진 사람들이 낮은 LDL을 가진
사람들과 식이, 운동, 사회경제적 지위, 다른 약물 등 여러 면에서 다르다는 것입니다.
이것들이 교란요인(confounders)입니다. 측정된 교란요인을 보정한 후에도 미측정
교란요인이 남아 있을 수 있습니다.

### 해결책: 자연 실험으로서의 유전학

수정 시, 각 사람은 각 유전자좌에서 부모로부터 하나의 대립유전자를 무작위로
물려받습니다. 이를 **멘델식 분리(Mendelian segregation)** 라고 합니다.
만약 유전 변이가 노출에 영향을 미친다면(예: LDLR 유전자의 변이가 LDL 콜레스테롤을
높인다면), 그 변이를 보유한 사람들은 평균적으로 평생 동안 더 높은 LDL에 노출됩니다
— 식이나 생활 방식 때문이 아니라 유전자형 때문입니다.

이러한 무작위 배정이 **멘델식 무작위화(MR)** 의 기초입니다. 유전 변이를
**도구 변수(instrumental variables)** 로 사용함으로써, MR은 미측정 교란을
우회하여 노출이 결과에 미치는 인과 효과를 추정합니다.

### MR의 주요 가정

유전적 도구가 유효하려면 세 가지 조건이 충족되어야 합니다:

1. **관련성(Relevance)**: 유전 변이가 노출과 연관되어야 합니다.
2. **독립성(Independence)**: 유전 변이가 교란요인과 연관되지 않아야 합니다
   (잘 설계된 연구에서 멘델식 분리에 따라 충족됨).
3. **배제 제한(Exclusion restriction)**: 유전 변이가 노출을 통해서만 결과에
   영향을 미쳐야 하며, 다른 경로를 통해서는 영향을 미치지 않아야 합니다.

### MR의 데이터 유형

| 데이터 유형 | 설명 | 예시 |
|---|---|---|
| 개인 수준 | 유전자형, 노출, 결과가 있는 1인 1행 데이터 | SNP 데이터, LDL 측정값, 병원 기록이 있는 10,000명 코호트 |
| 요약 수준 | GWAS의 SNP별 효과 추정치 | "SNP rs12345는 LDL에 대해 효과 0.05 (SE 0.01)를 가짐" |
| 다유전자 위험 점수(PRS) | 여러 SNP의 가중 합 | 100개 SNP로 계산된 LDL 콜레스테롤에 대한 PRS |

`mrwin`은 **PRS 도구**를 사용한 **개인 수준 데이터**를 사용합니다. PRS는
노출에 대한 누적 유전적 소인을 포착하는 연속형 도구 역할을 합니다.

---

## 승리 통계량: 목표 추정량

### 승리 통계량이란?

승리 통계량은 연구의 모든 개인 쌍을 비교합니다. 각 쌍 (A, B)에 대해 비교는
임상적 계층을 따릅니다:

1. 먼저, 가장 중요한 결과(예: 사망)에 대해 비교합니다. A가 사망하고 B가
   사망하지 않았다면, B가 승리합니다.
2. 첫 번째 결과에서 동점이면 다음 결과(예: 입원)로 이동합니다.
3. 한 사람이 승리하거나 모든 결과가 동점이 될 때까지 계층을 따라 내려갑니다.

**승리비(win ratio)** 는 더 높은 노출군의 승리 횟수를 더 낮은 노출군의
승리 횟수로 나눈 값입니다.

### 승리 통계량이 중요한 이유

| 방법 | 사건을 동등하게 취급? | 임상적 순서 존중? | 경쟁 위험 처리? |
|---|---|---|---|
| 첫 사건까지의 시간 | 예 | 아니오 | 아니오 |
| 복합 평가변수 분석 | 예 | 아니오 | 부분적 |
| 승리비 / 승리 통계량 | 아니오 | 예 | 예 |

승리 통계량은 Pocock 등(2012)에 의해 도입되었으며, 2022년에서 2024년 사이
36건 이상의 무작위 임상시험에서 채택되었습니다.

### 인과 승리비(cCWR)

표준 승리 통계량은 무작위 할당 또는 강한 무시 가능성 가정을 필요로 합니다.
관찰 역학에서는 이 가정이 방어되기 어렵습니다.

`mrwin`은 **연속형 인과 승리비(cCWR)** 를 유전적 도구의 분위수에 걸친
기울기 추정량으로 정의합니다:

- 모집단을 PRS로 정렬된 계층(예: 십분위)으로 나눕니다.
- 각 인접 계층 쌍 내에서 승리비를 계산합니다.
- 도구-표준화 기울기와 GLS 메타분석을 사용하여 계층 전체를 통합합니다.

결과는 강한 무시 가능성 가정을 필요로 하지 않는 **주변, 모집단 수준의 인과 효과**입니다.

---

## 핵심 공식: 평이한 언어로 설명

### 1. 계층적 비교 커널

두 개인 *i*와 *j*에 대해, 가장 높은 우선순위에서 가장 낮은 우선순위까지 결과를
비교합니다:

- 우선순위 *k*에서 *j*는 사건이 있었지만 *i*는 없었거나(또는 *i*가 더 오래 생존):
  *i*가 승리 (+1).
- 우선순위 *k*에서 *i*는 사건이 있었지만 *j*는 없었거나(또는 *j*가 더 오래 생존):
  *i*가 패배 (-1).
- 동점: 다음 우선순위로 이동.

**쉽게 말하면**: "임상적 중증도 순서를 존중했을 때, 이 사람이 상대방보다 더 나은
결과 프로파일을 가졌는가?"

### 2. 계층별 승리비

한 계층 내에서(예: 6번째 PRS 십분위 vs. 5번째):

```
theta_d = (더 높은 계층의 승리 횟수) / (패배 횟수)
log_theta_d = log(theta_d)
```

**쉽게 말하면**: "노출에 대한 유전적 소인이 약간 더 높은 사람들 중에서, 그들이
패배한 것보다 승리한 경향이 더 많았는가?"

### 3. 도구-표준화 기울기(ISG)

```
Delta_X_d = mean(계층 d의 노출) - mean(계층 d-1의 노출)
ISG_d = log_theta_d / Delta_X_d
```

**쉽게 말하면**: "유전학에 의해 예측된 대로, 노출의 단위당 증가에 따라 로그
승리비가 얼마나 변하는가?"

### 4. GLS 통합

모든 계층 쌍에 걸친 ISG 값은 축소 공분산 추정량을 사용한 일반화 최소제곱법으로
통합됩니다:

```
delta_GLS = 전체 계층에 걸친 통합 ISG
DS-CWR = exp(delta_GLS)
```

**쉽게 말하면**: "우선순위화된 결과 프로파일에 대한 노출의 전반적 인과 효과는
무엇인가?"

- **DS-CWR > 1**: 노출이 더 나은 결과 프로파일과 연관됨.
- **DS-CWR < 1**: 노출이 더 나쁜 결과 프로파일과 연관됨.
- **DS-CWR = 1**: 인과 효과의 증거 없음.

### 5. 승수 부트스트랩

불확실성은 다음과 같이 추정됩니다:
1. GWAS 가중치 섭동(외부 불확실성).
2. 무작위 승수 가중치 추출(내부 불확실성).
3. 각 부트스트랩 반복에 대해 전체 파이프라인 재계산.

이로써 두 가지 불확실성 원천을 모두 고려한 신뢰구간을 얻습니다.

### 6. 단계적 다면발현 진단(SDPD)

SDPD는 유전적 도구가 노출 이외의 경로를 통해 최우선순위 결과에 영향을 미치는지
(직접 다면발현) 검정합니다. 이는 SNP별 노출 및 결과 효과에 대해 MR-Egger
회귀를 적합하여 수행됩니다.

**쉽게 말하면**: "유전적 도구가 유효한가, 아니면 다른 경로를 통해 결과에
영향을 미치는가?"

---

## 임상 예시

### 예시 1: 심혈관 질환

**연구 질문**: 유전적으로 예측된 LDL 콜레스테롤이 전반적 심신질환 궤적을
악화시키는가?

**평가변수 계층** (가장 높은 우선순위에서 가장 낮은 순서로):
1. 사망
2. 심부전 입원
3. 신장 기능 저하

**해석**: DS-CWR = 0.82 (95% CI: 0.70 ~ 0.96)인 경우, 이는 유전적으로
예측된 LDL의 단위당 증가에 따라 우선순위화된 심신질환 결과 프로파일이 18%
악화됨을 의미합니다. CI가 1을 교차하지 않으므로 통계적으로 유의한 유해 효과를
시사합니다.

**임상적 의미**: 높은 LDL은 단순히 심장마비를 일으키는 것이 아니라, 전체 질병
궤적을 더 이른 사망, 더 많은 입원, 더 빠른 신장 기능 저하로 이동시키며, 이는
해당 임상적 순서로 나타납니다.

### 예시 2: 치매

**연구 질문**: 유전적으로 예측된 수면 장애가 사람들을 더 이른 중증 신경인지
결과로 이동시키는가?

**평가변수 계층**:
1. 사망
2. 치매 진단
3. 요양원 입소

**해석**: DS-CWR = 0.91 (95% CI: 0.78 ~ 1.06)인 경우, 점 추정치는
유해함을 시사하지만 CI가 1을 교차합니다. 인과 효과에 대한 강력한 증거는 없습니다.

**이것이 중요한 이유**: 수면이 치매에 미치는 영향에 대한 표준 MR에서, 치매가
발생하기 전에 사망한 사람들은 분석에서 손실됩니다. 이는 생존자 편향의 한
형태입니다: 유해한 수면 변이가 보호적으로 보일 수 있는 이유는 그 보유자들이
진단받을 만큼 충분히 오래 살지 못하기 때문입니다. `mrwin`은 사망을 최우선순위
사건으로 취급하여 이를 회피하므로, 차별적 생존이 조건화되어 제거되는 대신
추정량에 흡수됩니다.

### 예시 3: 암 역학

**연구 질문**: 유전적으로 예측된 BMI가 우선순위화된 암 결과를 악화시키는가?

**평가변수 계층**:
1. 암 사망
2. 진행 또는 전이
3. 재발

**해석**: DS-CWR = 0.75 (95% CI: 0.60 ~ 0.94)인 경우, 이는 유전적으로
예측된 BMI가 더 높을수록 암 궤적이 계층의 모든 수준에서 더 나쁜 결과로
이동함을 의미합니다.

**임상적 의미**: BMI는 단순히 암 위험을 증가시키는 것이 아니라, 재발부터 사망까지
전체 질병 경과를 악화시킵니다.

---

## 설치

```r
# GitHub에서 설치 (개발 버전)
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## 빠른 시작

### 시뮬레이션 데이터 (실제 코호트 불필요)

```r
library(mrwin)

# 심신질환 데이터셋 시뮬레이션
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# 모델 적합
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# 결과 검토
print(fit)
summary(fit)
plot(fit)
```

### 출력 읽기

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

| 항목 | 의미 | 확인 사항 |
|---|---|---|
| `delta_GLS` | 로그 척도 통합 효과 | 음수 = 유해, 양수 = 유익 |
| `DS-CWR` | 지수화된 효과 (승리비) | < 1 = 유해, > 1 = 유익, 1 = 효과 없음 |
| `95% CI` | 부트스트랩 신뢰구간 | 1을 교차하는가? |
| `Q` | 이질성 통계량 | 작은 p-값 = 효과가 계층 간에 변동 |
| `Warnings` | 구조화된 주의사항 | 해석 전에 항상 확인 |

### 결과 해석

DS-CWR = 0.82 (95% CI: 0.70 ~ 0.96)인 경우:

- 노출은 18% 더 나쁜 우선순위화된 결과 프로파일과 연관됩니다.
- CI가 1을 교차하지 않으므로 효과는 통계적으로 유의합니다.
- Q p-값 확인: 유의하면 효과가 노출 범위에 걸쳐 변동할 수 있습니다.
- 경고 확인: `sdpd_rejected`이면 다면발현이 추정치를 무효화할 수 있습니다.

---

## 실제 코호트 템플릿

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

## 출력 형식

### Print

```r
print(fit)
```

주요 결과 표시: 표본 크기, SNP, 계층, 부트스트랩 유효성, DS-CWR, CI,
Q 통계량 및 경고.

### Summary

```r
summary(fit)
```

다음을 포함한 형식화된 보고서 생성:
- 주요 DS-CWR 추정치, 주요(Fieller) CI, p-값, 선택적 다면발현 경계 CI
- 인접 계층 기울기 (log-theta, Delta-X, ISG, CWR)
- 이질성 Q 통계량
- 라벨이 붙은 참조로서의 델타법 구간 (Fieller 구간이 주요 구간임; 아래 *신뢰구간* 참조)
- SDPD 진단 (활성화된 경우)
- 구조화된 경고

### Tidy (테이블 및 추가 분석용)

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

### Plots

```r
plot(fit, type = "isg")       # 인접 ISG 기울기 (기본값)
plot(fit, type = "forest")    # 인접 CWR의 Forest plot
plot(fit, type = "bootstrap") # log-theta의 부트스트랩 분포
```

---

## 공변량 보정

교란요인이 측정된 경우(연령, 성별, 주성분), 순서형 IPTW 보정을 사용합니다:

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

IPTW 보정은 PRS 계층 간 공변량 균형을 맞추어 정밀도를 향상시킵니다.
SDPD는 직접 다면발현을 검정합니다.

---

## 장점과 한계

### 장점

- **미측정 교란 하의 인과 추론**: 관찰된 노출만이 아닌 유전적 도구를 사용.
- **임상적 계층 존중**: 사망이 입원보다 우선시되고, 입원이 생체지표 저하보다
  우선시됨.
- **경쟁 위험을 자연스럽게 처리**: 차별적 생존이 조건화되어 제거되지 않고
  추정량에 흡수됨.
- **단일 요약 측도**: 하나의 DS-CWR이 여러 구성요소별 분석을 대체.
- **진단 프레임워크**: SDPD가 배제 제한을 검정; Q 통계량이 이질성을 검정;
  구조화된 경고가 문제를 표시.
- **재현 가능**: 부트스트랩은 시드로 제어됨; 모든 출력은 주어진 시드에 대해
  결정론적.

### 한계

- **개인 수준 데이터 필요**: 현재 패키지는 유전자형, 노출, 결과가 포함된
  1인 1행 데이터가 필요합니다. 요약 데이터 전용 MR은 지원되지 않습니다.
- **비-붕괴성**: 승리비는 붕괴 가능하지 않습니다. 주변 DS-CWR은 조건부 승리비와
  다릅니다. 이는 추정량의 속성이지 버그가 아닙니다.
- **선형 PRS 도구**: 패키지는 선형 다유전자 위험 점수를 가정합니다.
  비선형 유전자-유전자 상호작용은 모델링되지 않습니다.
- **대각 GWAS 공분산**: 현재 부트스트랩은 전체 LD 공분산 행렬이 아닌
  SNP별 표준오차를 사용합니다. 이는 SNP가 연관 불균형 상태에 있을 때
  불확실성을 과소평가할 수 있습니다.
- **다면발현 취약성**: cCWR은 계층-오염 다면발현에 민감합니다.
  gamma = 0.05 정도의 작은 사망률 수준 다면발현도 포함률을 12%까지
  붕괴시킬 수 있습니다.
- **표본 크기 요구사항**: 바이오뱅크 규모의 표본(N > 100,000)은
  신뢰할 수 있는 추론을 위한 엄격한 통계적 전제조건입니다.
- **계층 수 `D`**: 십분위 방식 추정량은 계층 수를 선택해야 합니다. 연속적이고
  대역폭으로 제어되는 기울기 추정량(십분위 추정량은 그 정확한 박스카 특수 사례)이
  구현되고 검증되었으나, 외부 검사 결과 이는 십분위보다 잡음이 많고 `D` 민감도를
  줄이지 못하는 것으로 나타났습니다(도구-표준화 기울기는 비율이며, 더 미세한
  평활화는 그 분모를 축소시킵니다). 권장되는 방법은 연속 재매개변수화가 아니라
  `D`에 걸친 민감도 분석을 동반한 이산 추정량입니다(`inst/spec/wp15-continuous-isg.md`
  참조).

### 성능 및 추론 백엔드(선택적)

기본 백엔드는 변경되지 않았으나, 2단계 작업(`inst/spec/acceleration-roadmap.md`)은
`mrwin_controls()`를 통해 검증된 선택적 대안을 추가했습니다:

- `backend = "fast"` — 준2차의, 컴파일된(Rcpp) win/loss 커널(`K`개 우선순위
  수준에 대해 `Theta(N log^{K-1} N)`)로, 밀집 백엔드와 비트 단위로 동일하며
  추정량을 바이오뱅크 규모에 도달시킵니다(예: `K = 3`, `N = 80,000`을 1초 미만).
- `inference = "analytic"` — 승수 부트스트랩을 재현하는 닫힌 형식의 영향함수 분산
  (GWAS 가중치 불확실성에 대한 정확한 몬테카를로 항 포함)으로, 표본추출 성분에
  대한 부트스트랩 루프를 제거합니다.
- `stratification = "doubly_ranked"` — Tian/Burgess의 이중 순위 계층으로, PRS 순위
  구간화에 대한 더 약한 가정의 대안입니다(기본값은 `"prs_rank"`로 유지됨).

#### 신뢰구간

**주요하게 보고되는 구간은 Fieller 구간입니다.** 외부 보정 시뮬레이션 결과, 원래의
델타법 비율 구간은 너무 넓은(명목 제1종 오류에 도달하지 못함) 반면, 동일한 공분산
위의 Fieller 구성은 올바르게 보정되는 것으로 나타났습니다; 델타법 구간은 라벨이
붙은 참조로만 유지됩니다. 약한 도구 하에서 Fieller 구간은 거짓으로 좁은 집합이
아니라 정직하게 무한으로 보고됩니다.

---

## 가정

`mrwin` 프레임워크는 다음 가정에 의존합니다:

1. **관련성(Relevance)**: PRS가 노출과 연관되어야 합니다. 약한 도구는
   불안정한 추정치를 생성합니다(`weak_instrument` 경고 확인).

2. **독립성(Independence)**: PRS가 교란요인과 연관되지 않아야 합니다.
   이는 동질적 모집단에서 멘델식 분리에 의해 충족됩니다. 모집단 계층화는
   이를 위반할 수 있습니다; 주성분을 공변량으로 사용하십시오.

3. **배제 제한(Exclusion restriction)**: PRS가 노출을 통해서만 결과에
   영향을 미쳐야 합니다. SDPD는 최우선순위 평가변수에 대해 이를 검정합니다.
   기각되면 추정치가 무효일 수 있습니다.

4. **단조성(Monotonicity)**: 도구가 모든 개인에 대해 동일한 방향으로
   노출을 이동시켜야 합니다. 위반은 추정치에 편향을 초래할 수 있습니다.

5. **생존에 의한 효과 수정 없음**: 인과 효과가 생존자와 비생존자 간에
   달라서는 안 됩니다. 이는 검정 불가능합니다.

6. **올바른 우선순위 순서**: 임상적 계층은 분석 전에 지정되어야 합니다.
   패키지는 데이터로부터 이를 추론하지 않습니다.

---

## 경고 코드

| 코드 | 의미 | 조치 사항 |
|---|---|---|
| `weak_instrument` | 표현형 이동이 0에 가깝거나 Fieller CI가 무한 | 추정치를 불안정한 것으로 간주 |
| `sdpd_rejected` | MR-Egger 절편이 다면발현 시사 | 민감도 논의 추가 |
| `sdpd_underpowered` | SDPD를 위한 SNP 수 부족 | 검정력 한계 보고 |
| `positivity_failure` | 계층 내 IPTW ESS가 너무 작음 | 공변량 중첩 검사 |
| `bridged_strata` | 양성 실패로 인해 계층 건너뜀 | 브리징 스키마 보고 |
| `discordant_components` | 구성요소 MR 방향이 충돌 | 전역 효과를 과장하지 말 것 |

---

## 시뮬레이션 엔진

방법 검증을 위해 시나리오 그리드를 실행합니다:

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

| 시나리오 | 의미 | 예상 동작 |
|---|---|---|
| A_null | 인과 효과 없음 | DS-CWR이 1 근처, 기각률이 5% 근처 |
| B_valid_IV | 유효한 인과 효과 | DS-CWR이 예상 방향 |
| C_pleiotropy | 직접 다면발현 | SDPD 기각률 증가 |
| D_hierarchy_discordant | 반대 방향 구성요소 효과 | 경고 발생 |

---

## 이 패키지 인용하기

연구에 `mrwin`을 사용하는 경우 다음을 인용해 주십시오:

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

**저자:**
- Nguyen Thien Minh, 호치민시 의약대학, 베트남 (minhnt@ump.edu.vn)
- N. Ahmad Aziz (교신저자), 독일 신경퇴행성질환 센터(DZNE) 및 본 대학교,
  독일 (Ahmad.Aziz@dzne.de)

---

## 관련 참고문헌

**승리 통계량:**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**멘델식 무작위화:**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**다유전자 위험 점수:**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**다중 상태 모형:**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## 패키지 구조

```
mrwin/
  R/
    api.R              # 상위 수준 mrwin() 워크플로우
    kernel.R           # 계층적 비교 커널
    estimate.R         # DS-CWR 추정량, GLS 통합
    bootstrap.R        # 승수 부트스트랩 추론
    adjustment.R       # 순서형-IPTW, ESS, 브리징
    sdpd.R             # 단계적 다면발현 진단
    simulate.R         # 시뮬레이션 데이터 생성 과정
    simulation_engine.R # 시나리오 그리드 실행기
    validation.R       # 입력 검증
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # 성능 벤치마킹
    backend_sparse.R   # 희소 백엔드
    config.R           # 시뮬레이션 구성
    strata.R           # PRS 계층 할당
  tests/
    testthat/          # 273개 단위 테스트
  vignettes/           # 3개 비네트
  inst/spec/           # 구현 명세
```

---

## 라이선스

MIT
