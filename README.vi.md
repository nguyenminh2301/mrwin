# mrwin: Causal Win Statistics cho Tiêu chí Tổng hợp Phân cấp

[English](README.md) | [Tiếng Việt](README.vi.md) | [中文](README.zh.md) | [日本語](README.ja.md) | [한국어](README.ko.md) | [Français](README.fr.md) | [Deutsch](README.de.md)

## Chức năng của Gói

`mrwin` là một gói R trả lời một câu hỏi cụ thể trong dịch tễ học:

**Liệu một yếu tố phơi nhiễm được dự đoán bởi di truyền có làm cho tổng thể
kết cục của con người trở nên tốt hơn hay xấu đi, khi các kết cục được xếp
hạng lâm sàng theo mức độ nghiêm trọng?**

Trong nhiều bệnh, bệnh nhân không chỉ gặp một biến cố duy nhất. Một người
mắc bệnh tim-thận có thể tử vong, nhập viện vì suy tim, hoặc bị suy giảm
chức năng thận. Các biến cố này không thể hoán đổi cho nhau: tử vong nghiêm
trọng hơn nhập viện, và nhập viện nghiêm trọng hơn bất thường xét nghiệm.
Các phân tích tiêu chuẩn đối xử tất cả biến cố như nhau hoặc phân tích từng
biến cố một, làm mất đi thứ tự lâm sàng vốn quan trọng nhất đối với bệnh
nhân và bác sĩ lâm sàng.

`mrwin` kết hợp hai khung thống kê để giải quyết vấn đề này:

1. **Win statistics** -- so sánh từng cặp tôn trọng hệ thống phân cấp mức độ
   nghiêm trọng lâm sàng (tử vong > nhập viện > suy giảm dấu ấn sinh học).
2. **Mendelian randomization** -- biến công cụ di truyền hỗ trợ các tuyên bố
   nhân quả ngay cả khi các yếu tố gây nhiễu không được đo lường.

Kết quả là một **dose-standardized causal win ratio (DS-CWR)**: một con số duy
nhất tóm tắt liệu một yếu tố phơi nhiễm được dự đoán bởi di truyền có làm
cho hồ sơ kết cục được ưu tiên trở nên tốt hơn hay xấu đi.

---

## Mendelian Randomization: Giới thiệu Ngắn gọn

### Vấn đề: Nhiễu không Đo lường được

Trong dịch tễ học quan sát, chúng ta thường muốn biết liệu một yếu tố phơi
nhiễm (ví dụ: LDL cholesterol) có gây ra một kết cục (ví dụ: suy tim) hay
không. Khó khăn là những người có LDL cao cũng khác với những người có LDL
thấp ở nhiều khía cạnh khác -- chế độ ăn, tập thể dục, tình trạng kinh tế xã
hội, các loại thuốc khác. Đây là các yếu tố gây nhiễu. Ngay cả sau khi hiệu
chỉnh cho các yếu tố gây nhiễu đo lường được, các yếu tố gây nhiễu không đo
lường được vẫn có thể tồn tại.

### Giải pháp: Di truyền học như một Thí nghiệm Tự nhiên

Khi thụ thai, mỗi người nhận ngẫu nhiên một allele từ mỗi cha mẹ tại mỗi
locus di truyền. Điều này được gọi là **phân ly Mendel (Mendelian
segregation)**. Nếu một biến thể di truyền ảnh hưởng đến một yếu tố phơi
nhiễm (ví dụ: một biến thể trong gen LDLR làm tăng LDL cholesterol), thì
những người mang biến thể đó, trung bình, tiếp xúc với LDL cao hơn trong suốt
cuộc đời -- không phải vì chế độ ăn hay lối sống của họ, mà vì kiểu gen của
họ.

Sự phân công ngẫu nhiên này là cơ sở của **Mendelian randomization (MR)**.
Bằng cách sử dụng các biến thể di truyền như các **biến công cụ**, MR ước
lượng hiệu ứng nhân quả của yếu tố phơi nhiễm lên kết cục, vượt qua nhiễu
không đo lường được.

### Các Giả định Chính của MR

Để biến công cụ di truyền có giá trị, ba điều kiện phải được thỏa mãn:

1. **Tính liên quan (Relevance)**: Biến thể di truyền phải liên quan đến yếu
   tố phơi nhiễm.
2. **Tính độc lập (Independence)**: Biến thể di truyền không được liên quan
   đến các yếu tố gây nhiễu (điều này tuân theo phân ly Mendel trong các
   nghiên cứu được thiết kế tốt).
3. **Giới hạn loại trừ (Exclusion restriction)**: Biến thể di truyền chỉ
   được ảnh hưởng đến kết cục thông qua yếu tố phơi nhiễm, không thông qua
   các con đường khác.

### Các Loại Dữ liệu trong MR

| Loại dữ liệu | Là gì | Ví dụ |
|---|---|---|
| Cấp độ cá nhân | Một dòng cho mỗi người với kiểu gen, yếu tố phơi nhiễm, kết cục | Một đoàn hệ 10.000 người với dữ liệu SNP, đo LDL và hồ sơ bệnh viện |
| Cấp độ tóm tắt | Ước lượng hiệu ứng cho từng SNP từ GWAS | "SNP rs12345 có hiệu ứng 0,05 (SE 0,01) lên LDL" |
| Điểm nguy cơ đa gen (PRS) | Tổng có trọng số của nhiều SNP | PRS cho LDL cholesterol, tính từ 100 SNP |

`mrwin` sử dụng **dữ liệu cấp độ cá nhân** với một **công cụ PRS**. PRS đóng
vai trò như một công cụ liên tục nắm bắt khuynh hướng di truyền tích lũy đối
với yếu tố phơi nhiễm.

---

## Win Statistics: Estimand Mục tiêu

### Win Statistic là gì?

Một win statistic so sánh từng cặp cá nhân trong một nghiên cứu. Với mỗi cặp
(A, B), việc so sánh tuân theo một hệ thống phân cấp lâm sàng:

1. Trước tiên, so sánh trên kết cục quan trọng nhất (ví dụ: tử vong). Nếu A
   tử vong và B không, B thắng.
2. Nếu hòa trên kết cục đầu tiên, chuyển sang kết cục tiếp theo (ví dụ: nhập
   viện).
3. Tiếp tục xuống hệ thống phân cấp cho đến khi một người thắng hoặc tất cả
   các kết cục đều hòa.

**Win ratio** là số lần thắng của nhóm có yếu tố phơi nhiễm cao hơn chia cho
số lần thắng của nhóm có yếu tố phơi nhiễm thấp hơn.

### Tại sao Win Statistics Quan trọng

| Phương pháp | Đối xử biến cố như nhau? | Tôn trọng thứ tự lâm sàng? | Xử lý rủi ro cạnh tranh? |
|---|---|---|---|
| Thời gian đến biến cố đầu tiên | Có | Không | Không |
| Phân tích tiêu chí tổng hợp | Có | Không | Một phần |
| Win ratio / win statistic | Không | Có | Có |

Win statistics được giới thiệu bởi Pocock và cộng sự (2012) và đã được áp
dụng trong hơn 36 thử nghiệm lâm sàng ngẫu nhiên từ năm 2022 đến 2024.

### Causal Win Ratio (cCWR)

Các win statistics tiêu chuẩn yêu cầu phân công ngẫu nhiên hoặc giả định bỏ
qua mạnh (strong ignorability). Trong dịch tễ học quan sát, giả định này hiếm
khi có thể bảo vệ được.

`mrwin` định nghĩa **continuous Causal Win Ratio (cCWR)** như một estimand
gradient xuyên suốt các phân vị của công cụ di truyền:

- Chia quần thể thành các tầng có thứ tự theo PRS (ví dụ: decile).
- Trong mỗi cặp tầng liền kề, tính win ratio.
- Gộp chung xuyên suốt các tầng sử dụng gradient được chuẩn hóa theo công cụ
  và phân tích tổng hợp GLS.

Kết quả là một **hiệu ứng nhân quả biên, cấp độ quần thể** không yêu cầu giả
định bỏ qua mạnh.

---

## Công thức Cốt lõi: Giải thích bằng Ngôn ngữ Đơn giản

### 1. Kernel So sánh Phân cấp

Với hai cá nhân *i* và *j*, so sánh các kết cục từ mức ưu tiên cao nhất
xuống thấp nhất:

- Nếu *j* có biến cố tại mức ưu tiên *k* nhưng *i* không (hoặc *i* sống sót
  lâu hơn): *i* thắng (+1).
- Nếu *i* có biến cố tại mức ưu tiên *k* nhưng *j* không (hoặc *j* sống sót
  lâu hơn): *i* thua (-1).
- Nếu hòa: chuyển sang mức ưu tiên tiếp theo.

**Diễn đạt bằng lời**: "Người này có hồ sơ kết cục tốt hơn người kia không,
tôn trọng thứ tự mức độ nghiêm trọng lâm sàng?"

### 2. Win Ratio Theo Tầng

Trong một tầng (ví dụ: decile PRS thứ 6 so với thứ 5):

```
theta_d = (số lần thắng của tầng cao hơn) / (số lần thua)
log_theta_d = log(theta_d)
```

**Diễn đạt bằng lời**: "Trong số những người có khuynh hướng di truyền cao
hơn một chút đối với yếu tố phơi nhiễm, họ có xu hướng thắng thường xuyên
hơn thua không?"

### 3. Gradient được Chuẩn hóa theo Công cụ (ISG)

```
Delta_X_d = trung_bình(yếu_tố_phơi_nhiễm trong tầng d) - trung_bình(yếu_tố_phơi_nhiễm trong tầng d-1)
ISG_d = log_theta_d / Delta_X_d
```

**Diễn đạt bằng lời**: "Log win ratio thay đổi bao nhiêu cho mỗi đơn vị tăng
của yếu tố phơi nhiễm, như được dự đoán bởi di truyền?"

### 4. Gộp chung GLS

Các giá trị ISG trên tất cả các cặp tầng được gộp chung sử dụng Bình phương
Tối thiểu Tổng quát (Generalized Least Squares) với một ước lượng hiệp
phương sai co rút (shrinkage):

```
delta_GLS = ISG được gộp chung trên tất cả các tầng
DS-CWR = exp(delta_GLS)
```

**Diễn đạt bằng lời**: "Hiệu ứng nhân quả tổng thể của yếu tố phơi nhiễm lên
hồ sơ kết cục được ưu tiên là gì?"

- **DS-CWR > 1**: Yếu tố phơi nhiễm liên quan đến một hồ sơ kết cục tốt hơn.
- **DS-CWR < 1**: Yếu tố phơi nhiễm liên quan đến một hồ sơ kết cục xấu hơn.
- **DS-CWR = 1**: Không có bằng chứng về hiệu ứng nhân quả.

### 5. Multiplier Bootstrap

Độ không chắc chắn được ước lượng bằng cách:
1. Làm nhiễu các trọng số GWAS (độ không chắc chắn bên ngoài).
2. Rút các trọng số nhân ngẫu nhiên (độ không chắc chắn bên trong).
3. Tính toán lại toàn bộ pipeline cho mỗi lần lặp bootstrap.

Điều này cho ra một khoảng tin cậy có tính đến cả hai nguồn không chắc chắn.

### 6. Chẩn đoán Pleiotropy Giảm dần (SDPD)

SDPD kiểm tra liệu công cụ di truyền có ảnh hưởng đến kết cục ưu tiên cao
nhất thông qua các con đường khác ngoài yếu tố phơi nhiễm hay không
(pleiotropy trực tiếp). Việc này được thực hiện bằng cách khớp một hồi quy
MR-Egger trên các hiệu ứng trên yếu tố phơi nhiễm và kết cục cho từng SNP.

**Diễn đạt bằng lời**: "Công cụ di truyền có hợp lệ không, hay nó ảnh hưởng
đến kết cục thông qua các con đường khác?"

---

## Ví dụ Lâm sàng

### Ví dụ 1: Bệnh Tim mạch

**Câu hỏi nghiên cứu**: Liệu LDL cholesterol được dự đoán bởi di truyền có
làm xấu đi quỹ đạo tim-thận tổng thể?

**Hệ thống phân cấp tiêu chí** (từ ưu tiên cao nhất đến thấp nhất):
1. Tử vong
2. Nhập viện do suy tim
3. Suy giảm chức năng thận

**Diễn giải**: Nếu DS-CWR = 0,82 (95% CI: 0,70 đến 0,96), điều này có nghĩa
là cho mỗi đơn vị tăng LDL được dự đoán bởi di truyền, hồ sơ kết cục tim-thận
được ưu tiên xấu đi 18%. CI không chứa 1, gợi ý một hiệu ứng có hại có ý
nghĩa thống kê.

**Ý nghĩa lâm sàng**: LDL cao không chỉ gây đau tim -- nó dịch chuyển toàn bộ
quỹ đạo bệnh theo hướng tử vong sớm hơn, nhiều lần nhập viện hơn và suy giảm
chức năng thận nhanh hơn, theo đúng thứ tự lâm sàng đó.

### Ví dụ 2: Sa sút Trí tuệ

**Câu hỏi nghiên cứu**: Liệu rối loạn giấc ngủ được dự đoán bởi di truyền có
dịch chuyển con người theo hướng các kết cục thần kinh nhận thức nghiêm trọng
sớm hơn?

**Hệ thống phân cấp tiêu chí**:
1. Tử vong
2. Chẩn đoán sa sút trí tuệ
3. Nhập viện dưỡng lão

**Diễn giải**: Nếu DS-CWR = 0,91 (95% CI: 0,78 đến 1,06), ước lượng điểm gợi
ý tác hại, nhưng CI chứa 1. Không có bằng chứng mạnh mẽ về hiệu ứng nhân quả.

**Tại sao điều này quan trọng**: Trong MR tiêu chuẩn về giấc ngủ và sa sút
trí tuệ, những người tử vong trước khi phát triển sa sút trí tuệ bị mất khỏi
phân tích. Đây là một dạng thiên lệch người sống sót (survivor bias): một
biến thể giấc ngủ có hại có thể trông có vẻ bảo vệ vì những người mang nó
không bao giờ sống đủ lâu để được chẩn đoán. `mrwin` tránh điều này bằng cách
coi tử vong là biến cố ưu tiên cao nhất, do đó sự sống sót khác biệt được hấp
thụ vào estimand thay vì bị điều kiện hóa đi.

### Ví dụ 3: Dịch tễ học Ung thư

**Câu hỏi nghiên cứu**: Liệu BMI được dự đoán bởi di truyền có làm xấu đi các
kết cục ung thư được ưu tiên?

**Hệ thống phân cấp tiêu chí**:
1. Tử vong do ung thư
2. Tiến triển hoặc di căn
3. Tái phát

**Diễn giải**: Nếu DS-CWR = 0,75 (95% CI: 0,60 đến 0,94), điều này có nghĩa
BMI được dự đoán bởi di truyền cao hơn dịch chuyển quỹ đạo ung thư theo hướng
kết cục xấu hơn ở mọi cấp độ của hệ thống phân cấp.

**Ý nghĩa lâm sàng**: BMI không chỉ làm tăng nguy cơ ung thư -- nó làm xấu đi
toàn bộ diễn tiến bệnh từ tái phát đến tử vong.

---

## Cài đặt

```r
# Từ GitHub (phiên bản phát triển)
# install.packages("devtools")
devtools::install_github("nguyenminh2301/mrwin")
```

---

## Bắt đầu Nhanh

### Dữ liệu Mô phỏng (Không cần Đoàn hệ Thực)

```r
library(mrwin)

# Mô phỏng một tập dữ liệu tim-thận
cfg <- mrwin_config(n_outcome = 500, m_snps = 20, seed = 1)
dat <- mrwin_simulate(cfg)

# Khớp mô hình
fit <- mrwin(
  endpoint = mrwin_endpoint(dat$time, dat$status, c("death", "hf", "renal")),
  genotype = dat$G,
  exposure = dat$X,
  gwas = mrwin_gwas(dat$true_betas, rep(0.01, length(dat$true_betas))),
  controls = mrwin_controls(n_strata = 5, bootstrap = 100, seed = 2)
)

# Xem kết quả
print(fit)
summary(fit)
plot(fit)
```

### Đọc Kết quả Đầu ra

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

| Trường | Ý nghĩa | Cần chú ý gì |
|---|---|---|
| `delta_GLS` | Hiệu ứng gộp chung trên thang log | Âm = có hại, dương = có lợi |
| `DS-CWR` | Hiệu ứng lũy thừa hóa (win ratio) | < 1 = có hại, > 1 = có lợi, 1 = không hiệu ứng |
| `95% CI` | Khoảng tin cậy bootstrap | Có chứa 1 không? |
| `Q` | Thống kê không đồng nhất | Giá trị p nhỏ = hiệu ứng thay đổi giữa các tầng |
| `Warnings` | Cảnh báo có cấu trúc | Luôn kiểm tra trước khi diễn giải |

### Diễn giải Kết quả

Nếu DS-CWR = 0,82 (95% CI: 0,70 đến 0,96):

- Yếu tố phơi nhiễm liên quan đến hồ sơ kết cục được ưu tiên xấu đi 18%.
- CI không chứa 1, do đó hiệu ứng có ý nghĩa thống kê.
- Kiểm tra giá trị p của Q: nếu có ý nghĩa, hiệu ứng có thể thay đổi xuyên
  suốt phạm vi của yếu tố phơi nhiễm.
- Kiểm tra cảnh báo: nếu `sdpd_rejected`, pleiotropy có thể làm mất giá trị
  của ước lượng.

---

## Mẫu cho Đoàn hệ Thực

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

## Định dạng Đầu ra

### Print

```r
print(fit)
```

Hiển thị kết quả chính: cỡ mẫu, SNPs, số tầng, tính hợp lệ của bootstrap,
DS-CWR, CI, thống kê Q và cảnh báo.

### Summary

```r
summary(fit)
```

Tạo ra một báo cáo có định dạng với:
- Ước lượng DS-CWR chính, CI chính (Fieller), giá trị p và CI giới hạn
  pleiotropy tùy chọn
- Các gradient tầng liền kề (log-theta, Delta-X, ISG, CWR)
- Thống kê không đồng nhất Q
- CI delta-method được giữ lại như một tham chiếu có nhãn (CI Fieller là CI
  chính; xem mục *Khoảng tin cậy* bên dưới)
- Chẩn đoán SDPD (khi được kích hoạt)
- Cảnh báo có cấu trúc

### Tidy (cho bảng biểu và phân tích sâu hơn)

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

### Biểu đồ

```r
plot(fit, type = "isg")       # Các gradient ISG liền kề (mặc định)
plot(fit, type = "forest")    # Forest plot của CWR liền kề
plot(fit, type = "bootstrap") # Phân phối bootstrap của log-theta
```

---

## Hiệu chỉnh Hiệp biến

Khi các yếu tố gây nhiễu được đo lường (tuổi, giới tính, principal
components), sử dụng hiệu chỉnh ordinal IPTW:

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

Hiệu chỉnh IPTW cải thiện độ chính xác bằng cách cân bằng các hiệp biến xuyên
suốt các tầng PRS. SDPD kiểm tra pleiotropy trực tiếp.

---

## Điểm mạnh và Hạn chế

### Điểm mạnh

- **Suy luận nhân quả dưới nhiễu không đo lường được**: Sử dụng công cụ di
  truyền, không chỉ dựa vào yếu tố phơi nhiễm quan sát được.
- **Tôn trọng hệ thống phân cấp lâm sàng**: Tử vong được ưu tiên hơn nhập
  viện, và nhập viện được ưu tiên hơn suy giảm dấu ấn sinh học.
- **Xử lý rủi ro cạnh tranh một cách tự nhiên**: Sự sống sót khác biệt được
  hấp thụ vào estimand, không bị điều kiện hóa đi.
- **Một thước đo tóm tắt duy nhất**: Một DS-CWR thay thế nhiều phân tích
  thành phần riêng lẻ.
- **Khung chẩn đoán**: SDPD kiểm tra giới hạn loại trừ; thống kê Q kiểm tra
  không đồng nhất; cảnh báo có cấu trúc đánh dấu các vấn đề.
- **Có thể tái lập**: Bootstrap được kiểm soát bởi seed; tất cả các đầu ra
  là xác định với seed đã cho.

### Hạn chế

- **Yêu cầu dữ liệu cấp độ cá nhân**: Gói hiện tại yêu cầu dữ liệu một dòng
  cho một người tham gia với kiểu gen, yếu tố phơi nhiễm và kết cục. MR chỉ
  với dữ liệu tóm tắt không được hỗ trợ.
- **Không có tính collapsible**: Win ratio không có tính collapsible. DS-CWR
  biên khác với win ratio có điều kiện. Đây là một tính chất của estimand,
  không phải lỗi.
- **Công cụ PRS tuyến tính**: Gói giả định một điểm nguy cơ đa gen tuyến
  tính. Các tương tác gen-gen phi tuyến không được mô hình hóa.
- **Hiệp phương sai GWAS dạng chéo**: Bootstrap hiện tại sử dụng sai số
  chuẩn cho từng SNP, không phải ma trận hiệp phương sai LD đầy đủ.
  Điều này có thể đánh giá thấp độ không chắc chắn khi các SNP ở trạng thái
  linkage disequilibrium.
- **Tính dễ tổn thương với pleiotropy**: cCWR nhạy cảm với pleiotropy làm
  nhiễm hệ thống phân cấp. Pleiotropy ở mức tử vong nhỏ như gamma = 0,05 có
  thể làm sụp đổ độ phủ xuống còn 12%.
- **Yêu cầu về cỡ mẫu**: Mẫu quy mô biobank (N > 100.000) là điều kiện tiên
  quyết thống kê nghiêm ngặt cho suy luận đáng tin cậy.
- **Số tầng `D`**: ước lượng kiểu thập phân vị đòi hỏi chọn số tầng. Một ước
  lượng gradient liên tục, điều khiển bằng băng thông (mà ước lượng thập phân vị
  là trường hợp boxcar đặc biệt chính xác) đã được hiện thực và thẩm định, nhưng
  kiểm tra ngoại vi cho thấy nó nhiễu hơn ước lượng thập phân vị và không làm giảm
  độ nhạy theo `D` (gradient chuẩn-hoá-theo-công-cụ là một tỉ số, làm mịn càng
  nhỏ thì mẫu số càng co lại). Khuyến nghị là dùng ước lượng rời rạc kèm phân
  tích độ nhạy theo `D`, không tái tham số hoá liên tục
  (xem `inst/spec/wp15-continuous-isg.md`).

### Backend hiệu năng và suy luận (tuỳ chọn)

Backend mặc định không đổi, nhưng công việc Giai đoạn II
(`inst/spec/acceleration-roadmap.md`) đã bổ sung các lựa chọn tuỳ chọn đã thẩm
định qua `mrwin_controls()`:

- `backend = "fast"` — nhân win/loss dưới-bậc-hai, biên dịch (Rcpp)
  (`Theta(N log^{K-1} N)` cho `K` mức ưu tiên), giống hệt từng bit so với backend
  dense và đưa ước lượng tới quy mô biobank (ví dụ `K = 3`, `N = 80.000` dưới một
  giây).
- `inference = "analytic"` — phương sai hàm ảnh hưởng dạng đóng tái tạo lại
  multiplier bootstrap (kèm một số hạng Monte-Carlo chính xác cho độ không chắc
  chắn của trọng số GWAS), loại bỏ vòng lặp bootstrap cho thành phần lấy mẫu.

> **Phân tầng xếp-hạng-kép — không khuyến nghị (phát hiện phủ định).**
> `mrwin_doubly_ranked_strata()` (Tian/Burgess) đã được hiện thực, nhưng thẩm
> định hiệu chuẩn ngoại vi cho thấy `stratification = "doubly_ranked"` **không
> tương thích** với estimand DS-CWR giữa-các-tầng: nó cân bằng công cụ giữa các
> tầng, khiến tương phản tầng-liền-kề bị dẫn dắt bởi nhiễu gây nhiễu (sai số loại
> I ~1,0 khi có gây nhiễu). `mrwin()` vẫn chạy được nhưng phát cảnh báo
> `doubly_ranked_invalid`; mặc định là `"prs_rank"`. Xem
> `inst/spec/validation-findings.md`.

#### Khoảng tin cậy

**Khoảng tin cậy chính được báo cáo là khoảng Fieller.** Mô phỏng thẩm định hiệu
chuẩn cho thấy khoảng delta-method tỉ số ban đầu quá rộng (không đạt sai số loại
I danh nghĩa), trong khi cấu trúc Fieller trên cùng ma trận hiệp phương sai được
hiệu chuẩn đúng; khoảng delta-method chỉ được giữ lại như một tham chiếu có nhãn.
Khi công cụ yếu, khoảng Fieller được báo cáo trung thực là không giới hạn thay vì
một khoảng hẹp sai lệch.

---

## Các Giả định

Khung `mrwin` dựa trên các giả định sau:

1. **Tính liên quan (Relevance)**: PRS phải liên quan đến yếu tố phơi nhiễm.
   Công cụ yếu tạo ra các ước lượng không ổn định (kiểm tra cảnh báo
   `weak_instrument`).

2. **Tính độc lập (Independence)**: PRS không được liên quan đến các yếu tố
   gây nhiễu. Điều này được thỏa mãn bởi phân ly Mendel trong các quần thể
   đồng nhất. Phân tầng quần thể có thể vi phạm điều này; sử dụng principal
   components làm hiệp biến.

3. **Giới hạn loại trừ (Exclusion restriction)**: PRS chỉ được ảnh hưởng đến
   kết cục thông qua yếu tố phơi nhiễm. SDPD kiểm tra điều này cho tiêu chí
   ưu tiên cao nhất. Nếu bị bác bỏ, ước lượng có thể không có giá trị.

4. **Tính đơn điệu (Monotonicity)**: Công cụ phải dịch chuyển yếu tố phơi
   nhiễm theo cùng một hướng cho tất cả các cá nhân. Vi phạm có thể gây
   thiên lệch cho ước lượng.

5. **Không có điều chỉnh hiệu ứng bởi sự sống sót**: Hiệu ứng nhân quả không
   được thay đổi giữa những người sống sót và những người không sống sót.
   Điều này không thể kiểm tra được.

6. **Thứ tự ưu tiên chính xác**: Hệ thống phân cấp lâm sàng phải được chỉ
   định trước khi phân tích. Gói không tự suy luận nó từ dữ liệu.

---

## Mã Cảnh báo

| Mã | Ý nghĩa | Cần làm gì |
|---|---|---|
| `weak_instrument` | Dịch chuyển kiểu hình gần bằng 0 hoặc CI Fieller không giới hạn | Coi ước lượng là không ổn định |
| `sdpd_rejected` | Hệ số chặn MR-Egger gợi ý pleiotropy | Thêm thảo luận về độ nhạy |
| `sdpd_underpowered` | Quá ít SNP cho SDPD | Báo cáo hạn chế về lực thống kê |
| `positivity_failure` | IPTW ESS quá nhỏ trong một tầng | Kiểm tra sự chồng lấp hiệp biến |
| `bridged_strata` | Các tầng bị bỏ qua do thất bại positivity | Báo cáo sơ đồ bắc cầu |
| `discordant_components` | Hướng MR thành phần xung đột | Không tuyên bố quá mức hiệu ứng toàn cục |

---

## Công cụ Mô phỏng

Chạy lưới kịch bản để kiểm định phương pháp:

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

| Kịch bản | Ý nghĩa | Hành vi kỳ vọng |
|---|---|---|
| A_null | Không có hiệu ứng nhân quả | DS-CWR gần 1, tỷ lệ bác bỏ gần 5% |
| B_valid_IV | Hiệu ứng nhân quả hợp lệ | DS-CWR theo hướng kỳ vọng |
| C_pleiotropy | Pleiotropy trực tiếp | Tỷ lệ bác bỏ SDPD tăng |
| D_hierarchy_discordant | Hiệu ứng thành phần đối lập | Cảnh báo được phát ra |

---

## Trích dẫn Gói này

Nếu bạn sử dụng `mrwin` trong nghiên cứu của mình, vui lòng trích dẫn:

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

**Tác giả:**
- Nguyen Thien Minh, Đại học Y Dược Thành phố Hồ Chí Minh,
  Việt Nam (minhnt@ump.edu.vn)
- N. Ahmad Aziz (tác giả liên hệ), Trung tâm Bệnh Thoái hóa Thần kinh Đức
  (DZNE) và Đại học Bonn, Đức (Ahmad.Aziz@dzne.de)

---

## Tài liệu Tham khảo Liên quan

**Win statistics:**
- Pocock SJ, et al. The win ratio: a new approach to the analysis of
  composite endpoints in clinical trials. *Eur Heart J*. 2012;33(14):1744-1749.
- Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a
  composite endpoint based on prioritized components. *Biostatistics*.
  2016;17(1):178-191.
- Even Z, Josse A. Causal win ratio. *arXiv preprint*. 2025.

**Mendelian randomization:**
- Lawlor DA, et al. Mendelian randomization: using genes as instruments for
  making causal inferences in epidemiology. *Stat Med*. 2008;27(8):1133-1163.
- Davey Smith G, Hemani G. Mendelian randomization: genetic anchors for
  causal inference in epidemiological studies. *Hum Mol Genet*.
  2014;23(R1):R89-R98.
- Bowden J, et al. Mendelian randomization with invalid instruments: effect
  estimation and bias detection through Egger regression. *Int J Epidemiol*.
  2015;44(2):512-525.

**Polygenic risk scores:**
- Choi SW, et al. Tutorial: a guide to performing polygenic risk score
  analyses. *Nat Protoc*. 2020;15(9):2759-2772.

**Multi-state models:**
- Putter H, Fiocco M, Geskus RB. Tutorial in biostatistics: competing risks
  and multi-state models. *Stat Med*. 2007;26(11):2389-2430.

---

## Cấu trúc Gói

```
mrwin/
  R/
    api.R              # Quy trình mrwin() cấp cao
    kernel.R           # Kernel so sánh phân cấp
    estimate.R         # Ước lượng DS-CWR, gộp chung GLS
    bootstrap.R        # Suy luận multiplier bootstrap
    adjustment.R       # Ordinal-IPTW, ESS, bắc cầu
    sdpd.R             # Chẩn đoán Pleiotropy Giảm dần (SDPD)
    simulate.R         # Quy trình tạo dữ liệu mô phỏng
    simulation_engine.R # Chạy lưới kịch bản
    validation.R       # Kiểm tra đầu vào
    methods.R          # print, summary, plot, tidy, report
    benchmark.R        # Đo điểm chuẩn hiệu năng
    backend_sparse.R   # Backend thưa
    config.R           # Cấu hình mô phỏng
    strata.R           # Gán tầng PRS
  tests/
    testthat/          # 273 bài kiểm tra đơn vị
  vignettes/           # 3 vignette
  inst/spec/           # Đặc tả triển khai
```

---

## Giấy phép

MIT
