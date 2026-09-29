# Lộ trình hoàn thiện `mrwin` (Giai đoạn III): mục tiêu, kết quả kỳ vọng và kế hoạch thực hiện — Bản học thuật

| Mục | Nội dung |
|---|---|
| Phiên bản | 29/09/2026 |
| Trạng thái | **G0 hoàn thành 29/09/2026** (§6.1); D1 và D3 đã chốt; D2, D4, D5 còn mở (§8). Giai đoạn kế tiếp: G1 |
| Tài liệu gốc | `acceleration-roadmap.md`, `phase2-direction.md`, `wp18-theory-foundation.md`, `wp19-scalability-validation.md`, `validation-findings.md`, `project-checkpoint.md`, `algorithm-spec.md` |
| Bản phổ thông tương ứng | `completion-roadmap-plain.vi.md` (cùng mã KQ1–KQ5, G0–G6, D1–D5) |
| Ràng buộc vận hành | Không sử dụng hạ tầng điện toán đám mây cho tính toán hay lưu trữ dữ liệu (§6.8, §10) |

---

## Tóm tắt

**Bối cảnh.** `mrwin` hiện thực hoá *tỉ số thắng nhân quả chuẩn hoá theo liều* (DS-CWR). Đây
là một estimand biến công cụ, kết hợp so sánh cặp theo thứ bậc ưu tiên lâm sàng (win
statistics) với ngẫu nhiên hoá Mendel (MR) trên dữ liệu cá thể và công cụ là điểm nguy cơ đa gen
(PRS). Giai đoạn I (WP0–WP12) xây dựng đường ống ước lượng hoàn chỉnh. Giai đoạn II (WP13–WP19)
giải quyết nút thắt tính toán $\Theta(N^2/D)$, xây dựng phương sai giải tích dựa trên hàm ảnh
hưởng (influence function, IF), và phát hiện, sửa một khiếm khuyết hiệu chuẩn của khoảng tin cậy
chính.

**Mục tiêu.** Hoàn tất quá trình chuyển DS-CWR từ một bản hiện thực đã nhất quán nội tại thành
một **phương pháp đã được thẩm định ngoại tại, đã qua bình duyệt và được phát hành công khai**.

**Cách tiếp cận.** Bảy giai đoạn (G0–G6) được tổ chức quanh ba nguyên tắc: (i) không phát biểu
lý thuyết nào trước khi có kiểm thử tương ứng; (ii) ưu tiên thẩm định ngoại tại (hiệu chuẩn so
với chân lý đã biết) hơn nhất quán nội tại; (iii) định trước giao thức mô phỏng và kế hoạch phân
tích.

**Kết quả kỳ vọng.** KQ1: bài báo phương pháp đã bình duyệt. KQ2: `mrwin` 1.0.0 trên CRAN. KQ3:
nghiên cứu mô phỏng theo khung ADEMP, tái lập được. KQ4: ít nhất một ứng dụng trên dữ liệu thật,
báo cáo theo STROBE-MR. KQ5: tài liệu người dùng hoàn chỉnh. Thời gian ước tính đến lúc nộp bản
thảo là khoảng 16 tuần làm việc tập trung, với điều kiện dữ liệu thật được cấp đúng hạn.

---

## 1. Bối cảnh khoa học và khoảng trống nghiên cứu

### 1.1 Tiêu chí tổng hợp phân cấp và win statistics

Phân tích thời-gian-đến-biến-cố-đầu-tiên trên tiêu chí tổng hợp gán trọng số ngầm như nhau cho
các thành phần có mức độ nghiêm trọng khác nhau. Tỉ số thắng (win ratio) so sánh từng cặp theo
thứ bậc ưu tiên định trước [1]. Suy luận mẫu lớn dựa trên U-thống kê đa mẫu đã được xây dựng [2].
Estimand của win ratio phụ thuộc vào phân phối kiểm duyệt và độ dài theo dõi [3]. Vì vậy mọi
estimand dạng win ratio cần được định nghĩa kèm một **cửa sổ theo dõi** tường minh.

### 1.2 Ngẫu nhiên hoá Mendel

MR dùng biến thể di truyền làm biến công cụ để ước lượng hiệu ứng nhân quả khi có nhiễu không đo
lường được. Hai mối đe doạ chính là công cụ yếu, vốn gây chệch về phía liên kết quan sát [4], và
pleiotropy định hướng, được phát hiện qua hồi quy MR-Egger [5]. Các biến thể phi tuyến dùng
phân tầng (phương pháp phần dư, phương pháp xếp hạng kép) nhắm tới hiệu ứng *trong* tầng [6].

### 1.3 Khoảng trống

MR tiêu chuẩn nhắm tới một kết cục đơn lẻ hoặc thời gian đến biến cố đầu tiên. Khi kết cục không
tử vong bị cắt cụt bởi tử vong, phân tích có điều kiện trên sự sống sót gây chệch chọn lọc
(survivor bias). Hiện chưa có khung nào đồng thời đáp ứng bốn yêu cầu: (a) tôn trọng thứ bậc ưu
tiên lâm sàng, trong đó tử vong được hấp thụ vào estimand thay vì bị điều kiện hoá; (b) định danh
nhân quả nhờ công cụ di truyền; (c) có suy luận được hiệu chuẩn; (d) khả thi ở quy mô biobank.

### 1.4 Đóng góp dự kiến của công trình

- **C1 — Estimand và ước lượng:** DS-CWR, xây dựng từ các gradient chuẩn hoá theo công cụ (ISG)
  của các cặp tầng PRS liền kề, gộp bằng GLS.
- **C2 — Thuật toán:** tính số thắng/thua phân cấp có trọng số bằng đếm trội đa chiều, độ phức
  tạp $\Theta(N\log^{K-1}N)$ với $K\le 3$, cài đặt bằng Rcpp. Kết quả trùng khớp từng bit với
  kernel dày (dense).
- **C3 — Suy luận:** (i) phương sai IF dạng đóng cộng với số hạng bất định GWAS chính xác;
  (ii) phát hiện rằng khoảng delta-method cho tỉ số bị hiệu chuẩn sai, trong khi khoảng Fieller
  trên cùng ma trận hiệp phương sai được hiệu chuẩn đúng.
- **C4 — Phát hiện phủ định có giá trị phương pháp luận:** ISG liên tục (M1) và phân tầng xếp
  hạng kép (M2) đều vượt mọi kiểm tra nhất quán nội tại nhưng thất bại khi hiệu chuẩn ngoại tại.
  Riêng M2 không tương thích về mặt cấu trúc với tương phản *giữa* các tầng.

---

## 2. Estimand, ước lượng và giả định

### 2.1 Dữ liệu và kernel phân cấp

Với cá thể $i=1,\dots,N$ và $K$ thành phần kết cục xếp theo ưu tiên giảm dần
$k=1,\dots,K$, quan sát $(T_{ik},\Delta_{ik})$ gồm thời gian quan sát và chỉ báo biến cố;
$X_i$ là phơi nhiễm và $G_i\in\mathbb{R}^M$ là kiểu gen. Kernel phân cấp
$h(i,j)\in\{-1,0,+1\}$ phản đối xứng, $h(i,j)=-h(j,i)$. Kernel so sánh lần lượt từ $k=1$ và
dừng ở thành phần đầu tiên phân định được thắng/thua dưới kiểm duyệt. Kernel chỉ phụ thuộc
$(T,\Delta)$ của cặp, không phụ thuộc $\beta$ hay cách phân tầng. Tính bất biến này là nền tảng
của C2.

### 2.2 Công cụ và phân tầng

Công cụ được định nghĩa là $S_i=G_i^\top\hat\beta$, với $\hat\beta$ lấy từ GWAS ngoài. Các cá
thể được xếp hạng theo $S$ và chia thành $D$ tầng phân vị $\mathcal{S}_1,\dots,\mathcal{S}_D$
(phân tầng `prs_rank`, mặc định).

### 2.3 Ước lượng

Với $d=2,\dots,D$ và trọng số cá thể $w_i$ (bằng 1, hoặc là trọng số IPTW ổn định hoá):

$$
W_d=\sum_{i\in\mathcal{S}_d}\sum_{j\in\mathcal{S}_{d-1}} w_iw_j\,\mathbf 1\{h(i,j)=+1\},\qquad
L_d=\sum_{i\in\mathcal{S}_d}\sum_{j\in\mathcal{S}_{d-1}} w_iw_j\,\mathbf 1\{h(i,j)=-1\},
$$

$$
\hat\theta_d=\frac{W_d}{L_d},\qquad \Delta\hat X_d=\bar X_d-\bar X_{d-1},\qquad
\hat\delta_d=\frac{\log\hat\theta_d}{\Delta\hat X_d}\ \ (\text{ISG}).
$$

$$
\hat\delta_{\mathrm{GLS}}=\left(\mathbf 1^\top\tilde\Sigma^{-1}\mathbf 1\right)^{-1}\mathbf 1^\top\tilde\Sigma^{-1}\hat{\boldsymbol\delta},\qquad
\tilde\Sigma=(1-\rho)\hat\Sigma_{\mathrm{ISG}}+\rho\,\mu I,\qquad
\widehat{\mathrm{DS\text{-}CWR}}=\exp(\hat\delta_{\mathrm{GLS}}),
$$

$$
Q=(\hat{\boldsymbol\delta}-\mathbf 1\hat\delta_{\mathrm{GLS}})^\top\tilde\Sigma^{-1}(\hat{\boldsymbol\delta}-\mathbf 1\hat\delta_{\mathrm{GLS}})\ \overset{H_0}{\sim}\ \chi^2_{D-2}.
$$

### 2.4 Suy luận hiện hành

Ma trận hiệp phương sai chung $\mathrm{cov}_u$ của vectơ $(\log\hat\theta_d,\Delta\hat X_d)_d$
được ước lượng theo một trong hai cách:

- **Multiplier bootstrap:** $\xi_i\sim\mathrm{Exp}(1)$ và
  $\beta^\*\sim N(\hat\beta,\mathrm{diag}(\mathrm{se}^2))$.
- **Giải tích:** $\Sigma_{\mathrm{sampling}}=C^\top C$ dựa trên IF, cộng với
  $\Sigma_{\mathrm{gwas}}$ tính bằng tái chọn mẫu chỉ trên GWAS. Cách tái chọn mẫu này là số
  hạng chính xác, vì ước lượng điểm là hàm hằng từng khúc theo $\beta$.

Khoảng tin cậy chính là **khoảng Fieller**, cài đặt tại `R/bootstrap.R::.mrwin_fieller_ci`.
Với trọng số $\omega_d\propto 1/[\hat\Sigma_{\mathrm{ISG}}]_{dd}$, đặt
$\bar u_1=\sum_d\omega_d\log\hat\theta_d$ và $\bar u_2=\sum_d\omega_d\Delta\hat X_d$. Khoảng
Fieller là tập $\{\delta:(\bar u_1-\delta\bar u_2)^2\le z^2\,\mathrm{Var}(\bar u_1-\delta\bar u_2)\}$,
và có thể không bị chặn khi công cụ yếu. Giá trị $p$ tương ứng kiểm định
$H_0:\ \mathbb E\,\bar u_1=0$.

### 2.5 Giả định định danh

| Mã | Giả định | Kiểm tra/chẩn đoán trong gói |
|---|---|---|
| A1 | Liên quan: $S$ liên hệ với $X$ | Cảnh báo `weak_instrument`; khoảng Fieller không bị chặn |
| A2 | Độc lập: $S\perp$ nhiễu | Hiệp biến PC tổ tiên; IPTW thứ bậc |
| A3 | Loại trừ: $S$ chỉ tác động lên kết cục qua $X$ | SDPD (MR-Egger trên thành phần ưu tiên 1); CI có giới hạn pleiotropy |
| A4 | Đơn điệu: công cụ dịch chuyển $X$ cùng chiều ở mọi cá thể | Không kiểm định được; thảo luận |
| A5 | Không có biến đổi hiệu ứng theo tình trạng sống sót | Không kiểm định được; thảo luận |
| A6 | Thứ bậc ưu tiên được định trước | Quy trình định trước (G4) |

Đặc tính của estimand cần nêu tường minh trong bài báo: (i) estimand là **biên** (marginal) và
**không khả gộp** (non-collapsible); (ii) estimand phụ thuộc cửa sổ theo dõi $\tau$ và cơ chế
kiểm duyệt [3], do đó phải định nghĩa DS-CWR$(\tau)$; (iii) phân tầng phải tạo khác biệt về
**công cụ** giữa các tầng. Đây là lý do cấu trúc khiến phân tầng xếp hạng kép bị loại
(`validation-findings.md`).

### 2.6 Bốn vấn đề lý thuyết mở cần khép lại trước khi nộp bài

1. **Định lý định danh chính thức.** Sổ nghĩa vụ chứng minh (`wp18-theory-foundation.md` §2)
   có các dòng về thuật toán (L-C1, T-C2), phương sai (V-M3), hiệu chuẩn (C-R2) và phát hiện phủ
   định (D-M1, S-M2), nhưng **chưa có dòng cho định lý định danh** của DS-CWR dưới A1–A6. Đề
   xuất bổ sung dòng **I-E1**.
2. **Căn chỉnh ước lượng điểm và khoảng tin cậy.** $\hat\delta_{\mathrm{GLS}}$ là trung bình
   (GLS) của các tỉ số, trong khi khoảng Fieller có tâm tại tỉ số của các trung bình có trọng
   số, $\hat\delta_R=\bar u_1/\bar u_2$. Nếu gradient đồng nhất, tức $\delta_d\equiv\delta$ và
   $\log\theta_d=\delta\,\Delta X_d$ với mọi $d$, thì $\delta_R=\delta_{\mathrm{GLS}}=\delta$
   ở mức quần thể. Khi gradient **không đồng nhất**, hai đại lượng nói chung khác nhau. Khi đó
   khoảng chính không bảo đảm chứa ước lượng điểm được báo cáo, và độ phủ phải được đánh giá
   theo từng đại lượng đích.
3. **Hệ số ≈ 2 của delta-method.** Phương sai delta-method của tỉ số lớn gấp khoảng 2 lần phân
   tán thực nghiệm. Nguyên nhân chưa được xác định (`validation-findings.md`, mục "Optional
   follow-up").
4. **Tính bảo thủ khi $\sigma_\beta>0$.** Sai lầm loại I giảm xuống 0,007–0,020, gợi ý rằng bất
   định GWAS có thể đang được lan truyền quá mức.

---

## 3. Hiện trạng: bằng chứng đã có và hạng mục mở

### 3.1 Kết quả đã được thẩm định

| Kết quả | Bằng chứng | Nguồn |
|---|---|---|
| Kernel nhanh $K\in\{1,2,3\}$ trùng khớp kernel dày | Kiểm thử vi sai: 16k và 28k đoàn hệ ngẫu nhiên, 0 sai khác | `test-kernel-fast.R`, `tests/python/test_kernel_fast.py` |
| Khả năng mở rộng | $K=3$: nhanh hơn 72× ($n=2000$) và 300× ($n=5000$); $N=80\,000$ trong 0,78 s; số mũ thực nghiệm ≈ 1,1–1,2 | `benchmark-results.md` |
| Phương sai giải tích bằng bootstrap | Tỉ số se 0,997–1,009 (không hiệu chỉnh, $\sigma_\beta\in\{0,05;0,15;0,30\}$); 0,985–1,026 (IPTW) | `test-analytic-variance.R`, `wp17-analytic-variance.md` |
| Fieller được hiệu chuẩn đúng | Loại I 0,055, độ phủ 0,955 ($N=6000$, $M=200$); lưới WP19: loại I 0,007–0,060, độ phủ 0,960 | `test-fieller-primary.R`, `validation-findings.md` |
| Bộ kiểm thử | 93 nhóm, 0 lỗi (R 4.3.3, 06/2026) | `project-checkpoint.md` |

### 3.2 Phát hiện phủ định

| Hướng | Kết quả | Cơ chế |
|---|---|---|
| M1: ISG liên tục làm trơn bằng kernel | Nhiễu gấp 2,5–3,7 lần so với phân vị; không làm giảm độ nhạy theo $D$ | ISG là tỉ số; làm trơn mịn hơn làm mẫu số $\Delta X$ co lại |
| M2: phân tầng xếp hạng kép | Loại I ≈ 1,000; $\hat\delta$ trung bình −0,179 ± 0,017 dưới giả thuyết không | Cân bằng công cụ giữa các tầng nên tương phản giữa tầng bị nhiễu chi phối |
| Delta-method cho tỉ số | Loại I 0,000, độ phủ 0,995 | Phương sai bị ước lượng quá khoảng 2 lần (§2.6.3) |

### 3.3 Hạng mục mở

| Hạng mục | Nguồn | Ưu tiên | Giai đoạn |
|---|---|---|---|
| Định lý định danh (I-E1); căn chỉnh $\delta_{\mathrm{GLS}}$ và $\delta_R$ | §2.6 | P0 | G1 |
| Ledoit–Wolf chính xác hoặc được tài liệu hoá và thẩm định | `algorithm-spec.md` §4.6, §10 | P0 | G1 |
| Hệ số ≈ 2 của delta-method; tính bảo thủ khi $\sigma_\beta>0$ | `validation-findings.md` | P1 | G1 |
| $\Sigma_{\mathrm{GWAS}}$ đầy đủ (có LD), thay cho ma trận chéo | README, mục Hạn chế | P1 | G1 (D3) |
| AL-CWR; nội suy chệch theo Table 3 của v5; GPS | `algorithm-spec.md` §10 | P2 | G1 (D5, D3) |
| Parity gate `mrwin_verify_fast_dense_parity()` | `wp19-scalability-validation.md` T1 | P0 | G3 |
| Benchmark mở rộng trong R; kiểm định cận số mũ | WP19 T2 | P1 | G3 |
| Lưới thống kê đầy đủ (nhiều $N/D$, pleiotropy/SDPD, discordant, công cụ yếu, lực kiểm định) | WP19 T3 | P0 | G2 |
| Đổi mặc định `backend = "fast"` | WP19 T4 | P1 | G3 |
| `R CMD check --as-cran`, độ bao phủ ≥ 80% | WP19 T5; `release-checklist.md` | P0 | G3 |
| Bản thảo và nghĩa vụ T4 (truy vết mọi con số) | WP18 | P0 | G5 |

### 3.4 Sai lệch giữa tài liệu và hiện trạng (✅ đã xử lý trong G0, 29/09/2026)

Các mục dưới đây là trạng thái trước G0. Cách xử lý từng mục được ghi trong
`project-checkpoint.md`, mục "G0 record". Riêng CI: việc CI chỉ chạy trên `main` được giữ nguyên
có chủ đích, vì `main` là nhánh tích hợp và kiểm định cục bộ là chuẩn chính thức theo ràng buộc
không-đám-mây.

- `CITATION.cff`: phiên bản `0.1.0-arxiv-theory`; phần tóm tắt ghi "không chứa kết quả mô
  phỏng đã thẩm định". Nội dung này đã lỗi thời.
- `DESCRIPTION`: phiên bản `0.0.0.9000`.
- README ghi "273 bài kiểm tra" và "3 vignette", nhưng thư mục `vignettes/` **không tồn tại**
  trong kho. Quy tắc `*.Rmd` trong `.gitignore` đã loại các tệp này, trong khi `DESCRIPTION`
  vẫn khai báo `VignetteBuilder: knitr`.
- Trích dẫn Pocock et al. (2012) trong README ghi 33(14):1744–1749. Thông tin đúng là
  *Eur Heart J* 33(2):176–182 [1].
- `.github/workflows/ci.yml` chỉ kích hoạt trên `main`. `release-checklist.md` vẫn mang ngày
  2026-05-10.

---

## 4. Mục tiêu tối thượng và khung logic

### 4.1 Phát biểu mục tiêu tối thượng

> Xác lập DS-CWR như một **phương pháp suy luận nhân quả đã được thẩm định, bình duyệt và phổ
> biến**. Phương pháp cho phép kiểm định và ước lượng tác động nhân quả của một phơi nhiễm được
> dự đoán bởi di truyền lên **toàn bộ hồ sơ kết cục được xếp ưu tiên lâm sàng**, với suy luận
> được hiệu chuẩn, khả thi ở quy mô biobank, và có chẩn đoán vi phạm giả định.

### 4.2 Khung logic

| Cấp | Mô tả | Chỉ số kiểm chứng khách quan | Phương tiện kiểm chứng | Giả định và rủi ro |
|---|---|---|---|---|
| **Tác động** | DS-CWR được cộng đồng dịch tễ di truyền sử dụng cho các câu hỏi có tiêu chí tổng hợp | Trích dẫn, lượt tải CRAN, ứng dụng độc lập | Chỉ số trích dẫn, thống kê CRAN | Phương pháp được chấp nhận về khái niệm |
| **Kết quả** | KQ1–KQ5 (§5) | Như tiêu chí chấp nhận ở §5 | Thư chấp nhận, trang CRAN, báo cáo mô phỏng | Dữ liệu được cấp; không phát sinh khiếm khuyết cấu trúc mới |
| **Đầu ra** | Giao thức và kết quả mô phỏng; gói 1.0.0; bản thảo; phân tích ứng dụng | Tệp đã commit; parity gate xanh; bản thảo đã xây dựng | Kho mã (hash commit) | Tài nguyên tính toán cục bộ đủ dùng |
| **Hoạt động** | G0–G6 (§6) | Các mốc ở §7 | Sổ trạng thái trong `project-checkpoint.md` | Đủ nhân lực |
| **Đầu vào** | Nhân lực (2 tác giả chính), máy trạm/máy chủ cục bộ, dữ liệu đoàn hệ, khoảng 16 tuần | — | — | Được hội đồng đạo đức phê duyệt |

---

## 5. Kết quả kỳ vọng cuối cùng: định nghĩa vận hành và tiêu chí chấp nhận

| Mã | Định nghĩa vận hành | Tiêu chí chấp nhận |
|---|---|---|
| **KQ1** | Bài báo phương pháp (C1–C4 kèm ứng dụng) trên tạp chí có bình duyệt; arXiv v2 thay bản chỉ-lý-thuyết | Có thư chấp nhận. Mọi dòng trong sổ WP18 (kể cả I-E1) ở trạng thái *discharged*. Mọi con số truy được về một kiểm thử hoặc benchmark đã commit (WP18 T4) |
| **KQ2** | `mrwin` 1.0.0 trên CRAN, lưu trữ kèm DOI | `R CMD check --as-cran`: 0 ERROR, 0 WARNING. Độ bao phủ dòng ≥ 80%. Parity gate xanh. Mặc định: `backend="fast"` (tự lùi về dense khi $K\ge4$), Fieller, `prs_rank` |
| **KQ3** | Nghiên cứu mô phỏng ADEMP (§6.3), với giao thức commit trước khi chạy | Ở các ô công cụ không yếu và công cụ hợp lệ: loại I ≤ 0,064 và độ phủ ∈ [0,936; 0,964] với $n_{\mathrm{sim}}=1000$ (±2 MCSE). Tái lập được bằng một lệnh cục bộ |
| **KQ4** | Ít nhất một ứng dụng trên đoàn hệ thật, phân tích cục bộ, có kế hoạch phân tích định trước | Đáp ứng đủ danh mục STROBE-MR [7]. Có phân tích độ nhạy theo $D$, IPTW, SDPD và MR theo thành phần |
| **KQ5** | Tài liệu người dùng | ≥ 3 vignette xây dựng được trong `R CMD check`. README 7 ngôn ngữ đồng bộ về nội dung. Người dùng mới chạy xong ví dụ nhanh trong < 10 phút |

---

## 6. Kế hoạch thực hiện

Tuần 1 được tính từ đầu tháng 10/2026. Mỗi nhiệm vụ là một đơn vị kiểm thử độc lập (một
commit kèm kiểm thử), theo giao thức ở `acceleration-roadmap.md` §4.

### 6.1 G0 — Đồng bộ hiện trạng và môi trường cục bộ (tuần 1) — ✅ hoàn thành 29/09/2026

**Kết quả G0.**

- **G0.1:** đường cơ sở cục bộ (`tools/baseline/run-baseline.sh`) cho kết quả:
  - testthat: 887 phép kiểm, 0 lỗi;
  - `R CMD check --as-cran`: không có ERROR, không có WARNING do gói gây ra;
  - pytest: 45 đạt, 22 bỏ qua có chủ đích;
  - độ bao phủ 88,30%.

  Biên bản chi tiết ở `baseline-verification.md`.
- **G0.2:** đã sửa các sai lệch ở §3.4. Trong quá trình đó phát hiện thêm và sửa hai lỗi trích
  dẫn (Bebu & Lachin, Even & Josse), và thêm `out/` cùng các README dịch vào `.Rbuildignore`.
- **G0.3:** D3 đã chốt: bỏ GPS, hỗ trợ $\Sigma_{\mathrm{GWAS}}$ đầy đủ.
- **G0.4:** D1 đã chốt là dbGaP. Bộ hồ sơ nháp ở `tools/data-access/dbgap-dar-pack.md`; việc nộp
  do chủ nhiệm đề tài thực hiện.
- **Phát hiện chuyển sang G3:**
  - `print()` in thiếu khoảng trắng;
  - `tidy()` trả biên $e^{\pm50}$ và tên hàng `low` khi khoảng Fieller không bị chặn.

**Nhiệm vụ đã lên kế hoạch.**

- **G0.1** Dựng môi trường R ≥ 4.3 cục bộ; chạy `testthat`, `R CMD check --as-cran` và `pytest`
  để lập đường cơ sở, lưu kèm `sessionInfo()`.
- **G0.2** Sửa các sai lệch ở §3.4: thêm `!vignettes/*.Rmd` vào `.gitignore`, cập nhật
  `CITATION.cff` và phiên bản, sửa trích dẫn Pocock trong 7 README.
- **G0.3** Chốt phạm vi v1.0 (D3). Ghi quyết định vào `project-checkpoint.md`.
- **G0.4** Nộp hồ sơ xin dữ liệu (D1). Đây là việc nằm trên đường găng của G4.
- *Tiêu chí hoàn thành:* đường cơ sở xanh trên máy cục bộ; tài liệu khớp với hiện trạng.

### 6.2 G1 — Khép các vấn đề thống kê mở (tuần 2–5)

- **G1.1 Định danh (I-E1):** phát biểu và chứng minh DS-CWR$(\tau)$ dưới A1–A6, nêu rõ vai trò
  của khác biệt công cụ giữa các tầng. Hệ quả của chứng minh này giải thích vì sao M2 bị loại.
- **G1.2 Căn chỉnh điểm và khoảng:** so sánh ba phương án trên một DGM có gradient không đồng
  nhất (bổ sung vào G2):
  - (a) báo cáo $\hat\delta_R$ làm ước lượng điểm chính, nhất quán với Fieller, và giữ
    $\hat\delta_{\mathrm{GLS}}$ cùng $Q$ làm chẩn đoán đồng nhất;
  - (b) dùng trọng số GLS trong cấu trúc Fieller;
  - (c) giữ nguyên hiện trạng và chứng minh rằng sai lệch là không đáng kể.

  Tiêu chí chọn: độ phủ đối với đại lượng đích tương ứng, cộng với tính diễn giải.
- **G1.3 Ledoit–Wolf:** cài đặt hệ số co rút tối ưu theo Ledoit–Wolf [8], hoặc tài liệu hoá cách
  xấp xỉ hiện tại và thẩm định bằng độ phủ. Kiểm thử: đầu ra xác định dương, số điều kiện được
  cải thiện, đối chiếu với lời giải dạng đóng.
- **G1.4 Hệ số ≈ 2 và $\sigma_\beta>0$:** phân rã phương sai theo từng số hạng (Jacobian của tỉ
  số, co rút, trọng số GLS) để định vị nguồn thừa. Kiểm tra khả năng tính trùng
  $\Sigma_{\mathrm{gwas}}$.
- **G1.5 $\Sigma_{\mathrm{GWAS}}$ đầy đủ:** tuỳ chọn nhận ma trận hiệp phương sai (LD) trong
  `mrwin_gwas()`, với $\beta^\*\sim N(\hat\beta,\Sigma_{\mathrm{GWAS}})$. Nếu D3 loại hạng mục
  này thì ghi vào phần hạn chế.
- **G1.6 AL-CWR, Table 3, GPS:** làm hoặc loại khỏi phạm vi, kèm lý do (D5, D3). Nếu loại GPS
  thì gỡ giá trị `"gps"` khỏi `mrwin_controls()`.
- *Tiêu chí hoàn thành:* mỗi mục có kiểm thử, có dòng sổ WP18 tương ứng, và có ghi chú trong
  `algorithm-spec.md`.

### 6.3 G2 — Nghiên cứu mô phỏng theo khung ADEMP (tuần 4–10)

Giao thức theo Morris, White và Crowther [9]. Giao thức được commit thành
`inst/spec/simulation-protocol.md` **trước khi chạy**; hash của commit đóng vai trò dấu thời
gian định trước. Mọi sai lệch so với giao thức được ghi lại.

**Mục tiêu (A):**

- A1: hiệu chuẩn (sai lầm loại I).
- A2: độ chệch và độ phủ.
- A3: lực kiểm định, so sánh với các phương pháp đối chứng.
- A4: hiệu năng chẩn đoán (SDPD, $Q$, `weak_instrument`, `discordant_components`).
- A5: độ nhạy theo $D$.
- A6: khả năng mở rộng tính toán.

**Cơ chế sinh dữ liệu (D).** Dựa trên `mrwin_config()` và `mrwin_simulate()`. Phơi nhiễm
$X=\alpha_s S+\alpha_u U+\varepsilon$; hazard của thành phần $j$ là
$\alpha_{x,j}X+\nu_{u,j}U+\gamma_j S$, có frailty $\theta_f$; $K=3$.

| Yếu tố | Mức | Tham số |
|---|---|---|
| Cỡ mẫu $N$ | 5 000; 20 000; 100 000 | `n_outcome` |
| Số tầng $D$ | 5; 10; 20 | `n_strata` |
| Độ mạnh công cụ | $R^2_{X\sim S}\approx$ 0,1%; 1%; 5% | `alpha_s`, `m_snps` |
| Hiệu ứng nhân quả | $(0,0,0)$; $(-0{,}2)^{\times3}$; $(-0{,}4)^{\times3}$; bất đồng $(0{,}4;0;-0{,}4)$ | `alpha_x` |
| Pleiotropy trực tiếp ở ưu tiên 1 | $\gamma_1\in\{0;0{,}01;0{,}02;0{,}05\}$ | `gamma_direct` |
| Nhiễu | $\alpha_u\in\{0;0{,}8\}$ | `alpha_u` |
| Kiểm duyệt/theo dõi | Thấp; cao | `censoring_rate`, `max_follow_up` |
| Bất định GWAS | $\sigma_\beta\in\{0;0{,}05\}$ | `sigma_beta` |
| Gradient không đồng nhất | Tuyến tính; phi tuyến | **Mở rộng DGP mới** (phục vụ G1.2 và $Q$) |

Thiết kế gồm hai phần:

- **Khối lõi giai thừa:** $N\times D\times R^2\times\{\text{null},\text{valid}\}$, tức 54 ô.
- **Mở rộng từng-yếu-tố-một:** quanh ô tham chiếu ($N=20\,000$, $D=10$, $R^2=1\%$), cho
  $\gamma$, $\alpha_u$, kiểm duyệt, $\sigma_\beta$, bất đồng, phi tuyến, IPTW và engine suy luận.

Thiết kế giai thừa đầy đủ (> 2 500 ô) không cần thiết cho các mục tiêu A1–A6.

**Estimand (E).** Với giả thuyết không ($\alpha_x=\gamma=0$), kết cục độc lập với $S$, nên
$\delta^\*=0$ chính xác. Với các DGM khác, $\delta^\*_{\mathrm{GLS}}$ và $\delta^\*_R$ được xấp
xỉ bằng Monte Carlo cỡ lớn ($N_{\mathrm{truth}}=10^6$, khả thi nhờ backend nhanh). Số lần lặp
được chọn sao cho MCSE của chân lý nhỏ hơn một phần mười SE thực nghiệm của ô có $N$ nhỏ nhất.
Cách này thay cho quy trình xấp xỉ trước đây ($N=20\,000$, 8 lần lặp).

**Phương pháp (M).**

- M-a: DS-CWR kèm Fieller, bootstrap với $B=200$.
- M-b: DS-CWR kèm Fieller, giải tích.
- M-c: delta-method, chỉ làm tham chiếu.
- M-d: M-a và M-b có IPTW.
- Đối chứng: MR trên thời gian đến biến cố đầu tiên; MR theo từng thành phần
  (`mrwin_per_component_benchmark()`); win ratio quan sát theo tầng phơi nhiễm.
- Bootstrap chỉ chạy trên tập con các ô, vì chi phí tăng tỉ lệ với $B$.

**Thước đo hiệu năng (P).** Độ chệch, SE thực nghiệm, SE mô hình và sai số tương đối của SE mô
hình, độ phủ, sai lầm loại I, lực kiểm định, tỉ lệ Fieller không bị chặn, tỉ lệ bác bỏ của
SDPD, tỉ lệ phát cảnh báo, thời gian chạy và bộ nhớ. MCSE được tính theo [9]:
$\mathrm{MCSE}(\widehat{\text{bias}})=\widehat{\mathrm{SE}}_{\mathrm{emp}}/\sqrt{n_{\mathrm{sim}}}$ và
$\mathrm{MCSE}(\hat p)=\sqrt{\hat p(1-\hat p)/n_{\mathrm{sim}}}$. Với $n_{\mathrm{sim}}=1000$ ở
các ô hiệu chuẩn, MCSE tại $p=0{,}05$ là 0,0069. Với $n_{\mathrm{sim}}=500$ ở các ô mở rộng,
MCSE là 0,0097. Để so sánh, các ô WP19 trước đây chỉ có $M=150$–$200$ (MCSE 0,015–0,018),
chưa đủ độ chính xác cho một khẳng định công bố.

**Tiêu chí chấp nhận.**

- (i) Loại I ≤ 0,064 ở mọi ô không yếu và có công cụ hợp lệ. Đây là tiêu chí chính: không có
  lạm phát.
- (ii) Loại I < 0,036 được xếp là *bảo thủ*, phải báo cáo và giải thích (nối với G1.4).
- (iii) Độ phủ ∈ [0,936; 0,964] tại các ô không yếu, tính với đại lượng đích tương ứng (G1.2).
- (iv) Ô công cụ yếu và ô pleiotropy được báo cáo mô tả, không đặt ngưỡng đạt. Cần tái hiện và
  báo cáo trung thực độ phủ khoảng 12% tại $\gamma_1=0{,}05$ đã ghi trong README.

**Vận hành cục bộ.**

- Chạy thử (pilot) để đo $\hat t(N)$ cho mỗi lần khớp.
- Ngân sách CPU tính theo $C=\sum_{\text{ô}} n_{\mathrm{sim}}\,\hat t(N_{\text{ô}})$ (Phụ lục B).
- Chạy song song đa lõi (`parallel`) với luồng số ngẫu nhiên `L'Ecuyer-CMRG` riêng cho từng ô
  và từng lần lặp.
- Kết quả từng lần lặp được ghi ra đĩa để có thể chạy tiếp; tệp thô bị gitignore.
- Bảng tổng hợp commit tại `inst/spec/simulation-results.md`; mã chạy đặt tại
  `tools/simulation-study/`.

### 6.4 G3 — Kỹ thuật phần mềm và phát hành (tuần 4–9, song song với G2)

- **G3.1** `mrwin_verify_fast_dense_parity()`, mở rộng từ `mrwin_verify_sparse_dense_parity()`:
  so khớp nguyên chính xác khi không có trọng số, sai khác ≤ 1e-10 khi có trọng số, trên các
  ca ties và kiểm duyệt đối kháng. Chạy cục bộ và trong CI.
- **G3.2** Benchmark mở rộng trong `R/benchmark.R` với $N\in\{10^3,\dots,10^5\}$,
  $K\in\{1,3\}$, $D\in\{5,10\}$. Kiểm định cận: số mũ $p<1{,}3$ cho đường nhanh, $p\approx2$ cho
  đường tham chiếu.
- **G3.3** Đổi mặc định sang `backend = "fast"`, giữ `"dense"` để gỡ lỗi. Mặc định của
  `inference` được quyết định theo kết quả G2 (D4).
- **G3.4** Viết 3 vignette: mô phỏng nhanh, mẫu cho đoàn hệ thật, diễn giải chẩn đoán. Có thể
  thêm vignette thứ tư về lựa chọn backend và engine suy luận.
- **G3.5** Chuẩn bị phát hành: `NEWS.md`, `inst/CITATION`, phiên bản 1.0.0, `R CMD check
  --as-cran` trên máy cục bộ, `covr` ≥ 80%, cập nhật `release-checklist.md`, rồi nộp CRAN.
- *Tiêu chí hoàn thành:* KQ2.

### 6.5 G4 — Ứng dụng trên dữ liệu thật (song song; phụ thuộc thời điểm cấp dữ liệu)

- **Thiết kế:** MR một mẫu, dữ liệu cá thể. Ví dụ chính: LDL-C → hồ sơ tim–thận, với thứ bậc
  tử vong > nhập viện vì suy tim > suy giảm eGFR, định nghĩa ngưỡng trước.
- **Trọng số PRS:** lấy từ GWAS ngoài, không chồng lấn mẫu với đoàn hệ phân tích. Hiệu chỉnh
  cho tổ tiên bằng các PC.
- **Tính lực kiểm định trước khi xin dữ liệu:** dùng engine mô phỏng hiệu chỉnh theo tỉ lệ biến
  cố, thời gian theo dõi và $R^2$ của PRS trong đoàn hệ đích. README nêu rằng quy mô biobank
  ($N>10^5$) là điều kiện nghiêm ngặt cho suy luận đáng tin cậy. Vì vậy với các đoàn hệ cỡ
  $10^4$, cần đánh giá khả thi một cách định lượng.
- **Phương án gộp đa đoàn hệ** (câu hỏi nghiên cứu mới): gộp ở mức tử số và mẫu số
  $(\bar u_1,\bar u_2)$ rồi dựng khoảng Fieller, thay vì gộp các tỉ số, để tránh chệch của ước
  lượng tỉ số dưới công cụ yếu [4].
- **Kế hoạch phân tích định trước:** commit trong kho *trước khi truy cập dữ liệu*. Các phân tích
  độ nhạy gồm $D\in\{5,10,20\}$, IPTW, SDPD kèm CI giới hạn pleiotropy, và MR theo thành phần.
- **Báo cáo:** theo STROBE-MR [7].
- **Quản trị dữ liệu:** dữ liệu chỉ nằm trên máy chủ bảo mật của cơ sở. `.gitignore` đã chặn
  các định dạng dữ liệu. Chỉ commit mã và kết quả tổng hợp.

### 6.6 G5 — Bản thảo và công bố (tuần 8–16)

- Cấu trúc bản thảo theo `wp18-theory-foundation.md` §1. Bổ sung mục "Định danh" (I-E1) và mục
  "Căn chỉnh điểm và khoảng" (G1.2).
- Phần mô phỏng báo cáo theo ADEMP [9]. Phần ứng dụng báo cáo theo STROBE-MR [7]. Các phát hiện
  phủ định (C4) được trình bày như đóng góp phương pháp luận.
- Nghĩa vụ WP18 T4: đối chiếu từng con số trong bản thảo với tệp kiểm thử hoặc benchmark tương
  ứng.
- Đăng arXiv v2 → nộp tạp chí (D2) → phản hồi phản biện.
- Có thể viết thêm bài thứ hai dạng *Software Application Profile* sau khi gói lên CRAN.

### 6.7 G6 — Phát hành 1.0 và duy trì

- Lưu trữ bản phát hành kèm DOI. Đồng bộ README 7 ngôn ngữ. Lập quy trình tiếp nhận báo lỗi.
- Chương trình cho bài báo thứ hai:
  - chế độ dữ liệu tóm tắt;
  - $K\ge4$ với đếm trội $d$ chiều;
  - ước lượng LACE *trong* tầng, tận dụng `mrwin_doubly_ranked_strata()`;
  - GPS;
  - WP14, chỉ khi cần $B$ lớn mà không dùng được engine giải tích.

### 6.8 Ràng buộc "không điện toán đám mây": hệ quả thiết kế

- **Tính toán:** mọi tính toán (G1–G4) chạy trên máy trạm hoặc máy chủ của cơ sở. Cấu hình tham
  chiếu: 16–32 lõi, ≥ 64 GB RAM. Một ma trận kiểu gen $10^6\times100$ kiểu `double` chiếm khoảng
  0,8 GB.
- **Kiểm định phát hành:** chạy cục bộ là chuẩn chính thức. CI trên GitHub chỉ là lớp bổ sung,
  có thể tắt.
- **Dữ liệu:** theo chính sách truy cập hiện hành của UK Biobank (cần xác nhận lại khi nộp
  đơn), dữ liệu cá thể được phân tích trên nền tảng đám mây của họ. Dưới ràng buộc này, các
  lựa chọn ưu tiên là:
  - đoàn hệ trên dbGaP (ARIC, CHS, FHS, MESA), phân tích trên máy chủ bảo mật theo cam kết sử
    dụng dữ liệu;
  - Rhineland Study, qua đồng tác giả tại DZNE.

  Hệ quả là cỡ mẫu nhỏ hơn. Điều này làm tăng tầm quan trọng của bước tính lực kiểm định và của
  phương án gộp đa đoàn hệ (§6.5).

---

## 7. Tiến độ và đường găng

| Giai đoạn | T1 | T2–3 | T4–5 | T6–7 | T8–9 | T10–12 | T13–16 | > T16 |
|---|---|---|---|---|---|---|---|---|
| G0 | ■ | | | | | | | |
| G1 | | ■ | ■ | | | | | |
| G2 | | | ■ giao thức | ■ | ■ | ■ | | |
| G3 | | | ■ | ■ | ■ | ■ CRAN | | |
| G4 | hồ sơ | (chờ) | (chờ) | ■ SAP | ■ | ■ | ■ | |
| G5 | | | | | ■ | ■ | ■ nộp | phản biện |
| G6 | | | | | | | | ■ |

**Các mốc:** M1 (T1) có đường cơ sở. M2 (T5) G1 khép lại. M3 (T5) giao thức ADEMP được commit.
M4 (T10) có kết quả mô phỏng. M5 (T12) gói lên CRAN. M6 (T16) nộp bản thảo.

**Đường găng:** G0 → G1.1/G1.2 → G2 → G5.

**Dự phòng:** 2–3 tuần cho một phát hiện ngoại tại mới. Kinh nghiệm của dự án cho thấy điều này
đã xảy ra ba lần (R2, M1, M2).

Nếu G4 trễ, bản thảo được nộp với phần mô phỏng, còn ứng dụng được bổ sung ở vòng sửa. Phương án
này cần được kiểm tra lại theo chính sách của tạp chí đích.

---

## 8. Điểm quyết định

| Mã | Câu hỏi | Phương án | Tiêu chí | Khuyến nghị | Hạn chót |
|---|---|---|---|---|---|
| D1 | Đoàn hệ cho KQ4; có chấp nhận môi trường nghiên cứu tin cậy trên đám mây không | dbGaP / Rhineland / UK Biobank-RAP | Ràng buộc đám mây, cỡ mẫu, thời gian cấp | ✅ **Đã chốt (29/09/2026):** dbGaP (ARIC trước; MESA/CHS/FHS cho phân tích gộp), máy chủ cơ sở, không đám mây | T1 |
| D2 | Tạp chí đích | IJE / *Stat Med* / *Genet Epidemiol* / *Biostatistics* | Độc giả, định dạng, thời gian xử lý | IJE (định hướng sẵn trong `run_all.sh`) | T8 |
| D3 | Phạm vi v1.0 | Có/không GPS; có/không $\Sigma_{\mathrm{GWAS}}$ đầy đủ | Rủi ro phản biện so với chi phí | ✅ **Đã chốt (29/09/2026):** bỏ GPS (G1.6); làm $\Sigma_{\mathrm{GWAS}}$ đầy đủ (G1.5) | T1 |
| D4 | Engine suy luận mặc định | Bootstrap / giải tích | Kết quả hiệu chuẩn ở G2 | Quyết định sau M4 | T10 |
| D5 | AL-CWR và Table 3 trong bài 1 | Có / không | Có hoàn tất và thẩm định kịp trong G1 không | Chỉ đưa vào nếu đạt trước M2 | T5 |

---

## 9. Quản lý rủi ro

| Rủi ro | Khả năng | Tác động | Giảm thiểu | Dấu hiệu cảnh báo sớm |
|---|---|---|---|---|
| Cấp dữ liệu chậm | Cao | Trung bình | Nộp hồ sơ ở T1; phương án nộp bài trước khi có ứng dụng | Chưa có phản hồi sau 6 tuần |
| Đoàn hệ thiếu lực kiểm định (công cụ yếu, Fieller không bị chặn) | Trung bình–cao | Cao | Tính lực kiểm định trước; gộp đa đoàn hệ | $R^2$ PRS < 1% hoặc ít biến cố ưu tiên 1 |
| Phát hiện ngoại tại mới (tương tự R2/M1/M2) | Trung bình | Trung bình–cao | Dự phòng thời gian; coi là kết quả có giá trị | Loại I > 0,064 ở bất kỳ ô lõi nào |
| Căn chỉnh điểm/khoảng (G1.2) buộc đổi đầu ra chính | Trung bình | Trung bình | Chốt trước khi đổi mặc định và viết bản thảo | Độ phủ lệch giữa $\delta_R$ và $\delta_{\mathrm{GLS}}$ trong DGM phi tuyến |
| Nhạy cảm với pleiotropy | Đã biết | Cao | SDPD kèm CI giới hạn pleiotropy; nêu rõ trong phần Hạn chế | — |
| Phản biện về estimand (không khả gộp, phụ thuộc $\tau$) | Trung bình | Trung bình | Định lý I-E1; định nghĩa DS-CWR$(\tau)$; ví dụ diễn giải | — |
| Thiếu tài nguyên tính toán cục bộ | Thấp | Trung bình | Chạy pilot và lập ngân sách trước; ưu tiên khối lõi | Ngân sách ước tính vượt khả năng máy |

---

## 10. Nguyên tắc toàn vẹn khoa học và tái lập

1. **Không có phát biểu nào đi trước kiểm thử** (`acceleration-roadmap.md` §4.5.6).
2. **Thẩm định ngoại tại là tiêu chuẩn cuối cùng.** Nhất quán nội tại (fast = dense,
   giải tích = bootstrap) chỉ chứng minh rằng hai phép tính trùng nhau, không chứng minh con số
   đúng.
3. **Định trước:** giao thức mô phỏng và kế hoạch phân tích dữ liệu được commit trước khi chạy
   hoặc trước khi truy cập dữ liệu.
4. **Minh bạch phát hiện phủ định:** M1, M2 và delta-method được báo cáo đầy đủ.
5. **Tái lập:** kiểm soát seed; ghi `sessionInfo()`; mỗi bảng trong bài báo tái tạo được bằng
   một lệnh cục bộ.
6. **Quản trị dữ liệu và ràng buộc cục bộ:** không có dữ liệu cá thể nào trong kho mã; không
   dùng hạ tầng đám mây cho tính toán hay lưu trữ dữ liệu.

---

## 11. Tài liệu tham khảo

1. Pocock SJ, Ariti CA, Collier TJ, Wang D. The win ratio: a new approach to the analysis of
   composite endpoints in clinical trials based on clinical priorities. *Eur Heart J*.
   2012;33(2):176–182. doi:10.1093/eurheartj/ehr352
2. Bebu I, Lachin JM. Large sample inference for a win ratio analysis of a composite outcome
   based on prioritized components. *Biostatistics*. 2016;17(1):178–187.
   doi:10.1093/biostatistics/kxv032
3. Dong G, Huang B, Chang YW, Seifu Y, Song J, Hoaglin DC. The win ratio: impact of censoring and
   follow-up time and use with nonproportional hazards. *Pharm Stat*. 2020;19(3):168–177.
   doi:10.1002/pst.1977
4. Burgess S, Thompson SG. Avoiding bias from weak instruments in Mendelian randomization
   studies. *Int J Epidemiol*. 2011;40(3):755–764. doi:10.1093/ije/dyr036
5. Bowden J, Davey Smith G, Burgess S. Mendelian randomization with invalid instruments: effect
   estimation and bias detection through Egger regression. *Int J Epidemiol*.
   2015;44(2):512–525. doi:10.1093/ije/dyv080
6. Tian H, Mason AM, Liu C, Burgess S. Relaxing parametric assumptions for non-linear Mendelian
   randomization using a doubly-ranked stratification method. *PLoS Genet*.
   2023;19(6):e1010823. doi:10.1371/journal.pgen.1010823
7. Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the reporting of observational
   studies in epidemiology using Mendelian randomization: the STROBE-MR statement. *JAMA*.
   2021;326(16):1614–1621. doi:10.1001/jama.2021.18236
8. Ledoit O, Wolf M. A well-conditioned estimator for large-dimensional covariance matrices.
   *J Multivar Anal*. 2004;88(2):365–411.
9. Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods.
   *Stat Med*. 2019;38(11):2074–2102. doi:10.1002/sim.8086
10. Fieller EC. Some problems in interval estimation. *J R Stat Soc Series B*.
    1954;16(2):175–185.

---

## Phụ lục A. Ánh xạ giai đoạn ↔ gói công việc ↔ tệp

| Giai đoạn | WP / dòng sổ | Tệp chính |
|---|---|---|
| G0 | WP12 (tái kiểm), `release-checklist.md` | `DESCRIPTION`, `CITATION.cff`, `.gitignore`, `README*.md` |
| G1 | WP18 (I-E1 mới, C-R2), `algorithm-spec.md` §4.6–4.9 | `R/estimate.R`, `R/bootstrap.R`, `R/analytic_variance.R`, `R/api.R` |
| G2 | WP19 T3 | `R/simulate.R`, `R/simulation_engine.R`, `R/validate_calibration.R`, `tools/validation-scripts/` |
| G3 | WP19 T1, T2, T4, T5 | `R/backend_sparse.R`, `R/kernel_fast.R`, `src/fast_kernel.cpp`, `R/benchmark.R`, `.github/workflows/ci.yml` |
| G4 | WP18 §1 mục 7 (Application) | Mã phân tích (không kèm dữ liệu) |
| G5 | WP18 T2–T4 | `manuscript/` (gitignored), `wp18-theory-foundation.md` |
| G6 | — | `NEWS.md`, `inst/CITATION`, `README*.md` |

## Phụ lục B. Ước tính ngân sách tính toán cục bộ

Gọi $\hat t(N)$ là thời gian của một lần khớp đo trong pilot. Ngân sách được tính như sau:

$$
C_{\text{CPU-giờ}}=\frac{1}{3600}\sum_{\text{ô }c} n_{\mathrm{sim},c}\,\hat t(N_c)\,\kappa_c,
\qquad T_{\text{thực}}\approx C/n_{\text{lõi}},
$$

trong đó $\kappa_c$ là hệ số engine so với engine giải tích ở $\sigma_\beta=0$ ($\kappa=1$) và
được đo trực tiếp trong pilot. Giá trị tham chiếu từ `wp17-analytic-variance.md`: bootstrap chậm
hơn khoảng 7 lần; engine giải tích với $\sigma_\beta>0$ chỉ nhanh hơn bootstrap khoảng 1,5 lần,
vì vòng tái chọn mẫu GWAS vẫn còn. Giá trị $B$ dùng khi đo 7× chưa được ghi lại, nên $\kappa$ của
bootstrap phải được đo lại theo $B$ thực tế. Ví dụ minh hoạ (không phải số đo): nếu
khối lõi có 18 000 lần khớp ở $N=10^5$ với $\hat t\approx 20$ s trên engine giải tích, thì
$C\approx100$ CPU-giờ, tức khoảng 3–4 giờ trên 32 lõi. Con số chính thức được thay bằng kết quả
đo ở pilot (G2) trước khi khoá giao thức.
