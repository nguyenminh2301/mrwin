# Bộ hồ sơ xin dữ liệu dbGaP cho ứng dụng KQ4 (G0.4)

Soạn ngày: 29/09/2026. Quyết định liên quan: **D1 = dbGaP, phân tích trên máy chủ của cơ sở, không
dùng đám mây** (ghi trong `inst/spec/project-checkpoint.md`).

Tài liệu này chứa (i) danh mục việc chỉ chủ nhiệm đề tài và cơ sở làm được, (ii) bản nháp tiếng Anh
cho các ô của đơn xin dữ liệu (Data Access Request, DAR), và (iii) kế hoạch bảo mật dữ liệu cục bộ.
Mọi chi tiết đánh dấu **[xác minh]** phải được đối chiếu với trang dbGaP/NIH tại thời điểm nộp, vì quy
định và số phiên bản thay đổi theo thời gian.

Tệp nằm trong `tools/` nên không được đóng gói vào bản phát hành R.

---

## 1. Việc chỉ con người làm được (danh mục cho chủ nhiệm đề tài)

| # | Việc | Người phụ trách | Ghi chú |
|---|---|---|---|
| 1 | Bảo đảm cơ sở (UMP) có đăng ký trong eRA Commons và có Signing Official (SO) | Phòng Quản lý khoa học của cơ sở | Cơ sở ngoài Hoa Kỳ vẫn nộp được, nhưng phải có SO và tài khoản tổ chức **[xác minh]** |
| 2 | Chủ nhiệm đề tài có tài khoản eRA Commons vai trò PI | PI | Đăng nhập dbGaP bằng tài khoản này |
| 3 | Chỉ định IT Director của cơ sở | Cơ sở | Người này ký xác nhận kế hoạch bảo mật (mục 4) |
| 4 | Xin phê duyệt của Hội đồng đạo đức (IRB) của cơ sở | PI | Bắt buộc với các nhóm đồng thuận có cờ "IRB"; nên xin cho mọi nhóm |
| 5 | Chọn đúng nhóm đồng thuận (consent group) của từng nghiên cứu | PI | Xem mục 2.2; nghiên cứu phải khớp giới hạn sử dụng dữ liệu (DUL) |
| 6 | Cộng sự ở cơ sở khác (ví dụ DZNE) | PI và cộng sự | Cộng sự ngoài cơ sở của PI thường phải nộp DAR riêng; không được chuyển dữ liệu cá thể sang cơ sở khác nếu chưa được phép **[xác minh]** |
| 7 | Xác nhận tuân thủ yêu cầu bảo mật hiện hành của NIH | PI, IT Director | Từ 2025 NIH yêu cầu cơ sở nhận dữ liệu truy cập có kiểm soát tuân thủ NIST SP 800-171 **[xác minh]** |
| 8 | Nộp DAR, theo dõi, gia hạn hằng năm, nộp báo cáo đóng dự án | PI | Phê duyệt thường có hiệu lực một năm **[xác minh]** |

---

## 2. Nghiên cứu cần xin

### 2.1 Danh sách ưu tiên

| Ưu tiên | Nghiên cứu | Mã dbGaP **[xác minh]** | Lý do |
|---|---|---|---|
| 1 | ARIC (Atherosclerosis Risk in Communities) | phs000280 | Có đủ ba tầng kết cục: tử vong, nhập viện vì suy tim đã thẩm định, eGFR lặp lại |
| 2 | MESA (Multi-Ethnic Study of Atherosclerosis) | phs000209 | Đa sắc tộc; dùng cho phân tích độ nhạy theo tổ tiên |
| 3 | CHS (Cardiovascular Health Study) | phs000287 | Người cao tuổi, nhiều biến cố |
| 4 | FHS (Framingham Heart Study) | phs000007 | Theo dõi dài; cấu trúc gia đình cần xử lý riêng |

Khuyến nghị: xin **ARIC trước** để triển khai phân tích chính. Các đoàn hệ khác dùng cho phương án gộp
đa đoàn hệ (roadmap học thuật §6.5), vì một đoàn hệ cỡ khoảng 10^4 có thể thiếu lực kiểm định.

### 2.2 Nhóm đồng thuận và giới hạn sử dụng dữ liệu

- Ưu tiên các nhóm **HMB** (Health/Medical/Biomedical), vì phù hợp với câu hỏi tim–thận.
- Với nhóm **DS-CVD** (chỉ dùng cho bệnh tim mạch), cần lập luận rõ rằng kết cục thận được dùng như
  một thành phần của hồ sơ **tim–thận**. Nếu giới hạn không cho phép, loại thành phần thận khỏi phân
  tích trên nhóm đó **[xác minh DUL từng nhóm]**.
- Nhóm có cờ **NPU** (không dùng cho mục đích thương mại) không ảnh hưởng tới dự án học thuật này.

### 2.3 Biến cần có

| Nhóm | Biến |
|---|---|
| Di truyền | Kiểu gen đã imputation (các SNP của PRS LDL-C), thành phần chính (PC) tổ tiên, QC |
| Phơi nhiễm | LDL-C lúc ban đầu; tình trạng dùng thuốc hạ lipid (cho phân tích độ nhạy) |
| Hiệp biến | Tuổi, giới, trung tâm nghiên cứu, PC1–PC10 |
| Kết cục ưu tiên 1 | Tử vong do mọi nguyên nhân: có/không và ngày |
| Kết cục ưu tiên 2 | Nhập viện vì suy tim đã thẩm định: có/không và ngày |
| Kết cục ưu tiên 3 | eGFR các lần khám hoặc sự kiện thận (suy giảm eGFR ≥ 40% hoặc suy thận giai đoạn cuối): có/không và ngày |
| Theo dõi | Ngày vào nghiên cứu, ngày kết thúc theo dõi, lý do kiểm duyệt |

---

## 3. Bản nháp cho các ô của DAR (tiếng Anh)

### 3.1 Project title

> Causal win statistics for hierarchical cardio-renal outcomes using polygenic instruments

### 3.2 Research Use Statement

> **Objectives.** We will estimate the causal effect of genetically predicted low-density lipoprotein
> cholesterol (LDL-C) on a clinically prioritised cardio-renal outcome profile (all-cause death >
> hospitalisation for heart failure > kidney function decline) using the dose-standardised causal
> win ratio (DS-CWR), a Mendelian randomization estimand for hierarchical composite endpoints
> implemented in the open-source R package `mrwin`. Conventional Mendelian randomization analyses of
> non-fatal outcomes condition on survival, which can induce survivor bias. DS-CWR places death at
> the top of the hierarchy, so differential survival is absorbed into the estimand rather than
> conditioned away.
>
> **Study design.** One-sample Mendelian randomization with individual-level data. A polygenic risk
> score for LDL-C is constructed from published genome-wide association summary statistics that do
> not include the requested cohort, to avoid sample overlap. Participants are ranked by the score
> into ordered strata. Within each pair of adjacent strata we compute a pairwise hierarchical win
> ratio, standardise it by the between-stratum difference in LDL-C, and pool the standardised
> gradients by generalised least squares. Inference uses Fieller confidence intervals with
> multiplier-bootstrap or analytic influence-function variance that propagates GWAS-weight
> uncertainty.
>
> **Analysis plan.** The primary analysis, the outcome hierarchy, the definitions of each component
> and the sensitivity analyses (number of strata, inverse-probability weighting for measured
> covariates, pleiotropy diagnostics, component-wise Mendelian randomization and conventional
> time-to-first-event Mendelian randomization) will be pre-specified and time-stamped before data
> access. Before requesting individual-level data we will run a simulation-based power analysis
> calibrated to published cohort characteristics.
>
> **Why these data are required.** The method requires individual-level genotype, exposure and
> time-to-event data for all three prioritised outcomes in the same participants; summary
> statistics cannot support pairwise hierarchical comparisons. The requested cohort combines
> genotyping, measured LDL-C, adjudicated heart-failure hospitalisations, repeated kidney-function
> measures and mortality follow-up.
>
> **Dissemination.** Results will be reported in aggregate only, following the STROBE-MR guideline,
> in a peer-reviewed methods article and in the documentation of the `mrwin` package. No
> individual-level data will be shared or published.

### 3.3 Non-technical summary

> High cholesterol is known to raise the risk of heart attacks, but patients often experience
> several problems of different severity, such as death, hospital stays for heart failure and
> loss of kidney function. Most studies count only the first event, so a minor event can weigh as
> much as a death. We have developed a statistical method that compares people pair by pair,
> looking first at death, then at hospitalisation, then at kidney function, and uses inherited
> genetic differences as a natural experiment to estimate cause and effect. We will apply this
> method to the requested study to learn whether genetically higher cholesterol worsens the whole
> course of heart and kidney disease, and to show other researchers how to use the method.

### 3.4 Cloud computing

> No cloud computing will be used. All analyses will be performed on a dedicated, access-controlled
> server within the investigator's institution. Controlled-access data will not be transferred to
> any cloud service, personal device or third-party system.

### 3.5 Collaborators

> [Liệt kê cộng sự tại cơ sở của PI. Cộng sự tại cơ sở khác, ví dụ N. Ahmad Aziz (DZNE, Đức), cần
> được khai báo theo hướng dẫn hiện hành của dbGaP; nếu họ cần truy cập dữ liệu cá thể thì nộp DAR
> riêng **[xác minh]**.]

---

## 4. Kế hoạch bảo mật dữ liệu cục bộ (nộp kèm và để IT Director ký xác nhận)

1. **Nơi lưu trữ:** một máy chủ riêng trong mạng nội bộ của cơ sở, ổ đĩa mã hoá toàn phần, không
   mở ra Internet công cộng. Không lưu trên máy tính xách tay, ổ USB hay dịch vụ đám mây.
2. **Kiểm soát truy cập:** tài khoản cá nhân cho từng người được phê duyệt trong DAR; xác thực mạnh;
   nhật ký truy cập được lưu và rà soát định kỳ.
3. **Tách mã và dữ liệu:** kho mã `mrwin` chỉ chứa mã và kết quả tổng hợp. `.gitignore` đã chặn
   các định dạng dữ liệu (`data/`, `*.dta`, `*.sav`, `*.sas7bdat`, `*.parquet`, `*.feather`). Không
   commit bất kỳ tệp nào sinh ra từ dữ liệu cá thể.
4. **Kết quả xuất ra:** chỉ xuất số liệu tổng hợp (ước lượng, khoảng tin cậy, bảng đặc điểm). Không
   xuất ô có ít hơn ngưỡng tối thiểu người tham gia mà nghiên cứu gốc yêu cầu **[xác minh ngưỡng]**.
5. **Sao lưu:** chỉ sao lưu trong hạ tầng nội bộ, cũng được mã hoá.
6. **Kết thúc dự án:** xoá dữ liệu và các bản sao khi DAR hết hạn hoặc đóng dự án, và báo cáo việc
   xoá theo quy định.
7. **Sự cố:** báo cáo sự cố bảo mật cho dbGaP và cơ sở trong thời hạn quy định **[xác minh]**.
8. **Tuân thủ:** xác nhận tuân thủ yêu cầu bảo mật hiện hành của NIH cho dữ liệu truy cập có kiểm
   soát, được IT Director ký **[xác minh phiên bản yêu cầu]**.

---

## 5. Việc làm ngay sau khi nộp (không cần chờ dữ liệu)

- Viết kế hoạch phân tích định trước và commit vào kho **trước khi truy cập dữ liệu** (G4). Hash
  của commit là dấu thời gian.
- Chọn thống kê tóm tắt GWAS LDL-C **không chồng lấn mẫu** với đoàn hệ được xin. Các consortium lớn
  có thể đã bao gồm chính các đoàn hệ này, nên phải kiểm tra danh sách đoàn hệ thành phần.
- Chạy phân tích lực kiểm định bằng `mrwin_simulate()`, hiệu chỉnh theo cỡ mẫu, tỉ lệ biến cố và
  $R^2$ của PRS công bố cho đoàn hệ đích. Nếu lực không đủ, chuyển sang phương án gộp đa đoàn hệ.
