# Lộ trình hoàn thiện dự án `mrwin` — Bản dễ hiểu

Cập nhật: 29/09/2026. Người đọc: cộng sự, nghiên cứu sinh, bác sĩ lâm sàng, nhà tài trợ,
tức là những người cần hiểu dự án mà không phải đọc công thức.

Bản học thuật có cùng nội dung và cùng mã số mốc: `completion-roadmap-academic.vi.md`.

> **Trạng thái (29/09/2026):** ✅ **G0 đã hoàn thành.** Mọi kiểm tra chạy xanh trên máy cục bộ:
> 887 phép kiểm, 0 lỗi, độ bao phủ 88,30%. Tài liệu đã được sửa cho khớp thực tế. Nhóm đã chốt
> D1 (dbGaP, máy chủ cơ sở) và D3 (bỏ GPS; hỗ trợ tương quan gen LD). Chi tiết ở
> `baseline-verification.md` và `project-checkpoint.md`. **Việc tiếp theo: G1.**
Hai bản dùng chung các mã: **KQ1–KQ5** (kết quả cuối cùng), **G0–G6** (giai đoạn),
**D1–D5** (quyết định cần chốt).

---

## Tóm tắt trong một phút

- **Dự án làm gì?** `mrwin` là một phần mềm thống kê (gói R). Nó trả lời câu hỏi: *một yếu
  tố nguy cơ, ví dụ cholesterol LDL, có thật sự **gây ra** diễn tiến bệnh xấu hơn hay không,
  khi ta xét cùng lúc tử vong, nhập viện và suy giảm chức năng cơ quan theo đúng thứ tự
  nặng–nhẹ?*
- **Đã làm được gì?** Phần lõi đã chạy đúng. Thời gian tính cho dữ liệu hàng trăm nghìn người
  đã giảm từ mức không khả thi xuống dưới một giây. Nhóm cũng đã tự phát hiện và sửa một lỗi
  quan trọng về khoảng tin cậy.
- **Còn thiếu gì?** Một đợt kiểm nghiệm đầy đủ bằng mô phỏng, một ví dụ trên dữ liệu bệnh nhân
  thật, bài báo được bình duyệt, và bản phát hành chính thức.
- **Đích đến:** 5 kết quả cụ thể (KQ1–KQ5, mục 3). Nếu dữ liệu thật được cấp kịp thì cần khoảng
  **4 tháng** làm việc tập trung để nộp bài báo.
- **Ràng buộc:** không dùng điện toán đám mây. Mọi tính toán và mọi dữ liệu nằm trên máy của
  nhóm hoặc máy chủ của cơ sở (xem mục 8).

---

## 1. Dự án này để làm gì?

### Vấn đề

Bệnh nhân tim–thận có thể gặp nhiều biến cố: tử vong, nhập viện vì suy tim, suy giảm chức năng
thận. Các biến cố này **không nặng như nhau**. Nhiều phân tích hiện nay lại chỉ đếm "biến cố
đầu tiên", nên một lần xét nghiệm bất thường có thể được tính nặng như một ca tử vong. Thứ tự
quan trọng nhất với người bệnh và thầy thuốc bị mất đi.

### Ý tưởng 1: "giải đấu từng cặp" (win statistics)

Hãy tưởng tượng ta ghép từng cặp người và hỏi: *ai có diễn tiến tốt hơn?*

1. So sánh tử vong trước. Ai còn sống lâu hơn thì thắng.
2. Nếu hoà, so sánh nhập viện.
3. Nếu vẫn hoà, so sánh chức năng thận.

Đếm số trận thắng và số trận thua của một nhóm. **Tỉ số thắng (win ratio)** bằng số thắng chia
số thua. Cách này đã được dùng trong nhiều thử nghiệm lâm sàng ngẫu nhiên.

### Ý tưởng 2: "xổ số di truyền" (Mendelian randomization)

Người có LDL cao thường khác người có LDL thấp ở nhiều mặt khác: ăn uống, vận động, thu nhập,
thuốc men. Vì thế khó biết LDL có thật sự là **nguyên nhân** hay không. Tuy nhiên, mỗi người
nhận gen từ cha mẹ một cách ngẫu nhiên, giống như rút thăm. Người "rút" được nhiều biến thể làm
tăng LDL sẽ có LDL cao hơn **suốt đời**, và điều đó không do lối sống quyết định. Nhờ vậy, di
truyền trở thành một "thí nghiệm tự nhiên".

### Kết hợp hai ý tưởng

`mrwin` xếp mọi người theo "điểm di truyền" (điểm nguy cơ đa gen), chia thành các nhóm liền
kề, rồi tổ chức giải đấu giữa hai nhóm cạnh nhau. Kết quả được quy về **một con số, DS-CWR**:

- **DS-CWR < 1**: yếu tố nguy cơ làm diễn tiến bệnh **xấu đi**.
- **DS-CWR > 1**: yếu tố đó làm diễn tiến **tốt lên**.
- **DS-CWR = 1**: không có bằng chứng về tác động.

### Vì sao điều này quan trọng?

Ví dụ về sa sút trí tuệ: người mang một yếu tố có hại có thể **tử vong trước** khi kịp được
chẩn đoán sa sút trí tuệ. Phân tích thông thường khi đó có thể thấy yếu tố này như đang "bảo
vệ", vì người mắc đã không sống đủ lâu để được chẩn đoán. `mrwin` đặt tử vong lên đầu bảng
xếp hạng, nên cái chết không bị bỏ qua. Nhờ đó tránh được kiểu sai lệch "chỉ nhìn người còn
sống" (survivor bias).

---

## 2. Mục tiêu tối thượng

> **Biến `mrwin` thành một phương pháp nghiên cứu nhân quả đã được kiểm chứng, đã được giới
> khoa học thẩm định, dễ sử dụng và được người khác thực sự sử dụng.**

Mục tiêu này đứng trên ba trụ cột và hướng tới một nhóm người hưởng lợi:

| Trụ cột | Nghĩa là |
|---|---|
| **Khoa học** | Có bài báo trên tạp chí bình duyệt giải thích phương pháp và chứng minh nó đúng |
| **Phần mềm** | Bất kỳ ai cũng cài được bằng một lệnh, kết quả đáng tin |
| **Bằng chứng** | Có mô phỏng đầy đủ và ít nhất một ví dụ trên dữ liệu bệnh nhân thật |
| **Người dùng** | Nhà dịch tễ học, bác sĩ, nghiên cứu sinh tự dùng được mà không cần hỏi tác giả |

---

## 3. Kết quả kỳ vọng cuối cùng

Đây là những thứ **cầm nắm được** khi dự án xong.

| Mã | Kết quả | Nói đơn giản | Biết là xong khi nào |
|---|---|---|---|
| **KQ1** | Bài báo khoa học | Bài báo trên tạp chí có phản biện, ví dụ *International Journal of Epidemiology* | Có thư chấp nhận đăng; bản preprint mới trên arXiv thay bản "chỉ lý thuyết" |
| **KQ2** | Phần mềm phát hành chính thức | `mrwin` phiên bản 1.0.0 trên CRAN, kho phần mềm chính thức của R | Cài được bằng `install.packages("mrwin")`, vượt mọi bước kiểm tra tự động |
| **KQ3** | Bằng chứng mô phỏng | Chứng minh bằng hàng nghìn bộ dữ liệu giả lập rằng phương pháp không "báo động giả" quá mức | Tỉ lệ báo động giả gần 5%, khoảng tin cậy chứa giá trị thật khoảng 95% số lần; chạy lại được bằng một lệnh |
| **KQ4** | Ví dụ trên dữ liệu thật | Áp dụng cho một đoàn hệ bệnh nhân thật, ví dụ LDL và bệnh tim–thận | Kết quả được báo cáo theo chuẩn quốc tế STROBE-MR, phân tích trên máy của cơ sở |
| **KQ5** | Hướng dẫn sử dụng | 3 bài hướng dẫn (vignette), README 7 ngôn ngữ thống nhất | Người mới chạy được ví dụ đầu tiên trong vòng 10 phút |

Quy tắc chung: **mọi con số trong bài báo đều truy được về một phép kiểm tra đã lưu trong kho
mã**. Không có con số nào "nói miệng".

---

## 4. Chúng ta đang ở đâu?

Hãy hình dung dự án như **xây một toà nhà**.

| Hạng mục | Tình trạng | Ẩn dụ |
|---|---|---|
| Phần lõi: tính toán, kiểm tra dữ liệu, hiệu chỉnh yếu tố nhiễu, chẩn đoán, báo cáo (giai đoạn I) | ✅ Xong | Móng và khung nhà |
| Tốc độ: 72–300 lần nhanh hơn; 80.000 người chỉ mất 0,78 giây | ✅ Xong | Lắp thang máy |
| Cách tính độ bất định nhanh, không cần lặp lại hàng trăm lần | ✅ Xong | Hệ thống điện |
| Phát hiện và sửa lỗi khoảng tin cậy (khoảng cũ rộng gấp khoảng 2 lần, dùng khoảng Fieller thay thế) | ✅ Xong | Kiểm định an toàn tìm ra một lỗi và đã sửa |
| Kiểm nghiệm mô phỏng đầy đủ | ⏳ Mới xong một phần | Nghiệm thu toàn diện |
| Ví dụ dữ liệu thật | ❌ Chưa có | Có người ở thử |
| Bài báo được bình duyệt | ⏳ Mới có bản nháp và preprint lý thuyết | Hội đồng thẩm định |
| Phát hành chính thức | ❌ Chưa | Bàn giao chìa khoá |

### Một bài học quý đã rút ra

Trong quá trình làm, nhóm gặp **ba lần** tình huống: các phép kiểm tra "nội bộ" đều đạt (hai
cách tính khác nhau ra cùng một con số), nhưng khi kiểm tra "bên ngoài" (con số đó có đúng với
sự thật đã biết không?) thì lại phát hiện vấn đề:

1. Khoảng tin cậy cũ rộng gấp khoảng 2 lần mức cần thiết. **Đã sửa.**
2. Một cách ước lượng "liên tục" nghe hấp dẫn nhưng thực tế nhiễu hơn. **Đã loại bỏ, có ghi
   lại lý do.**
3. Một cách chia nhóm của nhóm tác giả khác (doubly-ranked) cho kết quả sai gần như 100% trong
   bối cảnh của `mrwin`. **Đã cảnh báo trong phần mềm và ghi lại.**

Bài học: *hai máy tính cùng ra một đáp số chưa chắc đáp số đó đúng.* Vì vậy giai đoạn tới đặt
trọng tâm vào **kiểm nghiệm bên ngoài**.

### Những chỗ tài liệu từng "lệch" với thực tế (✅ đã sửa ở G0, 29/09/2026)

- Tệp trích dẫn (`CITATION.cff`) vẫn ghi phần mềm là "bộ khung, chưa kiểm chứng".
- README nói có 3 bài hướng dẫn, nhưng thư mục hướng dẫn **không có trong kho mã**, vì một quy
  tắc trong `.gitignore` đã vô tình chặn các tệp `.Rmd`.
- Số phiên bản vẫn là 0.0.0.9000.
- Trích dẫn bài Pocock 2012 trong README ghi sai số tập/trang. Đúng là *Eur Heart J*
  33(2):176–182.

---

## 5. Kế hoạch từng bước

Thời gian tính theo **tuần làm việc tập trung**, bắt đầu từ tuần 1 (đầu tháng 10/2026). Đây là
ước tính; các con số sẽ được điều chỉnh sau bước chạy thử ở G2.

### G0 — Dọn nhà và chuẩn bị (tuần 1) — ✅ Hoàn thành 29/09/2026

- **Kết quả:**
  - Có script `tools/baseline/run-baseline.sh` để chạy lại mọi kiểm tra bằng một lệnh.
  - Mọi kiểm tra đều đạt; độ bao phủ 88,30%.
  - Sửa trích dẫn sai trong 7 README.
  - Phát hiện các bài hướng dẫn chưa từng được lưu vào kho mã; sẽ viết ở G3.
  - Đã chốt D1 và D3.
  - Hồ sơ xin dữ liệu dbGaP đã soạn sẵn ở `tools/data-access/dbgap-dar-pack.md`; chủ nhiệm đề tài
    cần nộp.
- **Mục đích:** biết chính xác mình đang đứng ở đâu, trên chính máy của nhóm.
- **Việc cần làm:**
  - Cài R và các gói cần thiết trên máy trạm của nhóm. Chạy lại toàn bộ kiểm tra để có "ảnh
    chụp" ban đầu.
  - Sửa các chỗ tài liệu lệch nêu ở mục 4.
  - Chốt phạm vi phiên bản 1.0: cái gì làm ngay, cái gì để sau (xem D3).
  - **Nộp đơn xin dữ liệu thật ngay trong tuần này**, vì đây là việc chờ lâu nhất.
- **Xong khi:** mọi kiểm tra chạy xanh trên máy của nhóm và tài liệu khớp với thực tế.

### G1 — Khép lại các câu hỏi thống kê còn mở (tuần 2–5)

- **Mục đích:** không để lại "chỗ hổng" nào mà người phản biện có thể chỉ ra.
- **Việc cần làm:**
  - Hoàn thiện một bước "làm mượt" ma trận (Ledoit–Wolf): làm đúng công thức chuẩn, hoặc ghi
    rõ cách làm gần đúng và chứng minh nó không làm sai kết quả.
  - Giải thích vì sao khoảng tin cậy cũ rộng gấp 2 lần.
  - Bảo đảm ước lượng điểm và khoảng tin cậy chính **đo cùng một đại lượng**. Hiện hai thứ
    này chỉ trùng nhau khi tác động đồng đều giữa các nhóm.
  - Tìm hiểu vì sao phương pháp hơi "quá thận trọng" khi tính thêm độ bất định của dữ liệu di
    truyền.
  - Cho phép khai báo tương quan giữa các biến thể gen (LD), hoặc ghi rõ đây là hạn chế.
  - Quyết định giữ hay bỏ ba tính năng phụ còn treo (AL-CWR, bảng hiệu chỉnh sai lệch, GPS).
- **Xong khi:** mỗi câu hỏi có câu trả lời, có phép kiểm tra đi kèm và được ghi vào sổ theo
  dõi.

### G2 — Kiểm nghiệm lớn bằng mô phỏng (tuần 4–10)

- **Mục đích:** chứng minh phương pháp **đáng tin** trong nhiều tình huống.
- **Cách làm:** tạo hàng nghìn bộ dữ liệu giả mà ta **biết trước đáp án**, rồi xem phương pháp
  có tìm đúng không. Giống như cho học sinh làm đề đã có đáp án.
- **Các tình huống thử:** không có tác động; có tác động thật; gen tác động "đi đường vòng"
  (pleiotropy); công cụ di truyền yếu; các biến cố đi ngược chiều nhau; nhiều mức cỡ mẫu
  (5.000 / 20.000 / 100.000 người) và số nhóm (5 / 10 / 20).
- **So sánh với:** các cách phân tích thông thường, để cho thấy `mrwin` mang lại điều gì mới.
- **Quy tắc:** viết sẵn kế hoạch mô phỏng và lưu vào kho mã **trước khi chạy**, để không ai
  nghi ngờ nhóm "chọn kết quả đẹp".
- **Xong khi:** tỉ lệ báo động giả khoảng 5% và độ phủ khoảng 95% ở các tình huống hợp lệ; mọi
  kết quả chạy lại được bằng một lệnh.

### G3 — Hoàn thiện phần mềm để phát hành (tuần 4–9, làm song song với G2)

- **Mục đích:** biến phần mềm "chạy được" thành phần mềm "phát hành được".
- **Việc cần làm:**
  - Thêm phép kiểm tra tự động: cách tính nhanh luôn cho **đúng từng chữ số** như cách tính
    chậm.
  - Đặt cách tính nhanh làm mặc định.
  - Khôi phục và viết lại 3 bài hướng dẫn: mô phỏng nhanh, mẫu cho đoàn hệ thật, cách đọc các
    cảnh báo.
  - Viết nhật ký thay đổi (NEWS), nâng phiên bản, cập nhật tệp trích dẫn.
  - Chạy kiểm tra chuẩn CRAN; đo độ bao phủ kiểm thử (mục tiêu ≥ 80%).
  - Nộp lên CRAN.
- **Xong khi:** CRAN chấp nhận gói.

### G4 — Ví dụ trên dữ liệu thật (song song, phụ thuộc thời điểm được cấp dữ liệu)

- **Mục đích:** cho thấy phương pháp dùng được ngoài đời thực.
- **Việc cần làm:**
  - **Trước khi xin dữ liệu:** dùng bộ mô phỏng để ước tính cần bao nhiêu người thì phương
    pháp đủ "nhạy" (tính cỡ mẫu).
  - Viết sẵn kế hoạch phân tích: câu hỏi, thứ tự biến cố, biến phơi nhiễm, cách xây điểm di
    truyền. Việc này làm **trước khi nhìn dữ liệu**.
  - Phân tích trên máy chủ bảo mật của cơ sở. **Dữ liệu không bao giờ được đưa vào kho mã.**
  - Làm thêm các phân tích độ nhạy (đổi số nhóm, hiệu chỉnh yếu tố nhiễu, kiểm tra
    pleiotropy).
- **Xong khi:** có phần "Ứng dụng" hoàn chỉnh trong bài báo, báo cáo theo STROBE-MR.

### G5 — Viết và nộp bài báo (tuần 8–16)

- **Mục đích:** được giới khoa học công nhận.
- **Việc cần làm:** hoàn thiện bản thảo theo khung đã có. Đăng bản preprint mới lên arXiv.
  Nộp tạp chí. Trả lời phản biện.
- **Nguyên tắc vàng:** không viết điều gì mà chưa có phép kiểm tra chứng minh. Các **phát hiện
  "âm tính"** (những gì đã thử mà không hiệu quả) cũng được báo cáo trung thực, vì đó là đóng
  góp có giá trị.
- **Xong khi:** bài được chấp nhận.

### G6 — Bàn giao và duy trì (sau khi phát hành)

- Lưu trữ phiên bản 1.0 kèm mã định danh DOI. Đồng bộ README 7 ngôn ngữ. Theo dõi phản hồi của
  người dùng.
- Lên danh sách cho **bài báo thứ hai**: dùng dữ liệu tóm tắt (không cần dữ liệu từng người),
  nhiều hơn 3 mức ưu tiên, và phân tích trong từng nhóm.

---

## 6. Lịch trình tổng thể

| Tuần | 1 | 2–3 | 4–5 | 6–7 | 8–9 | 10–12 | 13–16 | Sau 16 |
|---|---|---|---|---|---|---|---|---|
| G0 Chuẩn bị | ■ | | | | | | | |
| G1 Câu hỏi thống kê | | ■ | ■ | | | | | |
| G2 Mô phỏng | | | ■ | ■ | ■ | ■ | | |
| G3 Phần mềm | | | ■ | ■ | ■ | | | |
| G4 Dữ liệu thật | đơn xin | chờ | chờ | ■ | ■ | ■ | ■ | |
| G5 Bài báo | | | | | ■ | ■ | ■ nộp | phản biện |
| G6 Bàn giao | | | | | | ■ CRAN | | ■ |

**Đường găng** (chậm một việc là chậm cả dự án): G0 → G1 → G2 → G5. Việc **dễ chậm nhất** là
G4, vì phải chờ cấp dữ liệu. Vì thế đơn xin phải nộp ngay tuần 1. Nếu dữ liệu đến muộn, có thể
nộp bài với phần mô phỏng trước và bổ sung ứng dụng trong vòng sửa bài.

---

## 7. Những quyết định nhóm cần chốt

| Mã | Câu hỏi | Gợi ý |
|---|---|---|
| **D1** | Dùng đoàn hệ nào cho ví dụ thật? Có chấp nhận nền tảng đám mây của biobank không? | ✅ **Đã chốt 29/09/2026:** dbGaP (ARIC trước, sau đó MESA/CHS/FHS), phân tích trên máy chủ của cơ sở, không dùng đám mây. |
| **D2** | Nộp tạp chí nào? | *International Journal of Epidemiology* (đã định hướng từ trước); dự phòng: *Statistics in Medicine*, *Genetic Epidemiology*. |
| **D3** | Phiên bản 1.0 gồm những gì? | ✅ **Đã chốt 29/09/2026:** bỏ GPS (chuyển sang bài báo thứ hai); **có** hỗ trợ tương quan gen (LD), làm trong G1. |
| **D4** | Mặc định tính độ bất định bằng cách nào? | Quyết định **sau** G2, dựa trên kết quả mô phỏng. |
| **D5** | Có đưa AL-CWR và bảng hiệu chỉnh sai lệch vào bài báo đầu? | Chỉ đưa vào nếu làm xong và kiểm chứng kịp trong G1. |

---

## 8. Nguyên tắc "không dùng đám mây"

- **Tính toán:** mô phỏng, đo tốc độ và kiểm tra phát hành đều chạy trên **máy trạm hoặc máy
  chủ của cơ sở**. Gợi ý cấu hình: 16–32 lõi CPU, 64 GB RAM. Chạy song song trên nhiều lõi và
  lưu kết quả dở dang ra đĩa để có thể chạy tiếp nếu bị ngắt.
- **Dữ liệu bệnh nhân:** chỉ nằm trên máy chủ bảo mật đã được hội đồng đạo đức và đơn vị cấp
  dữ liệu phê duyệt. Kho mã đã có quy tắc chặn mọi tệp dữ liệu.
- **Lưu ý quan trọng:** theo chính sách hiện hành của UK Biobank (cần xác nhận lại khi nộp
  đơn), dữ liệu cá nhân chỉ được phân tích trên nền tảng đám mây của họ. Nếu giữ nguyên tắc
  "không đám mây", nên ưu tiên:
  - các đoàn hệ trên dbGaP như ARIC, CHS, Framingham, MESA, được phân tích tại máy chủ bảo mật
    của cơ sở theo cam kết sử dụng dữ liệu;
  - Rhineland Study, thông qua đồng tác giả tại DZNE.
  Các đoàn hệ này nhỏ hơn UK Biobank, nên **phải tính cỡ mẫu trước** (G4).
- **Mã nguồn:** vẫn lưu trên GitHub như hiện nay. Đây chỉ là nơi lưu mã, không chứa dữ liệu và
  không dùng để tính toán. Các kiểm tra chạy trên máy của nhóm là chuẩn chính thức; kiểm tra tự
  động trên GitHub chỉ là lớp bổ sung và có thể tắt nếu nhóm muốn.

---

## 9. Rủi ro và cách phòng

| Rủi ro | Nếu xảy ra thì sao | Cách phòng |
|---|---|---|
| Chờ cấp dữ liệu quá lâu | Chậm phần ứng dụng | Nộp đơn tuần 1; nộp bài với mô phỏng trước, bổ sung ứng dụng khi sửa bài |
| Đoàn hệ quá nhỏ, tín hiệu di truyền yếu | Khoảng tin cậy quá rộng, không kết luận được | Tính cỡ mẫu trước; cân nhắc gộp nhiều đoàn hệ |
| Kiểm nghiệm lại phát hiện lỗi mới | Phải sửa, mất thời gian | Dự phòng 2–3 tuần; coi đây là điều tốt vì lỗi được tìm ra trước khi đăng |
| Gen tác động "đi đường vòng" (pleiotropy) | Kết quả có thể sai lệch | Luôn báo cáo chẩn đoán SDPD và khoảng tin cậy có giới hạn pleiotropy |
| Phản biện khó hiểu cách diễn giải DS-CWR | Bài bị yêu cầu sửa nhiều | Viết rõ phần diễn giải, có ví dụ lâm sàng, nêu rõ giới hạn |

---

## 10. Bảng thuật ngữ

| Thuật ngữ | Giải thích ngắn |
|---|---|
| Win ratio (tỉ số thắng) | Số lần "thắng" chia số lần "thua" khi so sánh từng cặp theo thứ tự nặng–nhẹ |
| Mendelian randomization (MR) | Dùng gen như một "xổ số tự nhiên" để suy ra quan hệ nhân quả |
| Điểm nguy cơ đa gen (PRS) | Tổng có trọng số của nhiều biến thể gen, đo "khuynh hướng di truyền" |
| DS-CWR | Con số tóm tắt của `mrwin`: < 1 là có hại, > 1 là có lợi |
| Khoảng tin cậy Fieller | Cách tính khoảng tin cậy cho một tỉ số; là khoảng chính mà `mrwin` báo cáo |
| Pleiotropy | Gen ảnh hưởng kết cục qua con đường khác, không qua yếu tố đang nghiên cứu |
| SDPD | Bộ chẩn đoán pleiotropy của `mrwin` |
| Báo động giả (sai lầm loại I) | Kết luận "có tác động" trong khi thực tế không có; mục tiêu khoảng 5% |
| Độ phủ | Tỉ lệ số lần khoảng tin cậy chứa giá trị thật; mục tiêu khoảng 95% |
| CRAN | Kho phần mềm chính thức của R |
| STROBE-MR | Chuẩn quốc tế về cách báo cáo nghiên cứu MR |
