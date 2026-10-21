# SND 稽核 V2.0.12（接續 V2.0.11 稽核；V2.0.11 與舊檔保留）

## 輸入
作者把 V2.0.11 上傳到 GitHub 後回報：公式雖然能顯示，但「看起來還是不太對」，某些地方符號會異常高起（字母上方的橫線、帽子、波浪線，以及帶波浪符號的量的平方）。作者在 GitHub 的 issue 預覽（不送出）貼上測試片段並附截圖，經數輪對照確認兩件事：(a) 100% 縮放下相同，不是瀏覽器縮放造成；(b) 用規範逐字的公式行貼進預覽也重現，不是規範頁面特有的問題。
本版**沒有改動規範的數學、程序或數字**。

## 評估與處置
### 一、三種現象的判定
1. **`\bar f`（字母上的橫線）：屬實，要改。** 預覽中 `\bar f_b` 只畫出一小段飄在字母左上方的短線；同一位置的 `\overline{f}_b` 畫出完整橫線，帶波浪符號的 `\overline{\widetilde f}_b` 與一般的 `\overline{w}U` 也都正常。
2. **帶重音符號的量加平方（例如 f̃ 的平方）：屬實，要改。** 預覽中 `\widetilde f_{b,i}^{\,2}`、`\widehat p_0\widehat f^2`、`\widehat\sigma_{\rm rep}^2` 的指數都被疊到重音符號之上；寫成 `(\widetilde f_{b,i})^{2}` 等括號形式後，指數回到正常高度。
3. **帽子與波浪符號蓋在高字母（如 f）上的位置：判定為正常排版，不改。** 預覽中 `\widehat f`、`\hat f`、`\hat{f}` 畫出的結果相同，所以換寫法沒有幫助；重音符號落在高字母的頂端之上是 TeX 與 MathJax 的標準位置。若作者仍覺得太高，剩下的做法是換記號（例如 `\overset{\wedge}{f}`，但字形會變成尖角符號），我不建議。

### 二、處置
1. 只改兩類寫法，其餘一字不動：`\bar X` → `\overline{X}`；「帶重音符號的量＋指數」→ 括號包住該量再加指數（下標留在括號內；`^{\,2}` 的細空白省略）。適用的重音符號：`\widehat`、`\widetilde`、`\hat`、`\tilde`、`\overline`。
2. 範圍：V2.0.12 規範（由 V2.0.11 產生）、V2.0.01–V2.0.11 全部歷史規範、V2.0.02–V2.0.04 稽核，共 15 個檔案、204 條公式被改寫（33 種不重複的改寫；`\bar` 137 處、帶重音符號加指數 102 處）。V2.0.05 之後的稽核與其他文件沒有這兩類寫法。
3. 工具：`code/scripts/render_fix.py`（`apply`／`diff`／`selftest`）；`code/scripts/check_math.py` 新增規則：公式若會被 `render_fix` 再次改寫，就報為問題，以免日後新寫的公式又用回這兩種寫法。
4. 為什麼歷史檔也改：作者在 V2.0.11 要求歷史檔案的公式都要能正確顯示；這兩類寫法在歷史檔裡有同樣的排版問題。代價是歷史檔不再與上傳原檔逐位元對應；V2.0.11 發行時的狀態在 git commit cb217e8（tag v2.0.11），上傳原檔在 commit 25498f6。

### 三、規範 V2.0.12 被改寫的 27 條公式（逐條前後對照）
- 第 22 行：`\sigma_f^2=\int\widehat p_0\widehat f^2` → `\sigma_f^2=\int\widehat p_0(\widehat f)^2`
- 第 27 行：`\bar f_b=L^{-1}\sum_if_{b,i}` → `\overline{f}_b=L^{-1}\sum_if_{b,i}`
- 第 27 行：`\widehat g_j=K^{-1}\sum_b(f_{b,j}-\bar f_b)` → `\widehat g_j=K^{-1}\sum_b(f_{b,j}-\overline{f}_b)`
- 第 28 行：`\widetilde f_{b,j}=f_{b,j}-\bar f_b-\widehat g_j` → `\widetilde f_{b,j}=f_{b,j}-\overline{f}_b-\widehat g_j`
- 第 28 行：`\bar{\widetilde f}_b` → `\overline{\widetilde f}_b`
- 第 33 行：`U_b=\sum_{i=1}^{L}\widetilde f_{b,i}\,e^{\mathrm{i}\theta_{b,i}}\in\mathbb C,\qquad S_b^2=\sum_{i=1}^{L}\widetilde f_{b,i}^{\,2},\qquad V_b=\dfrac{L}{2(L-1)}S_b^2` → `U_b=\sum_{i=1}^{L}\widetilde f_{b,i}\,e^{\mathrm{i}\theta_{b,i}}\in\mathbb C,\qquad S_b^2=\sum_{i=1}^{L}(\widetilde f_{b,i})^{2},\qquad V_b=\dfrac{L}{2(L-1)}S_b^2`
- 第 57 行：`\mathbb E\sum_j(f_j-\bar f)^2=L\sigma_f^2-L\,\mathrm{Var}(\bar f)=\sigma_f^2(L-1-2S_L)` → `\mathbb E\sum_j(f_j-\overline{f})^2=L\sigma_f^2-L\,\mathrm{Var}(\overline{f})=\sigma_f^2(L-1-2S_L)`
- 第 57 行：`L^2\mathrm{Var}(\bar f)=\sigma_f^2(L+2\sum_{k}(L-k)\rho_k)` → `L^2\mathrm{Var}(\overline{f})=\sigma_f^2(L+2\sum_{k}(L-k)\rho_k)`
- 第 60 行：`\mathrm{Var}(\bar f)\ge0\Rightarrow1+2S_L\ge0\Rightarrow\kappa\le1+1/(L-1)` → `\mathrm{Var}(\overline{f})\ge0\Rightarrow1+2S_L\ge0\Rightarrow\kappa\le1+1/(L-1)`
- 第 96 行：`\widehat x^2\sim` → `(\widehat x)^2\sim`
- 第 117 行：`V^{\rm row}=\sum_b\frac{L}{2(L-1)}\sum_j(f_{b,j}-\bar f_b)^2` → `V^{\rm row}=\sum_b\frac{L}{2(L-1)}\sum_j(f_{b,j}-\overline{f}_b)^2`
- 第 126 行：`\sigma_f^2=\mathbb E[\widehat f^2]` → `\sigma_f^2=\mathbb E[(\widehat f)^2]`
- 第 131 行：`\widehat\sigma_f^2` → `(\widehat\sigma_f)^2`
- 第 177 行：`|\widehat\mu_{\rm rep}-\widehat\mu_{\rm orig}|^2/(\widehat\sigma_{\rm rep}^2+\widehat\sigma_{\rm orig}^2)\le\chi^2_{2,0.95}=5.991` → `|\widehat\mu_{\rm rep}-\widehat\mu_{\rm orig}|^2/((\widehat\sigma_{\rm rep})^2+(\widehat\sigma_{\rm orig})^2)\le\chi^2_{2,0.95}=5.991`
- 第 231 行：`\bar Q_c=n_c^{-1}\sum_{i\in c}e^{\mathrm i\theta_i}` → `\overline{Q}_c=n_c^{-1}\sum_{i\in c}e^{\mathrm i\theta_i}`
- 第 232 行：`U_{\rm mon}=\sum_{c}F_c\,\bar Q_c .` → `U_{\rm mon}=\sum_{c}F_c\,\overline{Q}_c .`
- 第 236 行：`\mathbb E[U_{\rm mon}\mid G]=\sum_cF_c\bar q_c=:B_{\rm mon}` → `\mathbb E[U_{\rm mon}\mid G]=\sum_cF_c\overline{q}_c=:B_{\rm mon}`
- 第 236 行：`\bar q_c` → `\overline{q}_c`
- 第 236 行：`B_{\rm res}=\sum_i\widetilde f_i(q_i-\bar q_{c(i)})` → `B_{\rm res}=\sum_i\widetilde f_i(q_i-\overline{q}_{c(i)})`
- 第 238 行：`\eta^2=\frac{\sum_cF_c^2/n_c}{\sum_i\widetilde f_i^{\,2}}\in[0,1]` → `\eta^2=\frac{\sum_cF_c^2/n_c}{\sum_i(\widetilde f_i)^{2}}\in[0,1]`
- 第 306 行：`\mathbb E_0[(\mathrm{Re}\,\bar wU)^2]=V` → `\mathbb E_0[(\mathrm{Re}\,\overline{w}U)^2]=V`
- 第 307 行：`\mathbb E|U_b|^2=\sum_j\widetilde f_j^2-\tfrac1{L-1}\bigl[(\sum_j\widetilde f_j)^2-\sum_j\widetilde f_j^2\bigr]=\tfrac L{L-1}\sum_j\widetilde f_j^2=2V_b` → `\mathbb E|U_b|^2=\sum_j(\widetilde f_j)^2-\tfrac1{L-1}\bigl[(\sum_j\widetilde f_j)^2-\sum_j(\widetilde f_j)^2\bigr]=\tfrac L{L-1}\sum_j(\widetilde f_j)^2=2V_b`
- 第 307 行：`(\mathrm{Re}\,\bar wU)^2=\tfrac12|U|^2+\tfrac12\mathrm{Re}(\bar w^2U^2)` → `(\mathrm{Re}\,\overline{w}U)^2=\tfrac12|U|^2+\tfrac12\mathrm{Re}((\overline{w})^2U^2)`
- 第 308 行：`\ell-1\propto\mathrm{Re}(\bar wU)` → `\ell-1\propto\mathrm{Re}(\overline{w}U)`
- 第 309 行：`|B|=\mathbb E_0[(\ell-1)\mathrm{Re}(\bar wU)]\le\sqrt{\mathbb E_0(\ell-1)^2}\sqrt{\mathbb E_0(\mathrm{Re}\,\bar wU)^2}=\sqrt{D_GV}` → `|B|=\mathbb E_0[(\ell-1)\mathrm{Re}(\overline{w}U)]\le\sqrt{\mathbb E_0(\ell-1)^2}\sqrt{\mathbb E_0(\mathrm{Re}\,\overline{w}U)^2}=\sqrt{D_GV}`
- 第 318 行：`\mathrm{Re}(\bar wU(S))` → `\mathrm{Re}(\overline{w}U(S))`
- 第 343 行：`\mathbb E[(\mathrm{Re}\,\bar wU)^2]/V` → `\mathbb E[(\mathrm{Re}\,\overline{w}U)^2]/V`

意義核對：`\bar X` 與 `\overline{X}` 都是「X 的上橫線」（樣本平均、複共軛的記號，與先前定義相同）；`A^s` 與 `(A)^s` 在 A 是單一符號加重音（含下標）時是同一個數學物件（先加重音、再取冪）。這 27 條只有這兩種變化，沒有別的。另 6 種只出現在舊版的改寫（`\sum_j\widetilde f^2`、`\mathbb E[\widehat f^2]`、`\sum\widetilde Z^2` 等）屬同樣兩類。

## 自查紀錄
- **【錯誤】自測期望值寫錯：** `render_fix.py` 第一版自測以為 `\bar x^2` 應得 `\overline{x}^2`，實際（正確）輸出是 `(\overline{x})^2`——橫線之上同樣會把指數頂高 → 改正自測的期望值，腳本邏輯不變。
- **【誤報】臨時檢查：** 為了獨立確認「只有公式被改」，另寫了一支逐行比較的臨時腳本，它把 ` ```math ` 區塊內的裸公式行誤判成文字差異（3 行）→ 這是臨時腳本只遮蔽行內公式的緣故；正式的比對（建置腳本內、含 ` ```math ` 區塊）全部通過。
- **【環境】** 本環境沒有安裝 pytest 於預設 python，測試改用專案自帶的無依賴執行器（`python tests/test_snd.py`，見 README）。

## 驗證狀態
| 項目 | 狀態 | 依據 |
|---|---|---|
| 公式分隔符無殘留、括號成對、無游離錢字號、無兩類排版隱患 | 已驗證 | `check_math.py`：15 個檔案、8298 條公式、0 個問題 |
| 公式可被 KaTeX 解析 | 已驗證 | KaTeX 0.16.47：8298 條（797 種不重複）、0 失敗 |
| 除公式外文字未改 | 已驗證 | 15 個檔案遮蔽公式後逐位元相同、行數相同；改寫冪等（再跑一次沒有任何改動） |
| 改寫規則正確 | 已驗證 | `render_fix.py selftest`：20 個公式案例（改寫、不該改、冪等）加 Markdown 遮蔽（雙反引號範例、程式碼區塊、` ```math ` 區塊）；規範的 27 條與舊版 6 種改寫已逐條核對意義 |
| 新規則確實會失敗 | 已驗證 | 人為樣本：`\bar`、行內帶重音符號加指數、顯示式帶重音符號加指數、括號不成對，全部 exit 1；正常樣本 exit 0 |
| 程式與數值 | 已驗證 | 18 項測試全過；40 項數值檢查全過 |
| 在 github.com 上「改」比「原」正常 | 已驗證（作者實測） | 作者在 issue 預覽對 6 組原／改的截圖：3 組橫線、3 組重音符號加平方，「改」全部較正常 |
| 全部檔案在 github.com 實際渲染 | **未驗證** | 沒有逐頁看過；請作者開規範第 28、33、307 行附近看一眼 |
| KaTeX 通過 = GitHub 通過 | **未驗證** | GitHub 用 MathJax，與 KaTeX 不是同一個解析器；通過是強證據，不是證明 |
| 帽子／波浪符號在高字母上的位置 | **未改（屬標準排版）** | 預覽中 `\widehat`、`\hat`、`\hat{f}` 三者相同；是否太高屬主觀，作者可決定要不要換記號 |
| 授權涵蓋 | 已驗證 | 新檔 `render_fix.py` 屬 `code/`（MIT）、本稽核檔屬 `docs/`（CC BY 4.0）；`check_licenses.py`：106 個檔案各屬一種授權（CC BY 4.0 41、MIT 65），無遺漏、無重複，0 個問題 |

「所有公式沒有錯誤」在這裡的意思是：**語法與結構沒有錯誤、能被解析、已知的兩類排版隱患已消除**。公式在數學上是否正確，仍由先前各版稽核的「已證明／數值支持／未驗證」標籤負責，本版沒有重新證明任何公式。

## 沒有改動／仍未解
V2.0.10 稽核列出的未解項目全部原樣保留（E2 物理假設、理想化 (I)(H)、§12.5.1 程序的實作、跨語言獨立 E1 實作、§12.5.5 θ 的物理施加、閘門是否保留等）；V2.0.11 稽核列出的授權待確認事項也原樣保留（`results/` 歸 CC BY 4.0、`.github/` 與 `.gitignore` 歸 MIT 的預設；登記前把 `LICENSES/CC-BY-4.0.txt` 對 creativecommons.org 官方文字核對一次；GitHub 側欄預期顯示「其他」）。

## 作者需要決定／確認
1. 開 GitHub 上的規範第 28、33、307 行附近看一眼：橫線是否完整、平方的 2 是否在正常位置。
2. 帽子與波浪符號蓋在高字母上的位置是否可以接受（我判定為標準排版，沒有改）。
3. 歷史檔（V2.0.01–V2.0.11）也套用了這兩類改寫；若寧可讓歷史檔保持 V2.0.11 的寫法，請告訴我，我會只保留 V2.0.12 規範與工具的改動。
4. 推送 V2.0.12 會讓倉庫的自動流程（`.github/workflows/release.yml`，依 `CITATION.cff` 的版本）建立公開的 tag 與 release `v2.0.12`。

## 標籤
已證明：同 V2.0.10（本版沒有新增或修改任何數學陳述）。已驗證：上表標示者。未驗證：上表標示者。
