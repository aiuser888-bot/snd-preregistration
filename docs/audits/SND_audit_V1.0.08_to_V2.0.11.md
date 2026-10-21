# SND 稽核 V2.0.11（接續 V2.0.10 稽核；V2.0.10 與舊檔保留）

## 輸入
作者把 V2.0.10 上傳到 GitHub 後回報兩件事：(1) 公式顯示成殘缺文字，要求做新版本、把歷史檔案全部修正並確認所有公式沒有錯誤；(2) 授權條款有三個問題：授權檔裡留著草稿提示、涵蓋範圍有缺口、GitHub 只認得 `LICENSE` 因而把整個倉庫標成 MIT。
本版**沒有改動規範的數學、程序或數字**。

## 評估與處置
### 一、公式渲染
1. **原因：屬實。** 原檔用 `\(` `\)` 與 `\[` `\]` 包公式；在 Markdown 裡 `\(` 是跳脫括號，GitHub 不渲染，只會吃掉反斜線，留下殘缺文字。這是格式問題，不是公式錯誤。
2. **處置：** 行內公式改為 GitHub 官方文件建議的寫法（對含有 Markdown 特殊字元的公式用 ``$`…`$`` 形式），獨立一行的公式改為 ` ```math ` 區塊。只改分隔符，其間的 LaTeX 逐位元不動。
3. **範圍：** V2.0.01–V2.0.10 全部規範（V2.0.10 移入 `protocol/history/`）、新的 V2.0.11 規範、V2.0.02–V2.0.04 稽核：14 個檔案、7235 條公式（顯示式 14 條）。V2.0.05–V2.0.10 的稽核與其他文件沒有舊式公式。

### 二、授權條款（三個問題都屬實）
1. **草稿提示：** `LICENSE-DOCS.txt` 只有 319 位元組，內容是 CC BY 的網址加一行「Choose another licence, e.g. CC0, if you prefer…」，沒有授權全文。
2. **涵蓋缺口：** MIT 只寫 `code/` 與 `legacy/`，CC BY 只寫 `protocol/` 與 `docs/`。`results/`、`metadata/`、根目錄的 README／CHANGELOG／VERSIONING／OSF 指南、CITATION.cff、.zenodo.json、SHA256SUMS.txt 與 `.github/` 沒有任何授權涵蓋。
3. **GitHub 標示：** GitHub 只認得 `LICENSE`（MIT 全文），整個倉庫被標成 MIT，讀者看不到規範文件是 CC BY 4.0。
4. **處置：** 刪除 `LICENSE-DOCS.txt`；`LICENSE` 改為說明文件（中英對照，列出兩種授權各涵蓋的路徑）；全文放入 `LICENSES/MIT.txt` 與 `LICENSES/CC-BY-4.0.txt`；新增逐檔對照 `REUSE.toml`（程式與自動化＝MIT；規範、文件、結果、中繼資料＝CC BY 4.0）；新增 `code/scripts/check_licenses.py`；README 開頭與授權段改為對照表；CITATION.cff 加註、兩份 Zenodo 中繼資料的 `notes` 改寫。
5. **GitHub 的限制：** 一個倉庫只能顯示一種授權，所以無法顯示雙授權；預期它改顯示「其他」（View license）。這是誠實的結果，不是缺陷；README 與 `LICENSE` 負責說清楚。

## 自查紀錄
- **【錯誤】檢查腳本自己的缺陷：** `check_math.py` 第一版把自己的 JS 檔當資料讀（`process.argv` 索引錯誤），KaTeX 步驟被誤報為「略過」而 exit 0 → 修正，並規定「node 與 katex 都在、但執行器失敗」算失敗。
- **【錯誤】數字：** 建置腳本一度把 CHANGELOG 的公式數寫成 13 個檔案、6172 條，漏算移入 history 的 V2.0.10 → 以 `check_math.py` 的實測輸出為準，更正為 14 個檔案、7235 條（顯示式 14 條）。
- **【錯誤】測試脆弱：** `test_theta_is_pure_function_of_seed_and_counter` 用字面值 "2.0.11" 當「另一個版本」，版本升級後與預設值相同而失敗 → 改用 "0.0.0-other"，並新增斷言「`seed_material` 預設版本標籤等於 `snd.__version__`」。
- **【錯誤】建置腳本：** Zenodo 中繼資料的字串取代順序有誤（會把 "V2.0.01-V2.0.10" 誤改成 "V2.0.01-V2.0.11"）→ 還原整個工作樹後重跑。
- **【缺陷】授權檢查：** 初版把「草稿提示」檢查也套到 CC BY 4.0 官方全文（逐字文本不該被改）→ 限於我們自己寫的檔案；另發現未加入 git 的新檔不在 `git ls-files` 內（檢查腳本以追蹤檔案為準，須先 `git add`）。

## 驗證狀態
| 項目 | 狀態 | 依據 |
|---|---|---|
| 公式分隔符無殘留、括號與 `\left`/`\right` 成對、無游離錢字號 | 已驗證 | `check_math.py`：14 個檔案、7235 條公式、0 個問題 |
| 公式可被 KaTeX 解析 | 已驗證 | KaTeX 0.16.47：7235 條（797 種不重複）、0 失敗 |
| 內容未改 | 已驗證 | 14 個檔案轉回舊寫法後與上傳原檔逐位元相同（V2.0.11 僅標題行不同）；去除公式後的正文也相同；`convert_math.py selftest` 通過 |
| 檢查腳本確實會失敗 | 已驗證 | `check_math.py` 6 種人為錯誤、`check_licenses.py` 8 種人為錯誤，全部 exit 1；正常樣本 exit 0 |
| 授權涵蓋 | 已驗證 | `check_licenses.py`：每個追蹤檔案恰屬一種授權，無遺漏、無重複；兩份全文齊全；CITATION／Zenodo 授權欄位一致 |
| 程式與數值 | 已驗證 | 18 項測試全過；40 項數值檢查全過 |
| 公式在 github.com 實際渲染 | **未驗證** | 本環境連不上 GitHub 的渲染服務；作者需在頁面上看一眼（例如規範第 33 行附近） |
| KaTeX 通過 = GitHub 通過 | **未驗證** | GitHub 用 MathJax，與 KaTeX 不是同一個解析器；通過是強證據，不是證明 |
| 其他檢視器（如 VS Code 內建預覽）顯示 | **未驗證，預期不渲染** | 這種寫法是 GitHub 專用 |
| GitHub 授權偵測改顯示「其他」 | **未驗證** | 依其偵測規則推論，未實測 |
| CC BY 4.0 全文與官方逐字相同 | **未驗證** | 本環境連不上 creativecommons.org；全文取自 TeX Live 的 doclicense 套件（官方純文字版），所有 8 節齊全、無 TeX 殘留；SHA-256 `9ba9550ad48438d0836ddab3da480b3b69ffa0aac7b7878b5a0039e7ab429411`。建議登記前對官方文字核對一次 |
| 官方 REUSE 工具 `reuse lint` | **未執行** | 無法安裝；`check_licenses.py` 是自寫的等價檢查 |

「所有公式沒有錯誤」在這裡的意思是：**語法與結構沒有錯誤、能被解析**。公式在數學上是否正確，仍由先前各版稽核的「已證明／數值支持／未驗證」標籤負責，本版沒有重新證明任何公式。

## 沒有改動／仍未解
V2.0.10 稽核列出的未解項目全部原樣保留（E2 物理假設、理想化 (I)(H)、§12.5.1 程序的實作、跨語言獨立 E1 實作、§12.5.5 θ 的物理施加、閘門是否保留等）。

## 作者需要決定／確認
1. 授權對應的預設：`results/` 歸 CC BY 4.0；`.github/` 與 `.gitignore` 歸 MIT。
2. 接受 GitHub 側欄顯示「其他」；若想讓側欄顯示單一授權，就得選一個主授權，但會讓另一類檔案的讀者被誤導。
3. 登記（OSF／Zenodo）前，把 `LICENSES/CC-BY-4.0.txt` 對 creativecommons.org 的官方文字核對一次。
4. 著作權人目前是筆名 Pietersite；若要加入真名或 ORCID，須同步改 `LICENSES/MIT.txt`、`REUSE.toml`、`CITATION.cff`、`.zenodo.json`。
5. 推送 V2.0.11 會讓倉庫的自動流程（`.github/workflows/release.yml`，依 `CITATION.cff` 的版本）建立公開的 tag 與 release `v2.0.11`。

## 標籤
已證明：同 V2.0.10（本版沒有新增或修改任何數學陳述）。已驗證：上表標示者。未驗證：上表標示者。
