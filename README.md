# Project SND — Pre-registration Protocol V2.0.12, Reference Implementation and Validation Record

**Status: pre-registration draft with reduced-scale simulation validation. No real experimental data exist or are included.**

**Licences:** code = MIT; protocol, documents, results and metadata = CC BY 4.0 (not MIT). See `LICENSE`.

## What is being tested
Whether a feature `x` that is **sealed first** is statistically dependent on a basis angle `θ` that is **assigned later by a fair random mechanism**. Under the null (θ independent of x) nothing is expected. The protocol specifies (i) a block-randomization test whose validity follows from the randomization design, and (ii) a one-sided amplitude upper confidence bound (Feldman–Cousins ordering of a Rice distribution). The protocol makes **no claim** about physical mechanisms; its interpretation layers (protocol §11) explicitly say that a rejection first calls for checks of the random number generator, assignment code, timestamps and analysis pipeline.

## Contents
| Path | What it is |
|---|---|
| `protocol/SND_protocol_V2.0.12.md` | The pre-registration text (current). Same content as V2.0.10; only the formula markup was changed so that github.com renders and typesets it properly (V2.0.11: delimiters; V2.0.12: `\bar` and accented squares). `protocol/history/` has V2.0.01–V2.0.11 |
| `docs/audits/SND_audit_V1.0.08_to_V2.0.12.md` (+ earlier) | Claim-by-claim audit: each claim labelled proved / numerically supported / failed / unverified, with seeds and failure cases |
| `CHANGELOG.md` | Per-version record of errors/contradictions found, fixes, and changed verification methods (replaces publishing AI conversation logs; AI chats are not published) |
| `docs/SND_evolution_for_everyone.md`, `docs/SND_evolution_summary.md` | Plain-language and technical development history |
| `code/snd/` | Reference implementation: block angle generation (Philox, rejection-sampled permutations), two-way centred statistic, randomization test, FC–Rice table and UCB, exact κ, leakage sensitivity (`leak.py`: type-I curve e_I(β), coverage under bias, environment-projected test U_mon), reduced-scale acceptance simulation |
| `code/tests/test_snd.py` | 18 tests (all pass, see below) |
| `code/scripts/` | `build_table.py` (rebuilds FC–Rice table and checks the protocol's checkpoints), `reproduce_all.py` (headline numbers), `leak_tables.py` (analytic leakage curves → `results/leak_curves.json`), `check_numbers.py` (recomputes and asserts the closed-form numbers quoted in protocol §6/§12) |
| `results/` | `fc_rice_table.npz` (SHA-256 below), `reproduce.json/.log`, `environment.txt` |
| `LICENSE`, `LICENSES/`, `REUSE.toml` | Licence notice (which licence covers which path), full licence texts (MIT, CC BY 4.0) and the machine-readable per-file mapping |
| `legacy/` | Original exploratory scripts with **known-issue notes** (not the reference implementation) |

## Reproduce
```
pip install numpy scipy
cd code
PYTHONPATH=. python scripts/build_table.py ../results/fc_rice_table.npz   # exit 0 = checkpoints and acceptance mass OK
PYTHONPATH=. python tests/test_snd.py                                    # 18 PASS
PYTHONPATH=. python scripts/leak_tables.py                               # leakage curves (instant)
PYTHONPATH=. python scripts/check_numbers.py                             # exit 0 = protocol's closed-form numbers consistent
PYTHONPATH=. python scripts/check_math.py                                # formula markup hygiene and github.com rendering hazards in all .md files (+ KaTeX parse if node and katex are installed)
PYTHONPATH=. python scripts/render_fix.py selftest                        # self-test of the two formula rewrite rules used since V2.0.12
python scripts/check_licenses.py                                         # every file covered by exactly one licence (REUSE.toml), licence texts present
PYTHONPATH=. python scripts/reproduce_all.py 1000                        # ~3 min, writes results/reproduce.json
```
FC–Rice table SHA-256: `ee982331fd46e2adaac987e00d944c6fb7c652f9b5752ee2126f22c7db8f5f27` (the `.npz` file may differ byte-wise across numpy versions; the checkpoint values and acceptance-mass range are the verifiable content).

## What has and has not been established
| Claim | Label |
|---|---|
| Randomization test has exact level α (given independent random assignment with the stated mechanism; preprocessing and stopping rule use x only) | Proved (protocol §4, Lemma 1) |
| Var(U\|x) = ΣV_b exactly under H₀; cross-block covariance 0 | Proved (Lemma 2 (ii),(iii),(vi) for H₀) |
| κ ≤ 1 + 1/(L−1) | Proved |
| 95% coverage of the one-sided UCB | **Numerically supported at reduced scale** (28 cells in the audit, lowest Wilson 95% lower bound 0.957; 15 cells re-run with this package: lowest 0.960, mean z 0.18); **not** proved; prediction formula is first-order |
| Full-scale CLT for the alternative (K = 10⁶, L = 10⁴) | Unverified (H₀ block-score kurtosis checked at L = 10⁴) |
| Bias from hardware states correlated with block-varying drift ("class C") | **Counter-example; not removable by centring.** V2.0.08 (protocol §12) *quantifies* its consequence (type-I curve e_I(β), coverage under bias), gives an exact randomization test and a claim gate for the part explained by **monitored** environment variables, and states the unmonitored remainder as an **assumption** (β* = 0.25 by default), not as a result |
| Unconditional exact level 0.05 | **Not claimed.** Exactness holds given assumption A1 (assignment independent of x and device state), which data cannot prove |
| Achievability of E[f̂²] ≥ 0.25, GPU budget, θ-zero alignment across labs | Unverified / hardware requirements |

Reduced-scale results transfer only for dimensionless quantities; compute cost, duration, bias thresholds and full-scale CLT quality do not (protocol §7).

## AI assistance and acknowledgements
Multiple AI assistants (Claude, ChatGPT, Gemini) were used and made to cross-check each other's drafts and derivations; this is **not** independent verification. Conversation logs are not published; `CHANGELOG.md` and `docs/audits/` record what was found and corrected. A friend provided anonymous discussion support.

## Limitations of this package
* The reference implementation covers the protocol's computational core; it is validated only at reduced scale and has not been independently reviewed.
* Randomization-test code re-creates one Philox generator per (block, redraw) — clear but slow; a production implementation for N ≈ 10¹⁰ needs streaming/GPU code (budget in protocol §7(E), CPU estimate only).
* `legacy/` scripts include known errors (see `legacy/README.md`).

## Licence and citation
This repository is multi-licensed. Every file is under exactly one licence (see `LICENSE` and `REUSE.toml`; full texts in `LICENSES/`):

| Licence | Applies to |
|---|---|
| MIT | `code/`, `legacy/`, `.github/`, `.gitignore` |
| CC BY 4.0 | `protocol/`, `docs/`, `results/`, `metadata/`, `CITATION.cff`, `.zenodo.json`, `SHA256SUMS.txt`, and the top-level `README.md`, `CHANGELOG.md`, `VERSIONING.md`, `OSF_REGISTRATION_GUIDE.md` |

GitHub's licence indicator can show only one licence per repository, so it cannot display this split; the protocol and the documents are **not** under the MIT licence. Cite via `CITATION.cff` (GitHub: "Cite this repository").

---
## 中文摘要
本套件包含 SND 預註冊規範 V2.0.12（內容同 V2.0.10，只修正公式在 GitHub 上的顯示與排版）、稽核紀錄（每項主張標示已證明／數值支持／未通過／未驗證）、經測試的參考實作與縮小規模驗證結果。**沒有任何真實實驗資料**；所有數值皆為模擬。核心保證：區塊隨機化檢定的第一類錯誤率在設計假設下精確成立（已證明）；上置信界的 95% 覆蓋率僅有縮小規模的模擬支持，尚未證明；已知反例（C 類硬體偏差）的後果由 §12 量化，可由已監測環境變數解釋的部分以精確檢定與閘門處理，未監測的殘餘只能以假設與實作稽核處理；其他未驗證項目列於上表與稽核檔。
