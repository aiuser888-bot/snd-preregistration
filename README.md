# Project SND — Pre-registration Protocol V2.0.13, Reference Implementation and Validation Record

**Status: pre-registration protocol with a reference implementation and reduced-scale simulation validation. No real experimental data exist or are included.**

**Licences:** code = MIT; protocol, documents, results and metadata = CC BY 4.0 (not MIT). See `LICENSE`.

## What is being tested
Whether a feature `x` that is **sealed first** is statistically dependent on a basis angle `θ` that is **assigned later by a fair random mechanism**. Under the null (θ independent of x) nothing is expected. The protocol specifies (i) a block-randomization test whose validity follows from the randomization design, and (ii) a one-sided amplitude upper confidence bound (Feldman–Cousins ordering of a Rice distribution). The protocol makes **no claim** about physical mechanisms; its interpretation layers (protocol §11) explicitly say that a rejection first calls for checks of the random number generator, assignment code, timestamps and analysis pipeline.

## Contents
| Path | What it is |
|---|---|
| `protocol/SND_protocol_V2.0.13.md` | The pre-registration text (current). Since V2.0.13 §3.1 defines the block angles θ bit for bit (Philox4x64-10, key and counter layout, rejection-sampled Fisher–Yates, serialisation and hash), §4 states how re-draws are computed independently and how a result that floating-point rounding could change is settled by a second tier in exact integer arithmetic (the data and seed stay untouched), §2 seals the centred array f̃ as the bit-level input, and §7(D) fixes the summation order and gives data-driven error bounds. `protocol/history/` has V2.0.01–V2.0.12 |
| `docs/audits/SND_audit_V1.0.08_to_V2.0.13.md` (+ earlier) | Claim-by-claim audit: each claim labelled proved / numerically supported / failed / unverified, with seeds and failure cases |
| `CHANGELOG.md` | Per-version record of errors/contradictions found, fixes, and changed verification methods (replaces publishing AI conversation logs; AI chats are not published) |
| `docs/SND_evolution_for_everyone.md`, `docs/SND_evolution_summary.md` | Plain-language and technical development history |
| `docs/TRACEABILITY.md` | Which protocol section, lemma or step each function implements and which test covers it (generated from the code's `Implements:` / `Verified by:` docstring lines by `code/scripts/check_trace.py`; the protocol has no equation numbers, so sections, lemma items and numbered steps are the reference points) |
| `code/snd/` | Reference implementation: block angle generation (`design.py`, numpy; `philox_ref.py`, the same algorithm with Python integers only, which is the executable form of protocol §3.1), fixed-order sums, error bounds and the cos/sin library check (`numerics.py`), the integer-only reference value of U and U_mon with exact comparisons (`highprec.py`), certified threshold comparisons of the decision rule (`decide.py`), two-way centred statistic, randomization test (`randtest.py`: one generator per (block, re-draw), no global random state, same result for any number of workers, second tier for near ties), FC–Rice table and UCB, exact κ, leakage sensitivity (`leak.py`: type-I curve e_I(β), coverage under bias, environment-projected test U_mon), reduced-scale acceptance simulation |
| `code/tests/test_snd.py` | 40 tests (all pass, see below) |
| `code/crosscheck/philox_theta.c` | A second implementation of the whole protocol §3.1 angle procedure in C (unsigned 64-bit integers, one 128-bit product, binary64), written from the protocol text and sharing no code with the Python implementation; `test_c_implementation_matches_python` compares it bit for bit with the known-answer vectors and with the numpy path (skipped on a machine without `gcc`). Same author as the reference, so it does not replace the third-party independent implementation that protocol 12.5.1 step 5 asks for |
| `code/tests/vectors/theta_kat.json`, `reference_kat.json` | Known-answer vectors for protocol §3.1 (cipher, stream, seed → key, block angles, angle hash, and one value from the ISO C++ `std::philox_engine` specification) and for §7(D)/§2 (fixed-point π, cos/sin, enclosures of U and U_mon, digest of f̃); regenerated and compared by `code/scripts/make_kat.py` |
| `code/scripts/` | `build_table.py` (rebuilds FC–Rice table and checks the protocol's checkpoints), `reproduce_all.py` (headline numbers), `leak_tables.py` (analytic leakage curves → `results/leak_curves.json`), `check_numbers.py` (recomputes and asserts the closed-form numbers quoted in protocol §6/§12), `check_numerics.py` (summation-error measurements against the exact reference value, worst-case bounds and the cost of the reference value, protocol §7(D) → `results/numerics_check.json`), `make_kat.py` (known-answer vectors), `check_trace.py` (protocol ↔ code references, protocol constants and quoted measured numbers, wording hygiene, version label and test count consistent across files), `make_checksums.py` (writes and verifies `SHA256SUMS.txt`) |
| `results/` | `fc_rice_table.npz` (SHA-256 below), `reproduce.json/.log`, `environment.txt`, `numerics_check.json` |
| `LICENSE`, `LICENSES/`, `REUSE.toml` | Licence notice (which licence covers which path), full licence texts (MIT, CC BY 4.0) and the machine-readable per-file mapping |
| `legacy/` | Original exploratory scripts with **known-issue notes** (not the reference implementation) |

## Reproduce
```
pip install numpy scipy
cd code
PYTHONPATH=. python scripts/build_table.py ../results/fc_rice_table.npz   # exit 0 = checkpoints and acceptance mass OK
PYTHONPATH=. python tests/test_snd.py                                    # 40 PASS
PYTHONPATH=. python scripts/leak_tables.py                               # leakage curves (instant)
PYTHONPATH=. python scripts/check_numbers.py                             # exit 0 = protocol's closed-form numbers consistent
PYTHONPATH=. python scripts/check_trace.py --check                       # protocol <-> code references, protocol constants = code constants, version labels and test count consistent, docs/TRACEABILITY.md is current
PYTHONPATH=. python scripts/make_kat.py --check                          # known-answer vectors are current
PYTHONPATH=. python scripts/check_numerics.py                            # exit 0 = summation-error checks and bounds of protocol 7(D) hold
python scripts/make_checksums.py --check                                 # SHA256SUMS.txt lists exactly the files present and every hash matches (also works on an exported tree: --root DIR)
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
| Block angles θ are defined bit for bit (protocol §3.1) and the numpy implementation reproduces the definition | **Verified (tests):** the integer-only reference `philox_ref.py` equals the numpy implementation on 361 (key, block, re-draw, L) cases; the cipher matches numpy's own Philox test data (when installed with numpy; both files were present in the build environment), numpy's generator for 44 further keys and counters including counter carries, and the value fixed by the ISO C++ `std::philox_engine` specification; a runtime check (`stream_conformance`) refuses to run on a numpy whose stream differs; the C implementation `code/crosscheck/philox_theta.c` (written from the protocol text) reproduces the known-answer vectors and 24 random (key, block, re-draw, L) cases of the numpy path bit for bit. **Not done:** the official Random123 known-answer file could not be fetched from the build environment, and a *third-party* independent implementation of the §3.1 procedure (protocol 12.5.1 step 5 asks for one) does not exist yet: the C implementation has the same author as the reference |
| Re-draws do not depend on global state, worker count or order | **Verified (tests):** all fields of the result are identical for 1 and 2 worker processes, for uneven, reordered and one-re-draw-per-chunk partitions of the re-draws, and after poisoning numpy's and Python's global random state; re-draws computed in another process give identical floats |
| U and V are summed in a fixed order; the summation error is far below the tolerance 10⁻⁶√V | Reduced-scale **measurement** (N ≤ 10⁶: ≤ 7.4×10⁻¹⁶ for U, ≤ 1.2×10⁻¹⁶ for V; ≤ 6.5×10⁻¹⁶√V against the exact reference value) and an **analytic worst-case bound** at registration scale (about 1.3×10⁻⁸, assuming the maths library's cos/sin err by at most 2⁻⁵⁰ absolute, which a fixed sample of 3088 angles checks but does not prove); **not measured** at N = 1.2×10¹⁰; cos/sin are library calls and are not claimed to be bit-identical across platforms |
| A result that rounding could change is settled in exact arithmetic (protocol §4, §7(D)) | **Verified (tests) at reduced scale:** the integer reference engine agrees with an independent decimal implementation and GNU `bc` to within 0.55 unit of 2⁻¹²⁸ (measured 0.50) and is pinned by known-answer vectors; genuine near ties on both sides of 0.05 are resolved and give the same result for any worker count; a deliberately corrupted floating-point value is corrected; a comparison still open at 1024 bits is reported as uncertain. **Not measured:** the cost at N = 1.2×10¹⁰ is an extrapolation (about 1.4–2.3 core-days per evaluation, pure Python); the chance that the second tier is triggered under H₀ (about 5×10⁻⁷ per analysis) is an analytic estimate under a Rayleigh idealisation, checked by simulation of the formula, not a property of real data |
| Achievability of E[f̂²] ≥ 0.25, GPU budget, θ-zero alignment across labs | Unverified / hardware requirements |

Reduced-scale results transfer only for dimensionless quantities; compute cost, duration, bias thresholds and full-scale CLT quality do not (protocol §7).

## AI assistance and acknowledgements
Multiple AI assistants (Claude, ChatGPT, Gemini) were used and made to cross-check each other's drafts and derivations; this is **not** independent verification. Conversation logs are not published; `CHANGELOG.md` and `docs/audits/` record what was found and corrected. A friend provided anonymous discussion support.

## Limitations of this package
* The reference implementation covers the protocol's computational core; it is validated only at reduced scale and has not been independently reviewed.
* Randomization-test code re-creates one Philox generator per (block, redraw) — clear but slow; a production implementation for N ≈ 10¹⁰ needs streaming/GPU code (budget in protocol §7(E), CPU estimate only). Such an implementation must reproduce the known-answer vectors of §3.1 and §7(D) and the summation order of §7(D); the second tier of §4 is implemented only in pure Python at reduced scale.
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
本套件包含 SND 預註冊規範 V2.0.13（相對 V2.0.12：把區塊角度 θ 的產生方式定義到位元層級並附已知答案向量與純整數參考實作，重抽的獨立計算、數值近平手的二級裁決（不改種子與資料，以整數精確算術重算而決定）、封存 f̃ 作為位元級輸入，以及求和順序與誤差界的固定；統計方法與數字不變）、稽核紀錄（每項主張標示已證明／數值支持／未通過／未驗證）、經測試的參考實作與縮小規模驗證結果。**沒有任何真實實驗資料**；所有數值皆為模擬。核心保證：區塊隨機化檢定的第一類錯誤率在設計假設下精確成立（已證明）；上置信界的 95% 覆蓋率僅有縮小規模的模擬支持，尚未證明；已知反例（C 類硬體偏差）的後果由 §12 量化，可由已監測環境變數解釋的部分以精確檢定與閘門處理，未監測的殘餘只能以假設與實作稽核處理；其他未驗證項目列於上表與稽核檔。
