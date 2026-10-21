# legacy/ — original exploratory scripts (NOT the reference implementation)

These are the working scripts used while developing the protocol (Oct 14–15, 2026). They are kept **for transparency and provenance only**. Use `code/` for the reference implementation. Several scripts contain **known errors**; do not reuse them.

| Script | Purpose | Status |
|---|---|---|
| fc.py | First FC–Rice table (nu ≤ 13) | Superseded: grid distorts u(x) for x > 11 (fixed in r3.py / `code/`) |
| r3.py | Extended FC table (current grid) + early bias test | Table: valid (reproduced by `code/scripts/build_table.py`). Bias test: **invalid** (leakage too weak) |
| a.py, sim.py, chk.py | Early counter-example / size / coverage checks (reduced scale) | Exploratory |
| r2.py | Bias × drift test | **Invalid** (leakage uniform/too weak) |
| r4.py | Bias × drift test | **Bug**: last line re-centres rows in all modes, so "raw" is actually row-centred |
| r5.py, r6.py | Coverage / phase tests | r6 inadequate (row-centring only, one nu); r5 exploratory |
| r7.py, r7b.py | Phase dependence of power and coverage (32 cells) | Valid (non-detection of phase dependence; not a proof) |
| r8.py | Centring × bias × drift factorial (A–F) | Valid; documents the C-class counter-example |
| r9.py, r10.py, r10b.py | Follow-up on a type-I excess (noise); Lemma 1 conditional variance; heavy-tail recheck | Valid |
| r11.py, r12.py | Reproduction-criterion success probabilities | Valid (Gaussian large-K model) |
| r13.py | H1 under serial correlation | **Superseded**: original AR(1) start-up was inconsistent (corr(Z0,Z1)=0.668 for rho=0.9); also used the old nu_eff convention. Not cited in v15 |
| r14.py | Closed-form kappa vs simulation | Valid, but compares with estimated sigma^2 (differences up to 0.005); see r16.py for true sigma^2 |
| r15.py | Lemma 2 direct variance / cross-block covariance check | Valid |
| r16.py | AR(1) start-up effect; exact kappa vs simulation (incl. negative rho) | Valid |
| r17.py (+ r17.out) | Coverage grid, rho in {-0.9..0.9}, nu in {2,3,4} | Valid after fixing variable shadowing; its printed "pred" column uses the **old** nu_eff convention (recomputed with the corrected form in v15) |
| r18.py, r19.py, r20.py | Small K/L; one conditional-tilt model; cross-block shared factors; CLT kurtosis at L=10^4 | Valid; "pred" columns use the old convention where printed |

| r22.py (+ r22.out) | V2.0.09: net (gated) exclusion probability/power, reproduction success probabilities and duration for V_target = 1.375e9 vs 1.5e9 | Valid (Gaussian large-K model; reproduces the earlier 0.888/0.956 with new seeds) |
| r23.py (+ r23.out) | V2.0.10: exact checks of Lemma 5 (design moments, β ≤ √D tight for an aligned tilt, TV bound, GF(2) hash identity) | Valid; superseded by `code/tests/test_snd.py::test_lemma5_*` |
| r21.py, r21b.py, r21c.py, r21d.py (+ .out files) | V2.0.08 leakage study: e_I(β) vs observed, coverage under bias, U_mon size/power, gated claims, η formula | Valid. Printed labels "delta" = β (normalized bias) in V2.0.08 notation. `r21c_part1_*.out` were produced with an extra (since removed) "ub≤0.5" column that is unattainable (see CHANGELOG V2.0.08); its other columns are valid |

Seeds are inside each script (`np.random.default_rng([...])`). Outputs were not all saved to files; `r17.out`, `r18a.out`, `r18b.out` are the saved ones. The numbers reported in the protocol's evidence section (v15) were recomputed from these results with the corrected prediction formula.
