# 規範與程式碼對照表（SND_protocol 與 code/）

本表由 `code/scripts/check_trace.py --write` 自 docstring 產生，請勿手改；`--check` 會比對本檔。程式碼的註解以 `Implements:` 標出所實作的規範章節（§ 為規範標題編號；`7(D)` 為 §7 的 (D) 項；`Lemma 5(c)` 為引理 5 的 (c) 項），`Verified by:` 標出測試。規範本身沒有公式編號。

## 規範項目 → 程式碼

| 規範項目 | 程式碼 | 測試 |
|---|---|---|
| 0 | `snd/kappa.py::rho_x_bernoulli_copula` | — |
| 2 | `snd/numerics.py::u_terms` | — |
| 2 | `snd/numerics.py::u_stat` | — |
| 2 | `snd/numerics.py::v_stat` | — |
| 2 | `snd/stats.py` | test_lemma1_variance, test_ucb95_and_scores, test_ftilde_digest_and_seal |
| 2 | `snd/stats.py::center_two_way` | — |
| 2 | `snd/stats.py::block_V` | — |
| 2 | `snd/stats.py::block_scores` | — |
| 2 | `snd/stats.py::lam` | — |
| 2 | `snd/stats.py::ftilde_digest` | — |
| 3.1 | `snd/design.py` | test_theta_matches_pure_python_reference |
| 3.1 | `snd/design.py::seed_material` | — |
| 3.1 | `snd/design.py::seed_to_key` | — |
| 3.1 | `snd/design.py::make_block_rng` | — |
| 3.1 | `snd/design.py::block_angles` | — |
| 3.1 | `snd/design.py::theta_rounding_bound` | — |
| 3.1 | `snd/philox_ref.py` | test_philox_cipher_known_answers, test_theta_matches_pure_python_reference, test_known_answer_vectors_file_is_current |
| 3.1 | `snd/philox_ref.py::philox4x64_10` | — |
| 3.1 | `snd/philox_ref.py::PhiloxStream` | — |
| 3.1 | `snd/philox_ref.py::key_from_digest` | — |
| 3.1 | `snd/philox_ref.py::block_theta` | — |
| 3.1 | `snd/philox_ref.py::theta_bytes` | — |
| 3.1 | `snd/randtest.py::draw_theta` | — |
| 3.1 | `snd/stats.py::ftilde_digest` | — |
| 3.1 | `scripts/check_trace.py` | test_trace_labels_parse |
| 3.1 | `scripts/check_trace.py::version_checks` | — |
| 3.1 | `scripts/make_kat.py` | test_known_answer_vectors_file_is_current, test_reference_vectors_file_is_current |
| 3.2 | `snd/design.py::design_check` | — |
| 4 | `snd/highprec.py::modulus_sq_bounds` | — |
| 4 | `snd/highprec.py::compare_moduli` | — |
| 4 | `snd/highprec.py::classify` | — |
| 4 | `snd/randtest.py` | test_randtest_size_under_H0, test_randtest_independent_of_workers_and_order, test_near_tie_resolution |
| 4 | `snd/randtest.py::abs_U_of_draw` | — |
| 4 | `snd/randtest.py::tail_counts` | — |
| 4 | `snd/randtest.py::split_ms` | — |
| 4 | `snd/randtest.py::randomization_test` | — |
| 5.1 | `snd/fcrice.py` | test_table_checkpoints, test_fc_rice_table_rebuilds_identically, test_ucb95_and_scores |
| 5.1 | `snd/fcrice.py::build_table` | — |
| 5.1 | `snd/fcrice.py::load_table` | — |
| 5.1 | `snd/fcrice.py::u_of_x` | — |
| 5.2 | `snd/fcrice.py::ucb95` | — |
| 7(B) | `snd/accept.py` | test_coverage_iid_and_serial |
| 7(B) | `snd/accept.py::wilson` | — |
| 7(B) | `snd/accept.py::simulate` | — |
| 7(D) | `snd/highprec.py` | test_highprec_sincos_against_independent_references, test_highprec_reference_sum_and_comparison, test_near_tie_resolution |
| 7(D) | `snd/highprec.py::pi_fixed` | — |
| 7(D) | `snd/highprec.py::sincos_fixed` | — |
| 7(D) | `snd/highprec.py::exact_weights` | — |
| 7(D) | `snd/highprec.py::reference_U` | — |
| 7(D) | `snd/numerics.py` | test_numerics_fixed_order_and_sum_checks, test_error_bound_and_library_qualification, test_near_tie_resolution |
| 7(D) | `snd/numerics.py::tree_sum` | — |
| 7(D) | `snd/numerics.py::kahan_sum` | — |
| 7(D) | `snd/numerics.py::dd_sum` | — |
| 7(D) | `snd/numerics.py::exact_sum` | — |
| 7(D) | `snd/numerics.py::sum_check` | — |
| 7(D) | `snd/numerics.py::class_plan` | — |
| 7(D) | `snd/numerics.py::error_bound_formula` | — |
| 7(D) | `snd/numerics.py::u_error_bound` | — |
| 7(D) | `snd/numerics.py::libm_conformance` | — |
| 7(D) | `snd/stats.py::block_scores` | — |
| 7(D) | `scripts/check_numerics.py` | test_numerics_fixed_order_and_sum_checks, test_error_bound_and_library_qualification |
| 10 | `snd/decide.py` | test_decision_comparisons_certified |
| 10 | `snd/decide.py::exclusion` | — |
| 10 | `snd/decide.py::reproduction_ii` | — |
| 12.1 | `snd/leak.py` | test_leak_curves, test_leak_cauchy_schwarz_and_exact_bias, test_u_mon_exact_size_and_unbiasedness, test_lemma5_* |
| 12.1 | `snd/leak.py::type1_of_beta` | — |
| 12.1 | `snd/leak.py::beta_of_B` | — |
| 12.1 | `snd/leak.py::beta_max` | — |
| 12.1 | `snd/leak.py::q_rms_for_beta` | — |
| 12.1 | `snd/leak.py::coverage_worst_phase` | — |
| 12.1 | `snd/leak.py::coverage_min_over_nu` | — |
| 12.2 | `snd/leak.py::u_mon` | — |
| 12.2 | `snd/numerics.py::u_mon_stat` | — |
| 12.2 | `snd/randtest.py::randomization_test` | — |
| 12.5.1 | `snd/design.py::theta_digest` | — |
| 12.5.1 | `snd/design.py::stream_conformance` | — |
| 12.5.1 | `snd/leak.py::seed_entropy_ok` | — |
| 12.5.3 | `snd/leak.py::chi2_divergence` | — |
| 12.5.3 | `snd/leak.py::beta_ub_from_D` | — |
| 12.5.3 | `snd/leak.py::size_ub_from_D` | — |
| 12.5.3 | `snd/leak.py::devices_D` | — |
| 12.5.3 | `snd/leak.py::hash_D_bound` | — |
| 12.5.3 | `snd/leak.py::best_of_k_size` | — |
| 12.5.3 | `snd/leak.py::retry_size_ub` | — |
| Lemma 2 | `snd/kappa.py` | test_kappa_values, test_coverage_iid_and_serial |
| Lemma 2 | `snd/kappa.py::rho_x_bernoulli_copula` | — |
| Lemma 2(ii) | `snd/kappa.py::kappa_exact` | — |
| Lemma 2(iv) | `snd/kappa.py::coverage_prediction` | — |
| Lemma 3 | `snd/leak.py` | test_leak_curves, test_leak_cauchy_schwarz_and_exact_bias, test_u_mon_exact_size_and_unbiasedness, test_lemma5_* |
| Lemma 3(a) | `snd/leak.py::beta_of_B` | — |
| Lemma 3(a) | `snd/leak.py::beta_max` | — |
| Lemma 3(b) | `snd/leak.py::type1_of_beta` | — |
| Lemma 3(c) | `snd/leak.py::coverage_worst_phase` | — |
| Lemma 3(c) | `snd/leak.py::coverage_min_over_nu` | — |
| Lemma 4 | `snd/leak.py` | test_leak_curves, test_leak_cauchy_schwarz_and_exact_bias, test_u_mon_exact_size_and_unbiasedness, test_lemma5_* |
| Lemma 4 | `snd/leak.py::u_mon` | — |
| Lemma 5 | `snd/design.py::seed_to_key` | — |
| Lemma 5 | `snd/leak.py` | test_leak_curves, test_leak_cauchy_schwarz_and_exact_bias, test_u_mon_exact_size_and_unbiasedness, test_lemma5_* |
| Lemma 5(b) | `snd/leak.py::beta_ub_from_D` | — |
| Lemma 5(c) | `snd/leak.py::size_ub_from_D` | — |
| Lemma 5(d) | `snd/leak.py::best_of_k_size` | — |
| Lemma 5(d′) | `snd/leak.py::retry_size_ub` | — |
| Lemma 5(e) | `snd/leak.py::devices_D` | — |
| Lemma 5(e)(f) | `snd/leak.py::hash_D_bound` | — |
| Lemma 5(f) | `snd/leak.py::seed_entropy_ok` | — |

## 沒有程式碼對應的規範章節

下列章節沒有被任何 `Implements:` 引用（敘述性、程序性或硬體章節，或其內容由上列相鄰項目涵蓋）：§1、§3、§3.3、§5、§6、§8、§9、§11、§12、§12.3、§12.4、§12.5、§12.5.2、§12.5.4、§12.5.5、§12.6、§12.7。

