# Uploading to OSF and Zenodo — checklist and pre-filled text

## Which platform for what
* **Zenodo**: permanent archive with a DOI; best for the whole package (protocol + code + results). Deposits are immutable once published (you can publish new *versions*).
* **OSF Registries ("OSF Preregistration" template)**: a time-stamped, frozen registration of the *plan*. Registrations cannot be edited after submission (only withdrawn). Best for the protocol text itself.
* Typical combination: Zenodo for the package (get the DOI) → OSF project containing the same files → OSF registration whose "Other" field cites the Zenodo DOI.

## Before you upload (things only you can fill in)
1. Author: currently the pen name `Pietersite` in `CITATION.cff`, `.zenodo.json`, `LICENSES/MIT.txt`, `REUSE.toml`. Optionally add a real name, ORCID and affiliation there before a Zenodo/OSF deposit.
2. Confirm licences (MIT for code; CC BY 4.0 for protocol, documents, results and metadata — see `LICENSE` and `REUSE.toml`) or choose others. If you change them, update `LICENSE`, `LICENSES/`, `REUSE.toml`, `CITATION.cff` and the Zenodo `license` field together, then run `code/scripts/check_licenses.py`.
3. Decide **embargo**: Zenodo allows an embargo date; OSF registrations can be embargoed (check the current maximum on OSF).
4. Decide whether to publish V2.0.11 as is or after the open items are addressed. A registration is permanent; an honest version history is possible via new Zenodo versions.
5. Confirm the §12.5 seed-acquisition choices: option A (one seed per device, drawn once after the data freeze T_F), raw QRNG output ≥ 512 bits with a conservative min-entropy estimate × length ≥ 160 bits, the external time-stamp provider for H_F, and that θ is an analysis-time label (§12.5.5). If θ must be physically set during data collection, §12.5.1 does not apply.

## Zenodo (manual upload)
1. Log in at zenodo.org → *New upload*.
2. Upload `SND_release_V2.0.11.zip` (or the files individually). Zenodo shows each file's MD5; compare with `SHA256SUMS.txt`/the MD5 list if you wish.
3. Metadata (copy from `metadata/zenodo.json`): Resource type *Publication → Report* (choose another type if you prefer), title, authors, description, keywords, licence, language, version `2.0.11`.
4. Related identifiers (optional): later add OSF project URL as "is supplemented by".
5. *Save* → check preview → *Publish* (a DOI is assigned; publishing cannot be undone, files can only be superseded by a new version).
Alternative: put `snd_release/` in a GitHub repository, enable the repository in Zenodo's GitHub integration, and create a release; Zenodo reads `.zenodo.json`/`CITATION.cff`.

## OSF Preregistration — suggested answers (English, adapt as needed)
**Title**: Project SND: block-randomization test and one-sided amplitude upper confidence bound for dependence between sealed features and later random basis angles (protocol V2.0.11).

**Description / hypotheses**: H0: the later randomly assigned angle θ is statistically independent of the sealed feature x. Alternative (modelling device only): p1(x|θ;a,b) = p0(x)[1 + f(x)(a cosθ + b sinθ)], r = |a+ib| ≤ 1e-3. Primary decision: randomization p-value ≤ 0.05. Secondary: 95% one-sided upper confidence bound on r; exclusion of r ≥ 1e-4 if the bound is below 1e-4 and p > 0.05.

**Design**: Blocks of L = 10^4 units; angles are a randomly rotated equispaced set, uniformly permuted within each block (counter-based Philox, QRNG supplies seed only). Angles are generated after features are sealed and centred (two-way centring using x only).

**Sampling / stopping**: Collect blocks until the information V ≥ V_target = 1.5e9 (about N ≈ 1.2e10 units if E[f̂²] = 0.25; about 11.4 months with 4 devices at 10 ms/unit). Stopping depends on x only. Feasibility gate: do not start if predicted N > 3e10.

**Variables**: x (sealed scalar feature from the frozen operator T), f̂ = frozen feature perturbation function (sup|f̂| ≤ 1), θ (angle). No outcome other than U, V defined in protocol §2.

**Analysis plan**: Statistic Λ = |U|²/V with U = Σ f̃ e^{iθ}; randomization test with M = 1999 redraws (level exactly 0.05); UCB95 = V^{-1/2}·u(|U|/√V) with the Feldman–Cousins–Rice table (checkpoints in §5.1). Invalid-analysis rules: audit-chain failure, incomplete angle records, design-check or uniformity-test failure, failed pre-validation (§7), numerical-accuracy failure.

**Exclusion / quality**: whole-block exclusion for hardware flags (decided before angles are generated).

**Other**: Reproduction criterion (protocol §10). Known limitations: see README and audit (class-C bias counter-example; coverage proved only numerically at reduced scale). Package DOI: <Zenodo DOI after publishing>.

**Important honesty note for the "Data" questions**: no data have been collected; no one has accessed real data; the only numbers are simulations.

## After publishing
Record the DOI and OSF URL in `README.md` of a new version; do not change published files — create a new version instead.
