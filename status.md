# Where we are, and how to pick this up

First written 17 October 2026, last updated 20 October 2026. This is the short hand-off. The long working record is `02_current_status.md`
(evidence base, extraction detail), and the design decisions are in `phase3_decisions.md`.

The project asks whether adding host genetic variants to a validated clinical score (Stanzani 2019) for invasive
mould disease (IMD) would change clinical decisions. It is a simulation study, with decision-curve net benefit at
the 5% threshold as the primary endpoint.

## Where the code lives, and the website

The repository is **lewis-id-lab/PROSE** (public since 19 October 2026), set as `origin`; the earlier personal copy, Russlewisbo/PROSE, is the remote `personal` and is no
longer updated. Push both `main` and the pre-specification tag (`git push origin main prespec-lock`) from this computer if either has changed.
`git status` should show no spreadsheet before you commit. The per-draw simulation files of the genetic layer
(`output/genetic_*_draws.csv`, 2 to 9 MB each) are gitignored and regenerate exactly from scripts 16 to 18.

**The website** is the Quarto site in `docs/`, live at https://lewis-id-lab.github.io/PROSE/ and linked from https://lewis-id-lab.github.io/. To update it: edit the `.qmd` pages (or the decision record and the status log,
which are included in the Pre-specified decisions and Status log pages), run `quarto render` from the project root, and commit `docs/` and `_freeze/` with the source. The workflow `.github/workflows/pages.yml`
publishes `docs/` to GitHub Pages on every push to `main` that changes it. Pages is enabled with the workflow build; to switch the repository back to private, which would take the site down on the
free plan, change its visibility in Settings.

## Starting on another computer

```bash
git clone https://github.com/lewis-id-lab/PROSE.git
cd PROSE
git fetch --tags        # needs the tag pushed first
```

You need R 4.6.1 and Quarto 1.8 or later (built with 1.10.18). Package versions are pinned in `renv.lock`
(80 CRAN packages). Start R in the project root. The `.Rprofile` activates renv, which installs itself on first
use if it is missing. Then run:

```r
renv::restore()
```

On Linux several packages compile from source (`mvtnorm` needs a Fortran compiler, and some others need system
libraries), which Posit Public Package Manager binaries avoid. After installing or updating any package, run
`renv::snapshot()` and commit `renv.lock`.

Check everything works, from the project root:

```r
source("tests/run_tests.R")     # the full suite, all passing at the last run
                                # (some tests skip if the gitignored data files are absent)
```

```bash
quarto render                    # rebuilds the site into docs/
```

All fitted results are committed in `output/`, so nothing needs re-running to continue. If you do re-run the
scripts, run them from the project root, in this order:

| script | what it does | time on 15 cores |
|---|---|---|
| `scripts/01_fit_correlation.R` | fits three correlation structures, writes the fitted matrix, profile table and comparison | about 30 min |
| `scripts/02_correlation_ensemble.R` | six-chain sampler for near-optimal matrices (150 members) | about 3 min |
| `scripts/03_decision_curves.R` | exact decision curves, ensemble bands, figures | under 1 min |
| `scripts/04_haplotype_frequencies.R` | phased 1000G haplotypes, rescaled to gnomAD NFE | seconds (genotypes cached; set `REFRESH <- TRUE` to refetch, needs network) |
| `scripts/05_check_bs_data.R` | checks `bs.xlsx` against the published counts and coefficients; needs the gitignored file | seconds |
| `scripts/06_pseudonymise_boscore.R` | builds the pseudonymised extract of the full BOSCORE file into `data/restricted/` (gitignored) | seconds |
| `scripts/08_reconcile_boscore.R` | reconciles `boscore.xlsx` with `bs.xlsx` (predictor mapping, inclusion rule, the 49 twin pairs), stops if any check fails, writes the pseudonymised base cohort | seconds |
| `scripts/09_baseline_cohort.R` | baseline decision curve and MCID on the base cohort (4,078 unique admissions) with a patient-level bootstrap, plus sensitivity analyses | about 10 s |
| `scripts/15_map_cd209_il10_ifng_evidence.R` | CD209 (G at rs4804800), IL10 and IFNG mapping and pooling, with a per-term hazard-ratio fallback; run after 13 | seconds |
| `scripts/16_genetic_layer.R` | first pass of the genetic layer: four terms on the refit baseline, 2,000 posterior draws under two readings of the truth, gain over the baseline against the MCID, check against simulation; run after 11, 12, 14 and 15 | about 70 s |
| `scripts/17_genetic_layer_extensions.R` | genetic layer extensions: the 150-member correlation ensemble, Shapley attribution to the four terms, sampling variability of a 4,078-admission study; run after 16 (builds the cached ensemble profile distributions on first run) | about 2.5 min on 15 cores |
| `scripts/18_genetic_ensemble_sampling.R` | sampling variability of a study with the clinical mix drawn from a copula member, plus controls; run after 17 | about 35 s |
| `scripts/19_genetic_sensitivity_cohorts.R` | the genetic layer (exact, attribution, sampling) on the published equation (4,078 and 4,127 rows) and the refit on 4,127 rows; set `SENS_ONLY` to run one setting at a time | about 1.5 min each |
| `scripts/20_genetic_donor_stratum.R` | donor genotype in the allo-HSCT stratum (D1 secondary analysis): the stratum from the cohort's allo flag (and, as a sensitivity analysis, from Table 1 margins), donor effects, gain against the stratum and whole-cohort MCIDs, sampling variability; needs the donor draws from scripts 12 and 14 | about 1 min |
| `scripts/21_aim3b_neutropenia.R` | Aim 3b: the genetic layer split by neutropenia status (six effect-modification scenarios, five scores), Brunel's stratified effects, and the interaction test in simulated studies; run after 16; with `AIM3B_SECOND_CLEC7A <- TRUE` it runs with five terms (second CLEC7A parameter) and writes `_second_clec7a` files | about 3.5 min |
| `scripts/22_clec7a_second_parameter.R` | the second CLEC7A parameter (effect of the A-G haplotype): rs7309123 rows mapped through the three-haplotype model using the strand resolution in `data/rs7309123_strand_resolution.csv`, pooled, and the layer with two CLEC7A parameters; run after 14 and 16 | about 1 min |
| `scripts/23_aim4_sweep.R` | Aim 4 sweep: 1 to 20 loci x OR 1.1 to 3.0 x risk-allele frequency 0.05 to 0.50, exact convolution, gain against the MCID at 5% and 20%, frontiers, published-equation and correlation-ensemble robustness, two figures; run after 11 and 16 | under 1 min |
| `scripts/24_attribution_second_clec7a.R` | Shapley attribution of the layer's gain to five terms (the two CLEC7A parameters separately); run after 16 and 22 | about 3 min |
| `scripts/25_aim4_realised_sweep.R` | Aim 4 sweep with a realised score: estimated weights with per-locus error sigma 0.05 to 0.40, 12 replicate discovery studies per architecture, frontiers and a figure; run after 23 | about 13 min |
| `scripts/26_five_term_robustness.R` | the five-term Aim 3b layer and Shapley attribution on the sensitivity cohorts and the correlation ensemble (set `ROBUST_ONLY` to run one setting at a time); run after 16, 22 and 24 | about 5 min per setting |
| `scripts/27_writeup_figures.R` | figures for the write-up: the headline gain against the MCID, the Aim 3b gain and power, and the attribution; run after 16, 21, 22 and 24 | seconds |
| `scripts/28_east_asian_haplotypes.R` | phased 1000 Genomes East Asian CLEC7A haplotype frequencies and the rs3901533 linkage cross-tab (genotypes cached; set `REFRESH <- TRUE` to refetch, needs network) | seconds |
| `scripts/29_nri_comparison.R` | tertiary NRI comparison (addendum 22): the haplotype panel and the PTX3 h2 term against the predecessor's three-SNP independent model, with net benefit and aROC under the same truth, and an effect grid; run after 16 | about 1 min |
| `scripts/30_prs_benchmark.R` | polygenic score benchmark (addendum 23): how strong a polygenic score, and how large its discovery study, would have to be to reach the MCID; run after 08 and 11 | not timed |
| `scripts/31_baseline_calibration_check.R` | calibration check (addendum 24): the panels, the Aim 4 exact sweep and the polygenic score benchmark under three comparisons (current, marginal, profile); stops unless the current comparison reproduces the committed results (the polygenic accuracy needed to 0.01 only; addendum 23); run after 16, 22, 23 and 30 | about 10 min | about 6 min |
| `scripts/32_calibration_check_realised.R` | the realised-score Aim 4 sweep under the three comparisons, with the oracle frontier; stops unless the current comparison reproduces script 25; run after 25 and 31 | about 23 min |
| `scripts/33_calibration_check_nri.R` | the NRI comparison and effect grid under the three comparisons; stops unless the current comparison reproduces script 29; run after 29 | about 1 min |
| `scripts/34_calibration_check_figures.R` | the headline and Aim 4 figures under the profile-preserving truth; run after 31 | seconds |
| *Scripts 17 to 21, 24, 26 and 27* | *each takes `PROSE_SCHEME=profile` (set with `Sys.setenv()` or `withr::with_envvar()` before sourcing) to run under the profile-preserving truth (addendum 25); they then write `_profile` files and leave the global ones alone. Times as in the global runs: 17 about 3 min, 18 about 2 min, 19 about 6 min, 20 about 3 min, 21 about 3 min for each of its two versions, 24 about 5 min, 26 about 6 min per setting* | |
| `scripts/37_aim3b_scheme_comparison_extended.R` | script 36 for the realised score and the five-term panel, with the shared gain functions in `R/scheme_comparison.R`; stops unless the four-term ceiling rows reproduce script 36 (addendum 28); run after 22 and 36 | under 1 min |
| `scripts/36_aim3b_scheme_comparison.R` | Aim 3b under the global and the profile-preserving truth at effects up to 4.9 times the pooled one, by where the effect sits and for the six scenarios (addendum 27); reads the pooled estimates and the committed draws | under 1 min |
| `scripts/35_calibration_check_headline.R` | the compact comparison table (current, marginal, profile) that the manuscript's `tbl-calibration` reads; run after 31, 32 and 33 | seconds |
| `scripts/14_map_clec7a_tlr4_evidence.R` | CLEC7A (Y238X haplotype) and TLR4 (G-T haplotype, carrier model) mapping and pooling; run after 13 | seconds |
| `scripts/13_verify_ptx3_rows.R` | checks rows against their sources (Herrero-Sanchez, Brunel), writes the corrections file and the verification log; run before 12 and 14 | seconds |
| `scripts/12_map_ptx3_evidence.R` | PTX3 pilot of the evidence-to-parameter mapping: rows to per-copy h2 effects, one effect per study, Bayesian pooling | seconds |
| `scripts/11_refit_baseline.R` | refits the baseline on the base cohort, with robust SE, patient bootstrap, optimism and the MCID; writes the fixed baseline parameters | about 1 min |
| `scripts/10_outcome_risk_sensitivity.R` | MCID under four readings of the outcome risk (published, refit, refit rule, observed), refit repeated inside a patient bootstrap | about 40 s |
| `scripts/07_empirical_baseline.R` | baseline decision curve and MCID on the empirical cohort, with bootstrap intervals and the comparison with the copula | under 1 min |

Script 01 scales with core count and would take hours on 4 cores. Cut `n_starts` if you need to.

## Where the project stands

| Phase | State |
|---|---|
| 0 Set-up, 1 Literature search | complete (831 unique records screened) |
| 2 Screening and extraction | round 1 complete; the user reports the human round-2 checks confirmed the extraction (record not in the repo); a direction check of the text-parsed rows is still to do |
| 3 Parameter derivation | **done provisionally** for PTX3, CLEC7A (two parameters), TLR4 and CD209 (addenda 7 to 9 and 16); IL10 and IFNG dropped (D14) |
| 4 Baseline | **done**; the base dataset is `boscore.xlsx` with the 49 twin admissions removed (D9), and the baseline model is a refit on it (D10, addendum 6): 4,078 admissions in 1,945 patients |
| 5 Genetic layer simulation | **first results in** (addenda 11 to 22 and 24 to 28): exact layer, sampling variability, correlation ensemble, attribution, sensitivity cohorts, donor genotype in the allo-HSCT stratum, Aim 3b, the second CLEC7A parameter, the East Asian check (addendum 21), the tertiary NRI comparison (addendum 22) and the calibration check (addendum 24, which made a profile-preserving truth primary) and the secondary analyses and the all-ancestry CLEC7A variant rerun under it (addenda 25 and 26), so every analysis is on it |
| 6 Sensitivity and architecture sweep | **Aim 4 sweep done** (addenda 16 and 17); sensitivity analyses done. Not done: the pathogen-specificity scenario, and the Aim 3b interaction test on the sensitivity cohorts |
| 7 Write-up | **full first draft** in `manuscript.qmd` (the Draft manuscript page): abstract, introduction, methods, results, discussion; references in `references.bib`, Vancouver style (`vancouver.csl`). Still to add: authors, ethics, data availability, funding, and two references (IMD burden; Feys et al. 2026, not yet in Crossref) |

### Evidence base (provisional: the human round-2 checks are reported complete but not recorded here; one study's direction was corrected from its paper on 19 October)

- 81 records included, 63 extracted, 2,008 variant rows. 970 are usable (823 estimates, 115 complete 2x3 genotype
  tables, 32 complete 2x2 tables).
- 13 variants reach three or more meta-analysable studies. Six reach four or more independent cohorts:
  rs16910526, rs2305619, rs7309123, rs3816527, rs4986790, rs4986791.
- Allele orientation: the 10 strand-flip cells are handled; the strand of rs7309123 is resolved for six of ten studies (grades A and B, including Kalkanci 2020 and Chen 2019 after their full texts were read, addendum 19), weakly for three and undetermined for one (Zhao 2022)
  (`data/rs7309123_strand_resolution.csv`).

### Pre-specification, locked at tag `prespec-lock` (commit `90e3656`)

| Decision | Settled |
|---|---|
| Primary genotype | recipient; donor is a secondary analysis in the allogeneic HSCT stratum only |
| Panel | four haplotype terms (PTX3, CLEC7A, TLR4, CD209); IL10 rs1800896 and IFNG rs2069705 dropped on 19 October 2026 (addendum 10, D14) |
| Dropped / excluded | TLR9 rs352140 dropped; Herrero-Sanchez 2018 excluded from the primary analysis |
| TLR4 pooling | carrier model with a continuity correction |
| Genetic term | recentred so prevalence stays at 3.3% |
| Baseline aROC target | 0.85 |
| Haplotype phase | coupling in every pair |
| Baseline cohort (D9, addendum 3) | **the 4,078 unique admissions of `boscore.xlsx`** (twins removed, event rows kept), patients resampled in the bootstrap; the 4,127 rows and the version 2 copula (D7) are sensitivity analyses |
| Baseline model (D10, addendum 6) | **maximum-likelihood refit of the seven predictors on the 4,078 admissions** (`output/baseline_refit_parameters.csv`); the published equation is a sensitivity analysis |
| MCID (20% more true cases, no extra false positives) | **+0.00458 net benefit at 5% (primary)**, +0.00160 at 20% (secondary), on the refit baseline. Sensitivity thresholds: +0.00404 (published equation) and +0.00445 (4,127 rows) |

`phase3_decisions.md` has a dated addendum written after the lock. Further changes there should be dated addenda,
not edits.

### Baseline model (Aim 2) and decision curves

- The Stanzani equation is fixed and never refitted. The joint distribution of its seven predictors is not
  published, so a Gaussian-copula correlation matrix is a fitted assumption.
- (Now a sensitivity analysis, see D8.) Version 2 fits all 16 published targets (chi-square 0.36, largest miss 0.4 standard errors). The fit is
  saturated: about nine independent constraints for nine parameters, so the good fit does not validate the
  structure. Malignancy comes out nearly uncorrelated with the rest of its block (0.04).
- Baseline net benefit at 5% is 0.01478 (14.7% flagged). The ensemble 90% range for the 5% MCID is 0.00397 to
  0.00460.
- Decision curves are exact functions of the 128 clinical profiles, so they are step functions. The 5% MCID sits on
  a plateau from about 4.2% to 7.8%.

### Admission-level data (`bs.xlsx`, gitignored patient-level data)

4,127 admissions and the seven predictors plus outcome (133 cases): the derivation dataset. `scripts/05_check_bs_data.R`
checks it. It reproduces Table 2 except two non-case counts that are one lower than published (high-risk
chemotherapy 776 against 777, prolonged neutropenia 1,097 against 1,098), a refit recovers the published
coefficients to within 0.006, and the apparent aROC is 0.853.

- `id` is the row number of `boscore.xlsx` (next section), which carries the patient link, so admissions per patient
  are recoverable.
- Only 75 of the 128 predictor profiles occur (35 with fewer than 5 admissions). The fitted copula is still about 10%
  of probability mass away from the real joint, and some of its correlations are wrong in sign (for example
  malignancy by CMV is -0.15 in the data against +0.23 fitted), so the saturated fit did not validate the structure.
- The MCID on the real joint is +0.00445 at 5% (bootstrap 90% interval 0.00412 to 0.00478) and +0.00142 at 20%
  (0.00116 to 0.00170). The values locked at the tag (+0.00428, +0.00140) are inside those intervals.
- **Decision made earlier (D8, addendum 2; superseded by D9):** the baseline cohort for the genetic layer is the empirical
  128-profile distribution, with the copula kept as a sensitivity analysis. Read it with `read_bs_data()` and
  `empirical_profile_counts()`; the profile counts themselves are not saved, because several profiles have fewer than 5 admissions.

### `boscore.xlsx`: the parent of `bs.xlsx` (gitignored, contains direct identifiers)

4,694 rows by 119 columns, ids 1 to 4,694, admissions 2007 to 2016, 133 EORTC proven or probable cases, and the
identifiers (`IDRICO`, `IDIND`, name, birth date). Its `id` is the `id` in `bs.xlsx`, so it is the key to patients. This is
the file that reconciles everything:

- **The seven `bs.xlsx` predictors are the `BOSCORE_*` columns, exactly** (all 4,127 rows agree): `BOSCORE_4` prior IA, `_6`
  steroids, `_8` active malignancy, `_9` high-risk chemotherapy, `_10` neutropenia, `_11` lymphopenia, `_15` CMV. The
  outcome `EORTC_PROV_PROB` matches too, with 133 cases in both. `pred` is the published equation (largest difference
  0.006), and `[pred >= 0.05]` flags 661 of the 4,127 (16.0%).
- **The inclusion rule:** the 4,127 are exactly the rows with complete data on the seven `BOSCORE_*` items. All 567
  excluded rows have missing items, and `pred` is missing for the same rows.
- **Patients:** 1,945 in the 4,127 (the paper says 1,944), 2.122 admissions per patient, 985 patients (50.6%) with one
  admission, at most 11. The 133 events fall in 128 patients (5 with two), 52 at a first admission and 81 at a later one.
  Among the 960 patients with more than one admission, 79% have a different predictor profile or outcome at different
  admissions, so covariates are admission-level while genotype will be patient-level.
- **49 admission numbers appear twice, and the pairs look like an event row plus a non-event twin.** In every pair the
  patient and predictor profile are identical and exactly one row is an event. The 4,127 rows are therefore 4,078 unique
  admissions, and 133 of 4,078 is 3.26%, which is the paper's stated 3.3% (133 of 4,127 is 3.22%). The published
  equation reproduces from the rows with the twins in, so it appears to have been fitted on them. Refit on the 4,078
  unique admissions, the coefficients are larger (CMV 1.42 against 1.14, prior mould 1.99 against 1.64, intercept
  -5.59 against -5.45).
- **Consequence for the baseline cohort:** on the 4,078 unique admissions the published equation predicts 3.06% against
  3.26% observed (aROC 0.859), and the MCID is +0.00404 at 5% and +0.00113 at 20%, against +0.00445 and +0.00142 on the
  4,127 rows (D8). Resampling patients instead of admissions widens the bootstrap intervals by 17 to 28% (5% MCID:
  0.00405 to 0.00489 on the rows, 0.00371 to 0.00439 on the unique admissions).
- **Why `boscore_tidy.xlsx` and the `.dta` did not reconcile:** they carry the older `SCORE*` items, not the
  `BOSCORE_*` (2017) definitions, hence different predictors and 108 events.
- **Refit baseline (D10, addendum 6):** the baseline is now the refit, with aROC 0.858 apparent and 0.853 corrected, calibration slope 0.97, and a
  primary MCID of +0.00458 (bootstrap 90% 0.00387 to 0.00566) and a secondary of +0.00160. The published equation is a sensitivity analysis.
- **Outcome-risk sensitivity (addendum 4):** the adopted MCID, which treats the published equation as truth, is the lowest of four
  readings at every threshold. Using the refit as truth gives +0.00452 at 5% (+12%) and +0.00137 at 20% (+21%); the observed
  outcomes give +0.00441 and +0.00142. The intervals overlap, and the published coefficients are inside the refit's intervals. The
  4,127-row MCID (+0.00445, +0.00142) sits close to these, so reporting against it brackets them.
- **Decided (D9, addendum 3 of `phase3_decisions.md`):** the baseline cohort is the 4,078 unique admissions and the bootstrap resamples
  patients. The MCID is +0.00404 at 5% and +0.00113 at 20%. The 49 pairs are confirmed as a recording artefact (addendum 5), so the
  4,127-row cohort is a sensitivity analysis only and the revert condition is retired.

### Full BOSCORE dataset (`BOSCORE_recoded_2019.xlsx`, gitignored, contains direct identifiers)

8,660 rows by 94 columns: the wider BOSCORE database, with names, birth dates and admission, death and follow-up
dates. It has a patient ID (`IDIND`, 2,545 patients) and an admission ID (`IDRICO`, 5,784 admissions, 2006 to 2018),
and several rows per admission (fever episodes). Only counts and coded values have been looked at; no identifier
value has been printed.

- It is **not** the 4,127-row analysis file and cannot be linked to it: there is no shared key (`bs.xlsx` `id` runs 2
  to 4,694, `IDRICO` 7,496 to 14,455). `boscore.xlsx` is the file that carries both keys.
- I did not find an EORTC proven/probable IMD outcome column, and the predictor definitions differ (for example
  neutropenia over 10 days is 40% of rows against 29% for prolonged neutropenia in `bs.xlsx`).
- A date filter does not reproduce the paper's cohort: 2007 to 2016 gives 4,799 admissions in 2,151 patients, against
  4,127 in 1,944 (2007 to 2015 gives 4,293 and 1,952). So admissions per patient for the analysis set is still not
  exactly recoverable.
- As an approximation, 2.27 admissions per patient overall: 1,244 patients (49%) have one admission, 496 have two,
  312 have three, up to a maximum of 16. This is usable for a clustering sensitivity analysis, not as the exact answer.

**Pseudonymised extract.** `scripts/06_pseudonymise_boscore.R` writes `data/restricted/boscore_pseudonymised.csv` (8,660
rows, 88 columns; gitignored). It drops names and birth dates, replaces patient, admission, episode and ward IDs with
random codes, generalises the admission date to the year, converts every other date to days from admission, and
top-codes age at 90. Patient and admission structure is preserved one to one, and a separate check found no original
name or ID anywhere in it. `data/boscore_pseudonymisation_manifest.csv` (tracked) lists what was done to each column.
It is still personal data, because anyone holding the original file can re-link it: keep it out of git and email. The
codes come from a seeded permutation (seed 6283), so they regenerate from the original file and no key file exists.
### `boscore_tidy.xlsx` (gitignored, contains direct identifiers)

8,037 rows by 412 columns, one row per admission with the fever episodes spread across numbered column blocks (`_2` to
`_12`), plus Stata survival fields (`_st`, `_d`, `_t`, `_t0`) and an outcome column (`EORTC_PROV_PROB`). It spans 2006 to
2021: 8,031 admissions in 3,571 patients. Read it with `col_types = "text"`, because `readxl` guesses types from the
first 1,000 rows and silently turns later values in columns 135 onward into missing.

- **It is a consistent superset of the earlier export.** All 5,784 admissions of `BOSCORE_recoded_2019.xlsx` are in
  it, and on those the clinical variables agree to 99.9% or better. It adds 2,247 admissions (1,582 in 2019 to 2021 and
  560 in 2007 to 2016), so the earlier export was incomplete.
- **It does not reproduce the paper's cohort.** For 2007 to 2016 it has 5,359 admissions in 2,377 patients, against
  4,127 in 1,944. The outcome flag gives 83 cases in that window, and cases coded as proven or probable mould in any
  episode number 156 (111 within 60 days of admission), against 133. No year, age or ward filter matches, and the
  score items that resemble the model's predictors differ from Table 2 by up to 4 percentage points even in the best
  window. Only corticosteroids (11.4% against 11.3%) and aspergillus colonisation (2.6% against 2.6%) line up.
- **Data quality to fix before relying on it:** the outcome flag misses 91 admissions that have a proven or probable
  episode; 44 episode dates in 2007 to 2016 fall before the admission date (up to 304 days); 29 shared admissions have
  a different patient ID than in the earlier export, and 113 a different age; one age is 1,069; 646 admissions have no
  score items; 5 admissions appear on more than one row; the last column has no name and two columns share a name.
- **What would reconcile it:** the Stata `.do` file (or the `.dta`) that built the paper's analysis set. The
  `_st` and `_d` fields show the analysis was run in Stata, and the `id` in `bs.xlsx` (2 to 4,694) looks like a Stata
  row number. This is no longer needed: `boscore.xlsx` (next section) supplies the rule.

**`boscore_tidy.dta` (gitignored) is the Stata source of that spreadsheet, and it does not reconcile either.** 8,038 rows by
73 columns, with no patient ID or name (only the hospital admission number `idrico`, admission and event dates, age and
death). It matches the spreadsheet exactly on all 8,031 shared admissions (age, length of stay and the seven analysis
variables), has one extra admission, and leaves the outcome missing for 6,400 admissions with no fever episode, which
the spreadsheet silently filled with 0. The analysis is a competing-risks set-up (`_st`, `_d`, `comptrisk`): 108 IMD
events, 398 competing events and 7,532 censored, with all rows used. Its analysis variables carry the same names as
`bs.xlsx` (`prior_IA`, `hd_corticoster`, `uncont_malign`, `highrisk_chemo`, `neutropenia_10d`, `lymphopenia_50`,
`cmv`), but they are not the same population:

| | `bs.xlsx` (the paper's 4,127) | `boscore_tidy.dta` (7,391 with all 7 predictors) |
|---|---|---|
| IMD events | 133 (3.2%) | 108 (1.5%) |
| Mean predicted risk of the published equation | 3.26% | 3.61% |
| aROC of the published equation | 0.853 | 0.747 |
| Refit, high-risk chemotherapy / CMV | 0.79 / 1.14 (published) | 0.10 / -0.25 |

So the `.dta` is a different, later dataset (2006 to 2021, different outcome adjudication), not the derivation cohort. It
cannot give admissions per patient (no patient ID), and the spreadsheet's `IDIND` is the only linkage. It may be
useful as a rough external check on the baseline model, but only after the outcome definitions are reconciled: the
equation over-predicts here (3.6% against 1.5%) and discriminates less well. `haven` reads it; it is not in
`renv.lock` because no script uses it yet.

### Evidence mapping (addenda 7 to 9, 16, 19 and 21; provisional)

`h2` is the PTX3 chromosome class carrying A at rs2305619 (A-C-A, A-A-A, A-A-G: 45.6% in 1000G EUR, 48.2% gnomAD NFE); h1 is
G-A-G. rs1840680 A is the same class; rs3816527 C marks only A-C-A. Each study's contrast is converted to a per-copy effect by
solving a haplotype model, one effect is taken per study, and the effects are pooled with a Bayesian random-effects model. The
pipeline now takes any haplotype block: CLEC7A uses the Y238X haplotype (C-G, 7.7%), TLR4 the G-T haplotype (6.1%) through carriers
against non-carriers (D3c).

| term | primary analysis | OR per copy (90% interval) | tau |
|---|---|---|---|
| PTX3 h2 | 5 studies | 1.35 (0.99 to 1.81) | 0.17 |
| CLEC7A Y238X | 4 studies | 0.84 (0.48 to 1.48), after the Fischer 2016 correction (addendum 16) | 0.25 |
| CLEC7A A-G (second parameter, secondary) | 3 strand-resolved European studies | 1.86 (1.26 to 2.63); with the two East Asian studies 1.28 (0.79 to 1.95), a sensitivity analysis (D25) | 0.17 |
| TLR4 G-T | 4 studies | 1.44 (0.77 to 2.69) | 0.27 |
| CD209 G at rs4804800 | 3 studies | 1.80 (1.14 to 2.73) | 0.21 |
| IL10 rs1800896 C | 1 study, HR fallback | 2.59 (1.01 to 5.02), prior-dominated; **dropped from the panel (D14)** | 0.33 |
| IFNG rs2069705 G | 2 studies, HR fallback | 0.82 (0.49 to 1.69), prior-dominated; **dropped from the panel (D14)** | 0.33 |

**Source checks (addendum 8).** Herrero-Sanchez 2018 was read in full: the reference label is reversed in every table of the paper,
so the extraction looked as if the study pointed the other way. That was an artefact, 30 rows are corrected, and the study's effect is
positive (the statement in addendum 7 is retracted). Brunel 2018 was read from a summary of the page: the design is confirmed but the
per-SNP counts are only in a supplementary figure, and a Hardy-Weinberg check suggests its rs3816527 "homozygous" group may be the
major-allele homozygotes, so its sign is not secure. **Since then (addenda 16 to 21):** the direction check of every text-parsed row is
closed (`data/direction_check_log.csv`; Fischer 2016 and one Sainz 2012 row corrected, 43 direction corrections in all in
`data/extraction_corrections.csv`, plus one comparison-text fix for Tanpaibule 2020); the rs7309123 strand is resolved for nine of ten
studies (`data/rs7309123_strand_resolution.csv`, only Zhao 2022 undetermined), which gives the second CLEC7A parameter; and the East
Asian studies point the opposite way, which the rs3901533 linkage explains as an ancestry-specific tag (addendum 21). The user reports the
human round-2 checks confirmed the extraction, but the record is not in the repository. The priors and most mapping rules are proposals.

### Haplotypes

Derived from phased 1000G EUR (503 samples, 1,006 chromosomes) and rescaled to gnomAD NFE margins. The PTX3 trio
has **three** common haplotypes: G-A-G 54%, A-C-A 40%, A-A-A 6% (alleles ordered rs2305619, rs3816527,
rs1840680). **h2 is settled (D11)** as the rs2305619/rs1840680 A class (A-C-A plus A-A-A, 48.2% in gnomAD NFE); the
three-haplotype split is an unrun sensitivity analysis. East Asian CLEC7A haplotype frequencies (1000 Genomes EAS) are in
`data/haplotype_frequencies_eas.csv` (addendum 21).

## What to do next

1. ~~Email Bologna.~~ Nothing left to ask. The 49 twin pairs are confirmed as an artefact (addendum 5), and the two
   off-by-one Table 2 counts and the question of whether the published equation was fitted with the twins are set
   aside by decision.
2. ~~Decide on the empirical joint.~~ Done: adopted as D8 (addendum 2).
3. **Round-2 verification:** add the record of the human checks to `data/`. The **direction check of the text-parsed rows is closed** (addenda 18 to 20, `data/direction_check_log.csv`): Tang, Tanpaibule, Chen and Fisher are correct, one Sainz row (V0261) and all of Fischer 2016 and Herrero-Sanchez 2018's reversed rows are corrected;
   only Lupianez's two IFNG rows (a term no longer in the panel) were not checked.
4. **Evidence-to-parameter mapping (all terms done provisionally, addenda 7 to 9 and 16):** the user reports (19 October) that the human round-2
   checks confirmed the extraction against the PDFs and changed none of the flagged items; the verification record is not in the repo, so
   add it to `data/`. IL10 and IFNG were dropped on 19 October (addendum 10), so the panel is PTX3, CLEC7A, TLR4 and CD209. **An error was found on 19 October (addendum 16, section 0):**
   Fischer 2016 (S073) had its direction reversed (the table lists the reference group first, and the pipeline read the first group as the numerator); 12 rows are corrected and every script that uses CLEC7A was rerun
   (the Y238X pooled OR moved from 0.98 to 0.84; the layer results by at most 0.00006 in gain). **The direction check of the other text-parsed rows is closed (addenda 18 to 20).**
   The second CLEC7A parameter (OR 1.9 per A-G copy, 1.26 to 2.63, from three strand-resolved studies after Kalkanci 2020 was added, addenda 16 and 19) is estimated and roughly doubles the CLEC7A-containing layer's gain; only Zhao 2022's strand remains undetermined.
   Open questions from addendum 8 that the checks did not settle: Herrero-Sanchez Table 5, Brunel, hazard ratios for the terms with enough odds-ratio studies.
5. **Genetic layer (addenda 11 to 17):** `R/genetic_layer.R`, `R/genetic_helpers.R` and scripts 16 to 25 add the panel to the refit baseline on `load_base_cohort()` (genotype fixed per patient, exact over 128 profiles by the dosage combinations,
   recentred to hold prevalence; since addendum 24 the primary truth has an intercept per clinical profile instead) and compare the gain with the MCID. Four-term panel (profile-preserving truth): the median gain is about +0.0006 at 5% against an MCID of +0.00458 (reached in at most 0.6% of draws) and +0.0011 at 20% against +0.00160 (reached in 17%; +0.0010 and 16% under the global-intercept truth).
   The correlation ensemble and the sensitivity cohorts change this little; a study of this size has a sampling SD of the observed gain of 0.0006 to 0.0012 at 5%, as large as the gain. Shapley: CD209 about 51%, PTX3 23%, TLR4 19%, CLEC7A 7%.
   **With the second CLEC7A parameter (OR 1.9 per A-G copy, three strand-resolved studies) the gain roughly doubles and CLEC7A supplies about 51% of it** (CD209 27%, PTX3 12%, TLR4 10%); this is conditional on three studies.
   Donor genotype in the allo-HSCT stratum (6.7% of admissions, admission-level flag from addendum 14) changes the true cases by -0.1 to +1.3 per 1,000 allo admissions at 5% against 30 needed, because the baseline already flags 77% of that stratum (donor evidence: one to two studies per term, none for CD209).
   **CLEC7A's lead holds** on the sensitivity cohorts and the 150-member correlation ensemble (addendum 18): about 50% of the gain (44.8% to 53.3% over every setting, reading, score and threshold), CD209 25% to 32%; it is the largest single contributor in 54% to 70% of draws.
   To do: confirm the reading of "realised" (addendum 11).
   **Calibration check (addendum 24, scripts 31 to 34).** One global intercept left the clinical score slightly miscalibrated within profiles (slope 0.98), but the best clinical-only rule under that truth beat it by at most +0.00027 in net benefit in any draw, so the gains were not a correction of that miscalibration. The profile-preserving truth was adopted as primary after the global results were known (disclosed in the manuscript). It lowers the panels' gains at 5% by 5% to 6%, moves the 5% architecture frontier out by two to three loci, and changes the continuous NRI of single-term scores; conclusions are unchanged. **Addendum 25 reran the sampling, ensemble, attribution, sensitivity-cohort, donor, Aim 3b and robustness analyses under the profile-preserving truth** (same seeds and draws; scripts 17 to 21, 24, 26 and 27 with `PROSE_SCHEME=profile`). Every conclusion is unchanged. The one notable difference is Aim 3b with all of the effect in non-neutropenic admissions: the oracle score reaches the 5% MCID in 2% of draws (four terms) and 6% (five), against 18% and 42% under the global-intercept truth. The reason is the construction, not a correction of the comparator: that scenario puts an effect about five times the pooled one into a low-risk stratum, and one global intercept flattens the clinical gradient there (addendum 26, which also reran the all-ancestry CLEC7A variant, so every analysis is now on the profile-preserving truth). Addendum 27 (script 36) checked the other scenarios at effects as large: the profile truth lowers their median gain by 11% to 28% (more where the effect is in the neutropenic stratum), but none has such an effect, so they differ by only 3% to 8%. Addendum 28 (script 37) repeated it for the five-term panel and the realised score: the same pattern, slightly stronger for five terms (31% against 26% lower when the large effect is in the neutropenic stratum), and one sign reversal where the realised score's gain is negligible under both truths; nothing in the paper changes.
6. **Aim 3b (addenda 15 and 17):** the expected direction was committed before the run (`dc60a4e`). The Feys hypothesis (effects larger in neutropenic admissions) and Brunel's stratified data (larger in non-neutropenic, P = 0.07 and 0.049) point opposite ways.
   Neutropenic admissions supply 80% to 95% of the whole-cohort gain under homogeneous effects. A common score with pooled estimates is the safer choice (a score built on the Feys reading is the most harmful when wrong). A cohort-sized study would detect a Feys-sized modification
   (58% power with four terms, 89% with five) and a Brunel-sized one (94% to 100%), better than predicted. The five-term gain roughly doubles in every scenario, and the five-term layer repeated on the sensitivity cohorts and the correlation ensemble (addendum 18) gives the same picture: the neutropenic stratum supplies 77% to 83% of the gain at 5% and 96% to 97% at 20% (all under the profile-preserving truth). To do: a genotyped neutropenia-stratified analysis would be the real test; the interaction test was not repeated on those settings.
7. **Aim 4 sweep (addenda 16, 17 and 24):** under the profile-preserving truth 40% of the 4,000 architectures reach the 5% MCID and 77% the 20% MCID (46% and 77% under the global-intercept truth). At OR 1.5 or below no number of loci up to 20 reaches the 5% MCID (under the global-intercept truth one architecture did, by +0.000007); at OR 2.0 it takes 9 to 12 loci at common frequencies (7 to 9 under the global-intercept truth);
   the four-term panel is equivalent to about 7 loci of OR 1.3. With estimated weights, per-locus precision of 0.1 or better costs almost nothing; at 0.2 to 0.4, the precision of the present evidence, 12% to 36% of the architectures that reach the MCID with exact weights no longer do. The frontier is slightly non-monotone in frequency near 0.5 (the gain peaks at about 0.45). Not done: the pathogen-specificity scenario.
   **Polygenic score benchmark (addendum 23, script 30):** with true weights a polygenic score needs an odds ratio per SD of 3.6 at 5% (4.2 under the profile-preserving truth) and 1.8 at 20%; with estimated weights at 5% and an odds ratio of 6, an accuracy of 0.74 and about 7,400 balanced discovery cases for 10,000 effective loci (16,000 under the profile-preserving truth; Fisher 2017 had 483 cases). Optimistic throughout; script 30's truth holds prevalence only to 0.2% to 1.0%.
8. **Strand and overlap:** only Zhao 2022's rs7309123 strand remains undetermined; the Spanish registry and Jena overlaps are only partly adjudicated.
9. **Tertiary NRI comparison (addenda 22 and 24): done.** Under the profile-preserving truth the four-term haplotype panel has a continuous NRI of about 34% (31% under the global-intercept truth) with a net-benefit gain an eighth of the MCID; the predecessor's
   three-SNP independent structure gives 25% at the pooled PTX3 effect (close to its reported 22.5%) against 34% for the h2 haplotype (13% under the global-intercept truth), so the continuous NRI of a single-term score depends on how the truth is built while its net-benefit gain does not; the net-benefit gain is four times larger for the independent structure (+0.00033 against +0.00008). The three-SNP score on real chromosomes is
   harmful in net benefit (-0.0004) but has a positive continuous NRI (29%). The predecessor's own effects and NRI form are not on record.
10. **Write-up:** full first draft on the Draft manuscript page (`manuscript.qmd`), revised on 20 October to make the profile-preserving truth primary. To do: authors, ethics and data-access statements, funding, the two missing references, and a
    review of the discussion.

Loose ends for write-up: the PubMed search date is unrecorded (Embase was 8 October 2026); the automated
screening decision file for all 831 records is not assembled; Skonieczna 2017's supplement is missing; whether an
Italian thesis title is required (F4).

## Repository map

| Path | Contents |
|---|---|
| `R/` | `baseline.R` (equation, profile metrics, net benefit), `copula.R` (correlation structure, fitting, exact profile probabilities), `decision_curve.R` (exact curves, comparison, plots), `haplotypes.R` (IPF rescaling, LD), `bs_data.R` (checks and empirical cohort from the admission-level data), `pseudonymise.R` (BOSCORE extract), `boscore_data.R` (base cohort, twin correction, patient bootstrap, reconciliation) |
| `scripts/` | the scripts above |
| `bs.xlsx` | admission-level derivation data. **Gitignored and not in the repository**; copy it over separately |
| `tests/testthat/` | the test suite, run with `tests/run_tests.R` |
| `data/` | baseline inputs, panel tables, allele frequencies, phased haplotype outputs; `data/restricted/` (gitignored) holds the pseudonymised extract |
| `output/` | fitted results; `output/v1/` is the archived first fit; `output/figures/` has PNGs |
| `renv.lock`, `renv/`, `.Rprofile` | pinned package versions; the private library in `renv/library/` is not tracked |
| `*.qmd`, `_quarto.yml`, `_setup.R` | the report site; the render list in `_quarto.yml` is explicit, so add new pages there. `manuscript.qmd` is the draft paper, with `references.bib` and `vancouver.csl` |
| `docs/` | rendered site (committed). `_freeze/` caches executed pages |
| `phase3_decisions.md` | the locked pre-specification record |
| `02_current_status.md` | the long working status document |
| `IMD_genetics_analysis_bundle (1).zip` | the original bundle: full extraction tables, round-2 lists, project description, README |

## Gotchas

- Positron does not allow forked R workers. The parallel code uses fresh R processes (`run_parallel()` in
  `R/copula.R`), which source `R/*.R` from the working directory, so run from the project root.
- Quarto pages are frozen (`freeze: auto`). The test table on the baseline page only refreshes after
  `rm -rf _freeze/baseline`.
- The prevalence target is 3.3% as pre-specified, but 133/4127 is 3.22%. The difference is about 0.3 standard
  errors and does not affect the fit.
- The README, project description and pre-analysis checklist inside the zip pre-date the targeted-recovery pass and
  supplementary files, so their counts are out of date. `02_current_status.md` outranks them. The round-2
  verification list is 907 rows in the bundle file but "1,380" in the text, and this is unresolved.
- The repository was private from 18 to 19 October 2026 and is **public again since 19 October** (lewis-id-lab/PROSE), with
  GitHub Pages publishing `docs/`. The history contains no spreadsheets or patient-level data, but the bundle zip, the decision
  record and the data summaries are visible. Switching back to private would take the site down on the free plan.
- `.posit/` (assistant settings) is ignored and untracked.
- `bs.xlsx` and `BOSCORE_recoded_2019.xlsx` are patient-level data, and the second has direct identifiers. All
  spreadsheets and `.sav`/`.dta` files are gitignored, so they will not reach the other computer through git. Move
  them by a secure route, and note that the two `bs.xlsx` tests skip without the file. Further copies sit in
  `~/Downloads` and `~/Projects`, and `~/riskscore/` (a Shiny app with an `rsconnect` folder) holds a copy of
  `bs.xlsx`, so check whether that app was deployed with the file.
- `boscore.xlsx` is the base dataset. `read_boscore()` skips names, birth dates and every column it does not need when it reads the
  file, and `scripts/08_reconcile_boscore.R` writes `data/restricted/base_cohort_pseudonymised.csv` (gitignored), which
  `load_base_cohort()` prefers, so analysis does not need the identified file. Without either file the `boscore` tests skip.
- Before using the identified file with any AI tool, check that the data-sharing agreement and ethics approval allow
  it. Analyse the pseudonymised extract (`data/restricted/`) and not the identified file.
