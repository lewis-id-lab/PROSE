# Where we are, and how to pick this up

Written 17 October 2026. This is the short hand-off. The long working record is `02_current_status.md`
(evidence base, extraction detail), and the design decisions are in `phase3_decisions.md`.

The project asks whether adding host genetic variants to a validated clinical score (Stanzani 2019) for invasive
mould disease (IMD) would change clinical decisions. It is a simulation study, with decision-curve net benefit at
the 5% threshold as the primary endpoint.

## Where the code lives, and the website

The repository is **lewis-id-lab/PROSE** (private), set as `origin`; the earlier personal copy, Russlewisbo/PROSE, is the remote `personal` and is no
longer updated. Push both `main` and the pre-specification tag (`git push origin main prespec-lock`) from this computer if either has changed.
`git status` should show no spreadsheet before you commit. The per-draw simulation files of the genetic layer
(`output/genetic_*_draws.csv`, 2 to 9 MB each) are gitignored and regenerate exactly from scripts 16 to 18.

**The website** is the Quarto site in `docs/`. To update it: edit the `.qmd` pages (or the decision record, which is included in the Pre-specified decisions page),
run `quarto render` from the project root, and commit `docs/` and `_freeze/` with the source. The workflow `.github/workflows/pages.yml` publishes `docs/` to GitHub Pages
on every push to `main`, but only when the repository is public, because the organisation is on the free GitHub plan and Pages is not available there for private
repositories. While the repository is private the job is skipped. To go live: make the repository public (this exposes the unpublished analysis, the decision record,
the results and the extraction bundle zip; the repository was private from 18 October for that reason), then enable Pages once with
`gh api -X POST repos/lewis-id-lab/PROSE/pages -f build_type=workflow` (or Settings, Pages, Source: GitHub Actions). The site would be at https://lewis-id-lab.github.io/PROSE/.

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
| `scripts/20_genetic_donor_stratum.R` | donor genotype in the allo-HSCT stratum (D1 secondary analysis): stratum built from Table 1 margins, donor effects, gain against the stratum and whole-cohort MCIDs, sampling variability; needs the donor draws from scripts 12 and 14 | about 1 min |
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
| 2 Screening and extraction | round 1 complete; **round 2 (human verification) not started** |
| 3 Parameter derivation | design decisions settled and locked; parameter derivation not started, gated on round 2 |
| 4 Baseline | **done**; the base dataset is `boscore.xlsx` with the 49 twin admissions removed (D9), and the baseline model is a refit on it (D10, addendum 6): 4,078 admissions in 1,945 patients |
| 5 Genetic layer simulation | **in progress** (addenda 11 to 13): exact layer on the base cohort with the four-term panel, correlation ensemble, per-term attribution, sampling variability, sensitivity cohorts and donor genotype in the allo-HSCT stratum; Aim 3b and the CLEC7A second parameter remain |
| 6 Sensitivity and architecture sweep | not started |
| 7 Write-up | not started |

### Evidence base (provisional: no row has been checked against a source PDF by a human)

- 81 records included, 63 extracted, 2,008 variant rows. 970 are usable (823 estimates, 115 complete 2x3 genotype
  tables, 32 complete 2x2 tables).
- 13 variants reach three or more meta-analysable studies. Six reach four or more independent cohorts:
  rs16910526, rs2305619, rs7309123, rs3816527, rs4986790, rs4986791.
- Allele orientation is unresolved for 10 study-by-variant cells (strand flips) and 7 rs7309123 cells (C/G
  palindrome).

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

### Evidence mapping (addenda 7 and 8; provisional)

`h2` is the PTX3 chromosome class carrying A at rs2305619 (A-C-A, A-A-A, A-A-G: 45.6% in 1000G EUR, 48.2% gnomAD NFE); h1 is
G-A-G. rs1840680 A is the same class; rs3816527 C marks only A-C-A. Each study's contrast is converted to a per-copy effect by
solving a haplotype model, one effect is taken per study, and the effects are pooled with a Bayesian random-effects model. The
pipeline now takes any haplotype block: CLEC7A uses the Y238X haplotype (C-G, 7.7%), TLR4 the G-T haplotype (6.1%) through carriers
against non-carriers (D3c).

| term | primary analysis | OR per copy (90% interval) | tau |
|---|---|---|---|
| PTX3 h2 | 5 studies | 1.35 (0.99 to 1.81) | 0.17 |
| CLEC7A Y238X | 4 studies | 0.98 (0.56 to 1.72) | 0.26 |
| TLR4 G-T | 4 studies | 1.44 (0.77 to 2.69) | 0.27 |
| CD209 G at rs4804800 | 3 studies | 1.80 (1.14 to 2.73) | 0.21 |
| IL10 rs1800896 C | 1 study, HR fallback | 2.59 (1.01 to 5.02), prior-dominated; **dropped from the panel (D14)** | 0.33 |
| IFNG rs2069705 G | 2 studies, HR fallback | 0.82 (0.49 to 1.69), prior-dominated; **dropped from the panel (D14)** | 0.33 |

**Source checks (addendum 8).** Herrero-Sanchez 2018 was read in full: the reference label is reversed in every table of the paper,
so the extraction looked as if the study pointed the other way. That was an artefact, 30 rows are corrected, and the study's effect is
positive (the statement in addendum 7 is retracted). Brunel 2018 was read from a summary of the page: the design is confirmed but the
per-SNP counts are only in a supplementary figure, and a Hardy-Weinberg check suggests its rs3816527 "homozygous" group may be the
major-allele homozygotes, so its sign is not secure. The other PTX3 studies are unchecked. The CLEC7A rs7309123 rows (26, from 8
studies) are not used because the SNP is a C/G palindrome with no control frequency recorded. The rows are otherwise unverified
(round 2 pending), and the priors and most mapping rules are proposals. Four questions are open for you at the end of addendum 8.

### Haplotypes

Derived from phased 1000G EUR (503 samples, 1,006 chromosomes) and rescaled to gnomAD NFE margins. The PTX3 trio
has **three** common haplotypes: G-A-G 54%, A-C-A 40%, A-A-A 6% (alleles ordered rs2305619, rs3816527,
rs1840680). **"h2" in the settled h1/h2 term needs an explicit definition.** I recommend the rs2305619/rs1840680
A class (A-C-A plus A-A-A), with the three-haplotype split as a sensitivity analysis.

## What to do next

1. ~~Email Bologna.~~ Nothing left to ask. The 49 twin pairs are confirmed as an artefact (addendum 5), and the two
   off-by-one Table 2 counts and the question of whether the published equation was fitted with the twins are set
   aside by decision.
2. ~~Decide on the empirical joint.~~ Done: adopted as D8 (addendum 2).
3. **Round-2 verification of the 299 panel rows** (`panel_verification_worklist.csv`, inside the zip): the 10
   strand-flip cells and 7 rs7309123 cells first, then the 33 agent-transcribed rows. For rs7309123, compare the
   reported allele's control frequency with gnomAD (G 0.43, C 0.57) in European cohorts.
4. **Evidence-to-parameter mapping (all terms done provisionally, addenda 7 to 9):** the user reports (19 October) that the human round-2
   checks confirmed the extraction against the PDFs and changed none of the flagged items; the verification record is not in the repo, so
   add it to `data/`. IL10 and IFNG were dropped on 19 October (addendum 10), so the panel is PTX3, CLEC7A, TLR4 and CD209. Still to do: the second CLEC7A parameter
   (needs the strand of rs7309123 per study; held at zero), and the open questions at the end of addendum 8 that the checks did not settle
   (Herrero-Sanchez Table 5, Brunel, hazard ratios for the terms with enough odds-ratio studies).
5. **Genetic layer (addenda 11 to 13):** `R/genetic_layer.R`, `R/genetic_helpers.R` and scripts 16 to 20 add the four terms to the refit baseline
   on `load_base_cohort()` (genotype fixed per patient, exact over 128 profiles by 81 dosage combinations, recentred to hold prevalence)
   and compare the gain with the MCID. Results: at 5% the median gain is +0.0004 to +0.0008 against an MCID of +0.00458, reached in
   at most about 1% of draws; at 20% (secondary) +0.0008 to +0.0013 against +0.00160, reached in 3.5% to 36%. Neither the 150-member correlation
   ensemble nor the sensitivity cohorts (published equation on 4,078 or 4,127 rows, refit on 4,127) change this materially. A study of this size has a
   sampling SD of the observed gain of 0.0006 to 0.0013 at 5%, the same size as the expected gain or larger, and it reaches the 5% MCID in at most
   5% of repetitions in any run. CD209 supplies half or more of the gain (Shapley), PTX3 and TLR4 most of the rest, CLEC7A almost none. Donor genotype
   in the allo-HSCT stratum (7.1% of admissions) adds 3 to 4 cases per 1,000 allo admissions at 5% against 21 needed (secondary analysis; the stratum
   is built from Table 1 margins because the base cohort has no allo flag, and donor evidence is one to two studies per term, none for CD209).
   To do: confirm the reading of "realised" (addendum 11); add an admission-level allo flag to the pseudonymised base cohort if the analysis file has one;
   the CLEC7A second parameter; Aim 3b.
6. **Aim 3b:** note the design tension. The Feys logic predicts a genotype effect in neutropenic patients, but the
   only direct stratified evidence (Brunel 2018, derived by subtraction and flagged ambiguous) shows it in
   non-neutropenic patients. State the expected direction before running it.
7. **Aim 4 sweep:** 1 to 20 loci by OR 1.1 to 3.0 by risk-allele frequency 0.05 to 0.50, using the exact
   convolution approach.

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
| `*.qmd`, `_quarto.yml`, `_setup.R` | the report site; the render list in `_quarto.yml` is explicit, so add new pages there |
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
- The repository was set to **private** on 18 October 2026. It was public until then. The history contains no
  spreadsheets or patient-level data, but the bundle zip, the decision record and the data summaries were visible,
  and any clones or forks made in that time cannot be recalled. GitHub Pages is not enabled, and a private
  repository needs a plan that supports private Pages before the site can be published.
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
