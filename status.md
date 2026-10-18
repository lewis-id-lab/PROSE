# Where we are, and how to pick this up

Written 17 October 2026. This is the short hand-off. The long working record is `02_current_status.md`
(evidence base, extraction detail), and the design decisions are in `phase3_decisions.md`.

The project asks whether adding host genetic variants to a validated clinical score (Stanzani 2019) for invasive
mould disease (IMD) would change clinical decisions. It is a simulation study, with decision-curve net benefit at
the 5% threshold as the primary endpoint.

## Before you leave this computer

Local `main` is one commit ahead of GitHub, and the pre-specification tag exists only locally.

```bash
git add -A        # bs.xlsx is gitignored; check `git status` shows no spreadsheet
git commit -m "Add hand-off status file, renv lockfile and bs.xlsx checks"
git push origin main
git push origin prespec-lock
```

## Starting on another computer

```bash
git clone https://github.com/Russlewisbo/PROSE.git
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
source("tests/run_tests.R")     # 62 tests, 264 expectations, all passing at the last run
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
| `scripts/07_empirical_baseline.R` | baseline decision curve and MCID on the empirical cohort, with bootstrap intervals and the comparison with the copula | under 1 min |

Script 01 scales with core count and would take hours on 4 cores. Cut `n_starts` if you need to.

## Where the project stands

| Phase | State |
|---|---|
| 0 Set-up, 1 Literature search | complete (831 unique records screened) |
| 2 Screening and extraction | round 1 complete; **round 2 (human verification) not started** |
| 3 Parameter derivation | design decisions settled and locked; parameter derivation not started, gated on round 2 |
| 4 Baseline reproduction | **done**; the real joint distribution is now available (`bs.xlsx`), and whether to switch to it is pending; admissions per patient still needs Bologna |
| 5 Genetic layer simulation | not started |
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
| Panel | four haplotype terms (PTX3, CLEC7A, TLR4, CD209); IL10 rs1800896 and IFNG rs2069705 optional |
| Dropped / excluded | TLR9 rs352140 dropped; Herrero-Sanchez 2018 excluded from the primary analysis |
| TLR4 pooling | carrier model with a continuity correction |
| Genetic term | recentred so prevalence stays at 3.3% |
| Baseline aROC target | 0.85 |
| Haplotype phase | coupling in every pair |
| Baseline cohort (D8, addendum 2) | **the empirical 128-profile distribution from `bs.xlsx`**; the version 2 copula (D7) is a sensitivity analysis |
| MCID (20% more true cases, no extra false positives) | **+0.00445 net benefit at 5% (primary)**, +0.00142 at 20% (secondary), on the empirical cohort. The values at the tag were +0.00428 and +0.00140 (copula) |

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

- `id` is unique per row, so it is an admission number. **Admissions per patient cannot be recovered**, and that
  still needs Bologna.
- Only 75 of the 128 predictor profiles occur (35 with fewer than 5 admissions). The fitted copula is still about 10%
  of probability mass away from the real joint, and some of its correlations are wrong in sign (for example
  malignancy by CMV is -0.15 in the data against +0.23 fitted), so the saturated fit did not validate the structure.
- The MCID on the real joint is +0.00445 at 5% (bootstrap 90% interval 0.00412 to 0.00478) and +0.00142 at 20%
  (0.00116 to 0.00170). The values locked at the tag (+0.00428, +0.00140) are inside those intervals.
- **Decision made (D8, addendum 2 of `phase3_decisions.md`):** the baseline cohort for the genetic layer is the empirical
  128-profile distribution, with the copula kept as a sensitivity analysis. Read it with `read_bs_data()` and
  `empirical_profile_counts()`; the profile counts themselves are not saved, because several profiles have fewer than 5 admissions.

### Full BOSCORE dataset (`BOSCORE_recoded_2019.xlsx`, gitignored, contains direct identifiers)

8,660 rows by 94 columns: the wider BOSCORE database, with names, birth dates and admission, death and follow-up
dates. It has a patient ID (`IDIND`, 2,545 patients) and an admission ID (`IDRICO`, 5,784 admissions, 2006 to 2018),
and several rows per admission (fever episodes). Only counts and coded values have been looked at; no identifier
value has been printed.

- It is **not** the 4,127-row analysis file and cannot be linked to it: there is no shared key (`bs.xlsx` `id` runs 2
  to 4,694, `IDRICO` 7,496 to 14,455).
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
### Haplotypes

Derived from phased 1000G EUR (503 samples, 1,006 chromosomes) and rescaled to gnomAD NFE margins. The PTX3 trio
has **three** common haplotypes: G-A-G 54%, A-C-A 40%, A-A-A 6% (alleles ordered rs2305619, rs3816527,
rs1840680). **"h2" in the settled h1/h2 term needs an explicit definition.** I recommend the rs2305619/rs1840680
A class (A-C-A plus A-A-A), with the three-haplotype split as a sensitivity analysis.

## What to do next

1. **Email Bologna** (you), now a short request: an anonymised patient code (or, at minimum, the number of patients
   with 1, 2, 3, ... admissions and how many of the 133 events were at a first versus a later admission), and
   whether the paper or the file is right on the two counts that differ by one. The three cross-tabulations are no
   longer needed, since `bs.xlsx` gives them directly.
2. ~~Decide on the empirical joint.~~ Done: adopted as D8 (addendum 2).
3. **Round-2 verification of the 299 panel rows** (`panel_verification_worklist.csv`, inside the zip): the 10
   strand-flip cells and 7 rs7309123 cells first, then the 33 agent-transcribed rows. For rs7309123, compare the
   reported allele's control frequency with gnomAD (G 0.43, C 0.57) in European cohorts.
4. **Evidence-to-parameter mapping:** one effect per study per term, as a log odds ratio per risk allele
   with nucleotide-pinned alleles; collapse overlapping cohorts (Spanish registry, Jena AML); map SNP effects to
   haplotype effects; Bayesian random-effects pooling with prediction intervals.
5. **Genetic layer:** draw effects from the pooled posterior, recentre, and report the ceiling (true effects in the
   score) and the realised value (published estimates in the score, a shrunk truth). Carry the correlation
   ensemble through. Decide how to handle admissions per patient (about 2.1 per patient, genotype fixed per
   patient).
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
| `R/` | `baseline.R` (equation, profile metrics, net benefit), `copula.R` (correlation structure, fitting, exact profile probabilities), `decision_curve.R` (exact curves, comparison, plots), `haplotypes.R` (IPF rescaling, LD), `bs_data.R` (checks and empirical cohort from the admission-level data), `pseudonymise.R` (BOSCORE extract) |
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
- Before using the identified file with any AI tool, check that the data-sharing agreement and ethics approval allow
  it. Analyse the pseudonymised extract (`data/restricted/`) and not the identified file.
