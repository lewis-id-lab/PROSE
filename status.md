# Where we are, and how to pick this up

Written 17 October 2026. This is the short hand-off. The long working record is `02_current_status.md`
(evidence base, extraction detail), and the design decisions are in `phase3_decisions.md`.

The project asks whether adding host genetic variants to a validated clinical score (Stanzani 2019) for invasive
mould disease (IMD) would change clinical decisions. It is a simulation study, with decision-curve net benefit at
the 5% threshold as the primary endpoint.

## Before you leave this computer

Local `main` is one commit ahead of GitHub, and the pre-specification tag exists only locally.

```bash
git add status.md .Rprofile renv.lock renv && git commit -m "Add hand-off status file and renv lockfile"
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
source("tests/run_tests.R")     # 45 tests, 209 expectations, all passing at the last run
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

Script 01 scales with core count and would take hours on 4 cores. Cut `n_starts` if you need to.

## Where the project stands

| Phase | State |
|---|---|
| 0 Set-up, 1 Literature search | complete (831 unique records screened) |
| 2 Screening and extraction | round 1 complete; **round 2 (human verification) not started** |
| 3 Parameter derivation | design decisions settled and locked; parameter derivation not started, gated on round 2 |
| 4 Baseline reproduction | **done** except for data only Bologna can supply |
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
| Predictor correlation | version 2, nine parameters |
| MCID (20% more true cases, no extra false positives) | **+0.00428 net benefit at 5% (primary)**, +0.00140 at 20% (secondary) |

`phase3_decisions.md` has a dated addendum written after the lock. Further changes there should be dated addenda,
not edits.

### Baseline model (Aim 2) and decision curves

- The Stanzani equation is fixed and never refitted. The joint distribution of its seven predictors is not
  published, so a Gaussian-copula correlation matrix is a fitted assumption.
- Version 2 fits all 16 published targets (chi-square 0.36, largest miss 0.4 standard errors). The fit is
  saturated: about nine independent constraints for nine parameters, so the good fit does not validate the
  structure. Malignancy comes out nearly uncorrelated with the rest of its block (0.04).
- Baseline net benefit at 5% is 0.01478 (14.7% flagged). The ensemble 90% range for the 5% MCID is 0.00397 to
  0.00460.
- Decision curves are exact functions of the 128 clinical profiles, so they are step functions. The 5% MCID sits on
  a plateau from about 4.2% to 7.8%.

### Haplotypes

Derived from phased 1000G EUR (503 samples, 1,006 chromosomes) and rescaled to gnomAD NFE margins. The PTX3 trio
has **three** common haplotypes: G-A-G 54%, A-C-A 40%, A-A-A 6% (alleles ordered rs2305619, rs3816527,
rs1840680). **"h2" in the settled h1/h2 term needs an explicit definition.** I recommend the rs2305619/rs1840680
A class (A-C-A plus A-A-A), with the three-haplotype split as a sensitivity analysis.

## What to do next

1. **Email Bologna** (you): the three 2x2 cross-tabulations (chemotherapy by neutropenia, neutropenia by
   lymphopenia, malignancy by chemotherapy) and the admissions-per-patient distribution. Malignancy by
   chemotherapy is the most informative, because the fit's near-zero malignancy correlation is the weakest part
   of the structure.
2. **Round-2 verification of the 299 panel rows** (`panel_verification_worklist.csv`, inside the zip): the 10
   strand-flip cells and 7 rs7309123 cells first, then the 33 agent-transcribed rows. For rs7309123, compare the
   reported allele's control frequency with gnomAD (G 0.43, C 0.57) in European cohorts.
3. **Evidence-to-parameter mapping:** one effect per study per term, as a log odds ratio per risk allele
   with nucleotide-pinned alleles; collapse overlapping cohorts (Spanish registry, Jena AML); map SNP effects to
   haplotype effects; Bayesian random-effects pooling with prediction intervals.
4. **Genetic layer:** draw effects from the pooled posterior, recentre, and report the ceiling (true effects in the
   score) and the realised value (published estimates in the score, a shrunk truth). Carry the correlation
   ensemble through. Decide how to handle admissions per patient (about 2.1 per patient, genotype fixed per
   patient).
5. **Aim 3b:** note the design tension. The Feys logic predicts a genotype effect in neutropenic patients, but the
   only direct stratified evidence (Brunel 2018, derived by subtraction and flagged ambiguous) shows it in
   non-neutropenic patients. State the expected direction before running it.
6. **Aim 4 sweep:** 1 to 20 loci by OR 1.1 to 3.0 by risk-allele frequency 0.05 to 0.50, using the exact
   convolution approach.

Loose ends for write-up: the PubMed search date is unrecorded (Embase was 8 October 2026); the automated
screening decision file for all 831 records is not assembled; Skonieczna 2017's supplement is missing; whether an
Italian thesis title is required (F4).

## Repository map

| Path | Contents |
|---|---|
| `R/` | `baseline.R` (equation, profile metrics, net benefit), `copula.R` (correlation structure, fitting, exact profile probabilities), `decision_curve.R` (exact curves, comparison, plots), `haplotypes.R` (IPF rescaling, LD) |
| `scripts/` | the four scripts above |
| `tests/testthat/` | the test suite, run with `tests/run_tests.R` |
| `data/` | baseline inputs, panel tables, allele frequencies, phased haplotype outputs |
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
- The repository is **public** and tracks the bundle zip. GitHub Pages is not enabled. The site would be public if
  you turn it on (Settings, Pages, branch `main`, folder `/docs`). Decide whether the zip and the unpublished
  extraction data should be in a public repository.
- `.posit/` (assistant settings) is ignored and untracked.
