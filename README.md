# PROSE: host genetics and decisions about invasive mould disease

A simulation study of whether adding host genetic variants (PTX3, CLEC7A, TLR4, CD209) to a validated clinical score for invasive mould disease
(Stanzani 2019) would change clinical decisions. The primary endpoint is decision-curve net benefit at the 5% threshold, judged against a
pre-specified minimum clinically important difference.

**Status: provisional.** The evidence mapping and the genetic layer are being built, and no result here is final. `status.md` is the short hand-off
(where the project stands, how to continue, known problems). `phase3_decisions.md` is the dated record of the design decisions. The pre-specification is
locked at the git tag `prespec-lock`.

## Layout

| path | contents |
|---|---|
| `R/`, `scripts/`, `tests/` | the model code, the numbered analysis scripts (run from the project root, in order), and the test suite (`source("tests/run_tests.R")`) |
| `output/` | fitted results and summaries; per-draw simulation files are gitignored and regenerate from scripts 16 to 18 |
| `data/` | small baseline inputs, panel tables, allele and haplotype frequencies |
| `*.qmd`, `_quarto.yml` | the report site; `docs/` is its rendered output |
| `renv.lock` | pinned package versions (`renv::restore()`) |

## Data

No patient-level data are in this repository. The admission-level files (`*.xlsx`, `*.dta`, `*.sav`) and the pseudonymised extracts in `data/restricted/`
are gitignored; they are personal data and must be copied between computers by a secure route. Tests that need them skip when they are absent.

## Updating the website

Edit the `.qmd` pages (the decision record is included in the Pre-specified decisions page), run `quarto render`, and commit `docs/` and `_freeze/`
with the source. `.github/workflows/pages.yml` publishes `docs/` to GitHub Pages on each push to `main` once the repository is public; see `status.md`.
