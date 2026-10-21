# PROSE: notes for AI assistants and new contributors

A simulation study of whether host genotype adds decision value to a clinical risk score for invasive mould disease (IMD): decision-curve net
benefit at a 5% threshold against a pre-specified minimum clinically important difference (MCID). R code in `R/`, runnable steps in `scripts/`,
tests in `tests/testthat/`, results in `output/`, a Quarto site in `docs/` (published to GitHub Pages). This file is a map and a list of pitfalls;
the decisions themselves are in `phase3_decisions.md` and the position is in `status.md` (short) and `02_current_status.md` (long).

## Where things are

- `phase3_decisions.md`: the decision record. Decisions D1 to D27 and addenda 1 to 28. Append a dated addendum for every change after the
  pre-specification lock (tag `prespec-lock`) and say what had been seen when the change was made. Do not rewrite earlier addenda; add a pointer.
- `status.md` (hand-off, script run-order table), `02_current_status.md` (working record), `index.qmd` (overview table): keep their phase tables and
  numbers in step with the addenda. `manuscript.qmd` carries its numbers as prose, not inline code, so a rerun means editing the prose.
- `_quarto.yml` lists the pages that are rendered; `AGENTS.md` and `README.md` are not on the site.

## The one modelling choice to know

The genetic truth is built with an intercept **per clinical profile** (D27, profile-preserving): the clinical score stays calibrated within every
profile, so a gain cannot be the correction of induced miscalibration. The first version used one global intercept (D5); it is kept as the
sensitivity analysis. `PROSE_SCHEME=profile` selects the primary construction in the genetic-layer code (`genetic_scheme()`); scripts 17 to 22, 24, 26
and 27 then write files ending in `_profile`, and the global files are left alone. Scripts 31 to 37 compare the two. Every analysis reported in the
manuscript is on the profile-preserving truth.

## Running things

- Run from the project root. Packages are pinned with renv. Use the base pipe `|>`.
- Tests: `testthat::test_dir("tests/testthat")`. All pass, and none skip when the data files below are present.
- Parallel work goes through `start_workers()` (`R/genetic_helpers.R`), which sets one BLAS thread per worker. Without that, workers oversubscribe the
  machine. After cancelling a run, check for orphaned workers (`pgrep -f workRSOCK`) and stop them.
- Do not leave helper functions in the global environment with the names of base functions: a `sub()` defined in a session made seven tests fail.

## Data (never commit)

Patient-level data are gitignored: `boscore.xlsx`, `bs.xlsx` and `BOSCORE_recoded_2019.xlsx` in the project root, and `data/restricted/`. Eight tests
need them. `data/restricted/` is rebuilt by script 08 (needs `boscore.xlsx` and `bs.xlsx`, checks the reconciliation) and script 06 (needs
`BOSCORE_recoded_2019.xlsx`). The per-draw simulation files are gitignored too (regenerated from fixed seeds). The repository is public.

## Git

- Commit in two parts: source, scripts, tests, outputs and figures first; then the regenerated `docs/` and `_freeze/` (from `quarto render`, no hand edits).
- Stage by explicit path, not `git add -A` at the root. Check what is staged before committing.
- Push only when asked. Before a push, scan the commits for spreadsheets, `data/restricted/`, per-draw files and identifier-like text.
- Run the full test suite and `quarto render` (check for unresolved cross-references) before committing.

## Standing caveats

The evidence is provisional (three to five studies per term; the human round-2 verification record is not in the repository). Continuous NRI is
unstable for single-term scores and depends on how the truth is built; report net benefit and categorical NRI. Script 30's polygenic-score truth
holds prevalence only to 0.2% to 1.0% (addendum 23); read the comparison between truths from script 31.
