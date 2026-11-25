# Extraction checks

Tools for checking the round-1 extraction (`extraction_variants_round1.csv` in the bundle, 2,008 rows) against its source papers. Nothing here edits the bundle or `data/extraction_corrections.csv`. Every check adds flags or proposals; a person with the PDF open makes each change.

**Internal consistency is not verification.** Two extractions that agree can both be wrong. A clean result lowers the chance of an error. It does not prove a value.

## What is and is not documented about round 1

The round-1 extraction script, its prompts, the model name and the extraction dates are not in this repository or the bundle. The `extraction_method` column records the pass for each row (round1, round1b, jobB, supplement, agent_read_not_model_extracted). The model name in the manuscript has to come from your own records. The re-extraction below is fully logged, so the verification pass can be described exactly.

## The three layers

1. **Arithmetic QC** (`qc_arithmetic.py`, no API, runs anywhere). Recomputes odds ratios from genotype counts under six contrast forms, checks each estimate against its interval, tests Hardy-Weinberg in controls and compares control allele frequencies with gnomAD. Writes `output/qc_*.csv` and, with `--workbook PATH`, a reviewer workbook. It can only see rows with both counts and an estimate.
2. **Blind re-extraction** (`run_agents.py extract`). An agent reads each source PDF and extracts the listed variants without seeing any round-1 value. Prompt: `prompts/extract_system.md`. Allowed values come from the bundle's `codebook.csv` (`schema.py`). Output is schema-constrained JSON.
3. **Judge** (`compare.py`, then `run_agents.py judge`). Rows are paired with the corrected round-1 table. Alleles named the other way round are reconciled before comparing. Every difference goes to a second agent, which reads the PDF, is shown the two versions as unlabelled A and B, and says which one matches the paper, quoting the text it read. `compare.py --workbook` then unblinds the verdicts and writes the reviewer workbook and `proposed_corrections.csv` (corrections-log format, `checked_on` and `checked_by` left empty for the human check).

## Running it (on the computer that holds the PDFs)

```sh
pip install anthropic pypdf openpyxl
export ANTHROPIC_API_KEY=...            # or: ant auth login
python3 extraction/run_agents.py extract --pdf-dir ~/path/to/pdfs --scope panel --run pilot --dry-run
python3 extraction/run_agents.py extract --pdf-dir ~/path/to/pdfs --scope panel --run pilot
python3 extraction/compare.py --run pilot
python3 extraction/run_agents.py judge --pdf-dir ~/path/to/pdfs --run pilot
python3 extraction/compare.py --run pilot --workbook
```

- `--scope panel` covers the 13 panel variants (the pilot). `--scope all` covers every variant in every study with a PDF.
- `sources.csv` maps studies to PDF file names, starting from the bundle's `pdf_source_audit.csv`. A `pages` range picks one abstract out of a conference supplement (S018, S021, S023). Eleven rows carry a note: the file is shared with another study or round 1b used a corrected source. Confirm those file names before running. Add a second row for a study to include its supplement.
- The default model is `claude-opus-5-5` at effort `high`. Server-side refusal fallback is on. `run_log.jsonl` records the model that actually answered each request, the prompt hash, the token usage and the time.
- A run can be resumed: studies already done are skipped unless `--redo` is given.
- Rough cost (an estimate, not measured): $10 to $20 for the panel pilot, $40 to $60 for every study, both passes included.

The run folders (`runs/<name>/`) hold the extracted values and the judge's quotations. Commit them once reviewed, as the record of the verification pass. The PDFs themselves never go in the repository.
