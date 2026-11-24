# Journal submission (Journal of Infectious Diseases, Major Article)

A separate Quarto project (its own `_quarto.yml`), not part of the website. Everything here is built from the files in `output/`; nothing is typed in by hand except the prose.

| File | What it is |
|---|---|
| `PROSE_manuscript.docx` | Main manuscript: title page (running title, 40-word summary, keywords, word count), structured abstract, text, Notes, references, figure legends, Tables 1 and 2 on separate pages. Double-spaced, line and page numbers. |
| `PROSE_supplement.docx` | Supplementary material: full methods, results and further limitations, Tables S1 to S5, Figures S1 to S6 embedded. |
| `figures/Figure_<n>.{svg,pdf,png}` | Main-text Figures 1 to 5 as separate files (vector SVG and PDF, PNG at 300 dpi). |
| `figures/Figure_S<n>.{svg,pdf,png}` | Supplementary figures (also embedded in the supplement). |
| `figures/PROSE_figures_main.pdf` | Figures 1 to 5 in one PDF, one per page, for reviewers. |
| `figures/manifest.csv` | Source script and print size of every figure. |
| `manuscript.qmd`, `supplement.qmd` | Sources of the two documents. |
| `reference.docx`, `make_reference_docx.R` | Word style template (Times New Roman 12 pt, double spacing, line and page numbers) and the script that builds it. |

## Rebuild (from the project root)

```sh
Rscript scripts/54_submission_figures.R          # figures (after the figure scripts 27, 34, 44, 52)
Rscript submission/make_reference_docx.R         # only if the Word template changes
quarto render submission                         # both .docx files
```

`submission/.Rprofile` loads the project's renv library when Quarto starts R in this folder. Code chunks run from the project root (`knitr: opts_knit: root.dir`).

## Keeping it in step with the website draft

The prose of `manuscript.qmd` here is a shortened version of the website's `manuscript.qmd` (JID limit: 3,500 words including the abstract; the count on the title page is computed from the file). The supplement carries the website draft's detailed text. A rerun that changes a number means editing the prose in all three files. The PRISMA wording follows addendum 46 (622 automated exclusions, 209 screened, 208 full texts, 81 included, 60 in the evidence base).

## Still to do before submission

- Placeholders in square brackets: authors, affiliations, corresponding and alternate corresponding authors, ethics approval and consent, author contributions, funding, conflicts of interest, AI-use statement, archive DOI, data access terms, prior presentation, the extraction and round-2 verification paragraph (supplement).
- JID house style uses US spelling (mold, hematologic); the text keeps British spelling for the copy editors to convert, or change it before submission.
- Figures 1, 3 and 4 were drawn 10 to 11 inches wide; at JID's full-page width (about 7 inches) their text prints at about 7 to 8 pt. Redraw at 7 inches if the journal asks.
- A cover letter is not included.
