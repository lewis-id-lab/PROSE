# Current status

Last updated 20 October 2026. Sections 1, 6c, 7 and 8 describe the current position. Sections 2 to 5 and 6 describe the evidence base and the baseline as they stood on 16 to 18 October 2026 and are kept as the working record; later changes to them are in 6c and in the dated addenda of `phase3_decisions.md`.

---

## 1. Where the project stands

| Phase | State |
|---|---|
| 0 Set-up | complete |
| 1 Literature search | complete — PubMed and Embase, 831 unique records screened |
| 2 Screening and extraction | **round 1 complete including Jobs A and B and the supplementary-file pass.** The user reports (19 October) that the human round-2 checks confirmed the extraction against the source PDFs; the record is not in the repository. Direction errors in two studies (Fischer 2016, and one row of Sainz 2012) were found and corrected from their papers on 19 October; the direction check of the other text-parsed rows is closed (addenda 18 to 20, `data/direction_check_log.csv`): the full texts of Chen 2019 and Fisher 2017 confirmed their rows |
| 3 Parameter derivation | **done provisionally** for the four-term panel (PTX3, CLEC7A with two parameters, TLR4, CD209); IL10 and IFNG dropped (D14). One to five strand-resolved or odds-ratio studies per term, so every effect is provisional |
| 4 Baseline | **done**: the refit baseline on the 4,078-admission base cohort (D9, D10), MCID +0.00458 at 5% and +0.00160 at 20% |
| 5 Genetic layer simulation | **first results in** (addenda 11 to 22 and 24 to 26): exact layer, sampling variability, correlation ensemble, attribution, sensitivity cohorts, donor genotype in the allo-HSCT stratum, Aim 3b, the second CLEC7A parameter, the East Asian check and the tertiary NRI comparison. Addendum 24 (20 October) made a profile-preserving truth primary, with the panels, the Aim 4 sweeps and the NRI rerun under it, and addenda 25 and 26 reran the sampling, ensemble, attribution, Aim 3b, donor and robustness analyses and the all-ancestry CLEC7A variant under it too, so every analysis is on it |
| 6 Sensitivity and architecture sweep | **Aim 4 sweep done** (addenda 16 and 17, exact and with estimated weights); sensitivity analyses done for the cohort, the baseline equation, the correlation structure and the sampling variability |
| 7 Write-up | **full first draft** (19 October 2026, revised 20 October for the profile-preserving truth) on the Draft manuscript page (`manuscript.qmd`, references in `references.bib`, Vancouver style); authors, ethics, data availability, funding and two references still to add. The analysis site is live at https://lewis-id-lab.github.io/PROSE/ (repository lewis-id-lab/PROSE, public since 19 October 2026) |

## 2. Search and screening

PubMed targeted pool 267 records (264 unique); Embase (Ovid, Embase-only) 785 records yielding 586
Embase-unique; after title deduplication the combined screening set was **831 unique records** — 263 from
PubMed and 568 from Embase, of which 57 have no abstract and 175 are conference abstracts. Raw per-database
exports were uploaded to Covidence rather than pre-deduplicated files, so Covidence's own duplicate count
stays meaningful.

Two audit-trail items remain open: the date the PubMed searches were run is still unrecorded (Embase was
8 October 2026), and the complete automated screening decision file for all 831 records has not been
assembled.

## 3. Extraction — round 1, with the Job A full-text pass applied

A recall audit found that the original 40,000-character text cap had truncated 38 of the 63 extracted
studies. Fourteen eligible studies were re-extracted from the **full uncapped text** with chunked output
(Job A, 16 October 2026); their records replace the originals.

| | n |
|---|---|
| Included records | 81 |
| Extracted from a usable source | **63** |
| Contributing at least one variant row | 61 |
| Variant rows | **1449** (983 before Job A, 1,420 after it) |
| Distinct rsIDs | **237** (was 189) |
| Rows meta-analysable (estimate or counts, with alleles) | **757** |
| Rows with full six-cell genotype counts | **254** (was 139) |

Rows per re-extracted study, before and after:

| study | paper | rows before | rows after |
|---|---|---|---|
| S003 | Cunha, Cristina 2017 | 9 | 13 |
| S049 | Skonieczna, Katarzyna 2017 | 18 | 173 |
| S054 | Matzaraki, Vasiliki 2022 | 5 | 6 |
| S080 | Fisher, Cynthia E. 2017 | 38 | 51 |
| S027 | Zhong, Jinjin 2022 | 31 | 64 |
| S071 | Sainz, Juan 2012 | 14 | 43 |
| S034 | Konuma, Takaaki 2026 | 3 | 10 |
| S042 | White, P. Lewis 2018 | 9 | 9 |
| S073 | Fischer, Mike 2016 | 13 | 17 |
| S057 | Lupiañez, Carmen B. 2016 | 74 | 100 |
| S075 | Lupiañez, C. B. 2020 | 56 | 54 |
| S044 | Gonçalves, Samuel M. 2021 | 4 | 6 |
| S076 | Sánchez-Maldonado, Jose Manuel 2020 | 41 | 134 |
| S026 | Tang, Tiantian 2020 | 18 | 90 |

Field completeness across the 1420 rows: valid rsID 1139, effect estimate
893, both confidence limits 799, p value
1206, both alleles 886, full genotype counts
254, flagged primary analysis 685.

**Nothing has been verified against a source PDF by a human.** All 1420 variant rows and the 63 study
rows carry `second_checker = PENDING ROUND 2 VERIFICATION`.

### What Job A fixed, and what it did not

PTX3 is transformed. rs2305619 went from 15 poolable rows and 2 with genotype counts to 45 and 32;
rs3816527 from 20 and 8 to 39 and 31; rs1840680 from 6 and 1 to 36 and 30. Most of that comes from
Tang 2020, whose full genotype tables sat beyond the old character cap.

**TLR4 is not fixed and cannot be fixed by extraction.** rs4986790 now has 9 studies and 2 poolable rows;
rs4986791 has 9 studies and 1. Fisher 2017 analysed both but published no estimate or counts for either —
they are in its supplemental Table 1, which is not in the supplied PDF. Fischer 2016 did not genotype them.
Konuma 2026 notes both are close to monomorphic in Asian populations. The binding constraint is the
missing supplementary files, not the extraction, so **TLR4 needs supplementary retrieval or author contact**
before the haplotype block can be parameterised.

Four studies have supplementary tables that are absent from the supplied PDFs and hold genotype counts or
null estimates: S049 (Suppl. Tables III and IV), S080 (supplemental Table 1), S027 (Table S1), S071
(Table S1).

### Provisional rows

33 rows were transcribed by the extracting agents directly from the full text after their model-call
budgets were exhausted, rather than produced by a model pass: S075's 12 ARNT2/CX3CR1 haplotype rows read
from Table 3, S073's 7 reference-group rows whose counts are derived by subtraction from printed category
totals, S034's 6 TLR rows, S027's 3 Table 6 haplotype rows plus the rs7525979 excluded-variant row,
S057's 2 genotyped-but-excluded rows, and S080's 2 analysed-but-no-estimate rows. Each carries
`extraction_method = agent_read_not_model_extracted` and a note naming it as such. All 33 are flagged P1 in
the round-2 list and should be checked first.

### Three variant-identification errors resolved from source

Round 1's flagged "missing" rsIDs included several that are typographical errors in the papers themselves,
now confirmed and captured under the correct identifier: `rs7931529` is a Table 3 column-header typo in
Lupiañez 2020 for rs7631529; `rs6744330` is a Results-text typo in Gonçalves 2021 for rs674430;
`rs381627` is a table-footnote typo in Tang 2020 for rs3816527; `rs47933` in Zhong 2022 is a truncation of
NLRC4 rs479333. Separately, rs6535939 and rs61329579 are retired rs numbers for TLR2 rs4585282, stated as
such in Skonieczna 2017's own Discussion.

Cunha 2017's 25 apparently-missed rsIDs were correctly **not** extracted: they appear only in the
"SNPs tagged" and LD-proxy columns of its supplementary tables. The paper genotyped five variants and five
were captured. Round 1 was right there.


### Job B — targeted recovery of the panel variants (16 October 2026)

Job B went back to 19 studies for 66 specific study-by-variant cells where a candidate panel variant was
reported but lacked either per-genotype counts or an effect estimate with an orientable allele. The
remaining 28 cells sat in studies Job A had already rebuilt from full text.

| status | cells |
|---|---|
| `recovered_counts` | 28 |
| `recovered_estimate` | 6 |
| `recovered_both` | 1 |
| `not_reported_in_text` | 22 |
| `in_supplement_not_available` | 9 |

131 rows were added, replacing 102. The table now holds **1449 variant rows**, of which
**757 are meta-analysable** — meaning they carry either a complete effect estimate with an orientable
effect allele, or a complete six-cell genotype table with both alleles named.

**The important finding is that the TLR4 problem was partly a question of looking for the wrong field.**
Three studies (Ceesay 2016, Kalkanci 2020, Grube 2013) print complete 2x3 genotype tables for both TLR4
SNPs and no odds ratio at all — they report only a chi-squared or Fisher p. Two more (Aki 2015,
Tanpaibule 2020) add partial or complete counts. An odds ratio under any genetic model can be computed
from those counts, so these are not author-contact targets.

The structural catch is that **every one of these cohorts has zero minor homozygotes in both cases and
controls** for rs4986790 and rs4986791 — for example Grube 2013 rs4986790 is 0/6/35 in cases and 0/7/100
in controls. Only a dominant or carrier model is estimable, and a continuity correction will be needed.
Tanpaibule 2020's rs4986790 is degenerate in the other direction (cases 33/0/0), which is why the paper
prints a dash rather than an estimate.

### What the panel looks like after both jobs

Counting a study as contributing if it supplies either an orientable estimate or a complete genotype table:

| rsID | gene | studies reporting | via estimate | via counts | meta-analysable |
|---|---|---|---|---|---|
| rs16910526 | CLEC7A | 12 | 5 | 2 | **6** |
| rs7309123 | CLEC7A | 10 | 7 | 4 | **8** |
| rs2305619 | PTX3 | 9 | 5 | 4 | **6** |
| rs4986790 | TLR4 | 9 | 1 | 5 | **5** |
| rs4986791 | TLR4 | 9 | 1 | 4 | **4** |
| rs3816527 | PTX3 | 7 | 4 | 3 | **5** |
| rs7248637 | CD209 | 6 | 4 | 1 | **4** |
| rs4804800 | CD209 | 5 | 4 | 2 | **4** |
| rs5743708 | TLR2 | 5 | 1 | 2 | **3** |
| rs1840680 | PTX3 | 4 | 3 | 2 | **3** |
| rs1800896 | IL10 | 3 | 3 | 2 | **3** |
| rs2069705 | IFNG | 3 | 2 | 1 | **3** |
| rs5743810 | TLR6 | 3 | 2 | 2 | **3** |

13 variants have three or more meta-analysable studies and 8 have four or more. That is the realistic
panel, and it is a materially different list from the provisional one: **TLR9 rs352140 has only one**
usable study despite appearing in three, while TLR4 is recoverable after all, via counts rather than
published estimates.

### What neither job could fix

- **In a supplement not supplied (9 cells).** Schnetzke 2024 Table S2 holds the counts for six panel variants; Herrero-Sanchez 2018 Tables S1-S5 hold one; Brunel 2018's per-SNP PTX3 results are in Online Supplementary Figure S2, the paper printing only a combined `rs2305619 and/or rs3816527` hazard ratio of 5.06 (1.68-15.2) in the non-neutropenic subgroup. These are supplement-retrieval targets.
- **Estimate printed but no nucleotide allele anywhere (10 cells).** Schnetzke 2024 and Zhao 2022 name their SNPs only by amino-acid change ("TLR2 Arg753Gln", "Dectin-1 Y238X") or by bare rsID in the table. Inferring the nucleotide from a protein-level name needs outside knowledge, so these were deliberately left null. Author contact, or drop.
- **Reported as a null with no numbers (3 cells).** The Schnetzke 2017 abstract genotyped all three PTX3 SNPs and reports only that no association was found. Author contact is the only route.
- **Data exists only inside a figure (2 cells).** Koldehoff 2013's entire TLR4 result is in Figure 1; the bar-chart labels do not survive text extraction.

### Three provenance caveats carried forward

Chai 2011's genotype counts are partly derived rather than printed; Cunha 2010's rs16910526 G/T alleles
were taken from a reference the paper cites rather than from the paper itself; and the TLR4 alleles in
Kalkanci 2020 and Grube 2013 were read from the papers' own "896A>G" notation rather than from a genotype
column. All three are flagged P1 in the round-2 list.

### Supplementary files (16 October 2026)

Six of the seven requested supplements were supplied and extracted: Fisher 2017, Brunel 2018,
Schnetzke 2024, Sainz 2012, Herrero-Sanchez 2018 and Zhong 2022. Only Skonieczna 2017 is outstanding.
**559 rows** were added, taking the table to **2008 variant rows**.

| study | rows | distinct rsIDs | six-cell counts | orientable estimate | allele pair |
|---|---|---|---|---|---|
| S071 Sainz 2012 | 102 | 27 | 50 | 97 | 102 |
| S027 Zhong 2022 | 300 | 15 | 75 | 44 | 60 |
| S080 Fisher 2017 | 75 | 25 | 0 | 0 | 75 |
| S029 Brunel 2018 | 30 | 2 | 0 | 0 | 0 |
| S045 Schnetzke 2024 | 30 | 11 | 0 | 0 | 26 |
| S035 Herrero-Sanchez 2018 | 22 | 2 | 12 | 16 | 22 |

**Two expectations were wrong and are corrected here.**

Fisher 2017's supplemental Table 1 is **genotyping and imputation quality control, not effect estimates** —
alleles, minor allele frequency, call rate, Hardy-Weinberg p and pass/fail across three cohorts for 25
SNPs. It does not advance TLR4 on the estimate route. What it does deliver is 25 nucleotide allele pairs
with per-cohort frequencies, which makes it the review's single largest allele source, since many SNPs
carry effect estimates elsewhere in the table with no alleles attached. Two cautions: the Alleles column
carries no footnote defining major/minor order, and it disagrees with Schnetzke 2024's explicitly
major/minor listing for rs3816527 (C/A against A/C), so the order cannot be assumed. Five SNPs failed QC
in Cohort 1 only — rs4833095, rs767455, rs3901533, rs4804800, rs7252229.

Schnetzke 2024's Table S2 gives three-genotype distributions for the stratification (n=186) and validation
(n=138) cohorts but **does not split them by invasive fungal disease status**, so no 2x2 can be
reconstructed and the per-genotype case and control counts remain unavailable. `genotype_counts_available`
is corrected to `no` for that study. rs1840680 does not appear in Table S2 at all, only in a footnote
stating complete linkage disequilibrium with rs2305619.

**Brunel 2018 delivered what was expected**, and it bears directly on Aim 3b. Supplementary Figure 2
reports rs2305619 and rs3816527 *separately*, which the main paper did not, in three strata: all episodes,
no neutropenia at presentation, and neutropenia at presentation. The non-neutropenic result that the
combined hazard ratio of 5.06 came from decomposes into rs2305619 homozygous 8/33 versus other 7/100
(P=0.008) and rs3816527 homozygous 9/37 versus other 5/92 (P=0.002), against P=0.7 and P=0.6 in the
neutropenic stratum. No per-SNP hazard ratio is printed anywhere — all three Cox models in Supplementary
Table 1 use the combined "and/or" variable.

**Herrero-Sanchez 2018 carries a scope caveat recorded on every row**: none of its supplementary tables
uses proven or probable invasive mould disease as the endpoint. They use *Candida* and/or *Aspergillus*
detection, fungal colonisation, or *Candida* colonisation. Those rows may not be eligible for the primary
synthesis.

### A third usability tier

Brunel's rows are complete 2x2 tables under a stated recessive model, not three-genotype tables, so the
original test did not recognise them. Rows now carry `usable_tier`:

| tier | rows | meaning |
|---|---|---|
| `estimate` | 823 | effect estimate, both confidence limits and the effect allele |
| `counts_2x3` | 115 | all six genotype cells with both alleles named |
| `counts_2x2` | 32 | complete 2x2 under a stated dominant, recessive, carrier or allelic model, with the effect allele named or identified as the minor allele |

The 2x2 tier is weaker: where the comparison names only "the minor allele", the nucleotide has to be
resolved against a frequency reference before the row can be oriented. Treat it as a separate sensitivity
stratum rather than pooling it blind with the other two.

### The panel after the supplements

| rsID | gene | studies reporting | via estimate | via 2x3 counts | via 2x2 counts | meta-analysable |
|---|---|---|---|---|---|---|
| rs16910526 | CLEC7A | 12 | 6 | 3 | 1 | **7** |
| rs7309123 | CLEC7A | 10 | 7 | 1 | 0 | **8** |
| rs2305619 | PTX3 | 9 | 5 | 1 | 1 | **7** |
| rs4986790 | TLR4 | 9 | 1 | 4 | 0 | **5** |
| rs4986791 | TLR4 | 9 | 1 | 3 | 0 | **4** |
| rs3816527 | PTX3 | 7 | 4 | 2 | 1 | **6** |
| rs4804800 | CD209 | 6 | 4 | 0 | 0 | **4** |
| rs7248637 | CD209 | 6 | 4 | 0 | 0 | **4** |
| rs1840680 | PTX3 | 5 | 3 | 0 | 0 | **3** |
| rs5743708 | TLR2 | 5 | 1 | 2 | 0 | **3** |
| rs1800896 | IL10 | 3 | 3 | 0 | 0 | **3** |
| rs2069705 | IFNG | 3 | 2 | 1 | 0 | **3** |
| rs5743810 | TLR6 | 3 | 2 | 1 | 0 | **3** |

13 variants reach three or more meta-analysable studies and 8 reach four or more — unchanged in count
from before the supplements, but rs16910526 gains a study and the two PTX3 tags gain Brunel. TLR9 rs352140
still has one usable study of three and should be dropped.

Those counts are of reporting studies, not of independent cohorts. Collapsing same-group series
(`cohort_independence.csv`) costs rs7309123 two cohorts (8 to 6) and rs4804800 and rs7248637 one each
(4 to 3). **After that adjustment 13 variants still reach three or more independent cohorts, but only six
reach four or more**: rs16910526 (7), rs2305619 (7), rs7309123 (6), rs3816527 (6), rs4986790 (5) and
rs4986791 (4). Use the adjusted figure when judging how much weight a pooled estimate can carry.

### Flagged for round 2

27 Brunel rows are marked `DERIVED BY SUBTRACTION` (non-case counts computed as N at risk minus failures
from the figure panels) and are the only derived numbers in the file. Brunel's Supplementary Figure 1 fails
an arithmetic check — its B and C row totals sum to 212 against a stated column total of 167 — and all
three of its panels are flagged `AMBIGUOUS`; Figures 2 and 3 reconcile exactly. Schnetzke's Table S2 writes
the dectin-1 SNP as rs72309123 where Table S1 writes rs 7309123; both are recorded as printed rather than
silently merged. Herrero-Sanchez's Table S6 caption says rs2078178 while a stray body label says
rs16910526, and its Table S2 prints two different counts (89 and 105) for the same control group at the
same percentage. Zhong's Tables S2 and S3 print count triples with no genotype labels, so 240 of its 300
rows have no effect allele; the labels are recoverable by cross-referencing Table S1 but that inference was
not made.

## 4. Source provenance

PDFs were originally matched to records by filename similarity, which failed for 13 of 81 records —
clustered on conference abstracts, because an abstract's PDF is often the whole proceedings book (one file
is 8.4 million characters). Provenance was re-derived by contiguous-title matching and the extractions
corrected:

- **56 verified** — the record's title is on the title page of the document used.
- **7 re-extracted** from the correct source (S010, S014, S018, S021, S023, S050, S051).
- **6 voided**, rows deleted, `extraction_status = voided_wrong_source`: S011, S012, S015, S016, S020, S025.
- **12 with no PDF supplied**: S007, S008, S009, S013, S017, S019, S022, S024, S032, S047, S048, S068.

Several of those 18 look like conference abstracts of papers already included as full publications
(S022/S072, S019/S075, S023/S073, S011/S010, S012/S016). If confirmed they are duplicate reports of the
same cohort rather than missing data, which changes the PRISMA numbers.

## 5. Variant annotation

All 239 rsIDs were queried against dbSNP on 16 October 2026 and the reported alleles compared with the
forward-strand alleles.

| orientation_flag | rows |
|---|---|
| fwd_ok | 701 |
| **strand_flip** | **188** |
| palindromic_unresolvable (A/T or C/G) | 130 |
| single_allele_only | 234 |
| no alleles reported | 444 |
| mismatch (allele absent from dbSNP) | 4 |
| no valid rsID | 306 |

These cover all 2008 rows. The classifier was corrected on 16 October 2026: it had tested reported
alleles against the full dbSNP allele set, which marks every pair at a multi-allelic site as ambiguous, and
it had treated true palindromes as unambiguous because a C/G pair equals its own complement. Both are
fixed; the test now uses the alleles segregating in Europeans and flags palindromes explicitly. The
previous figures of 85 strand flips and 332 ambiguous rows were wrong in both directions.

Restricted to the 299 usable rows behind the 13 panel variants: fwd_ok 226, strand_flip
38 (in 10 study x variant cells), palindromic 25 (all
rs7309123, 7 cells), incomplete allele pair 9.


These counts cover all 2008 rows in the table, including the supplementary-file rows. Many of the
`no_alleles` rows are Fisher 2017 QC entries and Zhong 2022 count triples carrying no genotype label, so
the figure measures how much is not yet orientable rather than how much is wrong.


37 distinct variants have at least one row confirmed on the opposite strand from the European segregating pair. Within the 13-variant panel the affected variants are **rs16910526, rs1800896, rs2069705 and rs5743810** — 38 rows across 10 study x variant cells.

**PTX3 rs2305619 is not among them**, contrary to an earlier version of this document. Of its 70 rows, 59 are forward-strand, 10 report no allele, and one (Schnetzke 2024) is `strand_flip_dbsnp_only`: it reports A/C where the alleles segregating in Europeans are A/G. That pair is not a clean complement of the European pair, so the discrepancy is more likely a reporting error at a multi-allelic site than a strand flip. It needs checking against the source but it is not a confirmed orientation error.

**rs7309123 is resolved**: dbSNP places it at chr12:10,119,994 in CLEC7A, not FCN2. The one study calling
it an L-ficolin variant wrote it only in its concluding sentence while naming rs17514136 and rs17549193
throughout its Methods and Results — a typo, not a gene-assignment conflict.

## 6. Baseline model — verified, with one gap

Everything needed is in the Stanzani paper and has been re-read from it: the equation and intercept, each
predictor's prevalence separately in cases and non-cases, and the validation targets.

**The joint distribution of the seven predictors is not published, and it is load-bearing.** Sampling
200,000 admissions independently from the published marginals and applying the fixed equation gives an IMD
prevalence of **2.27% against 3.3%**, and an **aROC of 0.775 against 0.85**. High-risk chemotherapy and
prolonged neutropenia are close to the same clinical state and cannot be drawn independently.

A Gaussian copula with a single exchangeable latent correlation gives:

| latent rho | simulated prevalence | simulated aROC |
|---|---|---|
| 0.00 | 2.27% | 0.769 |
| 0.20 | 2.92% | 0.829 |
| **0.30** | **3.32%** | **0.852** |
| 0.40 | 3.75% | 0.872 |
| published | 3.3% | 0.85 |

rho = 0.30 reproduces both global targets almost exactly but misses the per-predictor case-group
prevalences — uncontrolled malignancy 74% simulated against 59% published, corticosteroids 38% against
23%, high-risk chemotherapy 54% against 64%. One exchangeable parameter spreads correlation evenly across
pairs that are clinically very differently related, so a structured matrix fitted to the 16 published
targets is needed, with clinical constraints to make it identifiable.

Worth requesting from Bologna: 2x2 cross-tabulations for chemotherapy x neutropenia, neutropenia x
lymphopenia and malignancy status x chemotherapy, plus the distribution of admissions per patient. That
would convert an assumption into data.

**Update, 17 October 2026 — correlation structure fitted.** A structured latent correlation matrix
(non-negative dependence, clinically grouped) was fitted to the 16 published targets. The first version had
eight parameters and reduced chi-square against the targets from 346 under independence and 47 under the
exchangeable rho = 0.30 to 8.6. Two problems remained: uncontrolled malignancy in cases was simulated at 69%
against 59% published (2.3 standard errors), and two parameters sat at their upper bound because the
chemotherapy-lymphopenia correlation was tied to the block baseline, which also set malignancy's links.
Version 2 adds a separate chemotherapy-lymphopenia parameter (nine in all) and lets the block extras reach
0.7. Chi-square falls to 0.36, the largest miss is 0.4 standard errors, and the malignancy residual is gone.
The fitted matrix has malignancy almost uncorrelated with the rest of the block (about 0.04) and
chemotherapy, neutropenia and lymphopenia at 0.69 to 0.74 with each other. `d_chemo_neut` sits at its 0.7
cap, but the loss is lowest there, so the cap is not costing fit. Positive definiteness is no longer
limiting (smallest eigenvalue 0.26).

Read the near-zero chi-square with care. The 16 targets collapse to about nine independent constraints, since
the non-case prevalences follow from the case prevalences and the overall prevalence, and the model now has
nine parameters. A saturated fit is expected to be this good and does not validate the structure. A 150-matrix
ensemble (six chains, R-hat at most 1.09) shows the parameters are wide: the chemotherapy-neutropenia extra
0.24 to 0.68, steroids-CMV 0.02 to 0.60, the malignancy-chemotherapy extra 0.01 to 0.42 (5th to 95th
percentiles). The
Bologna cross-tabulations, especially malignancy by chemotherapy, would turn these into data.

Decision curves for thresholds of 1 to 30% are now exact functions of the profile distribution, which comes
from numerical integration of the fitted copula (no patient-level sampling). On the fitted cohort the minimum
clinically important difference is **+0.00428** net benefit at the 5% threshold (ensemble 90% range 0.00397
to 0.00460) and **+0.00140** at 20% (0.00117 to 0.00167). Baseline net benefit at 5% is 0.01478, with 21.4
true and 125.7 false positives per 1,000 admissions and 14.7% flagged. The curves are step functions because
there are only 128 clinical profiles. See `scripts/01_fit_correlation.R`, `02_correlation_ensemble.R` and
`03_decision_curves.R`, and the baseline and decision-curve pages of the report website. The first-version
outputs are kept in `output/v1/`.

**Update, 18 October 2026 — baseline cohort switched to the empirical joint (D8).** The admission-level derivation data
(`bs.xlsx`, 4,127 admissions) reproduces Table 2 apart from two non-case counts that are one lower than published, and a
refit recovers the published coefficients. Its real joint distribution shows the fitted copula is 0.098 away in total
variation distance, with several correlations wrong in sign, so the copula is now a sensitivity analysis and the
baseline cohort is the observed distribution of the 128 profiles (75 observed). The MCID rule is unchanged; its values
on the empirical cohort are **+0.00445 at 5%** (bootstrap 90% 0.00412 to 0.00478) and **+0.00142 at 20%**, replacing
+0.00428 and +0.00140. Admissions per patient and the rule defining the 4,127 admissions are still needed from
Bologna. A full BOSCORE export with identifiers (8,660 rows) cannot be linked to the 4,127 admissions; a pseudonymised
extract of it is in `data/restricted/` (gitignored). See `phase3_decisions.md` addendum 2 and `status.md`.

**Update, 18 October 2026 — base dataset reconciled and baseline cohort corrected (D9).** `boscore.xlsx`, the parent of
`bs.xlsx`, reconciles exactly: the seven predictors are the `BOSCORE_*` columns, the 4,127 are the rows complete on them,
and `id` is the key. 49 admission numbers appear on two rows (an event row and a non-event twin with the same patient and
profile), so the 4,127 rows are 4,078 unique admissions, and 133 of 4,078 is the paper's 3.3%. The baseline cohort is now
the 4,078 admissions in 1,945 patients, with patients resampled in the bootstrap. The MCID is **+0.00404 at 5%** (patient
bootstrap 90% 0.00370 to 0.00438) and **+0.00113 at 20%**, against +0.00445 and +0.00142 on the 4,127 rows; results will be
reported against both. Admissions per patient: 2.10 on average, 1,000 patients with one admission. The 49 twin pairs are confirmed as an artefact
(addendum 5) and no request to Bologna remains open. See `phase3_decisions.md` addenda 3 to 5 and `status.md`.

**Update, 18 October 2026 — outcome-risk sensitivity (addendum 4).** Scoring the published equation's decisions against a refit
truth, a refit rule and the observed outcomes gives a primary MCID of +0.00452, +0.00458 and +0.00441, against the adopted
+0.00404, and a secondary MCID of +0.00137, +0.00160 and +0.00142 against +0.00113. The adopted value is the lowest of the
four, the intervals overlap, and the published coefficients lie inside the refit's intervals. No decision changed; results
will be reported against the 4,127-row MCID as well, which brackets these. See `phase3_decisions.md` addendum 4.

**Update, 18 October 2026 — the baseline model is refit on the base cohort (D10).** The clinical block is now a maximum-likelihood
refit of the seven predictors on the 4,078 unique admissions (patient-clustered SE, patient bootstrap), replacing the published
equation, which becomes a sensitivity analysis. aROC is 0.858 apparent and 0.853 optimism-corrected, the calibration slope is
0.97, and every published coefficient lies inside the refit's 90% interval. The MCID is **+0.00458 at 5%** (bootstrap 90% 0.00387
to 0.00566) and **+0.00160 at 20%**, against +0.00404 and +0.00113 under the published equation. Aim 2 becomes estimating the
baseline and comparing it with the published equation. See `phase3_decisions.md` addendum 6 and `status.md`.

**Update, 18 October 2026 — PTX3 defined and the evidence mapping started (D11, addendum 7).** `h2` is the PTX3 chromosome class
carrying A at rs2305619 (A-C-A, A-A-A, A-A-G; 48.2% in gnomAD NFE); h1 is G-A-G; rs1840680 A is the same class and rs3816527 C
marks only A-C-A. Each study's contrast is converted to a per-copy h2 effect by solving a haplotype model, one effect is taken per
study, and the effects are pooled with a Bayesian random-effects model. Provisional PTX3 result from five studies: OR 1.35 per
h2 copy (90% interval 0.99 to 1.81), tau 0.17. The rows are unverified and most mapping rules and the priors are proposals open
to change. See `phase3_decisions.md` addendum 7 and `status.md`.

**Update, 18 October 2026 — rows checked against sources, and the mapping extended to CLEC7A and TLR4 (addendum 8).** Herrero-Sanchez
2018, read in full, has its reference label reversed in every table; the extraction followed it, so the study looked as if it pointed
the other way, which was an artefact (30 rows corrected, and the addendum 7 statement retracted). Brunel 2018's design is confirmed
but its per-SNP counts cannot be checked from the text and may be the wrong homozygote group. Provisional primary results, as odds
ratios per copy of the risk haplotype: PTX3 h2 1.35 (0.99 to 1.81), CLEC7A Y238X 0.98 (0.56 to 1.72), TLR4 G-T 1.44 (0.77 to 2.69).
The second CLEC7A effect needs the strand of rs7309123 resolved for each study. See `phase3_decisions.md` addendum 8 and `status.md`.

## 6b. Phase 3 gate — settled 16 October 2026

Recipient genotype is the primary exposure, with donor genotype as a pre-specified secondary restricted to
the allogeneic HSCT stratum. The panel is four haplotype terms — PTX3 (rs2305619/rs3816527/rs1840680),
CLEC7A (rs16910526/rs7309123, 3 haplotypes), TLR4 (rs4986790/rs4986791) and CD209 (rs4804800/rs7248637,
4 haplotypes) — with IL10 rs1800896 and IFNG rs2069705 optional (both dropped on 19 October 2026, D14). The minimum clinically important
difference is 20% more true cases of invasive mould disease at no increase in false positives. It was
stated as +0.00420 net benefit at the 5% threshold and +0.00191 at 20% on an indicative cohort, and was
re-derived on the calibrated cohort on 17 October 2026 as +0.00428 at 5% (primary) and +0.00140 at 20%
(secondary). Full record in `phase3_decisions.md`. All sub-decisions are now settled (17 October 2026): 5% is the primary
MCID threshold; TLR9 rs352140 is dropped (D3a); Herrero-Sanchez 2018 is excluded from the primary analysis
and run as a sensitivity analysis (D3b); TLR4 is pooled under a carrier model with a continuity correction
(D3c); the genetic term is recentred to hold prevalence (D5); the baseline aROC target is 0.85 (D6); and the
haplotype phase is coupling, confirmed arithmetically (E1).

The CLEC7A and CD209 terms were collapsed from single SNPs to haplotypes after a linkage check:
rs16910526 x rs7309123 sit at D' = 1.000 with the Y238X allele entirely on one rs7309123 background, and
rs4804800 x rs7248637 at r2 = 0.413, D' = 0.671 with all four haplotypes present. Entering either pair as
independent additive terms would repeat the predecessor thesis's error within a different gene.

## 6c. The genetic layer, Aims 3b and 4 (18 to 19 October 2026)

Everything below is provisional: the effects come from the provisional evidence mapping, and the design choices are proposals recorded in addenda 9 to 17 of `phase3_decisions.md`. Numbers are after the correction of Fischer 2016 (below).

**Evidence mapping (addenda 7 to 9, 16).** Pooled odds ratios per copy of each term's risk haplotype class (primary analysis, 90% interval): PTX3 h2 1.35 (0.99 to 1.81, 5 studies), CLEC7A Y238X 0.84 (0.48 to 1.48, 4), TLR4 G-T 1.44 (0.77 to 2.69, 4),
CD209 G at rs4804800 1.80 (1.14 to 2.73, 3). IL10 (one study) and IFNG (two, disagreeing) were dropped from the panel. **The second CLEC7A parameter**, the A-G haplotype (the G allele of rs7309123 on the non-Y238X background), is estimated at **OR 1.9 per copy (1.26 to 2.63)** from three strand-resolved studies (Schnetzke 2024, Kalkanci 2020, Sainz 2012), and stays between 1.5 and 1.9 under every alternative strand assumption.
The strand of rs7309123 in each of ten studies was determined from control frequencies and the papers (`data/rs7309123_strand_resolution.csv`): resolved with medium to high confidence for five (Fisher, Schnetzke, Tanpaibule, and, from their full texts, Kalkanci 2020 and Chen 2019), low to medium for one, weakly for three, undetermined for one (Zhao 2022). Kalkanci's primers also confirm that its rs16910526 labels are on the opposite strand while its rs7309123 labels are forward.

**An error found and corrected (19 October, addendum 16 section 0).** Fischer 2016 lists the reference group first in its table; the pipeline read the first group as the numerator, so its Y238X and rs7309123 effects were reversed. 12 rows are corrected with their evidence, the Y238X pooled OR moved from 0.98 to 0.84, and every script using CLEC7A was rerun: the layer gains moved by at most 0.00006. 139 rows in the four terms' analyses were read from comparison text with no reference marker;
the five used in primary analyses were inspected and the studies behind them have since been checked against their tables. The direction check (addenda 18 to 20) found Tang 2020, Tanpaibule 2020, Chen 2019 and Fisher 2017 correct and one reversed row in Sainz 2012 (corrected); only Lupianez 2015's two IFNG rows were not checked. Two earlier errors of the same family: Herrero-Sanchez 2018's reversed reference label (addendum 8) and a first interaction model for Aim 3b that left out the neutropenia main effect (found and fixed during the run).

**The genetic layer (addenda 11 to 14).** On the base cohort with the refit baseline, genotype fixed per patient, each term's risk-class dosage Hardy-Weinberg in gnomAD NFE frequencies, the gain over the baseline computed exactly over 128 clinical profiles by the dosage combinations, for 2,000 draws of the true effects from the pooled posteriors. **Since addendum 24 (20 October) the primary truth is profile-preserving:** an intercept per clinical profile keeps the clinical score exactly calibrated, so a gain cannot be the correction of induced miscalibration. The first version (one global intercept that holds prevalence, D5) is kept as the sensitivity analysis, and its numbers are given where they differ. The panels, the Aim 4 sweeps and the NRI were rerun under it in addendum 24, and the sampling, ensemble, attribution, sensitivity-cohort, donor, Aim 3b and robustness analyses in addendum 25; the numbers below are the profile-preserving ones, with the global-intercept values where they differ:
- **Primary result (profile-preserving truth).** With the true effects in the score (the ceiling) the median gain is +0.0006 at 5% (+0.00056, 90% range 0.00015 to 0.00127; 1.1 extra true cases per 1,000 against 4.6 needed) and +0.0011 at 20%; with the pooled estimates in the score (realised), +0.0004 and +0.0008. The MCID is +0.00458 at 5% (reached in no draw under the pooled-mean reading and 0.6% under the predictive reading) and +0.00160 at 20% (reached in 17% of draws for the ceiling score, 4% for the realised score). Under the global-intercept truth the 5% gain was +0.00059 and the 20% MCID was reached in 16% and 3%.
- **Calibration check (addendum 24).** Under the global-intercept truth the clinical score was slightly miscalibrated within profiles (slope 0.98, largest error 0.023). The best clinical-only rule under that truth, the marginal risk, beat the clinical score by at most +0.00027 in net benefit in any draw, so the original gains were not a correction of that miscalibration. The profile-preserving truth lowers the panels' gains at 5% by 5% to 6% and leaves them within 2% at 20%; the 5% architecture frontier and the continuous NRI move most (below).
- **Sampling variability of a study of this size** (4,078 admissions, 133 events): the standard deviation of the observed gain is 0.0006 to 0.0012 at 5% (0.0006 to 0.0013 under the global-intercept truth), as large as the expected gain or larger; an observed study reaches the 5% MCID in at most 5% of repetitions in any run.
- **Robustness.** The 150-member correlation ensemble moves the gain by at most 0.00015 (between-member SD 0.0001 to 0.0003); the sensitivity cohorts (published equation on 4,078 or 4,127 rows, refit on 4,127 rows) change the gain little and change which MCID is reached at 20% (the published equation sets a lower bar).
- **Attribution (Shapley), four terms:** CD209 supplies about 51% of the gain, PTX3 23%, TLR4 19%, CLEC7A 7% (52%, 23%, 19%, 7% under the global-intercept truth).
- **Donor genotype in the allo-HSCT stratum (D1 secondary, addenda 13 and 14).** The base cohort now carries an admission-level allo flag (derived from the composite item SCORE9 and a transplant date; 292 admissions as in Table 1, 43 events against 44; 272 admissions and 43 events in the base cohort). Allo admissions are the high-risk chemotherapy item and only 23 of 128 profiles occur, so the baseline is nearly calibrated there (odds multiplier 1.15) and flags 77% of them at 5%.
  Donor genotype changes the true cases per 1,000 allo admissions by -0.1 to +1.3 at 5% against 30 needed (+0.4 to +1.4 under the global-intercept truth), and adds 14 to 27 at 20% against about 11; the donor evidence is one study each for PTX3 and CLEC7A, two for TLR4, none for CD209. A smoothed stratum gives the same conclusions.

**With the second CLEC7A parameter (addenda 16, 17 and 19).** CLEC7A becomes an active term. The gain roughly doubles (ceiling +0.0011 at 5%, +0.0018 at 20%, and the 20% MCID is reached in 64% of draws; +0.0012 and 63% under the global-intercept truth), and **CLEC7A supplies about 51% of the Shapley gain, almost all through the A-G haplotype, with CD209 falling to 27%, PTX3 to 12% and TLR4 to 10%** (50%, 28%, 13% and 10% under the global-intercept truth). This rests on three studies and is the least mature part of the evidence.

**CLEC7A's lead holds (addendum 18).** The five-term attribution and Aim 3b layer, rerun on the three sensitivity cohorts and on the 150-member correlation ensemble, give CLEC7A about 50% of the gain in every setting (44.8% to 53.3% over every reading, score and threshold) and CD209 25% to 32%; CLEC7A is the largest single contributor in 54% to 70% of draws. The neutropenic stratum supplies 77% to 83% of the gain at 5% and 96% to 97% at 20% in every setting (profile-preserving truth; 78% to 83% and 95% to 97% under the global-intercept truth).

**Aim 3b, the neutropenia-stratified falsification test (addendum 15; five-term version in addendum 17).** The expected direction was committed before the run. The Feys reading (effects larger in neutropenic admissions) and Brunel 2018's stratified data (larger in non-neutropenic, difference P = 0.07 and 0.049) point opposite ways.
Neutropenic admissions hold 80% of the events and supply 80% to 95% of the whole-cohort gain under homogeneous effects (81% to 96% under the global-intercept truth). A single score with the pooled estimates is the safe choice; a score built on the Feys reading is the most harmful when wrong. In simulated studies of the cohort's size an interaction test is better powered than predicted (58% to 89% for a Feys-sized modification, 94% to 100% for Brunel-sized ones, false positives 5%). With all of the effect in non-neutropenic admissions the oracle score's chance of reaching the 5% MCID falls from 18% to 2% (four terms) and from 42% to 6% (five) under the profile-preserving truth; no other Aim 3b result moves by more than a few per cent. That scenario puts an effect about five times the pooled one into a low-risk stratum; one global intercept flattens the clinical gradient there, which makes it a different population and not a correction of the comparator (addendum 26).
Of the five predictions, two held, one mostly held, one did not hold in general and one (that the test would be underpowered) was wrong.

**Aim 4, the genetic-architecture requirements sweep (addenda 16 and 17).** 1 to 20 identical independent loci x OR 1.1 to 3.0 x frequency 0.05 to 0.50, by exact convolution. Under the profile-preserving truth 40% of the 4,000 architectures reach the 5% MCID and 77% the 20% MCID (46% and 77% under the global-intercept truth; 241 verdicts at 5% are lost and none gained). At OR 1.5 or below no number of loci up to 20 reaches the 5% MCID (under the global-intercept truth one architecture, OR 1.5 at frequency 0.45 with 20 loci, did, by +0.000007); at OR 2.0 it takes 9 to 12 loci at common frequencies (7 to 9 under the global-intercept truth). The four-term panel is equivalent to about 7 loci of OR 1.3. The frontier is slightly non-monotone in frequency near 0.5, because the gain peaks at about 0.45.
With estimated weights (independent per-locus error, SD 0.05 to 0.40) precision of 0.1 or better costs almost nothing (a discovery sample of about 1,000 balanced subjects), while at 0.2 to 0.4, the precision of the present evidence, 12% to 36% of the architectures that reach the MCID with exact weights no longer do. The realised-score frontier also moves out under the profile-preserving truth (at error 0.40, OR 1.8, frequency 0.3 and the 5% threshold, 20 loci become none within 20).

**Polygenic score benchmark (addendum 23, written up 20 October).** How strong, and how accurately estimated, would a polygenic score (G ~ N(0, 1), independent of the clinical profile) have to be to reach the MCID? With true weights the odds ratio per SD needed is 3.6 at 5% (score-alone aROC 0.77) and 1.8 at 20%; under the profile-preserving truth it is 4.2 and 1.8. With estimated weights at 5%, an odds ratio of 6 per SD needs an accuracy of 0.74 and, for 10,000 effective loci, about 7,400 balanced discovery cases (16,000 under the profile-preserving truth; Fisher 2017 had 483 cases); an odds ratio of 4 is unreachable under it. At 20% the bar is much lower. Everything assumed is optimistic, and script 30's truth holds prevalence only to 0.2% to 1.0%, so its cases needed are accurate to a few per cent.

**Tertiary NRI comparison (addendum 22).** NRI computed exactly, with net benefit and aROC under the same truth; rerun under the profile-preserving truth in addendum 24 (global-intercept values in brackets). The four-term haplotype panel has a continuous NRI of about 34% (90% range 19% to 53%; 31%) while its net-benefit gain at 5% is an eighth of the MCID; the categorical NRI at 5% is 2.4% (2.6%). At the pooled PTX3 effect, the predecessor's structure (three SNPs simulated independently) gives a continuous NRI of 25% (25%), close to its reported 22.5%. The h2 haplotype gives 34% (13%), so the earlier finding that linkage halves the continuous NRI does not hold under the profile-preserving truth; the net-benefit finding does (+0.00033 for the independent structure against +0.00008). The three-SNP score applied to real, linked chromosomes is harmful (gain -0.0004, aROC -0.005) and has a positive continuous NRI of 29% (13%). The continuous NRI is unstable: for the h2 term it depends on how the truth is built (13% against 34%) with the same net-benefit gain, because heterozygotes sit almost exactly at their profile's mean risk, and its event component is -18% (-41%). The predecessor's own effects and the form of its NRI are not on record, so the match is consistency, not reproduction.

**Write-up (19 October 2026, revised 20 October).** A full first draft is on the Draft manuscript page: abstract, introduction, methods (ADEMP), results, discussion with limitations and the comparison with the 22.5% result, and Vancouver-style references. The 20 October revision made the profile-preserving truth primary (addendum 24, D27), rebuilt the headline table, the architecture section and the NRI table from it, and disclosed that the change was made after the global-intercept results were known.

**Where the work lives.** Scripts 01 to 35 in `scripts/` (scripts 17 to 21, 24, 26 and 27 take `PROSE_SCHEME=profile` and then write `_profile` files) (run order in `status.md`), the model code in `R/`, the test suite in `tests/testthat/`, results in `output/`, pinned packages in `renv.lock`. The repository is lewis-id-lab/PROSE (public); the Quarto site in `docs/` is published to GitHub Pages on every push to `main`. No patient-level data are in the repository; the identified files were read locally with the user's permission and only counts or a binary flag were kept (the pseudonymised base cohort is gitignored).

## 7. What is open

1. **Round-2 verification.** Add the record of the human checks to `data/` (the direction check is closed, addenda 18 to 20); decide the three questions left in addendum 8 (Herrero-Sanchez Table 5, Brunel's counts and homozygote group, hazard-ratio studies).
2. **rs7309123 strand.** Only Zhao 2022 is undetermined (no allele letters printed). East Asian haplotype frequencies would let the Thai and Chinese studies enter the second parameter.
3. **Cohort overlap.** The Spanish registry series (Sainz / Lupiañez / Sánchez-Maldonado) and the Jena AML cohorts (S018, S023, S073; S023 is treated as the same cohort as S073) are not independent studies and have only been partly adjudicated.
4. **Thin evidence.** One to five studies per term; donor genotype one to two studies per term and none for CD209; TLR4's supplementary retrieval (Fisher 2017's supplement) is still outstanding; the priors on the pooled effect and on tau are proposals.
5. **Choices to confirm.** The reading of the "realised" score (addendum 11), the independence of genotype from the clinical profile, the rule that the Brunel pattern applies to every term in Aim 3b, and the profile-preserving truth as primary (D27, adopted after the global-intercept results were known).
6. **Not done.** The interaction test of Aim 3b on the sensitivity cohorts and the ensemble; the pathogen-specificity scenario of Aim 4; locus-specific error, winner's curse and linkage in the sweep; a genotyped neutropenia-stratified analysis (the real test of Aim 3b). The write-up has a full first draft (authors, ethics, data availability, funding and two references to add).

## 8. Suggested order of work

1. Add the round-2 verification record to the repository (the direction check is closed).
2. Resolve what can be resolved of the three undetermined strands, and decide whether to obtain East Asian haplotype frequencies.
3. ~~Repeat Aim 3b and the attribution with five terms on the correlation ensemble and the sensitivity cohorts~~ Done (addendum 18): CLEC7A's lead and the neutropenia pattern hold.
4. ~~Decide the write-up structure from the results above~~ Done: full first draft on the Draft manuscript page, with the tertiary NRI comparison (addendum 22). Original item: (the headline is that a four-term panel falls well short of the 5% MCID, and what architecture would not), with the figures of addenda 16 and 17.
5. Plan the genotyped, neutropenia-stratified study that the Aim 3b design analysis shows is feasible at this cohort size.
