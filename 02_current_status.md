# Current status

As at 16 October 2026.

---

## 1. Where the project stands

| Phase | State |
|---|---|
| 0 Set-up | complete |
| 1 Literature search | complete — PubMed and Embase, 831 unique records screened |
| 2 Screening and extraction | **round 1 complete including Jobs A and B and the supplementary-file pass; round 2 not started** |
| 3 Parameter derivation | not started — gated on round 2 |
| 4 Baseline reproduction | **in progress** — baseline coded, test harness in place and correlation structure fitted (17 October 2026); decision curves across thresholds and the Bologna cross-tabulations outstanding |
| 5 Genetic layer simulation | blocked on phase 3 |
| 6 Sensitivity and architecture sweep | blocked on phase 5 |
| 7 Write-up | not started |

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

## 6b. Phase 3 gate — settled 16 October 2026

Recipient genotype is the primary exposure, with donor genotype as a pre-specified secondary restricted to
the allogeneic HSCT stratum. The panel is four haplotype terms — PTX3 (rs2305619/rs3816527/rs1840680),
CLEC7A (rs16910526/rs7309123, 3 haplotypes), TLR4 (rs4986790/rs4986791) and CD209 (rs4804800/rs7248637,
4 haplotypes) — with IL10 rs1800896 and IFNG rs2069705 optional. The minimum clinically important
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

## 7. What blocks the genetic layer

1. Round-2 verification of the extraction — not started.
2. Allele orientation: 188 rows across the whole table are reported on the opposite strand from the European segregating pair, and 130 are palindromic and unresolvable from alleles alone. For the panel this reduces to 10 strand-flip cells and the 7 rs7309123 cells — see `panel_verification_worklist.csv`.
3. Cohort-overlap adjudication: the Spanish registry series (Sainz / Lupiañez / Sánchez-Maldonado, ten records 2007–2020) and the Jena AML cohorts (S018 n=186, S023 n=164) are not independent studies.
4. Final panel derivation from `panel_candidates.csv` — 13 variants now have three or more meta-analysable studies. Four provisional-panel entries — rs5743836 (TLR9), rs1800450 (MBL2), rs1800795 (IL6), rs1800972 (DEFB1) — have no coverage in the included set; CLEC1A rs2306894 has one study.
5. Allele frequencies: 7 entries in `candidate_snp_seed_table.csv` still PENDING against gnomAD v4 NFE.
6. ~~Haplotype frequencies for the PTX3 and TLR4 blocks have not been derived.~~ Derived 17 October 2026 for all four blocks directly from phased 1000 Genomes European genotypes (503 samples, 1,006 chromosomes) and rescaled to gnomAD NFE margins by iterative proportional fitting: `data/haplotype_frequencies_phased.csv`, `data/haplotype_ld_phased.csv`, `scripts/04_haplotype_frequencies.R`. The phase is coupling in every pair, confirmed from the phased chromosomes. The PTX3 trio has three common haplotypes (G-A-G 54%, A-C-A 40%, A-A-A 6% in 1000G EUR), not two; see the addendum to `phase3_decisions.md`.
7. **TLR4 supplementary retrieval or author contact** — Fisher 2017's supplemental Table 1 and three other supplementary files hold the estimates and counts that would make the TLR4 block poolable.

## 8. Suggested order of work

1. Code the baseline equation and its verification harness — all inputs confirmed, no dependencies.
2. Fit the correlation structure and report reproduction against all 16 published targets.
3. Build the decision-curve machinery on the baseline alone, so net benefit at 5% and 20% is working before any genetic term exists.
4. In parallel, run round-2 verification starting with the 13 provenance records, then the 310 P1 rows of 1380 flagged in total, then the six panel variants with four or more independent cohorts (rs16910526, rs2305619, rs7309123, rs3816527, rs4986790, rs4986791), verified line by line against source tables rather than sampled.
5. Derive the panel, then the frequencies and haplotypes, then the genetic layer.
