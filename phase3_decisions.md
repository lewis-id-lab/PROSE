# Phase 3 gate: decisions

Recorded 16 October 2026. This is the pre-specification record — it is written **before** any genetic
result exists and should not be revised once results are in view.

---

## D1. Primary genotype — **settled: recipient primary, donor secondary in the allo-HSCT stratum**

Agreed 16 October 2026. Recipient genotype is the primary exposure. Donor genotype enters as a
pre-specified secondary analysis restricted to the allogeneic HSCT stratum, which carries 33% of all IMD
events on 7% of admissions, and which gives Aim 3b a second host-context contrast alongside the
neutropenia split. The reasoning that forced the revision is kept below for the record.

**Originally stated:** donor genotype as the primary exposure.

**The biology supports you, the cohort does not.** After allogeneic transplant the reconstituted innate
immune system is donor-derived, so for PTX3, dectin-1 and the TLR genes the donor is arguably the right
genome — that is the Cunha 2014 *NEJM* argument. But two facts block donor as the *primary* exposure:

1. **The baseline cohort is overwhelmingly not transplanted.** In Stanzani's Table 1, allogeneic HSCT is the reason for admission in 248 of 3,994 non-IMD admissions and 44 of 133 IMD admissions — **292 of 4,127 admissions, 7.1%**. Donor genotype is undefined for the other 93%. A donor-genotype term cannot enter a model whose denominator is all haematology admissions.
2. **The evidence base cannot parameterise it.** Of 970 usable variant rows, only 90 are donor-stratified. Per panel variant, the number of studies contributing a donor-specific usable estimate is **at most one**, and rs4804800 and rs1840680 have none. No panel variant reaches two donor studies, let alone three.

Donor-as-primary would therefore shrink the analysis to a 7% subgroup — 292 admissions, 44 events —
parameterised from single studies.

**Proposed instead:** recipient genotype as the primary exposure, with **donor genotype as a
pre-specified secondary analysis restricted to the allogeneic HSCT stratum**. That stratum carries 33% of
all IMD events despite being 7% of admissions, so it is where a genetic effect has the most room to show,
and it is where the donor-genome argument actually applies. This also gives Aim 3b a second natural
contrast alongside the neutropenia split.

☑ Revision confirmed.

## D2. PTX3 and TLR4 as haplotypes — **accepted**

Consistent with the evidence and with the thesis's central methodological correction.

- **PTX3** enters as a single h1/h2 haplotype term, not three SNPs. In 1000 Genomes EUR, rs2305619 and rs1840680 sit at r² = 0.996 and are interchangeable; rs3816527 is at r² ≈ 0.79 with both; D′ ≈ 1 throughout. The trio collapses to two haplotypes at roughly 0.40–0.46 (h2) and 0.54 (h1). Simulating the three independently is what inflated the predecessor's result.
- **TLR4** enters as a single rs4986790–rs4986791 haplotype term. D′ = 0.981, r² = 0.945, GT haplotype frequency 0.056 in EUR, with about 0.3% recombinant chromosomes.

Parameter cost: these two blocks contribute 2 terms rather than 5.

☐ Outstanding: confirm the haplotype phase sign against a published haplotype table before simulating. The derivation in `haplotype_frequencies.csv` assumes the minor alleles are in coupling; D′ near 1 fixes the magnitude of D, not its direction.

## D3. Final panel — **settled, with the within-gene structure corrected**

The panel chosen was: PTX3 haplotype, TLR4 haplotype, CLEC7A rs16910526, CLEC7A rs7309123, CD209
rs4804800 and/or rs7248637, with IL10 rs1800896 and IFNG rs2069705 as possibles.

**Two of those terms had to be collapsed.** The panel as listed carries two CLEC7A SNPs and two CD209 SNPs
as separate additive terms, which is structurally the same error the thesis exists to correct — the
predecessor's mistake was treating linked SNPs within one gene as independent. Pairwise linkage
disequilibrium in 1000 Genomes EUR:

| pair | gene | r² | D′ | haplotypes present |
|---|---|---|---|---|
| rs16910526 × rs7309123 | CLEC7A | 0.104 | 1.000 | 3 |
| rs4804800 × rs7248637 | CD209 | 0.413 | 0.671 | 4 |

CLEC7A: D′ = 1 means no recombination between the two sites, so the rare Y238X allele (rs16910526, MAF
0.071) sits **entirely on one rs7309123 background**. The low r² reflects the frequency difference, not
independence. Three haplotypes exist — CG 0.071, majG 0.352, majmaj 0.578 — and the Cmaj haplotype is
absent. Entering both SNPs as independent additive terms double-counts the 7% of chromosomes carrying
Y238X.

CD209: D′ = 0.671 means recombination has occurred and all four haplotypes are present — GA 0.082,
Gmaj 0.042, majA 0.033, majmaj 0.843. The two SNPs are correlated at r² = 0.41, so again not independent.

**Final panel: four core haplotype terms, with two optional single SNPs.**

| term | variants | contributing studies | independent cohorts |
|---|---|---|---|
| PTX3 haplotype (h1/h2) | rs2305619, rs3816527, rs1840680 | 7 | 7 |
| CLEC7A haplotype (3 haplotypes) | rs16910526, rs7309123 | 11 | 9 |
| TLR4 haplotype | rs4986790, rs4986791 | 5 | 5 |
| CD209 haplotype (4 haplotypes) | rs4804800, rs7248637 | 5 | 4 |
| IL10 rs1800896 *(optional)* | rs1800896 | 3 | 3 |
| IFNG rs2069705 *(optional)* | rs2069705 | 3 | 3 |

That is **four core terms, six with the optional pair** — comfortably inside the roughly ten-parameter
ceiling that 133 events supports, and it leaves room for the Aim 4 sweep to add synthetic loci without
the baseline panel already being overfitted.

CLEC7A is now the best-evidenced term in the panel at nine independent cohorts, ahead of PTX3 at seven.

Three sub-decisions remain:

☐ **D3a.** Drop TLR9 rs352140 — reported by three studies, usable in one. Recommended.
☐ **D3b.** Herrero-Sanchez 2018 in or out. Its supplementary rows use *Candida*/*Aspergillus* detection or colonisation as the endpoint, not proven or probable invasive mould disease. It contributes a cohort to the PTX3 and CLEC7A terms; excluding it takes PTX3 from 7 to 6 and CLEC7A from 9 to 8.
☐ **D3c.** TLR4 pooling model. Both SNPs are usable only through genotype counts, and every contributing cohort has zero minor homozygotes in cases and controls, so only a dominant or carrier model is estimable and a continuity correction is required.

## D4. Minimum clinically important difference — **settled: 20% more true cases, no extra false positives**

**Stated:** a 20% improvement in prediction.

Net benefit is measured in true positives per patient, not in percent, so "20%" needs saying *of what*
before it can be tested. Three defensible readings, with what each requires numerically. The figures come
from the synthetic cohort at an exchangeable latent correlation of 0.30 (N = 400,000, prevalence 3.25%,
aROC 0.852) and should be re-derived once the calibrated cohort of section A2 exists.

Baseline clinical model, net benefit:

| threshold | NB (model) | NB (treat all) | % flagged | true positives / 1000 | false positives / 1000 |
|---|---|---|---|---|---|
| 5% | 0.01514 | −0.01843 | 13.3% | 21.0 | 111.7 |
| 10% | 0.01007 | −0.07501 | 6.0% | 15.1 | 44.9 |
| 20% | 0.00571 | −0.20939 | 2.5% | 9.5 | 15.4 |

What the genetic layer would have to add:

| reading | at the 5% threshold | at the 20% threshold |
|---|---|---|
| **(a)** 20% relative increase in net benefit | +0.00303 | +0.00114 |
| **(b)** 20% more true cases detected at no increase in false positives | +0.00420 | +0.00191 |
| **(c)** 20% relative increase in standardised net benefit | +0.00303 | +0.00114 |

For scale: a net-benefit increment of 0.00100 is one extra true case found per 1,000 admissions, which
across the 4,127-admission Bologna registry is about four extra cases over the decade.

**Recommendation: reading (b).** It is the one a clinician can state in a sentence — "the genetic layer
must find 20% more cases of invasive mould disease without treating anyone extra unnecessarily" — and it
is the strictest of the three. Reading (a) is simplest to compute and is the conventional decision-curve
formulation. Avoid framing the MCID in terms of reclassification: a "20% improvement" expressed as NRI is
precisely the predecessor thesis's 22.5% figure, and moving off that metric is one of the reasons for
redoing the work.

☑ **Settled: reading (b).** The genetic layer must find 20% more cases of invasive mould disease without
treating anyone extra unnecessarily. In net-benefit terms that is **+0.00420 at the 5% threshold** and
**+0.00191 at the 20% threshold**, on the indicative cohort. Re-derive both on the calibrated cohort
before locking; the ratios will hold, the absolute values may move.

☐ Outstanding: confirm whether the MCID must be met at both thresholds or either.

---

## Record of what is now fixed

| decision | status |
|---|---|
| Primary genotype | **fixed** — recipient primary; donor secondary, restricted to the allo-HSCT stratum |
| PTX3 as h1/h2 haplotype | **fixed** |
| TLR4 as rs4986790–rs4986791 haplotype | **fixed** |
| CLEC7A as a 3-haplotype term | **fixed** — collapsed from two SNPs after LD check |
| CD209 as a 4-haplotype term | **fixed** — collapsed from two SNPs after LD check |
| IL10 rs1800896, IFNG rs2069705 | optional, include or drop at the analyst's discretion |
| MCID | **fixed** — 20% more true cases at no increase in false positives; +0.00420 at 5%, +0.00191 at 20% |
| D3a drop TLR9 rs352140 | outstanding |
| D3b Herrero-Sanchez in or out | outstanding |
| D3c TLR4 dominant or carrier model | outstanding |
| MCID at both thresholds or either | outstanding |
