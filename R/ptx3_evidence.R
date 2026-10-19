# Evidence-to-parameter mapping for one haplotype term (addenda 7 and 8): extraction rows -> per-copy effect of the
# term's risk haplotype class (log OR) -> one effect per study -> pooled effect. Written for PTX3 (h2), and reused
# unchanged for CLEC7A (Y238X) and TLR4 (the G-T haplotype) by passing that block's haplotypes and loci.
# Requires R/evidence_mapping.R.
#
# The extraction rows are round 1. Where a row has been checked against its source and corrected, the correction is
# applied first (R/source_verification.R); everything else is unverified, so everything derived here is provisional.

# Studies whose ancestry is outside the European structure used for the haplotype frequencies. Only the leading
# descriptor counts (the text before the first bracket, semicolon, comma or colon), because a European cohort's note
# such as "Caucasian (predominantly; 4% other, including Asian and African)" must not flag it. Whole words only, so
# "Asian" does not match inside "Caucasian".
NON_EUROPEAN_PATTERN <- "\\b(Chinese|Thai|Japanese|Korean|Asian|African|Indian)\\b"
non_european_study <- function(ancestry) {
  lead <- sub("[(;,:].*$", "", ancestry)
  !is.na(lead) & grepl(NON_EUROPEAN_PATTERN, lead, ignore.case = TRUE, perl = TRUE)
}

# Subgroups that are whole-cohort analyses (not outcome, exposure or site strata). "donor genotype; whole cohort" counts, as do
# Ceesay 2016's "patients tested for SNP" and "HSCT donors" and Fisher 2017's "recipient (patient) genotype" and "recipient (PT)
# genotype" (whole cohort, at different time horizons), and Lupianez 2015's "patient (recipient) genotype" for its phase 1 and
# overall combined cohorts (addendum 9).
WHOLE_COHORT_PATTERN <- "whole cohort|^all episodes|stratification cohort|validation cohort|patients tested|HSCT donors|recipient \\(patient\\) genotype|patient \\(recipient\\) genotype|recipient \\(PT\\) genotype"

# Minor allele at a locus in the European haplotype frequencies (used only where a study says "the minor allele")
minor_allele <- function(hap, locus) {
  al <- hap[[paste0("a", locus)]]
  f <- tapply(hap$freq, al, sum)
  names(f)[which.min(f)]
}

# Whose genotype was analysed in a row? The row text decides; otherwise the study's single source; otherwise unstated.
row_genotype_source <- function(subgroup, study_source) {
  ifelse(grepl("donor", subgroup, ignore.case = TRUE), "donor",
  ifelse(grepl("recipient|patient", subgroup, ignore.case = TRUE), "recipient",
  ifelse(study_source %in% c("recipient", "donor"), study_source, "unstated")))
}

# One extraction row -> a one-row data frame with the per-copy effect of the risk class where it can be mapped.
#   counts_model: "trend" maps a 2x3 table through its additive trend (PTX3); "carrier" compares carriers of the effect
#   allele with non-carriers, with 0.5 added to a zero cell (TLR4, decision D3c).
map_extraction_row <- function(row, hap, loci = PTX3_LOCI, counts_model = c("trend", "carrier")) {
  counts_model <- match.arg(counts_model)
  locus <- match(row$rsid, loci)
  base <- data.frame(variant_row_id = row$variant_row_id, study_id = row$study_id, rsid = row$rsid, locus = locus,
                     tier = row$usable_tier, model = row$genetic_model, measure = row$effect_measure,
                     subgroup = row$subgroup, primary_flag = row$is_primary_analysis,
                     y = NA_real_, se = NA_real_, beta = NA_real_, se_beta = NA_real_, method = NA_character_,
                     flag = "", stringsAsFactors = FALSE)
  if (is.na(locus)) { base$flag <- "no resolvable SNP"; return(base) }
  fail <- function(msg) { base$flag <- msg; base }
  num <- function(x) suppressWarnings(as.numeric(x))

  # A study reported on the opposite strand: complement the alleles, and later the genotypes written in its text
  flip <- identical(as.character(row$orientation_flag), "strand_flip")
  if (flip) {
    if (!is.na(row$effect_allele)) row$effect_allele <- complement_base(row$effect_allele)
    if (!is.na(row$other_allele)) row$other_allele <- complement_base(row$other_allele)
  }
  strand_note <- if (flip) "reported on the opposite strand, complemented" else ""

  # The reported alleles must be alleles at this locus. If the effect allele is but the other allele is not (a reporting
  # error, such as A/C at the A/G locus rs2305619), the other allele is taken to be the locus's remaining allele.
  loc_alleles <- unique(hap[[paste0("a", locus)]])
  note <- strand_note
  e0 <- row$effect_allele; o0 <- row$other_allele
  if (!is.na(e0) && !(e0 %in% loc_alleles)) return(fail(paste0("effect allele ", e0, " is not an allele at this locus (strand?)")))
  if (!is.na(e0) && (is.na(o0) || !(o0 %in% loc_alleles))) {
    rest <- setdiff(loc_alleles, e0)
    if (length(rest) == 1) {
      row$other_allele <- rest
      note <- paste(c(note, paste0("other allele reported as ", ifelse(is.na(o0), "missing", o0), ", replaced by ", rest))[c(nzchar(note), TRUE)], collapse = "; ")
    }
  }
  add_note <- function(flag) paste(c(note, flag)[c(nzchar(note), nzchar(flag))], collapse = "; ")

  if (row$usable_tier == "estimate") {
    lo <- log_or_from_estimate(row$effect_estimate, row$ci_lower, row$ci_upper)
    if (is.na(lo$y)) return(fail(add_note(lo$flag)))
    contrast <- parse_contrast(row$comparison_as_reported, row$genetic_model, row$effect_allele, row$other_allele, flip_text = flip)
    mapped <- map_to_h2(lo$y, lo$se, contrast, hap, locus)
    base[, c("y", "se")] <- c(lo$y, lo$se)
    base[, c("beta", "se_beta")] <- mapped[, c("beta", "se_beta")]
    base$method <- paste0("estimate, ", contrast$method)
    base$flag <- add_note(mapped$flag)
    return(base)
  }
  if (row$usable_tier == "counts_2x3") {
    if (is.na(row$effect_allele) || is.na(row$other_allele)) return(fail(add_note("effect allele not stated")))
    cs <- num(c(row$n_cases_hom_effect, row$n_cases_het, row$n_cases_hom_other))
    ct <- num(c(row$n_ctrl_hom_effect, row$n_ctrl_het, row$n_ctrl_hom_other))
    if (anyNA(c(cs, ct))) return(fail(add_note("incomplete 2x3 table")))
    e <- row$effect_allele; o <- row$other_allele
    if (counts_model == "carrier") {
      lo <- log_or_from_counts_2x2(cs[1] + cs[2], cs[3], ct[1] + ct[2], ct[3])
      contrast <- list(type = "genotype", numerator = c(paste0(e, e), paste(sort(c(e, o)), collapse = "")), reference = paste0(o, o), flag = "")
      method <- "counts 2x3, carriers against non-carriers"
      extra <- if (any(c(cs[1] + cs[2], cs[3], ct[1] + ct[2], ct[3]) == 0)) "0.5 added to a zero cell" else ""
    } else {
      lo <- log_or_from_counts_2x3(cs, ct)
      contrast <- list(type = "allele", allele = e, flag = "")
      method <- "counts 2x3, additive trend"
      extra <- "additive trend mapped as an allelic effect (approximation)"
    }
    if (is.na(lo$y)) return(fail(add_note(lo$flag)))
    mapped <- map_to_h2(lo$y, lo$se, contrast, hap, locus)
    base[, c("y", "se")] <- c(lo$y, lo$se)
    base[, c("beta", "se_beta")] <- mapped[, c("beta", "se_beta")]
    base$method <- method
    base$flag <- add_note(paste(c(extra, mapped$flag)[c(nzchar(extra), nzchar(mapped$flag))], collapse = "; "))
    return(base)
  }
  if (row$usable_tier == "counts_2x2") {
    # Layout (Brunel 2018): events and non-events among homozygotes for the effect allele, against all others
    homozygote_comparison <- identical(as.character(row$genetic_model), "recessive") || grepl("homozyg", row$comparison_as_reported, ignore.case = TRUE)
    if (!homozygote_comparison) return(fail(add_note("2x2 table is not a homozygote comparison")))
    e <- row$effect_allele
    inferred <- FALSE
    if (is.na(e)) { e <- minor_allele(hap, locus); inferred <- TRUE }
    o <- setdiff(unique(hap[[paste0("a", locus)]]), e)
    if (length(o) != 1) return(fail(add_note("alleles cannot be resolved")))
    lo <- log_or_from_counts_2x2(num(row$n_cases_hom_effect), num(row$n_cases_het), num(row$n_ctrl_hom_effect), num(row$n_ctrl_het))
    if (is.na(lo$y)) return(fail(add_note(lo$flag)))
    ee <- paste0(e, e); het <- paste(sort(c(e, o)), collapse = ""); oo <- paste0(o, o)
    mapped <- map_to_h2(lo$y, lo$se, list(type = "genotype", numerator = ee, reference = c(het, oo), flag = ""), hap, locus)
    base[, c("y", "se")] <- c(lo$y, lo$se)
    base[, c("beta", "se_beta")] <- mapped[, c("beta", "se_beta")]
    base$method <- "counts 2x2, homozygotes against the rest"
    base$flag <- add_note(paste(c(if (inferred) paste0("effect allele inferred as the European minor allele (", e, ")"), mapped$flag)[c(inferred, nzchar(mapped$flag))], collapse = "; "))
    return(base)
  }
  fail("tier not mapped")
}

# Priority of a row within a study: counts-based first, then allelic, dominant or carrier, recessive, other
model_priority <- function(tier, model) {
  ifelse(tier == "counts_2x3", 1L,
  ifelse(model %in% c("additive_per_allele", "allelic"), 2L,
  ifelse(model %in% c("dominant", "carrier"), 3L,
  ifelse(model == "recessive", 4L, 5L))))
}

# Map every usable row for a block's SNPs and attach the study-level attributes needed for the selection rules
ptx3_row_effects <- function(vr, st, hap, loci = PTX3_LOCI, counts_model = "trend", corrected_ids = character()) {
  d <- vr[vr$rsid %in% loci & !is.na(vr$meta_analysable) & vr$meta_analysable == "yes" &
            vr$usable_tier %in% c("estimate", "counts_2x3", "counts_2x2"), , drop = FALSE]
  out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) map_extraction_row(d[i, ], hap, loci, counts_model)))
  si <- st[, c("study_id", "first_author", "year", "genotype_source", "ancestry_as_reported", "population_directness")]
  out <- merge(out, si, by = "study_id", all.x = TRUE, sort = FALSE)
  out$genotype <- row_genotype_source(out$subgroup, out$genotype_source)
  out$donor <- out$genotype == "donor"
  out$whole_cohort <- grepl(WHOLE_COHORT_PATTERN, out$subgroup, ignore.case = TRUE)
  out$cohort <- ifelse(grepl("stratification", out$subgroup, ignore.case = TRUE), "stratification",
                ifelse(grepl("validation", out$subgroup, ignore.case = TRUE), "validation", "whole"))
  out$combined <- grepl("overall combined", out$subgroup, ignore.case = TRUE)
  out$non_european <- non_european_study(out$ancestry_as_reported)
  out$odds_ratio_scale <- !(out$measure %in% c("HR", "RR", "subdistribution_HR"))
  out$priority <- model_priority(out$tier, out$model)
  out$corrected <- out$variant_row_id %in% corrected_ids
  out[order(out$study_id, out$rsid, out$variant_row_id), ]
}

# Herrero-Sanchez 2018 reports several endpoints. Its Table 3 endpoint, detection of Candida or Aspergillus, includes
# colonisation; Table 5 is invasive aspergillosis in the whole cohort. This swaps which of them counts as the whole-cohort analysis.
use_s035_ia_endpoint <- function(re) {
  s <- re$study_id == "S035"
  ia <- grepl("invasive aspergillosis \\(IA\\) vs no IA", re$subgroup) & !grepl("prophylaxis", re$subgroup, ignore.case = TRUE)
  re$whole_cohort[s] <- ia[s]
  re
}

# One effect per study. Rules, in order:
#   1. a row must have a mapped effect, be a whole-cohort analysis, and use the chosen genotype source (recipient by
#      default; donor genotype is the secondary analysis of decision D1; rows whose source is unstated are left out);
#      hazard ratios are excluded unless allow_hr;
#   2. studies in `exclude_studies` are left out (Herrero-Sanchez 2018 in the primary analysis, decision D3b);
#   3. the SNP is the first of snp_order that the study has an eligible row for (rs2305619, then rs1840680 for h2);
#   4. within that SNP and cohort, the row with the best priority (counts, allelic, dominant or carrier, recessive,
#      other), then the row flagged primary, then the first row;
#   5. a study that reports several independent cohorts (Schnetzke 2024) has them combined by fixed-effect inverse variance.
ptx3_study_effects <- function(re, snp_order = c("rs2305619", "rs1840680"), exclude_studies = "S035",
                               allow_hr = FALSE, genotype = c("recipient", "donor"), allow_non_european = TRUE) {
  genotype <- match.arg(genotype)
  ok <- !is.na(re$beta) & re$whole_cohort & re$rsid %in% snp_order &
    re$genotype == genotype & (allow_hr | re$odds_ratio_scale) &
    (allow_non_european | !re$non_european) & !(re$study_id %in% exclude_studies)
  e <- re[ok, , drop = FALSE]
  if (!nrow(e)) return(data.frame())
  rows <- lapply(split(e, e$study_id), function(s) {
    snp <- snp_order[snp_order %in% s$rsid][1]
    s <- s[s$rsid == snp, ]
    per_cohort <- lapply(split(s, s$cohort), function(cc) {
      # An overall combined analysis contains the phase analyses it combines, so it is preferred to them
      combined <- if (is.null(cc$combined)) rep(FALSE, nrow(cc)) else cc$combined
      cc <- cc[order(!combined, cc$priority, cc$primary_flag != "yes", cc$variant_row_id), ]
      cc[1, ]
    })
    pc <- do.call(rbind, per_cohort)
    if (nrow(pc) > 1) {
      comb <- combine_fixed(pc$beta, pc$se_beta)
      best <- pc[1, ]
      best$beta <- comb$y; best$se_beta <- comb$se
      best$method <- paste0(best$method, "; ", nrow(pc), " cohorts combined by fixed-effect inverse variance")
      best$flag <- paste0("cohorts disagree (Q = ", sprintf("%.1f", comb$q), " on ", nrow(pc) - 1, " df)")
      best$cohort <- paste(pc$cohort, collapse = "+")
      return(best)
    }
    pc
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Bayesian random-effects pooling of a set of study effects, as a one-row summary plus the draws
pool_study_effects <- function(effects, seed) {
  if (!nrow(effects)) return(NULL)
  b <- bayes_meta_grid(effects$beta, effects$se_beta, prior_sd_mu = 1, prior_scale_tau = 0.5, seed = seed)
  list(summary = data.frame(k = nrow(effects), studies = paste(effects$study_id, collapse = " "),
                            mu_lo = b$mu[1], mu = b$mu[2], mu_hi = b$mu[3],
                            or = exp(b$mu[2]), or_lo = exp(b$mu[1]), or_hi = exp(b$mu[3]),
                            tau_lo = b$tau[1], tau = b$tau[2], tau_hi = b$tau[3],
                            pred_lo = b$prediction[1], pred_hi = b$prediction[3], prob_below_zero = b$prob_mu_below_zero),
       draws = b$draws)
}
