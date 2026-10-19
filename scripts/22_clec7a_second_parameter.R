# The second CLEC7A parameter (addendum 16). Run from the project root, after scripts 14 and 16.
#
# CLEC7A has three haplotypes (rs16910526 and rs7309123, gnomAD NFE): A-C 57.4% (the reference), A-G 34.9% and C-G 7.7% (Y238X). Y238X's effect (beta1) comes from
# the rs16910526 rows (script 14). The second parameter, beta2, is the per-copy effect of A-G (the G allele of rs7309123 on the non-Y238X background), and is
# identified by the rs7309123 rows once the strand of each study is known. rs7309123 is a C/G palindrome; the strand of each study was determined from its control
# frequencies and the papers (data/rs7309123_strand_resolution.csv). Studies with no usable evidence (grade D) are left out of the primary analysis; after the full texts of Kalkanci 2020 and Chen 2019 were read (addendum 19) only Zhao 2022 is grade D.
#
# East Asian studies (Thai, Chinese Han) are mapped through East Asian haplotype frequencies from phased 1000 Genomes (script 28, data/haplotype_frequencies_eas.csv): the Thai study with
# CDX (Dai in Xishuangbanna) and the Chinese Han studies with CHB + CHS. In those populations G is the common allele of rs7309123 (0.77 to 0.80), A-G is the common haplotype
# and A-C the reference minority (0.20 to 0.23), and Y238X is almost absent, so the same per-copy effect of A-G is estimated against a very different background.
#
# Each rs7309123 row is mapped to beta2 by solving the three-haplotype model with beta1 held at its pooled median (the C-G chromosomes carry 18% of the G alleles, so the
# effect of this plug-in is small). Non-European studies are left out of the primary analysis: the European haplotype frequencies do not describe them (G is 0.43 in
# Europeans and 0.81 in East Asians), which makes rule 7 of addendum 7 untenable for this SNP.
#
# Writes   output/clec7a_second_study_effects.csv, clec7a_second_pooled.csv, clec7a_second_posterior_draws.csv, clec7a_second_layer_summary.csv

source("R/baseline.R"); source("R/bs_data.R"); source("R/decision_curve.R"); source("R/boscore_data.R")
source("R/refit_baseline.R"); source("R/evidence_mapping.R"); source("R/ptx3_evidence.R"); source("R/source_verification.R")
source("R/genetic_layer.R"); source("R/genetic_helpers.R")

vr <- read_bundle_table("extraction_variants_round1.csv"); st <- read_bundle_table("extraction_studies_round1.csv")
vr <- apply_corrections(vr, read.csv("data/extraction_corrections.csv", stringsAsFactors = FALSE))$data
strand <- read.csv("data/rs7309123_strand_resolution.csv", stringsAsFactors = FALSE)
beta1 <- median(read.csv("output/clec7a_posterior_draws.csv")$mu)

hap <- block_haplotypes("CLEC7A", c("C", NA), source = "gnomAD_NFE")
hap2 <- hap
hap2$risk <- as.integer(hap$a1 == "A" & hap$a2 == "G")                      # the class whose effect is estimated: A-G
hap2$offset <- beta1 * as.integer(hap$a1 == "C" & hap$a2 == "G")            # Y238X at its pooled effect
message("haplotypes: ", paste(sprintf("%s %.3f", hap$haplotype, hap$freq), collapse = ", "), " | beta1 (held) = ", sprintf("%.3f", beta1))

rows_for <- function(v) ptx3_row_effects(v, st, hap2, loci = c("rs16910526", "rs7309123"), counts_model = "carrier")
eas_tab <- read.csv("data/haplotype_frequencies_eas.csv", stringsAsFactors = FALSE)
eas_hap <- function(pop) {
  d <- eas_tab[eas_tab$population == pop & eas_tab$freq > 0, ]
  h <- data.frame(haplotype = d$haplotype, freq = d$freq / sum(d$freq), a1 = d$allele_1, a2 = d$allele_2, stringsAsFactors = FALSE)
  h$risk <- as.integer(h$a1 == "A" & h$a2 == "G"); h$offset <- beta1 * as.integer(h$a1 == "C" & h$a2 == "G")
  h
}
eas_population <- c(S077 = "Thai proxy (CDX)", S074 = "Chinese Han (CHB+CHS)", S081 = "Chinese Han (CHB+CHS)")
# Rows of European studies with the European haplotypes and rows of the East Asian studies with the East Asian ones
rows_all <- function(v) {
  eu <- rows_for(v[!(v$study_id %in% names(eas_population)), ])
  ea <- do.call(rbind, lapply(names(eas_population), function(sid) { x <- v[v$study_id == sid, ]
    usable <- x$rsid %in% c("rs16910526", "rs7309123") & !is.na(x$meta_analysable) & x$meta_analysable == "yes" & x$usable_tier %in% c("estimate", "counts_2x3", "counts_2x2")
    if (!any(usable)) return(NULL); ptx3_row_effects(x, st, eas_hap(eas_population[[sid]]), loci = c("rs16910526", "rs7309123"), counts_model = "carrier") }))
  rbind(eu, ea)
}
re <- rows_for(vr)
re_all <- rows_all(vr)
vr_flip <- vr
unresolved <- strand$study_id[strand$grade == "D"]
vr_flip$orientation_flag[vr_flip$rsid %in% "rs7309123" & vr_flip$study_id %in% unresolved] <- "strand_flip"
re_flip <- rows_for(vr_flip)
re_flip_all <- rows_all(vr_flip)

all_studies <- unique(re$study_id)
excl <- function(grades) setdiff(all_studies, strand$study_id[strand$grade %in% grades])
se <- function(r, ...) ptx3_study_effects(r, snp_order = "rs7309123", ...)
scenarios <- list(
  "primary: strand resolved (grades A and B), European, odds ratios" = se(re, exclude_studies = excl(c("A", "B")), allow_non_european = FALSE),
  "adding the non-European strand-resolved studies (S077 Thai, S074 Chinese Han), with East Asian haplotype frequencies" = se(re_all, exclude_studies = excl(c("A", "B")), allow_non_european = TRUE),
  "adding hazard ratios (S080)" = se(re, exclude_studies = excl(c("A", "B")), allow_non_european = FALSE, allow_hr = TRUE),
  "adding the Jena cohort (grade C)" = se(re, exclude_studies = excl(c("A", "B", "C")), allow_non_european = FALSE),
  "all studies, unresolved ones taken as forward" = se(re_all, exclude_studies = character(), allow_non_european = TRUE),
  "unresolved studies taken as the opposite strand" = se(re_flip_all, exclude_studies = character(), allow_non_european = TRUE),
  "East Asian strand-resolved studies only (S077, S074), with East Asian haplotype frequencies" = se(re_all, exclude_studies = setdiff(all_studies, c("S077", "S074")), allow_non_european = TRUE),
  "non-European studies with the European haplotype frequencies (not valid, for comparison)" = se(re, exclude_studies = excl(c("A", "B")), allow_non_european = TRUE))
pooled <- Map(pool_study_effects, scenarios, seq_along(scenarios) + 2209L)
keep <- !vapply(pooled, is.null, logical(1))
tab <- do.call(rbind, Map(function(p, nm) cbind(scenario = nm, p$summary), pooled[keep], names(scenarios)[keep]))
primary <- scenarios[[1]]
study_tab <- primary[, c("study_id", "first_author", "year", "rsid", "cohort", "tier", "model", "measure", "method", "y", "se", "beta", "se_beta", "flag")]
study_tab$or_per_copy <- exp(study_tab$beta)
unres <- do.call(rbind, lapply(names(scenarios)[5:6], function(n) { x <- scenarios[[n]]; if (nrow(x)) cbind(scenario = n, x[, c("study_id", "ancestry_as_reported", "model", "measure", "y", "se", "beta", "se_beta")]) }))
message("\nPrimary study effects (per copy of A-G, beta1 held):"); show_df(study_tab[, c("study_id", "tier", "model", "measure", "y", "se", "beta", "se_beta", "or_per_copy")])
message("\nPooled:"); show_df(cbind(tab[, c("scenario", "k", "studies")], tab[, c("mu_lo", "mu", "mu_hi", "or", "or_lo", "or_hi", "tau", "prob_below_zero")]))

dir.create("output", showWarnings = FALSE)
write.csv(study_tab, "output/clec7a_second_study_effects.csv", row.names = FALSE)
write.csv(tab, "output/clec7a_second_pooled.csv", row.names = FALSE)
write.csv(pooled[[1]]$draws, "output/clec7a_second_posterior_draws.csv", row.names = FALSE)
write.csv(unres, "output/clec7a_second_unresolved_scenarios.csv", row.names = FALSE)
write.csv(pooled[[2]]$draws, "output/clec7a_second_posterior_draws_all_ancestries.csv", row.names = FALSE)
ancestry_effects <- do.call(rbind, lapply(c(2L, 7L), function(i) { x <- scenarios[[i]]; if (nrow(x)) cbind(scenario = names(scenarios)[i], x[, c("study_id", "ancestry_as_reported", "tier", "model", "measure", "y", "se", "beta", "se_beta", "flag")]) }))
write.csv(ancestry_effects, "output/clec7a_second_ancestry_study_effects.csv", row.names = FALSE)
# Do the two ancestries agree? Fixed-effect combination within each ancestry, then the difference
fx <- function(x) combine_fixed(x$beta, x$se_beta)
a_eu <- fx(scenarios[[1]]); a_ea <- fx(scenarios[[7]])
d_ea <- a_ea$y - a_eu$y; se_d <- sqrt(a_eu$se^2 + a_ea$se^2)
heterogeneity <- data.frame(group = c("European (S045, S056, S071)", "East Asian (S074, S077)"), studies = c(a_eu$k, a_ea$k), beta_per_copy = c(a_eu$y, a_ea$y), se = c(a_eu$se, a_ea$se),
                            q_within_group = c(a_eu$q, a_ea$q), difference_east_asian_minus_european = d_ea, se_difference = se_d, z = d_ea / se_d, p_two_sided = 2 * pnorm(-abs(d_ea / se_d)))
write.csv(heterogeneity, "output/clec7a_second_ancestry_heterogeneity.csv", row.names = FALSE)
message("\nAncestry heterogeneity of beta2:"); show_df(heterogeneity)
message("\nStudy effects with ancestry-specific haplotypes (all strand-resolved studies, and the East Asian ones):"); show_df(ancestry_effects[, c("scenario", "study_id", "model", "y", "se", "beta", "se_beta")] |> transform(scenario = substr(scenario, 1, 30)))

# ---- The layer with two CLEC7A parameters --------------------------------------------------------------------------------------------
at <- c(0.05, 0.20)
base <- select_cohort(load_base_cohort(), include_twins = FALSE)
lp <- profile_logit(load_refit_params())
q <- tabulate(profile_index(base$X), 2L^length(PREDICTORS)) / nrow(base$X)
mcid_tab <- read.csv("output/dca_baseline_refit_summary.csv")
mcid <- setNames(mcid_tab$mcid_b[match(at, mcid_tab$threshold)], as.character(at))

p <- genetic_term_frequencies()
terms5 <- c("PTX3", "CLEC7A", "CLEC7A_AG", "TLR4", "CD209")
grid <- joint_grid(list(dosage_grid(p["PTX3"]), haplotype_gene_part(hap, c(CLEC7A = "C-G", CLEC7A_AG = "A-G")), dosage_grid(p["TLR4"]), dosage_grid(p["CD209"])))
grid$G <- grid$G[, terms5]
stopifnot(abs(sum(grid$prob) - 1) < 1e-12)
cells <- genetic_cells(q, lp, grid)

effects_file <- read.csv("output/genetic_layer_effects.csv")
n_draw <- 2000L
b2 <- pooled[[1]]$draws
set.seed(6724)
b2_draw <- list(mean = b2$mu[seq_len(n_draw)], predictive = rnorm(n_draw, b2$mu[seq_len(n_draw)], b2$tau[seq_len(n_draw)]))
mk <- function(h, with_b2) { m <- effects_file[effects_file$truth == h & effects_file$draw <= n_draw, c("PTX3", "CLEC7A", "TLR4", "CD209")]
  cbind(PTX3 = m$PTX3, CLEC7A = m$CLEC7A, CLEC7A_AG = if (with_b2) b2_draw[[h]] else 0, TLR4 = m$TLR4, CD209 = m$CD209) }
beta_hat <- c(pooled_estimates(load_term_draws())[c("PTX3", "CLEC7A")], CLEC7A_AG = median(b2$mu), pooled_estimates(load_term_draws())[c("TLR4", "CD209")])[terms5]
beta_hat0 <- beta_hat; beta_hat0[["CLEC7A_AG"]] <- 0
message("\npooled estimates carried by the realised score: ", paste(names(beta_hat), sprintf("%.3f", beta_hat), collapse = ", "))

cl <- start_workers()
run <- function(with_b2) { eff <- list(mean = mk("mean", with_b2), predictive = mk("predictive", with_b2)); summarise_layer(exact_layer(cl, cells, eff, if (with_b2) beta_hat else beta_hat0, at, mcid)) }
lay <- rbind(cbind(second_parameter = "estimated", run(TRUE)), cbind(second_parameter = "held at zero", run(FALSE)))
# Sensitivity: beta2 estimated from all strand-resolved studies, each with its own ancestry's haplotype frequencies
run_alt <- function(b2a) {
  set.seed(6725)
  ba <- list(mean = b2a$mu[seq_len(n_draw)], predictive = rnorm(n_draw, b2a$mu[seq_len(n_draw)], b2a$tau[seq_len(n_draw)]))
  mk_a <- function(h) { m <- effects_file[effects_file$truth == h & effects_file$draw <= n_draw, c("PTX3", "CLEC7A", "TLR4", "CD209")]
    cbind(PTX3 = m$PTX3, CLEC7A = m$CLEC7A, CLEC7A_AG = ba[[h]], TLR4 = m$TLR4, CD209 = m$CD209) }
  bh <- beta_hat; bh[["CLEC7A_AG"]] <- median(b2a$mu)
  summarise_layer(exact_layer(cl, cells, list(mean = mk_a("mean"), predictive = mk_a("predictive")), bh, at, mcid))
}
lay <- rbind(lay, cbind(second_parameter = "estimated from all ancestries", run_alt(pooled[[2]]$draws)))
parallel::stopCluster(cl)
write.csv(lay, "output/clec7a_second_layer_summary.csv", row.names = FALSE)
message("\nGenetic layer with and without the second CLEC7A parameter (gain over the baseline):")
show_df(lay[, c("second_parameter", "truth", "scenario", "threshold", "gain", "lo", "hi", "extra_cases_per_1000", "p_meets_mcid")])
