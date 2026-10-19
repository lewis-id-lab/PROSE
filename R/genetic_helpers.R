# Helpers shared by the genetic-layer scripts 18 to 20: worker processes, chunking, grouped summaries.
# Requires R/genetic_layer.R.

# Fixed number of chunks, so results do not depend on the number of cores
chunk_ids <- function(n, k = 30L) split(seq_len(n), rep_len(seq_len(k), n))

# Fresh worker processes (forking is unavailable in Positron) that source the model code
start_workers <- function(cores = max(1L, min(15L, parallel::detectCores() - 1L)), root = getwd(),
                          files = c("R/baseline.R", "R/decision_curve.R", "R/genetic_layer.R")) {
  cl <- parallel::makeCluster(cores)
  parallel::clusterCall(cl, function(root, files) { setwd(root); for (f in files) source(f); invisible(NULL) }, root, files)
  cl
}

parallel_rows <- function(cl, chunks, FUN, ctx) do.call(rbind, parallel::parLapply(cl, chunks, FUN, ctx))

# Apply f to the subset of d for each combination of the key columns and stack the results
by_key <- function(d, keys, f) {
  k <- unique(d[keys])
  do.call(rbind, lapply(seq_len(nrow(k)), function(i) {
    s <- Reduce(`&`, lapply(keys, function(n) d[[n]] == k[[n]][i]))
    cbind(k[i, , drop = FALSE], f(d[s, , drop = FALSE]), row.names = NULL)
  }))
}

# Summary of observed gains over simulated studies
summarise_sampling <- function(d) data.frame(
  mean_gain = mean(d$delta_nb), median_gain = stats::median(d$delta_nb), sd_gain = stats::sd(d$delta_nb),
  lo = stats::quantile(d$delta_nb, 0.05, names = FALSE), hi = stats::quantile(d$delta_nb, 0.95, names = FALSE),
  p_gain_positive = mean(d$delta_nb > 0), p_reaches_mcid = mean(d$delta_nb >= d$mcid), p_reaches_sample_mcid = mean(d$delta_nb >= d$mcid_sample),
  mcid = stats::median(d$mcid), median_events = stats::median(d$events))

show_df <- function(x) print(cbind(x[vapply(x, is.character, logical(1))], round(x[vapply(x, is.numeric, logical(1))], 5)), row.names = FALSE)

# ---- Layer runners used by the sensitivity-cohort script -----------------------------------------------------------------------------

# Exact layer for named effect matrices (truth -> draws x terms); one row per truth, draw, score and threshold
exact_layer <- function(cl, cells, effects, beta_hat, at, mcid) {
  do.call(rbind, lapply(names(effects), function(h) {
    ctx <- list(cells = cells, eff = effects[[h]], bh = beta_hat, at = at, mcid = mcid)
    cbind(truth = h, parallel_rows(cl, chunk_ids(nrow(effects[[h]])), function(ch, ctx) {
      r <- run_genetic_layer(ctx$cells, ctx$eff[ch, , drop = FALSE], ctx$bh, ctx$at, ctx$mcid)
      r$draw <- ch[r$draw]
      r
    }, ctx))
  }))
}

summarise_layer <- function(res) {
  key <- paste(res$truth, res$draw, res$threshold)
  base <- res[res$scenario == "baseline", ]
  res$delta_tp <- res$tp - base$tp[match(key, paste(base$truth, base$draw, base$threshold))]
  by_key(res[res$scenario != "baseline", ], c("truth", "scenario", "threshold"), function(d) data.frame(
    gain = stats::median(d$delta_nb), lo = stats::quantile(d$delta_nb, 0.05, names = FALSE), hi = stats::quantile(d$delta_nb, 0.95, names = FALSE),
    extra_cases_per_1000 = 1000 * stats::median(d$delta_tp), p_meets_mcid = mean(d$meets_mcid), mcid = d$mcid[1],
    auc = stats::median(d$auc), flagged = stats::median(d$flagged)))
}

# Shapley attribution for named effect matrices, first n_attr draws, ceiling and realised scores
attribution_layer <- function(cl, cells, effects, beta_hat, at, n_attr = 1000L) {
  terms <- colnames(cells$G)
  do.call(rbind, lapply(names(effects), function(h) {
    ctx <- list(cells = cells, eff = effects[[h]][seq_len(n_attr), , drop = FALSE], bh = beta_hat, at = at, terms = terms)
    cbind(truth = h, parallel_rows(cl, chunk_ids(n_attr), function(ch, ctx) {
      do.call(rbind, lapply(ch, function(i) do.call(rbind, lapply(c("ceiling", "realised"), function(sc) {
        s <- shapley_gain(ctx$cells, ctx$eff[i, ], if (sc == "ceiling") ctx$eff[i, ] else ctx$bh, ctx$at)
        cbind(draw = i, score = sc, threshold = ctx$at, total = s$total,
              stats::setNames(as.data.frame(s$phi), paste0("phi_", ctx$terms)),
              stats::setNames(as.data.frame(s$alone), paste0("alone_", ctx$terms)),
              stats::setNames(as.data.frame(s$leave_out), paste0("loo_", ctx$terms)))
      }))))
    }, ctx))
  }))
}

summarise_attribution <- function(a, terms) by_key(a, c("truth", "score", "threshold"), function(d) {
  out <- data.frame(total = stats::median(d$total))
  for (t in terms) {
    out[[paste0("phi_", t)]] <- stats::median(d[[paste0("phi_", t)]])
    out[[paste0("share_", t)]] <- sum(d[[paste0("phi_", t)]]) / sum(d$total)
    out[[paste0("alone_", t)]] <- stats::median(d[[paste0("alone_", t)]])
    out[[paste0("loo_", t)]] <- stats::median(d[[paste0("loo_", t)]])
  }
  out
})

# Sampling variability by resampling the cohort's patients; `runs` is a named list of list(eff = draws x terms, seed = int)
sampling_layer <- function(cl, cg, lp, p, cells, runs, beta_hat, at, mcid) {
  do.call(rbind, lapply(names(runs), function(r) {
    ctx <- list(cg = cg, lp = lp, p = p, cells = cells, eff = runs[[r]]$eff, bh = beta_hat, at = at, mcid = mcid, seed = runs[[r]]$seed)
    cbind(analysis = r, parallel_rows(cl, chunk_ids(nrow(runs[[r]]$eff)), function(ch, ctx) {
      set.seed(ctx$seed + ch[1])
      do.call(rbind, lapply(ch, function(i) {
        sh <- genetic_scores(ctx$cells, ctx$eff[i, ], ctx$bh)
        s <- simulate_study(ctx$cg, ctx$lp, ctx$p, sh, ctx$eff[i, ], ctx$bh, ctx$at, ctx$mcid)
        cbind(draw = i, s[s$scenario != "baseline", ])
      }))
    }, ctx))
  }))
}
