# Haplotype frequencies from phased genotypes, and rescaling to new allele-frequency margins.

# Table of haplotype counts across loci. `chrom` has one row per phased chromosome and one
# column of nucleotide alleles per locus; `alleles` is a named list giving the two alleles
# analysed at each locus (the order fixes the dimension order of the table).
haplotype_table <- function(chrom, alleles) {
  loci <- names(alleles)
  ok <- Reduce(`&`, lapply(loci, function(l) chrom[[l]] %in% alleles[[l]]))
  used <- chrom[ok, loci, drop = FALSE]
  fac <- lapply(loci, function(l) factor(used[[l]], levels = alleles[[l]]))
  names(fac) <- loci
  tab <- table(fac)
  attr(tab, "n_used") <- nrow(used)
  attr(tab, "n_dropped") <- sum(!ok)
  tab
}

# Iterative proportional fitting of a haplotype table to target one-way margins. This keeps
# every odds-ratio (interaction) structure of the starting table and moves only the allele
# frequencies, and it extends to any number of loci. Cells that start at zero stay at zero.
ipf_haplotypes <- function(tab, margins, tol = 1e-12, maxit = 5000) {
  stopifnot(length(margins) == length(dim(tab)))
  p <- unclass(tab) / sum(tab)
  dn <- dimnames(tab)
  for (it in seq_len(maxit)) {
    old <- p
    for (k in seq_along(dim(p))) {
      cur <- apply(p, k, sum)
      f <- margins[[k]][dn[[k]]] / cur
      f[!is.finite(f)] <- 0
      p <- sweep(p, k, f, `*`)
    }
    if (max(abs(p - old)) < tol) break
  }
  err <- max(unlist(lapply(seq_along(dim(p)), function(k) abs(apply(p, k, sum) - margins[[k]][dn[[k]]]))))
  attr(p, "iterations") <- it
  attr(p, "max_margin_error") <- err
  p
}

# LD summary of a 2x2 table of haplotype frequencies. The first level of each dimension is
# taken as allele 1.
ld_stats <- function(tab2) {
  p <- unclass(tab2) / sum(tab2)
  pa <- sum(p[1, ]); pb <- sum(p[, 1])
  D <- p[1, 1] - pa * pb
  dmax <- if (D >= 0) min(pa * (1 - pb), (1 - pa) * pb) else min(pa * pb, (1 - pa) * (1 - pb))
  c(freq_1 = pa, freq_2 = pb, D = D, dprime = D / dmax,
    r2 = D^2 / (pa * (1 - pa) * pb * (1 - pb)))
}

# Alternative rescaling for a pair: keep D' (with its sign) and move the margins.
rescale_pair_dprime <- function(tab2, margin_a, margin_b) {
  s <- ld_stats(tab2)
  pa <- margin_a[[1]]; pb <- margin_b[[1]]
  dmax <- if (s[["D"]] >= 0) min(pa * (1 - pb), (1 - pa) * pb) else min(pa * pb, (1 - pa) * (1 - pb))
  D <- s[["dprime"]] * dmax
  f11 <- pa * pb + D
  out <- matrix(c(f11, pb - f11, pa - f11, 1 - pa - pb + f11), 2, 2,
                dimnames = dimnames(tab2))
  out
}
