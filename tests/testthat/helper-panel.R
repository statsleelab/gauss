## Shared access to the synthetic reference panel in inst/extdata, plus the R
## reference implementations the correlation kernels are checked against.
##
## The panel is 3 populations x 20 samples x 200 SNPs on chromosome 1, written
## by data-raw/make_test_fixture.R.  It uses the same file formats as the real
## 33KG panel, so every code path that reads a reference panel is exercisable
## from a fresh clone.

panel_pops     <- c("PPA", "PPB", "PPC")
panel_pop_size <- c(20L, 20L, 20L)
panel_sup_pops <- c("SPX", "SPX", "SPY")
panel_n_snp    <- 200L
panel_bp       <- seq_len(panel_n_snp) * 1000L

panel_file <- function(name) {
  path <- system.file("extdata", name, package = "gauss")
  if (!nzchar(path)) {
    skip(paste0("test fixture '", name, "' is not installed"))
  }
  path
}

## The file arguments every entry point takes.
panel <- function() {
  list(
    index       = panel_file("testpanel_index.gz"),
    data        = panel_file("testpanel_geno.gz"),
    pop_desc    = panel_file("testpanel_pop_desc.txt"),
    z           = panel_file("testpanel_zscore.txt"),
    z_unmatched = panel_file("testpanel_zscore_unmatched.txt"),
    af1         = panel_file("testpanel_af1.txt"),
    annot       = panel_file("testpanel_annotation.txt")
  )
}

## pop_wgt_df is read positionally: column 1 = population, column 2 = weight.
panel_wgt_df <- function(wgt = c(0.5, 0.3, 0.2), pops = panel_pops) {
  data.frame(pop = pops, wgt = wgt, stringsAsFactors = FALSE)
}

## The analyses print progress to stdout; swallow it so test output stays
## readable.  Errors propagate unchanged.
quietly <- function(expr) {
  out <- tempfile()
  sink(out)
  on.exit({
    sink()
    unlink(out)
  }, add = TRUE)
  force(expr)
}

## ------------------------------------------------------- panel contents ----
##
## BGZF is gzip with an extra header field, so R's own gzip reader can read the
## panel files sequentially.  That gives the tests an independent view of the
## genotypes the C++ side loads by seeking.

.panel_cache <- new.env(parent = emptyenv())

panel_records <- function() {
  if (is.null(.panel_cache$rec)) {
    n_pop <- length(panel_pops)
    fields <- strsplit(
      trimws(readLines(gzfile(panel_file("testpanel_geno.gz")))),
      "[[:space:]]+"
    )
    .panel_cache$rec <- list(
      index = utils::read.table(
        gzfile(panel_file("testpanel_index.gz")),
        header = FALSE, stringsAsFactors = FALSE,
        col.names = c("rsid", "chr", "bp", "a1", "a2", "af1ref", "fpos")
      ),
      geno = t(vapply(fields, function(f) f[seq_len(n_pop)], character(n_pop))),
      af1 = t(vapply(fields,
                     function(f) as.numeric(f[n_pop + seq_len(n_pop)]),
                     numeric(n_pop)))
    )
  }
  .panel_cache$rec
}

## Genotype strings of panel variant `i`, exactly as the C++ side sees them:
## one string per population, in panel order.
panel_geno_vec <- function(i, pops = seq_along(panel_pops)) {
  unname(panel_records()$geno[i, pops])
}

## The same variant as numeric dosages, one vector per population.
panel_dosage <- function(i, pops = seq_along(panel_pops)) {
  lapply(panel_geno_vec(i, pops),
         function(s) as.numeric(strsplit(s, "", fixed = TRUE)[[1]]))
}

## Row index of the panel variants at the given base-pair positions.
panel_row <- function(bp) match(bp, panel_records()$index$bp)

## --------------------------------------------- reference implementations ----

## Correlation over the pooled samples of the supplied populations.
ref_pooled_cor <- function(x, y) {
  stats::cor(unlist(x), unlist(y))
}

## Mixture correlation: moments are population averages weighted by wgt, i.e.
## E_w[f] = sum_k w_k mean_k(f).  This is what CalCor(x, y, pop_wgt) computes.
ref_wgt_cor <- function(x, y, wgt) {
  mean_k <- function(v, f) sum(wgt * vapply(v, f, numeric(1)))
  ex  <- mean_k(x, mean)
  ey  <- mean_k(y, mean)
  ex2 <- mean_k(x, function(a) mean(a^2))
  ey2 <- mean_k(y, function(a) mean(a^2))
  exy <- sum(wgt * mapply(function(a, b) mean(a * b), x, y))
  (exy - ex * ey) / sqrt((ex2 - ex^2) * (ey2 - ey^2))
}

## The mixture covariance, transcribed from Equation (2) of the DISTMIX paper
## (Lee et al., Bioinformatics 2015, doi:10.1093/bioinformatics/btv348, S2.1):
##
##   mu = sum_k w_k mu_k                                            (1)
##   C  = sum_k w_k C_k + sum_k w_k (mu_k - mu)(mu_k - mu)^T        (2)
##
## with C_k the unbiased variance-covariance matrix of population k.  Written in
## the paper's centred form on purpose: the expanded form the C++ code
## accumulates is equivalent only when sum_k w_k == 1, and keeping the centred
## version here means the test would catch it if that ever stopped holding.
ref_wgt_cov <- function(x, y, wgt) {
  mu_xk <- vapply(x, mean, numeric(1))
  mu_yk <- vapply(y, mean, numeric(1))
  mu_x <- sum(wgt * mu_xk)                                  # Eq. (1)
  mu_y <- sum(wgt * mu_yk)
  within <- sum(wgt * mapply(stats::cov, x, y))             # sum_k w_k C_k
  between <- sum(wgt * (mu_xk - mu_x) * (mu_yk - mu_y))     # centred outer product
  within + between
}

## The mixture correlation the paper's step 2 asks for: normalise C by the
## genotype standard deviations.  This is what distmix() and computeLD() build
## their LD matrices from.
ref_wgt_cov_cor <- function(x, y, wgt) {
  ref_wgt_cov(x, y, wgt) /
    sqrt(ref_wgt_cov(x, x, wgt) * ref_wgt_cov(y, y, wgt))
}

## Write `lines` to a temporary file and return its path.
temp_lines <- function(lines) {
  path <- tempfile(fileext = ".txt")
  writeLines(lines, path)
  path
}
