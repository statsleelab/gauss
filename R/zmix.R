#' Estimate ancestry proportions from Z-scores
#'
#' Runs the zmix preparation (prep_zmix5) and constrained quadratic programming
#' to estimate ancestry proportions at the population or superpopulation level.
#'
#' @param input_file File name of input data containing rsid, chr, bp, a1, a2, and z.
#' @param reference_index_file File name of reference panel index data.
#' @param reference_data_file File name of reference panel data.
#' @param reference_pop_desc_file File name of reference panel population description data.
#' @param percentile Percentile cutoff for normalized variance (see prep_zmix5).
#' @param interval Interval for selecting SNPs (see prep_zmix5).
#' @param level Resolution of ancestry proportions: "superpopulation" or "population".
#' @return A data.frame of ancestry proportions at the requested level.
#' @export
zmix <- function(input_file,
                 reference_index_file,
                 reference_data_file,
                 reference_pop_desc_file,
                 percentile = 0.99,
                 interval = 1,
                 level = c("superpopulation", "population")) {
  level <- match.arg(level)

  if (level == "superpopulation") {
    mat <- prep_zmix5_sup(
      input_file = input_file,
      reference_index_file = reference_index_file,
      reference_data_file = reference_data_file,
      reference_pop_desc_file = reference_pop_desc_file,
      percentile = percentile,
      interval = interval
    )
  } else {
    mat <- prep_zmix5(
      input_file = input_file,
      reference_index_file = reference_index_file,
      reference_data_file = reference_data_file,
      reference_pop_desc_file = reference_pop_desc_file,
      percentile = percentile,
      interval = interval
    )
  }
  mat <- as.matrix(mat)
  if (ncol(mat) < 2) {
    stop("zmix: prep_zmix5 returned an invalid matrix.")
  }

  keep <- is.finite(rowSums(mat))
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) {
    stop("zmix: no valid rows after filtering.")
  }

  y <- mat[, 1, drop = FALSE]
  x <- mat[, -1, drop = FALSE]

  pop_desc <- utils::read.table(
    reference_pop_desc_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("Population_Abbreviation", "Super_Population") %in% names(pop_desc))) {
    stop("zmix: reference_pop_desc_file must include Population_Abbreviation and Super_Population columns.")
  }
  pops <- pop_desc$Population_Abbreviation
  suppops <- pop_desc$Super_Population

  if (level == "superpopulation") {
    sup_order <- suppops[!duplicated(suppops)]
    if (ncol(x) != length(sup_order)) {
      stop("zmix: column count mismatch between prep_zmix5_sup output and superpopulation metadata.")
    }
    colnames(x) <- sup_order
  } else {
    if (ncol(x) != length(pops)) {
      stop("zmix: column count mismatch between prep_zmix5 output and population metadata.")
    }
    colnames(x) <- pops
  }

  dmat <- crossprod(x)
  dvec <- drop(crossprod(y, x))
  ncolx <- ncol(x)
  amat <- cbind(rep(1, ncolx), diag(ncolx), -diag(ncolx))
  bvec <- c(1, rep(0, ncolx), rep(-1, ncolx))

  sol <- quadprog::solve.QP(
    Dmat = dmat,
    dvec = dvec,
    Amat = amat,
    bvec = bvec,
    meq = 1
  )
  weights <- sol$solution
  weights <- weights / sum(weights)
  weights <- round(weights, 5)
  weights <- weights / sum(weights)

  if (level == "superpopulation") {
    data.frame(
      SuperPopulation = colnames(x),
      Weight = weights,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Population = pops,
      SuperPopulation = suppops,
      Weight = weights,
      stringsAsFactors = FALSE
    )
  }
}
