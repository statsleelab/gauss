#' PGC2 Schizophrenia GWAS Ancestry Proportion
#'
#' Estimated ancestry proportions of the PGC SCZ2 cohort over the 33KG
#' reference panel, in the shape every mixture analysis expects for
#' \code{pop_wgt_df}: column 1 the population abbreviation, column 2 the
#' mixing proportion.
#'
#' Estimated with \code{\link{afmix}} from the study allele frequencies in
#' \code{data/PGC2_Chr22_ilmn1M_AF1.txt}; 21 of the panel's 29 populations
#' carry a positive weight. The proportions are rescaled to sum to 1, which the
#' mixture covariance assumes; \code{afmix()} itself does not rescale its
#' output, so the values as first published summed to 1.061. Rescaling changes
#' only the scale, not any population's share relative to the others. See
#' \code{data-raw/make_PGC2_SCZ_ANC_Prop.R}.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 21 rows and 2 columns:
#' \describe{
#'   \item{pop}{33KG population abbreviation.}
#'   \item{wgt}{Mixing proportion; non-negative, summing to 1.}
#' }
#'
#' @usage data(PGC2_SCZ_ANC_Prop)
#'
#' @keywords datasets
#'
#' @references
#' Schizophrenia Working Group of Psychiatric Genomics Consortium (2014). 
#' Biological insights from 108 schizophrenia associated genetic loci. Nature
#' \href{https://www.nature.com/articles/nature13595}{Link to article}
#'
#' @examples
#' \dontrun{
#' data(PGC2_SCZ_ANC_Prop)
#' str(PGC2_SCZ_ANC_Prop)
#' }
#'
#' @name PGC2_SCZ_ANC_Prop
NULL