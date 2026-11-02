#' JEPEG Functional SNP Annotation, version 1.0
#'
#' The SNP-to-gene functional annotation that \code{\link{jepeg}} and
#' \code{\link{jepegmix}} read through their \code{annotation_file} argument.
#' Each row ties one variant to one gene under one functional category and
#' carries the weight that category's test applies to it.
#'
#' A variant appears once per (gene, category) pair it belongs to, so the
#' 509,570 rows cover 393,303 distinct variants over 21,458 genes. 1,651 rows
#' carry \code{.} in place of an rsid; those variants are matched on position
#' and alleles instead.
#'
#' The weight's scale and sign are category-specific — the JEPEG method uses a
#' different proxy measure for each, a deleteriousness score for protein
#' function and a free-energy change for microRNA binding, among others — so
#' weights are comparable within a category and not across categories.
#' \code{ReadAnnotation} recognises exactly the six category names below and
#' skips any row with another, reporting the count.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 509,570 rows and 8 columns:
#' \describe{
#'   \item{rsid}{dbSNP rs identifier, or \code{.} where none is recorded.}
#'   \item{chr}{Chromosome, \code{1}-\code{22}, \code{X} or \code{Y}.}
#'   \item{bp}{Base pair position.}
#'   \item{a1}{Reference allele.}
#'   \item{a2}{Alternative allele.}
#'   \item{geneid}{Gene symbol.}
#'   \item{categ}{Functional category: \code{PROTEIN} (260,518 rows),
#'     \code{TRANS_EQTL} (94,925), \code{CIS_EQTL} (75,118), \code{TFBS}
#'     (58,695), \code{WTH_TARGET} (19,660, within a microRNA target site) or
#'     \code{WTH_HAIR} (654, within a microRNA hairpin).}
#'   \item{wgt}{Signed weight for this variant in this category.}
#' }
#'
#' @usage data(JEPEG_SNP_Annotation.v1.0)
#'
#' @keywords datasets
#'
#' @references
#' Lee D, Williamson VS, Bigdeli TB, Riley BP, Webb BT, Fanous AH, Kendler KS,
#' Vladimirov VI, Bacanu SA (2015). JEPEG: a summary statistics based tool for
#' gene-level joint testing of functional variants. \emph{Bioinformatics}
#' 31(8):1176-1182. \doi{10.1093/bioinformatics/btu816}
#'
#' @examples
#' \dontrun{
#' data(JEPEG_SNP_Annotation.v1.0)
#' str(JEPEG_SNP_Annotation.v1.0)
#' }
#'
#' @name JEPEG_SNP_Annotation.v1.0
NULL


#' PGC2 Schizophrenia Z-scores, 3 Mb Window on Chromosome 10
#'
#' Association Z-scores for 721 variants in a 2.99 Mb window on chromosome 10
#' (103,008,802-105,997,104), from the Psychiatric Genomics Consortium's second
#' schizophrenia GWAS wave (PGC SCZ2).
#'
#' The small window makes this the example input for the functions whose cost
#' grows with the square of the number of variants in the region:
#' \code{\link{computeLD}} and \code{\link{dist}}.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 721 rows and 6 columns:
#' \describe{
#'   \item{rsid}{dbSNP rs identifier.}
#'   \item{chr}{Chromosome; \code{10} throughout.}
#'   \item{bp}{Base pair position.}
#'   \item{a1}{Reference allele.}
#'   \item{a2}{Alternative allele.}
#'   \item{z}{Association Z-score.}
#' }
#'
#' @usage data(PGC2_3Mb)
#'
#' @keywords datasets
#'
#' @references
#' Schizophrenia Working Group of the Psychiatric Genomics Consortium (2014).
#' Biological insights from 108 schizophrenia-associated genetic loci.
#' \emph{Nature} 511(7510):421-427. \doi{10.1038/nature13595}
#'
#' @examples
#' \dontrun{
#' data(PGC2_3Mb)
#' str(PGC2_3Mb)
#' }
#'
#' @name PGC2_3Mb
NULL


#' PGC2 Schizophrenia Allele Frequencies, Chromosome 22
#'
#' Study reference-allele frequencies for the 13,362 chromosome 22 variants of
#' the Psychiatric Genomics Consortium's second schizophrenia GWAS wave (PGC
#' SCZ2) that were measured on the Illumina 1M array.
#'
#' This is the example input for \code{\link{afmix}}, which estimates a
#' cohort's ancestry proportions by matching its allele frequencies against the
#' reference panel. It covers the same variants, in the same order, as
#' \code{\link{PGC2_Chr22_ilmn1M_Z}}, and the proportions estimated from it are
#' shipped as \code{\link{PGC2_SCZ_ANC_Prop}}.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 13,362 rows and 6 columns:
#' \describe{
#'   \item{rsid}{dbSNP rs identifier.}
#'   \item{chr}{Chromosome; \code{22} throughout.}
#'   \item{bp}{Base pair position.}
#'   \item{a1}{Reference allele.}
#'   \item{a2}{Alternative allele.}
#'   \item{af1}{Frequency of \code{a1} in the study cohort.}
#' }
#'
#' @usage data(PGC2_Chr22_ilmn1M_AF1)
#'
#' @keywords datasets
#'
#' @references
#' Schizophrenia Working Group of the Psychiatric Genomics Consortium (2014).
#' Biological insights from 108 schizophrenia-associated genetic loci.
#' \emph{Nature} 511(7510):421-427. \doi{10.1038/nature13595}
#'
#' @examples
#' \dontrun{
#' data(PGC2_Chr22_ilmn1M_AF1)
#' str(PGC2_Chr22_ilmn1M_AF1)
#' }
#'
#' @name PGC2_Chr22_ilmn1M_AF1
NULL


#' PGC2 Schizophrenia Z-scores, Chromosome 22
#'
#' Association Z-scores for the 13,362 chromosome 22 variants of the
#' Psychiatric Genomics Consortium's second schizophrenia GWAS wave (PGC SCZ2)
#' that were measured on the Illumina 1M array, spanning 16,054,667-51,213,613.
#'
#' This is the example input for \code{\link{jepeg}}. It covers the same
#' variants, in the same order, as \code{\link{PGC2_Chr22_ilmn1M_AF1}}.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 13,362 rows and 6 columns:
#' \describe{
#'   \item{rsid}{dbSNP rs identifier.}
#'   \item{chr}{Chromosome; \code{22} throughout.}
#'   \item{bp}{Base pair position.}
#'   \item{a1}{Reference allele.}
#'   \item{a2}{Alternative allele.}
#'   \item{z}{Association Z-score.}
#' }
#'
#' @usage data(PGC2_Chr22_ilmn1M_Z)
#'
#' @keywords datasets
#'
#' @references
#' Schizophrenia Working Group of the Psychiatric Genomics Consortium (2014).
#' Biological insights from 108 schizophrenia-associated genetic loci.
#' \emph{Nature} 511(7510):421-427. \doi{10.1038/nature13595}
#'
#' @examples
#' \dontrun{
#' data(PGC2_Chr22_ilmn1M_Z)
#' str(PGC2_Chr22_ilmn1M_Z)
#' }
#'
#' @name PGC2_Chr22_ilmn1M_Z
NULL


#' PGC3 Schizophrenia Z-scores, Genome-wide
#'
#' Association Z-scores for 791,425 variants across chromosomes 1-22, from the
#' Psychiatric Genomics Consortium's third schizophrenia GWAS wave (PGC SCZ3),
#' restricted to variants measured on the Illumina 1M array.
#'
#' This is the example input for \code{\link{fiqt}}, whose FDR inverse quantile
#' transformation corrects Winner's Curse bias and therefore needs a
#' genome-wide set of Z-scores rather than one region.
#'
#' @docType data
#'
#' @format A \code{data.frame} with 791,425 rows and 6 columns:
#' \describe{
#'   \item{rsid}{dbSNP rs identifier.}
#'   \item{chr}{Chromosome, \code{1}-\code{22}.}
#'   \item{bp}{Base pair position.}
#'   \item{a1}{Reference allele.}
#'   \item{a2}{Alternative allele.}
#'   \item{z}{Association Z-score.}
#' }
#'
#' @usage data(PGC3_SCZ_ilmn1M_Z)
#'
#' @keywords datasets
#'
#' @references
#' Trubetskoy V, Pardinas AF, Qi T, et al. (2022). Mapping genomic loci
#' implicates genes and synaptic biology in schizophrenia. \emph{Nature}
#' 604(7906):502-508. \doi{10.1038/s41586-022-04434-5}
#'
#' @examples
#' \dontrun{
#' data(PGC3_SCZ_ilmn1M_Z)
#' str(PGC3_SCZ_ilmn1M_Z)
#' }
#'
#' @name PGC3_SCZ_ilmn1M_Z
NULL
