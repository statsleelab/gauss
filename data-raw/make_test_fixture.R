## Build the synthetic reference panel used by the test suite.
##
##   Rscript data-raw/make_test_fixture.R
##
## The real 33KG panel is multi-gigabyte and is not distributed with the
## package, so nothing that touches a reference panel can be tested from a
## fresh clone.  This script writes a miniature panel with the same file
## formats -- 3 populations x 20 samples x 200 SNPs -- into inst/extdata/.
## Everything is deterministic given the seed below; the generated files are
## committed, so this script only needs to be re-run when the fixture itself
## changes.
##
## Files written to inst/extdata/:
##   testpanel_pop_desc.txt         population description (tab-delimited)
##   testpanel_index.gz             BGZF index: rsid chr bp a1 a2 af1ref fpos
##   testpanel_geno.gz              BGZF genotypes: 3 genotype strings + 3 af1
##   testpanel_zscore.txt           study Z-scores, every SNP present in panel
##   testpanel_zscore_unmatched.txt study Z-scores incl. SNPs absent from panel
##   testpanel_af1.txt              study allele frequencies (afmix / cpw2)
##   testpanel_annotation.txt       SNP annotation (jepeg / jepegmix)
##   testpanel_index_badfpos.gz     index with fpos = -1 at bp 100,000

set.seed(20260905)

out_dir <- file.path("inst", "extdata")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## ---------------------------------------------------------------- BGZF ----
##
## BGZF is gzip with a mandatory extra subfield ("BC") carrying the compressed
## block size, which is what makes random access possible.  A block is
##
##   16 bytes of fixed header | 2 bytes BSIZE-1 | deflate data | CRC32 | ISIZE
##
## R's gzfile() writes a standard 10-byte gzip header followed by exactly the
## "deflate data | CRC32 | ISIZE" tail we need, so we borrow that and prepend
## our own 18-byte BGZF header.

BGZF_HEADER <- as.raw(c(
  0x1f, 0x8b, 0x08, 0x04,       # ID1, ID2, CM=deflate, FLG=FEXTRA
  0x00, 0x00, 0x00, 0x00,       # MTIME = 0 (keeps output reproducible)
  0x00, 0xff,                   # XFL, OS = unknown
  0x06, 0x00,                   # XLEN = 6
  0x42, 0x43, 0x02, 0x00        # SI1='B', SI2='C', SLEN = 2
))

## The 28-byte empty block that marks a well-formed end of file.
BGZF_EOF <- as.raw(c(
  0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
  0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
  0x00, 0x00, 0x00, 0x00
))

uint16le <- function(x) {
  x <- as.integer(x)
  as.raw(c(bitwAnd(x, 255L), bitwAnd(bitwShiftR(x, 8L), 255L)))
}

## deflate stream + CRC32 + ISIZE, taken from a throw-away gzip file.
deflate_body <- function(bytes) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  con <- gzfile(tmp, "wb")
  writeBin(bytes, con)
  close(con)
  raw_gz <- readBin(tmp, "raw", n = file.size(tmp))
  raw_gz[-seq_len(10)]
}

bgzf_block <- function(bytes) {
  stopifnot(length(bytes) <= 65280L)   # leave room for the 256-byte overhead
  body <- deflate_body(bytes)
  bsize <- 18L + length(body)
  stopifnot(bsize <= 65536L)
  c(BGZF_HEADER, uint16le(bsize - 1L), body)
}

## Write `lines` as a BGZF file.  With one_line_per_block = TRUE every line
## starts a new block, so its BGZF virtual offset is simply
## (byte offset of the block) << 16 -- which is the `fpos` the index stores.
## Returns those virtual offsets.
write_bgzf <- function(path, lines, one_line_per_block = FALSE) {
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  pos <- 0                                   # byte offset of the next block
  offsets <- numeric(length(lines))
  if (one_line_per_block) {
    for (i in seq_along(lines)) {
      blk <- bgzf_block(charToRaw(paste0(lines[i], "\n")))
      offsets[i] <- pos * 65536
      writeBin(blk, con)
      pos <- pos + length(blk)
    }
  } else {
    ## pack the lines into ~16 KB chunks; the index is read sequentially, so
    ## the offsets are not meaningful here
    chunk <- character(0)
    nbytes <- 0
    flush_chunk <- function() {
      if (!length(chunk)) return(invisible(NULL))
      blk <- bgzf_block(charToRaw(paste0(paste(chunk, collapse = "\n"), "\n")))
      writeBin(blk, con)
      invisible(NULL)
    }
    for (ln in lines) {
      if (nbytes + nchar(ln) + 1 > 16384) {
        flush_chunk(); chunk <- character(0); nbytes <- 0
      }
      chunk <- c(chunk, ln)
      nbytes <- nbytes + nchar(ln) + 1
    }
    flush_chunk()
  }
  writeBin(BGZF_EOF, con)
  offsets
}

## ------------------------------------------------------------- panel ------

n_snp   <- 200L
pops    <- c("PPA", "PPB", "PPC")
pop_n   <- c(20L, 20L, 20L)
sup_pop <- c("SPX", "SPX", "SPY")
n_pop   <- length(pops)
chr     <- 1L

bp   <- seq_len(n_snp) * 1000L
rsid <- sprintf("rs%06d", 100000L + seq_len(n_snp))

allele_pairs <- list(c("A", "G"), c("C", "T"), c("A", "T"), c("C", "G"))
picked <- allele_pairs[sample.int(length(allele_pairs), n_snp, replace = TRUE)]
a1 <- vapply(picked, `[`, character(1), 1L)
a2 <- vapply(picked, `[`, character(1), 2L)

## Allele frequencies: a shared baseline plus population-specific drift, with
## every fifth SNP strongly differentiated so that the ancestry-informative
## marker selection in zmix() has something to find.
base_af <- runif(n_snp, 0.25, 0.75)
af <- matrix(0, nrow = n_snp, ncol = n_pop,
             dimnames = list(rsid, pops))
informative <- seq(1L, n_snp, by = 5L)
for (k in seq_len(n_pop)) {
  drift <- rnorm(n_snp, 0, 0.05)
  drift[informative] <- drift[informative] + (k - 2) * 0.18
  af[, k] <- pmin(pmax(base_af + drift, 0.20), 0.80)
}

## Haplotypes come from a latent AR(1) process along the chromosome, so LD
## decays with distance instead of every SNP being independent.
rho <- 0.6
sim_haplotypes <- function(n_hap, freq) {
  z <- matrix(0, nrow = n_snp, ncol = n_hap)
  z[1L, ] <- rnorm(n_hap)
  for (j in 2:n_snp) {
    z[j, ] <- rho * z[j - 1L, ] + sqrt(1 - rho^2) * rnorm(n_hap)
  }
  (z <= matrix(qnorm(freq), nrow = n_snp, ncol = n_hap)) * 1L
}

## dosage[[k]] is n_snp x pop_n[k], entries in {0,1,2}
dosage <- vector("list", n_pop)
for (k in seq_len(n_pop)) {
  h1 <- sim_haplotypes(pop_n[k], af[, k])
  h2 <- sim_haplotypes(pop_n[k], af[, k])
  dosage[[k]] <- h1 + h2
}

## A monomorphic SNP makes every correlation involving it 0/0, which poisons
## the whole LD matrix with NaN.  Resample any such variant.
for (k in seq_len(n_pop)) {
  for (j in seq_len(n_snp)) {
    tries <- 0L
    while (length(unique(dosage[[k]][j, ])) == 1L) {
      tries <- tries + 1L
      stopifnot(tries < 100L)
      dosage[[k]][j, ] <- rbinom(pop_n[k], 2, af[j, k])
    }
  }
}

geno_str <- matrix("", nrow = n_snp, ncol = n_pop)
for (k in seq_len(n_pop)) {
  geno_str[, k] <- apply(dosage[[k]], 1L, function(g) paste(g, collapse = ""))
}

## Per-population allele frequency, as stored in the genotype record.
af1_panel <- matrix(0, nrow = n_snp, ncol = n_pop)
for (k in seq_len(n_pop)) {
  af1_panel[, k] <- round(rowMeans(dosage[[k]]) / 2, 5)
}
## Pooled frequency, as stored in the index.
af1_ref <- round(rowSums(vapply(dosage, rowSums, numeric(n_snp))) /
                   (2 * sum(pop_n)), 5)

## ------------------------------------------------------------- writing ----

pop_desc_path <- file.path(out_dir, "testpanel_pop_desc.txt")
writeLines(c(
  paste("Population_Abbreviation", "Number_of_Subject", "Super_Population",
        "Population_Description", sep = "\t"),
  paste(pops, pop_n, sup_pop,
        c("Test population A", "Test population B", "Test population C"),
        sep = "\t")
), pop_desc_path)

geno_lines <- vapply(seq_len(n_snp), function(j) {
  paste(c(geno_str[j, ], formatC(af1_panel[j, ], format = "f", digits = 5)),
        collapse = " ")
}, character(1))
fpos <- write_bgzf(file.path(out_dir, "testpanel_geno.gz"), geno_lines,
                   one_line_per_block = TRUE)

index_lines <- sprintf("%s %d %d %s %s %s %s",
                       rsid, chr, bp, a1, a2,
                       formatC(af1_ref, format = "f", digits = 5),
                       format(fpos, scientific = FALSE, trim = TRUE))
invisible(write_bgzf(file.path(out_dir, "testpanel_index.gz"), index_lines))

## Study summary statistics: every second SNP is "measured".  A handful are
## written with the alleles the other way round, which exercises the
## orientation fix in ReadReferenceIndex().
measured <- seq(1L, n_snp, by = 2L)
z <- round(rnorm(length(measured)), 6)
swapped <- measured %in% sample(measured, 12L)
z_a1 <- ifelse(swapped, a2[measured], a1[measured])
z_a2 <- ifelse(swapped, a1[measured], a2[measured])

z_lines <- c(
  "rsid chr bp a1 a2 z",
  sprintf("%s %d %d %s %s %s", rsid[measured], chr, bp[measured],
          z_a1, z_a2, format(z, scientific = FALSE, trim = TRUE))
)
writeLines(z_lines, file.path(out_dir, "testpanel_zscore.txt"))

## The same file plus five variants that do not exist in the panel -- what a
## real GWAS input always looks like.  Their positions sit inside the analysis
## window used by the tests so that they actually reach the analysis code.
extra <- sprintf("rsx%05d %d %d %s %s %s", seq_len(5L), chr,
                 seq(60000L, 100000L, by = 10000L) + 1L, "A", "G",
                 format(round(rnorm(5L), 6), scientific = FALSE, trim = TRUE))
writeLines(c(z_lines, extra),
           file.path(out_dir, "testpanel_zscore_unmatched.txt"))

## Study allele frequencies for afmix()/cpw2(): an exact 50/30/20 mixture of
## the three panel populations, so the recovered weights have a known target.
mix_wgt <- c(0.5, 0.3, 0.2)
af1_study <- round(as.vector(af1_panel %*% mix_wgt), 5)
writeLines(c(
  "rsid chr bp a1 a2 af1",
  sprintf("%s %d %d %s %s %s", rsid, chr, bp, a1, a2,
          formatC(af1_study, format = "f", digits = 5))
), file.path(out_dir, "testpanel_af1.txt"))

## Annotation for jepeg()/jepegmix(): 8 genes of 25 consecutive SNPs each.
categs <- c("PROTEIN", "TFBS", "WTH_HAIR", "WTH_TARGET", "CIS_EQTL", "TRANS_EQTL")
gene_id <- sprintf("GENE%02d", rep(seq_len(8L), each = 25L))
annot_categ <- categs[(seq_len(n_snp) - 1L) %% length(categs) + 1L]
annot_wgt <- round(runif(n_snp, 0.5, 3), 3)
writeLines(c(
  "rsid chr bp a1 a2 geneid categ wgt",
  sprintf("%s %d %d %s %s %s %s %s", rsid, chr, bp, a1, a2,
          gene_id, annot_categ,
          format(annot_wgt, scientific = FALSE, trim = TRUE))
), file.path(out_dir, "testpanel_annotation.txt"))

## A copy of the index whose entry at bp = 100,000 points at no genotype record.
## Reading a variant's record starts with a seek to the `fpos` the index gives; a
## negative one can only mean the index and the genotype file do not belong
## together, and the failure used to go unchecked. The position is inside the
## window the tests analyse, so the bad entry is actually reached rather than
## filtered out beforehand. Used by the error-path tests.
bad_index <- index_lines
bad_row <- which(bp == 100000L)
stopifnot(length(bad_row) == 1L)
bad_index[bad_row] <- sub("[0-9]+$", "-1", bad_index[bad_row])
invisible(write_bgzf(file.path(out_dir, "testpanel_index_badfpos.gz"), bad_index))

cat("fixture written to ", out_dir, "\n", sep = "")
print(file.info(list.files(out_dir, full.names = TRUE))["size"])
