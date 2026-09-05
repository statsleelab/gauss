## The fixture panel is the foundation every other test stands on, so check
## that it still holds the invariants the C++ side assumes.

test_that("the index describes 200 variants on chromosome 1", {
  index <- panel_records()$index

  expect_equal(nrow(index), panel_n_snp)
  expect_true(all(index$chr == 1L))
  expect_equal(index$bp, panel_bp)
  expect_false(anyDuplicated(index$rsid) > 0)
  expect_true(all(index$a1 != index$a2))
  expect_true(all(index$a1 %in% c("A", "C", "G", "T")))
  expect_true(all(index$a2 %in% c("A", "C", "G", "T")))
})

test_that("every index entry points at a distinct genotype record", {
  fpos <- panel_records()$index$fpos

  expect_equal(length(unique(fpos)), panel_n_snp)
  expect_true(all(diff(fpos) > 0))
  ## fpos is a BGZF virtual offset; one record per block means the low 16 bits
  ## (the offset within the uncompressed block) are always zero.
  expect_true(all(fpos %% 65536 == 0))
})

test_that("genotype strings match the declared population sizes", {
  ## Nothing in the package verifies this, and a mismatch would silently
  ## misalign every population in the mixture kernels.
  desc <- utils::read.table(panel_file("testpanel_pop_desc.txt"),
                            header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  expect_equal(desc$Population_Abbreviation, panel_pops)
  expect_equal(desc$Number_of_Subject, panel_pop_size)
  expect_equal(desc$Super_Population, panel_sup_pops)

  geno <- panel_records()$geno
  for (k in seq_along(panel_pops)) {
    expect_true(all(nchar(geno[, k]) == panel_pop_size[k]),
                info = panel_pops[k])
  }
})

test_that("genotypes are reference-allele counts in {0,1,2}", {
  geno <- panel_records()$geno
  expect_false(any(grepl("[^012]", geno)))
})

test_that("no variant is monomorphic in any population", {
  ## A monomorphic variant divides by zero in every correlation kernel and the
  ## resulting NaN spreads through the whole LD matrix.
  for (i in seq_len(panel_n_snp)) {
    spread <- vapply(panel_dosage(i), function(d) length(unique(d)), integer(1))
    expect_true(all(spread > 1L), info = paste("variant", i))
  }
})

test_that("stored allele frequencies match the genotypes", {
  af1 <- panel_records()$af1
  for (i in seq_len(panel_n_snp)) {
    computed <- vapply(panel_dosage(i), function(d) mean(d) / 2, numeric(1))
    expect_equal(af1[i, ], round(computed, 5), tolerance = 1e-9,
                 info = paste("variant", i))
  }
  ## Everything stays clear of the default 0.01 / 0.05 frequency cutoffs.
  expect_true(all(af1 > 0.05 & af1 < 0.95))
})

test_that("the study Z-score file covers half the panel", {
  z <- utils::read.table(panel_file("testpanel_zscore.txt"), header = TRUE,
                         stringsAsFactors = FALSE)

  expect_equal(names(z), c("rsid", "chr", "bp", "a1", "a2", "z"))
  expect_equal(nrow(z), panel_n_snp / 2)
  expect_equal(z$bp, panel_bp[seq(1L, panel_n_snp, by = 2L)])
  expect_true(all(is.finite(z$z)))

  ## Some rows carry the alleles the other way round, which is what makes the
  ## orientation fix in ReadReferenceIndex() testable.
  index <- panel_records()$index[panel_row(z$bp), ]
  expect_gt(sum(z$a1 != index$a1), 0)
  expect_true(all(z$a1 == index$a1 | z$a1 == index$a2))
})

test_that("the unmatched-input file adds variants absent from the panel", {
  z <- utils::read.table(panel_file("testpanel_zscore_unmatched.txt"),
                         header = TRUE, stringsAsFactors = FALSE)
  extra <- z[grepl("^rsx", z$rsid), ]

  expect_equal(nrow(extra), 5L)
  expect_true(all(is.na(panel_row(extra$bp))))
})
