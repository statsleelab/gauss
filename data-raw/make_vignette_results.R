## Regenerate the analysis results the vignettes display.
##
## Every vignette that calls an analysis function needs the 33KG reference
## panel: 43 GB of genotypes that are not, and cannot be, distributed with the
## package.  The vignettes therefore show the real call but fall back to the
## results saved here, so that they build on a machine that has no panel --
## which is every machine except the maintainer's.
##
## Run this script whenever a change moves the numbers (the 1.1.0 mixture
## covariance fix did), then rebuild the vignettes and the pkgdown site.
##
##   GAUSS_REF_33KG=/path/to/ref/33KG Rscript data-raw/make_vignette_results.R
##
## It takes on the order of an hour: distmix() alone is ~21 minutes on the
## 1 Mb window below, and jepeg()/jepegmix() scan all of chromosome 22.
##
## The calls below must stay character-for-character identical to the ones the
## vignettes display, or the shipped results stop being the results of the
## code the reader is shown.

library(gauss)

panel_dir <- Sys.getenv("GAUSS_REF_33KG", unset = "ref/33KG")
if (!dir.exists(panel_dir)) {
  stop("reference panel not found at '", panel_dir,
       "'; set GAUSS_REF_33KG to the directory holding 33kg_index.gz")
}

reference_index_file    <- file.path(panel_dir, "33kg_index.gz")
reference_data_file     <- file.path(panel_dir, "33kg_geno.gz")
reference_pop_desc_file <- file.path(panel_dir, "33kg_pop_desc.txt")

for (f in c(reference_index_file, reference_data_file, reference_pop_desc_file)) {
  if (!file.exists(f)) stop("missing reference panel file: ", f)
}

## The same paths the vignettes use, so the results are the results of the
## documented code and not of some local copy of the inputs.
gauss_data <- function(name) {
  path <- system.file("data", name, package = "gauss")
  if (!nzchar(path)) stop("'", name, "' is not in the installed package")
  path
}

out_dir <- "inst/extdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_result <- function(object, name) {
  path <- file.path(out_dir, paste0("vignette_", name, ".rds"))
  saveRDS(object, path, version = 3, compress = "xz")
  cat(sprintf("  wrote %-34s %8.1f KB\n", path, file.size(path) / 1024))
  invisible(path)
}

run <- function(name, expr) {
  cat("running ", name, " ...\n", sep = "")
  ## The analyses report progress on stdout; keep the log readable.
  sink(tempfile())
  on.exit(sink(), add = TRUE)
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  sink()
  on.exit()
  cat(sprintf("  %.1f s\n", elapsed))
  save_result(value, name)
  invisible(value)
}

data("PGC2_SCZ_ANC_Prop", package = "gauss", envir = environment())
af1_cutoff <- 0.001

## ----------------------------------------------------------- afmix ---------
## vignettes/afmix_example.Rmd
input_file <- gauss_data("PGC2_Chr22_ilmn1M_AF1.txt")

run("afmix", afmix(input_file = input_file,
                   reference_index_file = reference_index_file,
                   reference_data_file = reference_data_file,
                   reference_pop_desc_file = reference_pop_desc_file))

## ------------------------------------------------- computeLD, dist ---------
## vignettes/computeLD_example.Rmd and vignettes/dist_example.Rmd
input_file <- gauss_data("PGC2_3Mb.txt")

run("computeLD", computeLD(chr = 10,
                           start_bp = 104000001,
                           end_bp = 105000000,
                           pop_wgt_df = PGC2_SCZ_ANC_Prop,
                           input_file = input_file,
                           reference_index_file = reference_index_file,
                           reference_data_file = reference_data_file,
                           reference_pop_desc_file = reference_pop_desc_file,
                           af1_cutoff = af1_cutoff))

run("dist", dist(chr = 10,
                 start_bp = 104000001,
                 end_bp = 105000000,
                 wing_size = 500000,
                 study_pop = "EUR",
                 input_file = input_file,
                 reference_index_file = reference_index_file,
                 reference_data_file = reference_data_file,
                 reference_pop_desc_file = reference_pop_desc_file,
                 af1_cutoff = af1_cutoff))

run("distmix", distmix(chr = 10,
                       start_bp = 104000001,
                       end_bp = 105000000,
                       wing_size = 500000,
                       pop_wgt_df = PGC2_SCZ_ANC_Prop,
                       input_file = input_file,
                       reference_index_file = reference_index_file,
                       reference_data_file = reference_data_file,
                       reference_pop_desc_file = reference_pop_desc_file,
                       af1_cutoff = af1_cutoff))

## ----------------------------------------------- jepeg, jepegmix -----------
## vignettes/jepeg_example.Rmd
input_file      <- gauss_data("PGC2_Chr22_ilmn1M_Z.txt")
annotation_file <- gauss_data("JEPEG_SNP_Annotation.v1.0.txt")

## jepeg() is called without af1_cutoff in the vignette; jepegmix() with it.
run("jepeg", jepeg(study_pop = "GBR",
                   input_file = input_file,
                   annotation_file = annotation_file,
                   reference_index_file = reference_index_file,
                   reference_data_file = reference_data_file,
                   reference_pop_desc_file = reference_pop_desc_file))

run("jepegmix", jepegmix(pop_wgt_df = PGC2_SCZ_ANC_Prop,
                         input_file = input_file,
                         annotation_file = annotation_file,
                         reference_index_file = reference_index_file,
                         reference_data_file = reference_data_file,
                         reference_pop_desc_file = reference_pop_desc_file,
                         af1_cutoff = af1_cutoff))

cat("done\n")
