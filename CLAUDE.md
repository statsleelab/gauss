# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`gauss` — an R package (Rcpp + RcppEigen) for analysing GWAS summary statistics from cohorts of
**mixed ancestry**, using an external 32,953-sample / 29-population reference panel (33KG). Nearly
all computation is C++; the R layer is a thin veneer over `.Call` wrappers.

Two companion documents already exist and are worth reading before non-trivial work:

* **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — full architecture: data flow, file formats,
  algorithms, generated files, caution areas.
* **[docs/REVIEW_PLAN.md](docs/REVIEW_PLAN.md)** — 56 catalogued defects with severities and
  fixes, referenced by ID (`C-1`, `PERF-3`, …). Findings below cite these IDs.

## Commands

```bash
# Regenerate Rcpp bindings — REQUIRED after changing any [[Rcpp::export]] signature
Rscript -e 'Rcpp::compileAttributes()'

# Regenerate man/*.Rd and NAMESPACE (run compileAttributes FIRST — see below)
Rscript -e 'devtools::document()'

# Compile + load into a session without installing
Rscript -e 'devtools::load_all()'

# Install
R CMD INSTALL --no-multiarch --with-keep.source .

# Full check (currently FAILS — see "Known-broken" below)
Rscript -e 'devtools::check()'

# Rebuild README.md from README.Rmd
Rscript -e 'devtools::build_readme()'

# Rebuild the pkgdown site (do NOT pass clean=TRUE — it would delete the two docs above)
Rscript -e 'pkgdown::build_site()'
```

### Compiling a single file with warnings on

The package builds warning-suppressed by default. To see real warnings for one translation unit:

```bash
RINC=$(Rscript -e 'cat(R.home("include"))')
RCPP=$(Rscript -e 'cat(system.file("include",package="Rcpp"))')
REIG=$(Rscript -e 'cat(system.file("include",package="RcppEigen"))')
cd src && clang++ -std=gnu++17 -fsyntax-only -Wall -Wextra -Wconditional-uninitialized \
  -I. -I"$RINC" -I"$RCPP" -I"$REIG" -include rcpp_eigen_wrap.h dist.cpp
```

There are ~90 pre-existing `-Wsign-compare` warnings (D-1); filter them out to see anything new.
Run `make clean` in `src/` if stale `.o`/`.so` files shadow your edits.

### Tests

**There is no test suite.** No `tests/`, no testthat, no CI. `inst/dev-tests/` holds 27 gitignored
ad-hoc driver scripts that hard-code absolute paths and need the multi-GB panel — they are not
runnable from a fresh clone and assert nothing.

Practical consequence: **numeric changes cannot be verified by any existing check.** Several open
findings (C-1, C-3, C-5, C-6) alter results without altering behaviour. If you touch the
linear algebra or the correlation kernels, build a fixture first (T-1 proposes a synthetic
3-population × 20-sample × 200-SNP BGZF panel, <200 KB).

## Architecture in one page

### The pipeline every entry point shares

There are 13 exported analyses in `src/<name>.cpp`, but they are variations on one sequence
implemented in `src/gauss.cpp`:

```
read_ref_desc                     panel metadata -> Arguments
init_pop_flag_vec                 [dist/qcat/jepeg]  resolve `study_pop` -> 0/1 flags
init_pop_flag_wgt_vec             [*mix/computeLD]   resolve pop_wgt_df  -> flags + weights
ReadInputZ / ReadInputAf          user file -> map<MapKey, Snp*>, all tagged type 2
ReadReferenceIndex[All]           attach rsid + `fpos`; reclassify to type 0/1; fix orientation
ReadAnnotation                    [jepeg only] attach geneid + (category, weight)
MakeSnpVec / MakeSnpVecMix        seek fpos, compute af1ref/af1mix, apply MAF cutoff
ReadGenotype                      seek fpos, load genotype strings for flagged populations
run_<analysis>                    Eigen linear algebra
FreeGenotype + manual delete of every Snp*
```

### Three concepts that explain most of the code

**`Snp::type_` (0/1/2)** — `0` = in panel but unmeasured (imputation target), `1` = in panel and
measured (predictor), `2` = measured but absent from panel (passed to output, excluded from all
maths). Compared as bare integers in ~20 places across 13 files; there is no enum.

**`Snp::fpos_`** — the BGZF *virtual* offset of that variant's genotype record. The panel is never
loaded; every access is `bgzf_seek(fpos)` + read one record. Type-2 SNPs never get one and keep
the default `-1`, and **no `bgzf_seek` return value is checked anywhere** (C-5, P-3).

**Allele orientation is settled once**, in `ReadReferenceIndex*`, by rewriting the *study* side to
match the panel (flip `z`, swap `a1`/`a2`, re-key the map). This is why the genotype-flipping code
in `ReadGenotype` is commented out and `Snp::flip_` is vestigial.

### Genotypes are ASCII strings

One `std::string` per population, one character per individual, each in `{'0','1','2'}` (reference
allele count). All correlation kernels in `src/util.cpp` stream over those characters. A 33KG
genotype record is ~33 KB.

## Repo-specific gotchas

**`compileAttributes()` must run before `document()`.** Roxygen for C++ functions lives as `//'`
comments in the `.cpp` files. `compileAttributes()` copies them into `R/RcppExports.R`, and only
then does roxygen see them. Editing a C++ doc block without this step silently changes nothing.

**Never hand-edit** `R/RcppExports.R`, `src/RcppExports.cpp`, `man/*.Rd`, `NAMESPACE`, `README.md`,
`src/Makevars*`, or `docs/**` other than `ARCHITECTURE.md` / `REVIEW_PLAN.md`. All are generated.

**`CalWgtCov()` (`src/util.cpp:100`) is on the critical path of every mixture analysis** —
`distmix`, `computeLD`, `jepegmix`, `prep_recessive_impute` all reduce to it. It currently returns
a covariance inflated by `m²` (C-1, the package's most consequential defect). Fixing it changes
every published mixture number, so it needs a test fixture and a deliberate re-validation decision
— not a drive-by patch.

**Error paths leak the entire SNP map.** `Snp` objects are bare `new`, freed only by an explicit
loop at the bottom of each entry point; `Rcpp::stop()` throws past it. This includes the routine
`"Not enough number of SNPs loaded"` message (M-1). Put new validation *before* allocation, or
convert to `unique_ptr` ownership first.

**Four copies of the same parser.** `ReadInputZ`/`ReadReferenceIndex` in `gauss.cpp` have private
near-twins in `zmix.cpp` (`read_input_zmix`, `read_ref_index_zmix`) and `cpw2.cpp`
(`read_input_cpw2`, `read_ref_index_cpw2`). Parser fixes must be applied to all four. Likewise
`dist`/`distmix`, `qcat`/`qcatmix` and `jepeg`/`jepegmix` are pairwise copy-edits — a change to
one almost always belongs in its twin.

**`pop_wgt_df` is read positionally** (`pop_wgt_df[0]`, `pop_wgt_df[1]`) in five entry points, so
it silently misinterprets any frame whose first two columns aren't `(pop, wgt)`. `afmix()` returns
*three* columns (`sup.pop, pop, wgt`), so its output does **not** feed `distmix()`/`jepegmix()`
directly despite what the README workflow implies (C-7). The canonical shape is
`data(PGC2_SCZ_ANC_Prop)`.

**Numerics fail silently rather than loudly.** `CholeskyMat` ignores `LLT::info()` (N-1);
`MakePosDef`/`RmvPC` no-op when the eigensolver fails (N-5); every correlation is an unguarded
division, so a monomorphic SNP yields `0/0` and the `NaN` spreads through the whole LD matrix
(N-3). Prefer failing loudly when touching this code.

**The reference panel is not in the repo.** `ref/` is gitignored and multi-gigabyte, so no panel
code path is exercisable from a fresh clone.

## Known-broken (do not treat as regressions you caused)

* `devtools::check()` fails — vignettes execute at build time but need undeclared packages
  (`tidyverse`, `data.table`, `kableExtra`, `dplyr`, `reshape2`) and read `../ref/33KG/`, which
  isn't distributed (R-5).
* `stats` and `utils` are used at runtime but absent from `Imports`/`NAMESPACE` (R-2).
* `data/` ships 63 MB of raw `.txt`, well past CRAN's 5 MB tarball limit (R-1).
* `gauss::dist()` masks `stats::dist()` (R-3).
* `src/etc/` and `src/zmixrcpp.cpp_` are stale pre-Rcpp/GSL duplicates, not live code (D-2).
* `NAMESPACE` uses `exportPattern("^[[:alpha:]]+")`, so experimental scaffolding
  (`cpw2`, `prep_zmix`–`prep_zmix4`) is exported and documented alongside the real API (R-4).
