# GAUSS — Architecture

**Package:** `gauss` (Genome Analysis Using Summary Statistics) · **Version:** 1.0
**Repository:** <https://github.com/statsleelab/gauss> · **Site:** <https://statsleelab.github.io/gauss>
**Document status:** descriptive reference for contributors. For known defects and remediation
priorities see [REVIEW_PLAN.md](REVIEW_PLAN.md); cross-references below use its finding IDs
(`C-1`, `PERF-3`, …).

---

## 1. Purpose

GAUSS analyses **GWAS summary statistics** — per-variant Z-scores or allele frequencies — without
requiring individual-level genotypes. Its distinguishing feature is that it works for **cohorts of
mixed ancestry**, by modelling the study sample as a weighted mixture of populations drawn from a
large external reference panel.

The reference panel (33KG) holds 32,953 samples across 29 populations: 20,281 European,
10,800 East Asian, 817 African, 533 Native American, 522 South Asian.

Five capabilities, in the order a user typically needs them:

| Goal | Function(s) | Needs |
|---|---|---|
| Estimate the cohort's ancestry mixture | `afmix()` (from allele frequencies), `zmix()` (from Z-scores only) | panel |
| Compute ancestry-informed LD | `computeLD()`, `simulateLD()` | panel + mixture weights |
| Impute Z-scores of unmeasured variants | `dist()` (homogeneous), `distmix()` (mixed) | panel (+ weights) |
| Gene-level / transcriptome-wide association | `jepeg()` (homogeneous), `jepegmix()` (mixed) | panel + SNP annotation (+ weights) |
| Correct Winner's Curse bias | `fiqt()` | nothing — pure R |

The decision the whole API hangs on: **is the cohort ethnically homogeneous?** If yes, use the
plain variants (`dist`, `qcat`, `jepeg`) with a single `study_pop` name. If not, estimate mixing
proportions first, then use the `*mix` variants.

Citation: Lee & Bacanu, *Bioinformatics* (2024), <https://doi.org/10.1093/bioinformatics/btae203>.

---

## 2. Repository layout

```
DESCRIPTION  NAMESPACE          package metadata; NAMESPACE is roxygen-generated
configure                       copies src/Makevars.in -> src/Makevars  (see §10, R-6)
configure.win                   copies src/Makevars.in -> src/Makevars.win

R/                              thin R layer
  gauss.R                       package-level roxygen block (useDynLib, exportPattern)
  zmix.R                        zmix(): the only substantial hand-written R function
  fiqt.R                        fiqt(): pure-R Winner's Curse correction
  PGC2_SCZ_ANC_Prop.R           dataset documentation stub
  RcppExports.R                 GENERATED — .Call wrappers + roxygen copied from src/

src/                            all computation
  gauss.{h,cpp}                 Arguments struct + shared I/O engine (~1,400 lines)
  snp.{h,cpp}                   Snp: per-variant state container
  gene.{h,cpp}                  Gene, Categ: JEPEG gene-level test
  util.{h,cpp}                  correlation kernels + Eigen matrix helpers
  bgzf.{c,h}, khash.h           VENDORED samtools BGZF reader
  rcpp_eigen_wrap.h             warning-suppressing wrapper around <RcppEigen.h>
  <entry-point>.cpp             one file per exported analysis (13 files)
  RcppExports.cpp               GENERATED — .Call registration table
  Makevars.in                   template; Makevars/Makevars.win are generated
  etc/                          STALE duplicates from the pre-Rcpp / GSL era (see D-2)

man/                            GENERATED — roxygen2 .Rd files
docs/                           GENERATED — pkgdown site (this file is a hand-written exception)
vignettes/                      7 .Rmd articles (currently unbuildable — see R-5)
data/                           example inputs, 63 MB of raw .txt (see R-1)
ref/33KG/                       reference panel — gitignored, user-supplied
inst/dev-tests/                 27 ad-hoc driver scripts — gitignored, not a test suite
```

---

## 3. Major R entry points

`NAMESPACE` uses `exportPattern("^[[:alpha:]]+")`, so **everything is exported** — including
experimental scaffolding (see R-4). The functions divide into three tiers.

### 3.1 User-facing API

| Function | Defined in | Returns | Notes |
|---|---|---|---|
| `afmix(input_file, ref×3, interval=1000)` | `src/afmix.cpp` | `data.frame(sup.pop, pop, wgt)` | ancestry from allele frequencies |
| `zmix(input_file, ref×3, percentile=0.9, interval=10, level=)` | `R/zmix.R` | `data.frame(Population, SuperPopulation, Weight)` | ancestry from Z-scores; the only R-side algorithm |
| `computeLD(chr, start_bp, end_bp, pop_wgt_df, input_file, ref×3, af1_cutoff)` | `src/computeLD.cpp` | `list(snplist, cormat)` | mixture LD for measured SNPs |
| `simulateLD(..., sim_size, ...)` | `src/simulateLD.cpp` | `list(snplist, cormat)` | LD via resampled pseudo-individuals |
| `dist(chr, start_bp, end_bp, wing_size, study_pop, input_file, ref×3, af1_cutoff)` | `src/dist.cpp` | `data.frame(rsid…z, pval, info, type)` | homogeneous-cohort imputation. **Masks `stats::dist` — see R-3** |
| `distmix(..., pop_wgt_df, ...)` | `src/distmix.cpp` | same, with `af1mix` | mixed-cohort imputation |
| `jepeg(study_pop, input_file, annotation_file, ref×3, af1_cutoff)` | `src/jepeg.cpp` | `data.frame(geneid, chisq, df, jepeg_pval, …)` | gene-level TWAS |
| `jepegmix(pop_wgt_df, …)` | `src/jepegmix.cpp` | same | mixed-cohort TWAS |
| `qcat(...)`, `qcatmix(...)` | `src/qcat.cpp`, `src/qcatmix.cpp` | `data.frame(…, qcat_m, qcat_t, qcat_chisq, qcat_pval)` | conditional association test; no vignette |
| `fiqt(z, min.p=1e-300)` | `R/fiqt.R` | numeric vector | pure R; no panel needed |

`ref×3` is shorthand for the three panel arguments that appear in every signature:
`reference_index_file`, `reference_data_file`, `reference_pop_desc_file`.

### 3.2 Preparation helpers (exported, documented, but really internal)

`prep_qcat()`, `prep_recessive_impute()` — return the raw ingredients (Z-vector, LD blocks) so an
analysis can be prototyped in R rather than C++. `prep_recessive_impute()` additionally returns
dominant- and recessive-coded correlation matrices.

### 3.3 Superseded experiments (exported — should not be)

`cpw2()` — an `afmix()` variant applying an arcsine-square-root variance-stabilising transform to
allele frequencies before the regression. Returns `(pop, wgt)`, a different schema from `afmix()`.

`prep_zmix()`, `prep_zmix2()`, `prep_zmix3()`, `prep_zmix4()` — successive attempts at SNP-pair
selection for the zmix regression, all superseded by `prep_zmix5()`/`prep_zmix5_sup()`, which are
what `R/zmix.R` actually calls. They differ only in which SNP pairs enter the design matrix:

| Variant | Pairing rule |
|---|---|
| `prep_zmix5` | all pairs among the ancestry-informative subset |
| `prep_zmix5_sup` | as above, correlations pooled to superpopulation level |
| `prep_zmix4` | fixed `offset`: SNP *i* paired with SNP *i+offset*, striped by `interval` |
| `prep_zmix3` | SNP *i* paired with the next `steps` SNPs |
| `prep_zmix2` | offset pairing, loading genotypes one SNP at a time (see M-3) |
| `prep_zmix` | original all-pairs form |

> Their `@return` documentation is wrong across the board — they return a numeric **design
> matrix**, not a weights data frame. See DOC-1.

---

## 4. C++/Rcpp components

### 4.1 `Arguments` (`src/gauss.h`)

A single mutable configuration struct threaded through every function by reference. It carries
user parameters (`chr`, `start_bp`, `end_bp`, `wing_size`, file paths, `af1_cutoff`), panel
metadata loaded at startup (`ref_pop_vec`, `ref_pop_size_vec`, `ref_sup_pop_vec`, `num_pops`), the
resolved population selection (`pop_flag_vec`, `pop_wgt_vec`, `pop_wgt_map`), and tuning constants
set in the constructor:

| Field | Default | Meaning |
|---|---|---|
| `lambda` | 0.1 | ridge term added to LD-matrix diagonals |
| `min_abs_eig` | 1e-5 | eigenvalue floor in `MakePosDef` |
| `eig_cutoff` | 0.01 | PC-retention threshold in `RmvPC`/`CountPC` |
| `min_num_measured_snp` / `min_num_unmeasured_snp` | 10 | minimum window occupancy |
| `total_num_categ` | 6 | JEPEG functional categories |
| `categ_cor_cutoff` | 0.8 | JEPEG collinearity pruning threshold |
| `denorm_norm_w` | 3 | JEPEG low-variance pruning denominator |

> The constructor does **not** initialise `chr`, `start_bp`, `end_bp`, `wing_size`, `af1_cutoff`,
> `study_pop`, `num_pops` or `num_samples`. Entry points that need them set them; entry points
> that do not (`afmix`, `cpw2`, `zmix`) leave them indeterminate. Safe today only because the
> code paths that read them are never reached from those functions.

### 4.2 `Snp` (`src/snp.{h,cpp}`)

Plain container, one per variant, always heap-allocated with bare `new` and freed by hand-written
loops at the end of each entry point (see M-1). Key state:

* Identity — `rsid_`, `chr_`, `bp_`, `a1_`, `a2_`
* Statistics — `z_`, `info_`, `af1study_` (from user input), `af1ref_` (panel, single-population),
  `af1mix_` (panel, weighted mixture)
* Location — `fpos_`: **the BGZF virtual offset of this SNP's genotype record**, the linchpin of
  the whole random-access design. Defaults to `-1` (see C-5)
* `type_` — the central three-state classifier:

  | `type_` | Meaning | Used for |
  |---|---|---|
  | `0` | in panel, **not** measured in the study | imputation targets |
  | `1` | in panel **and** measured | predictors; LD matrix rows |
  | `2` | measured but **absent** from the panel | passed through to output, excluded from all maths |

* JEPEG payload — `geneid_`, `categ_map_` (category → weight)
* `genotype_vec_` — one ASCII string per selected population, loaded lazily and freed eagerly

### 4.3 `Gene` / `Categ` (`src/gene.{h,cpp}`)

Encapsulates one gene's JEPEG test. `RunJepeg()` / `RunJepegmix()` tally how many SNPs carry each
of the six functional categories, build a `Categ` for each non-empty one, then delegate to
`CalJepegPval()` / `CalJepegmixPval()`. The two paths are near-identical; they differ only in how
the genotype correlation matrix `CorG` is built (`CalCor` vs. `CalWgtCov`-derived).

### 4.4 `util.{h,cpp}` — the numeric kernels

Two families:

**Genotype correlation over ASCII strings.** These are the hot loops of the package. Each walks
the digit characters accumulating `Σx, Σx², Σy, Σy², Σxy`:

* `CalCor(vector<string>&, vector<string>&)` — pools all selected populations; used by the
  non-`mix` functions
* `CalCor(string&, string&)` — single population; used by the `prep_zmix*` family
* `CalWgtCov(vector<string>&, vector<string>&, vector<double>& w)` — the mixture covariance;
  used by every `*mix` function. **This function is the package's most consequential defect —
  see C-1**
* `CalCorSup(...)` — superpopulation pooling for `prep_zmix5_sup`

**Eigen wrappers.** `MakePosDef` (eigenvalue clipping), `InvMat` (`fullPivLu().inverse()`),
`CholeskyMat`, `CnvrtCovToCor`, `RmvPC`/`CountPC` (PC retention), plus `MpMatMat`/`MpNumMat`
thin aliases retained from the GSL era. Several are dead (see D-2).

### 4.5 `bgzf.{c,h}`, `khash.h` — vendored

A decade-old samtools BGZF reader, providing `bgzf_open`, `bgzf_seek`, `bgzf_getc`, `bgzf_close`.
`util.cpp::BgzfGetLine()` builds lines from it character by character. Errors are recorded into
`fp->error` but no caller reads it, and the `fprintf(stderr)` diagnostics were commented out to
silence an `R CMD check` warning (see P-3).

---

## 5. Data flow

### 5.1 The common pipeline

Nearly every entry point follows the same nine phases. Understanding this sequence is most of
understanding the codebase.

```
 1. Parse arguments into `Arguments`; resolve af1_cutoff default (0.01)

 2. read_ref_desc(args)
      33kg_pop_desc.txt -> ref_pop_vec / ref_pop_size_vec / ref_sup_pop_vec, num_pops

 3. Population selection — one of:
      init_pop_flag_vec(args)       [dist, qcat, jepeg]
        match `study_pop` against population OR superpopulation names -> 0/1 flags
      init_pop_flag_wgt_vec(args)   [*mix, computeLD, simulateLD]
        match pop_wgt_df names against ref_pop_vec -> flags + parallel weight vector

 4. ReadInputZ(snp_map, args, All)  or  ReadInputAf(...)
      user file -> map<MapKey, Snp*>, every SNP tagged type = 2
      All=false: filter to [start_bp - wing_size, end_bp + wing_size]   (dist, distmix, qcat, …)
      All=true : keep everything, genome-wide                            (jepeg, jepegmix)

 5. ReadReferenceIndex(snp_map, args)      -- adds panel-only SNPs as type 0
    ReadReferenceIndexAll(snp_map, args)   -- updates existing entries only
      For each index line, look up BOTH allele orientations:
        found as (a1,a2)  -> type 1, record fpos
        found as (a2,a1)  -> type 1, record fpos, FLIP Z (z := -z), swap alleles, re-key the map
        found as neither  -> type 0 (ReadReferenceIndex) / ignored (…All)
        found as both     -> error: "input file contains duplicates"

 6. [jepeg only] ReadAnnotation(snp_map, args)
      attach geneid + (category, weight) pairs; same two-orientation matching

 7. MakeSnpVec(snp_vec, snp_map, args)      -- af1ref from pooled flagged populations
    MakeSnpVecMix(snp_vec, snp_map, args)   -- af1mix = Σ w_k · af1_k
      seek fpos, read the record, compute the frequency, keep if
        af1_cutoff < af1 < 1 - af1_cutoff

 8. ReadGenotype(snp_vec, args)
      seek fpos, split the record, retain only genotype strings where pop_flag_vec[k] == 1

 9. run_<analysis>(snp_vec, args)  ->  build R output  ->  FreeGenotype + delete every Snp*
```

**Allele orientation** is settled once, in phase 5, by rewriting the *study* side to match the
panel (flipping Z, swapping `a1`/`a2`). This is why the genotype-flipping code in `ReadGenotype`
is commented out and `Snp::flip_` is vestigial (D-2).

**Window vocabulary** used throughout `dist`/`qcat`/`prep_*`:

```
        wing_size            prediction window            wing_size
    |<------------>|<---------------------------->|<------------>|
    |                    extended window                         |
    start_bp - wing                                    end_bp + wing
```
Predictors (`type 1`) are drawn from the *extended* window; imputation targets (`type 0`) only
from the *prediction* window. Results are reported for the prediction window alone.

### 5.2 `zmix()` — the one flow that crosses the R/C++ boundary twice

```
R: zmix(level=)
     |
     +-> C++ prep_zmix5 / prep_zmix5_sup
     |     read input + index (private copies of the phase-4/5 parsers)
     |     keep type-1 SNPs, decimate by `interval`
     |     cal_af_norm_var(): normalised variance  var(af1) / (mean·(1-mean))  per SNP
     |     stats::quantile(..., percentile)  <- called back into R from C++
     |     keep SNPs above the cutoff = "ancestry-informative markers"
     |     for every pair (i,j):  row = [ z_i·z_j , r_1(i,j), …, r_p(i,j) ]
     |     return that design matrix                              <- O(n²) rows, see PERF-3
     |
     +-> R: drop non-finite rows; y = column 1, X = the rest
           label X's columns from the pop_desc file
           quadprog::solve.QP(Dmat = X'X, dvec = X'y,
                              subject to  Σw = 1,  0 ≤ w ≤ 1)
           round to 5 dp, renormalise, return
```

The statistical idea: `E[z_i z_j] ≈ Σ_k w_k · r_k(i,j)`, so the mixing proportions are the
non-negative, sum-to-one least-squares coefficients regressing observed Z-score products on
per-population LD. `afmix()` solves the analogous problem using allele frequencies instead of
Z-score products, in C++ via `CxxInv · Cxy` rather than a constrained solver.

---

## 6. Important file formats

All are **whitespace-delimited text**. The two panel data files are BGZF-compressed
(`.gz`, block-gzip — *not* plain gzip; random access depends on this).

### 6.1 Study inputs

**Z-score file** — header line required and skipped:
```
rsid chr bp a1 a2 z
rs1000427 22 36890105 A G -1.49697409010693
```

**Allele-frequency file** (`afmix`, `cpw2`) — identical but the 6th column is `af1`:
```
rsid chr bp a1 a2 af1
rs1000427 22 36890105 A G 0.11598
```

**SNP annotation file** (`jepeg`, `jepegmix`) — header required:
```
rsid chr bp a1 a2 geneid categ wgt
. 7 20180286 T C 7A5 WTH_TARGET 2.2
```
`categ` must be exactly one of `PROTEIN`, `TFBS`, `WTH_HAIR`, `WTH_TARGET`, `CIS_EQTL`,
`TRANS_EQTL`, mapping to internal category numbers 0–5 and output labels `PFS`, `TFB`, `STR`,
`TAR`, `CIS`, `TRN`. **Any other value silently corrupts the weight matrix — see C-3.**
`rsid` may be `.`; matching is by `(chr, bp, a1, a2)`, never by rsid.

**Mixing-proportion data frame** (`pop_wgt_df`) — read **positionally**: column 1 = population
abbreviation (upper-cased on ingest), column 2 = weight. The canonical shape is
`data(PGC2_SCZ_ANC_Prop)`: a 2-column `(pop, wgt)` frame. Note `afmix()` returns *three* columns,
so its output must have the first dropped before use (see C-7).

### 6.2 Reference panel (33KG)

**`33kg_pop_desc.txt`** — tab-delimited, header required, 29 data rows:
```
Population_Abbreviation  Number_of_Subject  Super_Population  Population_Description
ACB                      164                AFR               African Caribbeans in Barbados
```
`Number_of_Subject` must equal the corresponding genotype string length; nothing verifies this
(see V-5). Row order defines column order everywhere downstream — `R/zmix.R` re-reads this file
independently and relies on matching order.

**`33kg_index.gz`** — BGZF, **no header**, one line per panel variant:
```
rsid chr bp a1 a2 af1ref fpos
rs140052487 1 54353 A C 0.00173 33217
```
`fpos` is the BGZF *virtual* offset (`coffset << 16 | uoffset`) of this variant's record in the
genotype file. `af1ref` is parsed and then discarded — recomputing it from the genotype record is
the package's dominant I/O cost (see PERF-1).

**`33kg_geno.gz`** — BGZF, no header. Each record is `2 × num_pops` fields — for 33KG, 58:
```
<geno_1> <geno_2> … <geno_29>   <af1_1> <af1_2> … <af1_29>
```
Each `geno_k` is a string of `Number_of_Subject[k]` digit characters in `{0,1,2}` — the count of
the reference allele for each individual in population *k*. For 33KG a single record is roughly
33 KB of genotype text, and reading one is the unit cost of nearly every panel access.

---

## 7. External libraries

| Dependency | Kind | Used for |
|---|---|---|
| **Rcpp** (≥ 1.0.7) | `Imports` + `LinkingTo` | R↔C++ glue, `NumericMatrix`/`DataFrame`, `Rcpp::stop` |
| **RcppEigen** | `Imports` + `LinkingTo` | all dense linear algebra (`MatrixXd`, `LLT`, `SelfAdjointEigenSolver`, `FullPivLU`) |
| **quadprog** | `Imports` | `solve.QP` — the constrained solver behind `zmix()` |
| **R internals** | — | `R::pnorm5`, `R::pchisq` from `<Rmath.h>`; `stats::quantile` called from C++ via `Environment::namespace_env` |
| **zlib** | linked via R | required by the vendored BGZF reader |
| **BGZF / khash** | vendored source | `src/bgzf.{c,h}`, `src/khash.h`, from samtools |
| knitr, rmarkdown | `Suggests` | vignettes |

**GSL was removed** in commit `e9a77a2`; every `gsl_matrix` operation was ported to Eigen. Traces
remain in commented-out blocks and in `src/etc/RcppGSL.cpp` (D-2).

Undeclared but used at runtime: **`stats`** (`pnorm`, `qnorm`, `p.adjust`, and the C++ callback to
`quantile`) and **`utils`** (`read.table`) — see R-2. The vignettes additionally need `tidyverse`,
`data.table`, `kableExtra`, `dplyr`, `reshape2`, none of which are declared (R-5).

---

## 8. Major algorithms

### 8.1 DIST / DISTMIX — Gaussian imputation of Z-scores

Treats the vector of Z-scores as multivariate normal with covariance equal to the LD matrix.
For an unmeasured SNP with LD vector `b₂₁` to the measured set:

```
B₁₁        = LD among measured SNPs, with `lambda` added to the diagonal (ridge)
             -> MakePosDef (clip eigenvalues at min_abs_eig) -> InvMat
ẑ₂         = b₂₁ B₁₁⁻¹ z₁
info       = |b₂₁ B₁₁⁻¹ b₁₂|          (imputation quality, reported as `info`)
reported z = ẑ₂ / sqrt(info)           (variance-normalised)
```

`dist` builds `B₁₁` with `CalCor` over pooled genotypes of one population. `distmix` builds it as
`CalWgtCov(i,j) / (σᵢσⱼ)` with `σ = sqrt(CalWgtCov(i,i))` — the mixture path, and the one affected
by C-1.

Cost: `O(m²)` correlations plus two `O(m³)` factorisations, where `m` is the measured count in the
extended window.

### 8.2 QCAT / QCATMIX — conditional association testing

Decorrelates the Z-vector by a Cholesky factor rather than inverting:

```
num_eig = CountPC(B₁₁, eig_cutoff)      # effective rank
L       = chol(B₁₁)                     # NOTE: no MakePosDef first, and info() unchecked (N-1)
r       = cor( L⁻¹z₁ , L⁻¹b )           # b = the SNP's LD row
T       = sqrt(num_eig - 3) · r         # Fisher-style transform (NaN if num_eig < 3, see N-2)
chisq   = (num_eig - 3) · r²
```

### 8.3 JEPEG / JEPEGMIX — gene-level TWAS

Aggregates SNP Z-scores within a gene into per-functional-category statistics, then combines them:

```
W      = k × n weight matrix; W[categ, snp] = annotation weight × sqrt(info)
U      = W z                      # one score per functional category
CovU   = W · CorG · Wᵀ            # CorG = SNP-genotype LD within the gene
normU  = U / sqrt(diag(CovU))     # per-category Z, converted to a p-value

prune categories that are
  (a) collinear:      |CorU(i,j)| > categ_cor_cutoff      (C-10)
  (b) low-variance:   CovU(i,i) < WWt(i,i) / denorm_norm_w

X, CovX = U, CovU restricted to surviving categories
chisq   = Xᵀ CovX⁻¹ X              df = number of survivors
p       = pchisq(chisq, df, lower = FALSE)
```

Genes are formed by sorting `gene_snp_vec` by `geneid` and slicing contiguous runs
(`MakeGeneStartEndVec`). Genes whose categories are all pruned yield a placeholder row (C-11).

### 8.4 AFMIX — ancestry from allele frequencies

Splits SNPs into `interval` interleaved subsets (SNP *i* joins subset `i mod interval`), which
decorrelates the LD structure within each subset. For each subset it forms the covariance matrix
of `[af1_study, af1_pop1, …, af1_popP]`, then solves an unconstrained least-squares system:

```
W_i = Cxx⁻¹ · Cxy          # Cxx = pop×pop block, Cxy = pop×study block
W  += W_i / interval        # average across subsets
```
Negative weights are clipped to zero and the result rounded to 3 dp. **The result is not
normalised to sum to 1** (V-2). `cpw2()` is the same algorithm with an arcsine-square-root
transform applied to every frequency first.

### 8.5 ZMIX — ancestry from Z-scores alone

See §5.2. Two ideas make it work: (a) ancestry-informative markers are selected by normalised
allele-frequency variance across populations, `var(af1) / (mean·(1−mean))`, keeping the upper tail;
(b) the mixing weights are recovered from `E[z_i z_j] ≈ Σ_k w_k r_k(i,j)` by a **constrained**
quadratic program (`Σw = 1`, `w ≥ 0`), which is what distinguishes it from `afmix`'s unconstrained
solve-then-clip.

### 8.6 simulateLD — resampling-based LD

Draws `w_k · sim_size` pseudo-individuals with replacement from each population's genotype column,
assembles a synthetic genotype matrix, and takes its empirical correlation. Intended as a
validation cross-check on `computeLD()`. Currently affected by C-6 (unfilled columns) and P-1
(non-reproducible RNG).

### 8.7 FIQT — Winner's Curse correction (pure R)

FDR inverse-quantile transformation: convert Z to two-sided p-values, apply Benjamini–Hochberg,
map the adjusted p-values back through `qnorm`, and restore the sign. Z-scores beyond the
`min.p` boundary are passed through unshrunk.

### 8.8 prep_recessive_impute — non-additive coding

Recodes prediction-window genotypes as dominant (`0→0`, `1,2→1`) and recessive (`2→1`, else `0`),
then returns additive/dominant/recessive correlation matrices against the additive measured set —
letting a caller impute non-additive effects in R. Calls `UpdateSnpToMinorAllele` first so the
minor allele is always the reference (recessive coding is not orientation-symmetric).

---

## 9. Testing structure

**There is no automated test suite.** No `tests/` directory, no `testthat`, no CI workflow.

What exists instead:

* **`inst/dev-tests/`** — 27 hand-run driver scripts (`test_dist_33kg.R`, `test_zmix_pgc.R`,
  `test_computeLD_33kg.R`, …). They hard-code absolute paths, require the multi-GB 33KG panel,
  produce output for visual inspection rather than assertions, and are `.gitignore`d — so they are
  absent from a fresh clone.
* **Vignettes** — the closest thing to end-to-end coverage, but they read `../data/` and
  `../ref/33KG/` and depend on undeclared packages, so `R CMD build` cannot execute them (R-5).
* **Figure scripts** — `gauss_figure1_manhattan_plot.R` and
  `gauss_figure_S1_computeLD_comp_simulateLD_33kg.R` reproduce publication figures; the latter is
  effectively a manual `computeLD` vs `simulateLD` consistency check.

**Consequence for contributors:** no change to numeric code is currently verifiable. Several
review findings (C-1, C-3, C-5, C-6) alter results without altering behaviour, and nothing in the
repository would notice. See T-1 in the review for a proposed fixture-based suite — the key
enabler is a small synthetic panel (3 populations × 20 samples × 200 SNPs, BGZF-compressed,
<200 KB) that makes the whole pipeline runnable in seconds.

---

## 10. Generated files

Do not hand-edit these; regenerate them.

| Path | Generated by | Trigger |
|---|---|---|
| `src/RcppExports.cpp`, `R/RcppExports.R` | `Rcpp::compileAttributes()` | any change to an `// [[Rcpp::export]]` signature |
| `man/*.Rd`, `NAMESPACE` | `roxygen2::roxygenise()` / `devtools::document()` | any change to a roxygen block |
| `README.md` | `knitr::knit("README.Rmd")` / `devtools::build_readme()` | edits to `README.Rmd` |
| `docs/**` (except this file and `REVIEW_PLAN.md`) | `pkgdown::build_site()` | releases; config in `_pkgdown.yml` |
| `src/Makevars`, `src/Makevars.win` | `configure` / `configure.win` at install time | — |
| `src/*.o`, `src/gauss.so` | compiler | — |

Two consequences worth internalising:

1. **Roxygen for the C++ functions lives in the `.cpp` files**, as `//'` comments above the
   `[[Rcpp::export]]` tag. `compileAttributes()` copies them into `R/RcppExports.R`, and only then
   does roxygen see them. Editing a C++ function's documentation therefore requires
   `compileAttributes()` *before* `document()`, or the `.Rd` will not change.
2. `docs/` is `.Rbuildignore`d, so this file and `REVIEW_PLAN.md` ship with the repository but not
   with the installed package. A `pkgdown::build_site()` run will not delete them (pkgdown leaves
   unknown files alone), but a `clean = TRUE` rebuild would.

Untracked-but-present locally: `.Rproj.user/`, `.Rhistory`, `src/*.o`, `src/gauss.so`, `ref/`,
`inst/dev-tests/`. Tracked but shouldn't be: seven `.DS_Store` files (R-8).

---

## 11. Areas that require extra caution

### 11.1 `CalWgtCov()` is on the critical path of every mixture analysis

`src/util.cpp:100`. `distmix`, `computeLD`, `jepegmix` and `prep_recessive_impute` all reduce to
this one function. It currently returns a covariance inflated by `m²` (C-1), which re-weights
populations by sample size and numerically erases the between-population term. **Fixing it changes
every published mixture number**, so it needs the test fixture in place first, plus a deliberate
decision about re-validating prior results.

### 11.2 `fpos` and the `type` classifier are load-bearing

The random-access design rests on `fpos` being a valid BGZF virtual offset. Type-2 SNPs never get
one (`fpos_ == -1`), and no `bgzf_seek` return value is checked anywhere (C-5, P-3). A seek that
fails leaves the stream on an arbitrary record and the code reads it as if nothing happened. Any
new panel access must (a) skip type-2 SNPs and (b) check the seek.

Equally: the `type` values 0/1/2 are compared as bare integers in ~20 places across 13 files. There
is no enum and no single definition. Changing their meaning requires a full sweep.

### 11.3 Manual memory management with throwing error paths

`Snp` objects are `new`ed in the parsers and freed only by an explicit loop at the bottom of each
entry point. `Rcpp::stop()` throws past that loop, so **every error path leaks the entire SNP
map** — including the routine `"Not enough number of SNPs loaded"` message (M-1). When adding an
early return or a new validation check, either place it before allocation or convert the map to
`unique_ptr` ownership first.

### 11.4 Positional column access at the R/C++ boundary

`pop_wgt_df[0]` / `pop_wgt_df[1]` appears in five entry points. It silently misinterprets any data
frame whose first two columns are not `(pop, wgt)` — which is exactly what `afmix()` returns
(C-7). Adding a new consumer of `pop_wgt_df` propagates the trap.

### 11.5 Four copies of the same parser

`ReadInputZ`/`ReadReferenceIndex` in `gauss.cpp` have near-identical private twins in `zmix.cpp`
(`read_input_zmix`, `read_ref_index_zmix`) and `cpw2.cpp` (`read_input_cpw2`,
`read_ref_index_cpw2`). A fix to allele matching, duplicate detection or validation must be applied
to all of them. Likewise `dist`/`distmix`, `qcat`/`qcatmix` and `jepeg`/`jepegmix` are pairwise
copy-edits — a change to one almost always belongs in its twin.

### 11.6 Silent numerical degradation

Several failure modes produce plausible numbers rather than errors:

* `CholeskyMat` ignores `LLT::info()`, returning garbage for indefinite matrices (N-1)
* `MakePosDef` and `RmvPC` silently no-op when the eigensolver fails, usually because `NaN`
  already entered (N-5)
* every correlation is an unguarded division; a monomorphic SNP gives `0/0` and the `NaN`
  propagates through the whole LD matrix (N-3)
* recessive coding in `prep_recessive_impute` makes zero-variance columns common, not exceptional

When touching the linear algebra, prefer failing loudly.

### 11.7 Unbounded allocation from user parameters

`prep_zmix5` computes `int total_rows = n*(n-1)/2` — signed overflow above n ≈ 65,536, and ~12 GB
of R memory at n = 10,000, which the `R/zmix.R` defaults can reach on a genome-wide input
(PERF-3). `interval = 0` sends the SNP-selection loop into an infinite `push_back` (V-3). Validate
new numeric parameters at the boundary.

### 11.8 Regenerate, don't hand-edit

See §10. In particular, editing a `//'` roxygen block in a `.cpp` file has no effect on the
`.Rd` until `Rcpp::compileAttributes()` runs.

### 11.9 The reference panel is not in the repository

`ref/` is `.gitignore`d and multi-gigabyte. Nothing in a fresh clone can exercise the panel code
paths. Any local verification requires downloading 33KG separately — which is precisely why the
synthetic fixture panel (§9) is the highest-leverage infrastructure addition available.
