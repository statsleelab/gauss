# GAUSS — Comprehensive Code Review

**Repository:** `statsleelab/gauss` · **Version reviewed:** 1.0 (`main` @ `7622f7d`)
**Date:** 2026-09-05 · **Scope:** whole package (R, `src/`, build system, docs, packaging)
**Status:** Review only — no code was modified.

---

## 1. Project architecture and purpose

GAUSS (*Genome Analysis Using Summary Statistics*) is an R package for analysing GWAS
summary statistics from cohorts of mixed ancestry, backed by a large external reference
panel (33KG: 32,953 samples / 29 populations). Almost all computation is C++ behind Rcpp;
the R layer is a thin veneer.

### 1.1 Layers

| Layer | Files | Role |
|---|---|---|
| R API | `R/zmix.R`, `R/fiqt.R`, `R/RcppExports.R` | `zmix()` (prep + quadratic program), `fiqt()` (pure R), generated `.Call` wrappers |
| Rcpp entry points | `src/{dist,distmix,computeLD,simulateLD,afmix,cpw2,zmix,qcat,qcatmix,jepeg,jepegmix,prep_qcat,prep_qcatmix}.cpp` | One exported function per analysis |
| Shared engine | `src/gauss.cpp` / `gauss.h` | `Arguments` config struct, all reference-panel and input I/O, SNP-set construction |
| Domain objects | `src/snp.{h,cpp}`, `src/gene.{h,cpp}` | `Snp` (per-variant state), `Gene` (JEPEG gene-level test) |
| Numerics | `src/util.{h,cpp}` | Correlation/covariance kernels over genotype strings, Eigen matrix helpers |
| Vendored I/O | `src/bgzf.{c,h}`, `src/khash.h` | BGZF random-access reader (from samtools) |

### 1.2 Central data flow (shared by nearly every entry point)

```
read_ref_desc(args)                 # parse 33kg_pop_desc.txt -> ref_pop_vec / sizes / superpops
init_pop_flag_vec | init_pop_flag_wgt_vec   # which reference populations are in play, and weights
ReadInputZ / ReadInputAf            # user summary stats -> map<MapKey, Snp*>  (type 2)
ReadReferenceIndex[All]             # BGZF index -> attach rsid + byte offset fpos; align alleles
                                    #   (type 1 = measured & in panel, type 0 = panel-only)
MakeSnpVec | MakeSnpVecMix          # seek fpos, compute af1ref / af1mix, apply MAF cutoff
ReadGenotype                        # seek fpos, load per-population genotype strings ("012...")
run_<analysis>                      # Eigen linear algebra over LD matrices
FreeGenotype + manual delete of every Snp*
```

The key design decisions are (a) genotypes are kept as ASCII digit strings, one string per
population, and correlations are computed by streaming over those characters; (b) the
reference panel is not loaded — each SNP's genotype record is fetched by `bgzf_seek` to a
byte offset stored in the index; (c) all `Snp` objects are raw `new`/`delete` with manual
cleanup at the end of each entry point.

### 1.3 Notable characteristics

* **No tests.** There is no `tests/` directory and no `testthat` infrastructure. `inst/dev-tests/`
  holds ~27 ad-hoc driver scripts that require the multi-GB 33KG panel and are not runnable in CI.
* **Heavy duplication.** `afmix.cpp`, `cpw2.cpp` and `zmix.cpp` each carry their own private
  near-identical copies of `read_input_*` and `read_ref_index_*`; `dist`/`distmix`,
  `qcat`/`qcatmix`, `jepeg`/`jepegmix` are pairwise copy-edits.
* **Large experimental surface.** 20 functions are exported, of which `cpw2`, `prep_zmix`,
  `prep_zmix2`–`prep_zmix5`, `prep_zmix5_sup`, `prep_qcat` and `prep_recessive_impute` are
  internal/experimental scaffolding rather than user API.

---

## 2. Findings

Severity scale: **Critical** (wrong scientific results or crashes) · **High** (wrong results in
plausible cases, or blocks distribution) · **Medium** (robustness, performance, maintainability) ·
**Low** (hygiene).

Difficulty: **S** ≤ 1 h · **M** ≈ half a day · **L** ≈ multiple days.

---

### 2.1 Correctness bugs

#### C-1 · `CalWgtCov()` returns a covariance inflated by *m²*, silently destroying the mixture correction — **Critical** — difficulty **M**
**File:** [src/util.cpp:100-124](../src/util.cpp) (`CalWgtCov`)
**Callers:** `distmix`, `computeLD`, `jepegmix` (`gene.cpp:571-586`), `prep_recessive_impute`

The mixture covariance is documented (and intended) to be
`Σ_k w_k·Cov_k(X,Y) + Σ_k w_k·μ_xk·μ_yk − (Σ_k w_k·μ_xk)(Σ_k w_k·μ_yk)`.
The within-population term is computed as

```cpp
double factor = ((double)m) / (m - 1);
wsumcov += wgt_val*factor*(m*sumxy-sumx*sumy);   // line 117-118
```

`Cov_k = (m·Σxy − Σx·Σy) / (m·(m−1))`, so the code produces `w_k · m² · Cov_k` — a factor of
`m²` too large, where `m` is that population's sample size. Two consequences:

1. Populations are effectively re-weighted by `w_k·m_k²` instead of `w_k`. In the 33KG panel
   population sizes span roughly two orders of magnitude, so the user-supplied mixing
   proportions are largely overridden by panel sample sizes.
2. The between-population mean terms (`wsum_mi_mj − wsum_mi*wsum_mj`), which are *the* reason
   mixture LD differs from within-population LD, are ~4 orders of magnitude smaller than the
   inflated first term and are numerically discarded. DISTMIX/computeLD LD therefore behaves
   like a size-weighted within-population average, not a mixture.

The `cov/(std_i·std_j)` normalisation cancels only a *uniform* scale factor; because `m_k`
varies across populations, it does not cancel here.

**Consequence:** ancestry-informed LD, `distmix()` imputed Z-scores, `jepegmix()` gene p-values
and `prep_recessive_impute()` correlations are systematically biased — the headline claim of the
package. Note the pre-existing code used integer `m/(m-1)` (== 1), so this is a long-standing
bug that commit `e9a77a2` altered but did not fix.

**Fix:** divide by `m*(m-1)`:
```cpp
wsumcov += wgt_val * (m*sumxy - sumx*sumy) / (static_cast<double>(m) * (m - 1));
```
Guard `m < 2` (skip the population or return `NaN`). Then re-validate `computeLD()` against a
direct pooled-genotype LD computation on a small locus before releasing.

---

#### C-2 · `exit(EXIT_FAILURE)` kills the user's R session — **Critical** — difficulty **S**
**File:** [src/gauss.cpp:1289-1293](../src/gauss.cpp) (`ReadAnnotation`)

```cpp
if(!in_annotation){
  Rcpp::Rcout << "ERROR: can't open snp annotation data file '"<<annotation_file<<"'"<<std::endl;
  exit(EXIT_FAILURE);
}
```

**Consequence:** a mistyped annotation path in `jepeg()`/`jepegmix()` terminates the entire R
process — unsaved work is lost, and this is an outright CRAN policy violation.
**Fix:** `Rcpp::stop("ERROR: can't open snp annotation data file '" + annotation_file + "'");`
(matching every other I/O error path in the file).

---

#### C-3 · Uninitialised `categ_num` used for unknown annotation categories — **Critical** — difficulty **S**
**File:** [src/gauss.cpp:1298, 1318-1332](../src/gauss.cpp) (`ReadAnnotation`)
Confirmed by the compiler: `gauss.cpp:1336: warning: variable 'categ_num' may be uninitialized when used here [-Wconditional-uninitialized]`.

`categ_num` is declared without an initialiser outside the loop and assigned only inside an
`if/else if` chain covering six known category names. Any other value in the `categ` column
(typo, new category, extra whitespace, a truncated line) leaves `categ_num` holding the previous
row's value — or, on the first row, an indeterminate value.

**Consequence:** SNP weights are silently filed under the wrong functional category (or an
out-of-range one), corrupting the JEPEG weight matrix `W` and the resulting gene p-values.
Reading indeterminate memory is undefined behaviour.
**Fix:** initialise `int categ_num = -1;`, and `continue` (with an accumulated warning count) when
the category is unrecognised.

---

#### C-4 · `init_pop_flag_vec()` reads out of bounds when a name is both a population and a superpopulation — **High** — difficulty **S**
**File:** [src/gauss.cpp:1017-1063](../src/gauss.cpp)

`pop_vec` is assigned only in the two mutually exclusive branches
(`in_pop && !in_sup_pop`) and (`!in_pop && in_sup_pop`); the error branch covers
(`!in_pop && !in_sup_pop`). The fourth case — the study population name appearing in **both**
`ref_pop_vec` and `ref_sup_pop_vec` — leaves `pop_vec` empty, and the loop at line 1056 then
indexes `pop_vec[i]` for `i` up to `num_pops`.

**Consequence:** out-of-bounds `std::vector` read (UB): crash, or a silently wrong population
selection for `dist()`, `qcat()` and `jepeg()`. A reference panel whose superpopulation labels
overlap population abbreviations triggers it.
**Fix:** add an explicit `else` that calls `Rcpp::stop()` with an "ambiguous population name"
message, and prefer `pop_vec.at(i)` in the loop.

---

#### C-5 · `MakeSnpVec()` / `MakeSnpVecMix()` seek to `fpos == -1` for input-only SNPs — **High** — difficulty **S**
**Files:** [src/gauss.cpp:546-606](../src/gauss.cpp), [src/gauss.cpp:634-692](../src/gauss.cpp)

Both functions iterate over **every** entry of `snp_map`, including type-2 SNPs (present in the
user's input but absent from the reference panel), whose `fpos_` is still the constructor default
`-1` (`snp.cpp:31`). `bgzf_seek(fp, -1, SEEK_SET)` fails; the return value is never checked
(true of all 8 `bgzf_seek` call sites), so `BgzfGetLine` reads whatever record the stream happens
to sit on. The parsed allele counts then produce a meaningless `af1ref`/`af1mix`, which
(a) decides via the MAF cutoff whether the SNP enters `snp_vec` at all, and (b) is reported to
the user as the `af1ref`/`af1mix` column.

**Consequence:** the `af1ref`/`af1mix` column of `dist()`, `distmix()`, `qcat()` and
`prep_qcat()` output contains garbage for every type-2 row, and inclusion of those rows is
non-deterministic. Also a wasted full-record read (~1 MB per SNP) — see P-1.
**Fix:** `if ((it_sm->second)->GetType() == 2) continue;` at the top of both loops, decide
explicitly whether type-2 SNPs should be passed through to the output, and check
`bgzf_seek(...) < 0` everywhere.

---

#### C-6 · `simulateLD()` leaves unfilled genotype columns and correlates over them — **High** — difficulty **M**
**File:** [src/simulateLD.cpp:132-200](../src/simulateLD.cpp)

`geno_mat` is allocated as `num_measured × sim_size`, but only
`Σ_k (int)(w_k · sim_size)` columns are ever written. Integer truncation means this sum is
almost always strictly less than `sim_size` (up to `num_pops − 1` short, and much more if the
weights do not sum to 1). The trailing columns stay at their zero-initialised value and are
included in every `CalCor()` call.

**Consequence:** simulated LD is attenuated toward zero by a block of identical zero
"individuals"; the bias grows as `sim_size` shrinks or the number of populations grows.
**Fix:** compute `n_total = Σ pop_num_sim` and allocate/correlate over `n_total` columns
(or distribute the truncation remainder via largest-remainder rounding).

---

#### C-7 · `afmix()` output cannot be fed to `distmix()` / `computeLD()` / `jepegmix()` / `simulateLD()` — **High** — difficulty **S**
**Files:** [src/afmix.cpp:100-110](../src/afmix.cpp) vs. [src/distmix.cpp:47-55](../src/distmix.cpp), `computeLD.cpp:41-49`, `jepegmix.cpp:36-44`, `simulateLD.cpp:54-62`

`afmix()` returns `data.frame(sup.pop, pop, wgt)`. Every consumer reads the weight table
**positionally**:
```cpp
std::vector<std::string> pop_vec_in     = as<std::vector<std::string>>(pop_wgt_df[0]);
std::vector<double>      pop_wgt_vec_in = as<std::vector<double>>(pop_wgt_df[1]);
```
so it takes `sup.pop` as the population ID and tries to coerce the `pop` character column to
`double`. The workflow in `README.md` presents these as directly composable; the only hint that
a column must be dropped is a commented-out `wgt.df[,-1]` in `vignettes/afmix_example.Rmd:97`.

**Consequence:** the documented end-to-end workflow errors out (or, worse, if a user reorders
columns, silently uses superpopulation labels that match nothing — see V-2). `cpw2()` returns
only `(pop, wgt)`, so the two ancestry estimators have incompatible output shapes.
**Fix:** select columns by name (`pop_wgt_df["pop"]`, `pop_wgt_df["wgt"]`) with a validated
fallback to positional access, and make `afmix()`, `cpw2()` and `zmix()` return a single
consistent schema.

---

#### C-8 · Duplicate input rows silently leak and overwrite — **Medium** — difficulty **S**
**Files:** `ReadInputZ` [src/gauss.cpp:172-175](../src/gauss.cpp), `ReadInputAf` (:254-256),
`read_input_zmix` (`zmix.cpp:1106-1108`), `read_input_cpw2` (`cpw2.cpp:~300`)

```cpp
snp = new Snp();  ...  snp_map[mkey] = snp;
```
If two input rows share `(chr, bp, a1, a2)`, the second assignment overwrites the pointer and the
first `Snp` is never deleted. There is no warning, and the surviving record is whichever came last.

**Consequence:** silent data loss plus a per-duplicate leak. Note `ReadReferenceIndex` *does*
detect duplicates (case 4) — but only ones visible from the reference side.
**Fix:** use `snp_map.emplace(mkey, snp)` / `insert`, `delete` the rejected object, and report the
duplicate count to the user.

---

#### C-9 · `GetTopCateg()` can return a removed category — **Low** — difficulty **S**
**File:** [src/gene.cpp:884-895](../src/gene.cpp)

`top_index` starts at `0` and is only updated when `top_p > p && !rmv`. If category 0 was flagged
for removal (collinear or low-variance) and has the smallest p-value, it is returned anyway.
**Consequence:** `top_categ` / `top_categ_pval` in `jepeg()`/`jepegmix()` output can name a
category that was excluded from the reported `chisq`/`df`.
**Fix:** initialise `top_index` to the first non-removed index (and return a sentinel if all are
removed).

---

#### C-10 · Collinearity pruning depends on already-removed categories — **Low** — difficulty **S**
**File:** [src/gene.cpp:390-398](../src/gene.cpp) (and the `Jepegmix` twin at :662)

The pruning loop marks `j` for removal when `|Cor(i,j)| > cutoff` for any `i < j`, without
checking whether `i` itself is already marked. A chain of correlated categories over-prunes,
lowering `df_` and making the gene test conservative.
**Fix:** skip `i` when `categ_vec_[i].GetRmv()` is true.

---

#### C-11 · Genes with `df_ == 0` emit a placeholder row — **Low** — difficulty **S**
**File:** [src/gene.cpp:440-548](../src/gene.cpp)

`geneid_`, `chisq_`, `jepeg_pval_`, `top_*` are all assigned **inside** `if(df_){...}`. When every
category is pruned, the row keeps the constructor defaults: `geneid = "."`, `chisq = -1`,
`jepeg_pval = -1`, `top_snp = "."`.
**Consequence:** the returned data frame contains unlabelled rows with impossible p-values that
will pollute any downstream `p.adjust()` / QQ plot.
**Fix:** set `geneid_` unconditionally, use `NA_REAL` instead of `-1`, or skip the row entirely.

---

### 2.2 Memory and resource management

#### M-1 · Every `Rcpp::stop()` on a hot error path leaks the whole SNP map — **High** — difficulty **M**
**Files:** `run_dist` ([src/dist.cpp:141-147](../src/dist.cpp)), `run_distmix` (`distmix.cpp:152-158`),
`run_qcat` (`qcat.cpp:157-162`), `qcatmix.cpp`, `computeLD.cpp:88-92`, `simulateLD.cpp:118-122`,
`prep_qcat.cpp:93-99`, `prep_qcatmix.cpp:141-143`, plus every `Rcpp::stop` inside
`MakeSnpVec*`/`ReadGenotype`/`BgzfGetLine`.

`Snp` objects are allocated with bare `new` and deleted only by a manual loop at the very end of
each entry point. `Rcpp::stop()` throws, unwinding past that loop.

The most common message in normal use — `"Not enough number of SNPs loaded"` — is exactly such a
path. On a genome-wide `jepeg()` map this is hundreds of MB per failed call, and it is
**not** reclaimed until the R session exits. The `BGZF*` handles opened in `MakeSnpVec*`,
`ReadGenotype`, `cal_af_norm_var` and `afmix_vec` leak the same way.

**Fix:** store `std::unique_ptr<Snp>` (or values) in the map, and wrap `BGZF*` in a small RAII
handle with a `bgzf_close` deleter. This removes ~15 hand-written cleanup loops at the same time.

---

#### M-2 · Allele-swap re-keying leaks the displaced `Snp` — **Medium** — difficulty **S**
**Files:** [src/gauss.cpp:355-362](../src/gauss.cpp) (`ReadReferenceIndex`), `gauss.cpp:487-494`
(`ReadReferenceIndexAll`), `gauss.cpp:1350-1354` (`ReadAnnotation`), `zmix.cpp:1163-1166`,
`cpw2.cpp` equivalent

```cpp
MapKey new_key(chr, bp, a1, a2);
snp_map[new_key] = it2->second;   // if new_key already exists, its Snp* is dropped
snp_map.erase(it2);
```
If the input file contains both allele orientations for the same locus, the previous occupant of
`new_key` is overwritten and leaked — and one of the two records disappears without warning.
**Fix:** check `find(new_key)` first; on collision, report the conflict and `delete` the loser.

---

#### M-3 · `ReadGenotypeOne` / `FreeGenotypeOne` reopen the panel file per SNP pair — **Medium** — difficulty **S**
**File:** [src/zmix.cpp:752-766](../src/zmix.cpp) (`prep_zmix2`)

Inside the SNP-pair loop, `ReadGenotypeOne` is called twice per pair, and each call does a full
`bgzf_open` … `bgzf_close` (`gauss.cpp:811-873`). For `n` SNPs that is `O(n²)` opens of a
multi-GB gzipped file.
**Consequence:** `prep_zmix2` is effectively unusable; also risks exhausting file descriptors.
**Fix:** hoist a single `BGZF*` out of the loop, or (better) drop `prep_zmix2` — see D-1.

---

#### M-4 · `.o` / `.so` build artefacts are present in `src/` — **Low** — difficulty **S**
36 MB of `src/*.o` and `src/gauss.so` sit in the working tree. They are `.gitignore`d (so not
committed) but they *are* picked up by a local `R CMD build`/`install` and can mask a stale
source. `src/etc/` (a second, older copy of `dist.cpp`, `jepeg.cpp`, `qcat.cpp`, `bgzf.c` and a
`gauss_main.cpp`) and `src/zmixrcpp.cpp_` are also present.
**Fix:** `make clean` before building; add `src/etc` and `*.cpp_` to `.Rbuildignore`; delete
`src/etc/` and `zmixrcpp.cpp_` from the repository.

---

### 2.3 Input validation

#### V-1 · No stream-state checking on any parsed line — **High** — difficulty **M**
**Files:** all `ReadInput*` / `ReadReferenceIndex*` / `ReadAnnotation` / `read_ref_desc`
(`gauss.cpp`, `zmix.cpp`, `cpw2.cpp`)

Every parser is `buffer >> rsid >> chr >> bp >> a1 >> a2 >> z;` with no `if (!buffer) …`. A blank
line, a comment line, a file with the wrong column count, a comma-separated file, or a header the
user did not intend to skip all yield a fully-formed but bogus `Snp` (`chr = 0`, `bp = 0`,
`a1 = a2 = ""`, `z = 0`) inserted into the map. There is also no check that `z` is finite, that
`chr` is in range, or that `bp > 0`.

**Consequence:** silently wrong analyses on malformed input, with no diagnostic. Combined with
C-8, repeated blank lines collide on the same key and leak.
**Fix:** validate `buffer` after extraction, count and report skipped lines, and validate the
header against the documented column names before reading.

---

#### V-2 · Population weights are accepted without any validation — **High** — difficulty **S**
**File:** [src/gauss.cpp:1097-1116](../src/gauss.cpp) (`init_pop_flag_wgt_vec`), and the
`pop_wgt_df` block duplicated in `distmix.cpp:47`, `computeLD.cpp:41`, `jepegmix.cpp:36`,
`simulateLD.cpp:54`, `prep_qcatmix.cpp:57`

`init_pop_flag_wgt_vec` walks `ref_pop_vec` and picks up whatever it finds in `pop_wgt_map`.
Nothing checks that:
* every user-supplied population actually matched a reference panel population (a typo means the
  population is silently dropped and the remaining weights no longer sum to 1);
* weights are non-negative;
* weights sum to 1;
* at least one population matched (an all-typo table gives an empty `pop_wgt_vec`, after which
  `CalWgtCov` reads `pop_wgt_vec[i]` out of bounds).

`afmix()` itself returns **un-normalised** weights (`afmix.cpp:203-211` computes `sum_pop_wgt`
but never divides by it), while `zmix()` normalises — so even the package's own outputs disagree.

**Consequence:** silently mis-weighted mixture LD; out-of-bounds read in the worst case.
**Fix:** in `init_pop_flag_wgt_vec`, error on any unmatched population name, error on
negative weights, warn-and-renormalise when `|Σw − 1| > tol`, and error when nothing matched.

---

#### V-3 · `interval` and `percentile` are unvalidated; `interval <= 0` hangs — **Medium** — difficulty **S**
**Files:** [src/zmix.cpp:64-68, 110-118](../src/zmix.cpp) (and the same pattern in `prep_zmix2/3/4/5_sup`),
[src/afmix.cpp:43-47, 137-148](../src/afmix.cpp), `cpw2.cpp:44-48`

```cpp
for(int i=0;;i++){
  int index = i*args.interval;
  if(index < measured_snp_vec.size()){ snp_vec.push_back(...); } else break;
}
```
With `interval = 0`, `index` is always 0 and the loop pushes the same pointer forever until the
process is killed by the OOM killer. Negative values are also accepted.

In `afmix_vec`/`cpw2_vec` the complementary failure occurs: `for(int i=0; i<interval; i++)` with
`interval > n_snps` yields empty `snp_subvec`s, so `af1_mat` is 0×(p+1); `CalCov` then takes the
mean of an empty vector (`0/0 = NaN`), poisoning `W_mat` for every subsequent iteration.

`percentile` is likewise passed straight to `stats::quantile` with no `[0,1]` check.

**Fix:** validate at the R/Rcpp boundary — `interval >= 1`, `interval <= n_snps`,
`0 < percentile < 1` — with clear `Rcpp::stop` messages.

---

#### V-4 · `fiqt()` errors on `NA` input — **Medium** — difficulty **S**
**File:** [R/fiqt.R:7-14](../R/fiqt.R)

`pvals[pvals < min.p] <- min.p` uses a logical index that is `NA` wherever `z` is `NA`, and R
raises `NAs are not allowed in subscripted assignments`. GWAS Z-score vectors routinely contain
missing values.
`min.p` is also unvalidated (`min.p <= 0` makes `qnorm(min.p/2, lower=FALSE)` infinite, so the
`abs(z) > …` guard never fires).
**Fix:** handle `NA` explicitly (`idx <- !is.na(pvals) & pvals < min.p`), return `NA` for `NA`
inputs, and use `p.adjust(..., n = sum(!is.na(pvals)))`. Validate `0 < min.p < 1`.

---

#### V-5 · `simulateLD()` indexes genotype strings with an unchecked random index — **Medium** — difficulty **S**
**File:** [src/simulateLD.cpp:143-150, 168-178](../src/simulateLD.cpp)

`ran_index` is drawn from `[0, ref_pop_size_vec[k] − 1]`, taken from the population description
*text file*, then used as `geno_str[ran_index]` on a genotype string from the *data file*. If the
two disagree by even one sample (an easy mistake when swapping panels), this is an out-of-bounds
`std::string` read.
**Fix:** validate `geno_str.length() == ref_pop_size_vec[k]` once per SNP and `Rcpp::stop` on
mismatch, or clamp using `geno_str.size()`.

---

### 2.4 Numerical issues

#### N-1 · `CholeskyMat()` never checks `LLT::info()` — **High** — difficulty **S**
**File:** [src/util.cpp:271-274](../src/util.cpp); callers [src/qcat.cpp:204](../src/qcat.cpp), `qcatmix.cpp:232`

```cpp
Eigen::LLT<Eigen::MatrixXd> llt(m);
result = llt.matrixL();
```
If `m` is not positive definite Eigen sets `info() == NumericalIssue` and `matrixL()` returns
partially-computed garbage. In `run_qcat`/`run_qcatmix` `B11` is built from empirical
correlations with only `+lambda` on the diagonal and — unlike `dist`/`distmix` — **no**
`MakePosDef` call, so indefiniteness is entirely plausible for a dense LD block.
**Consequence:** silent nonsense QCAT T-scores and χ² values.
**Fix:** check `llt.info() != Eigen::Success` and either call `MakePosDef` first or `Rcpp::stop`.

---

#### N-2 · `sqrt(num_eig - 3)` produces `NaN` for small windows — **Medium** — difficulty **S**
**File:** [src/qcat.cpp:227, 240](../src/qcat.cpp); `qcatmix.cpp:250, 265`

`num_eig` comes from `CountPC()` and is not bounded below. With `num_eig < 3` the T-score is
`NaN` and the χ² is negative, producing `NaN`/1.0 p-values with no warning.
**Fix:** require `num_eig >= 4` and `Rcpp::stop`/warn otherwise.

---

#### N-3 · Unguarded division by zero in the correlation kernels — **Medium** — difficulty **M**
**Files:** [src/util.cpp:49-70, 126-150, 152-168](../src/util.cpp) (`CalCor` overloads),
`zmix.cpp:1221-1246` (`CalCorSup`), all `cov/(stdi*stdj)` sites in `distmix.cpp`, `computeLD.cpp`,
`gene.cpp`, `prep_qcatmix.cpp`

Every correlation is `numer/denor` with no check that `denor > 0`. A monomorphic SNP within the
selected populations gives `0/0 = NaN`, which then propagates into `B11`, through `MakePosDef`
(where `SelfAdjointEigenSolver` fails and — see N-5 — the failure is swallowed) and into every
imputed Z-score in the window.

The risk is highest in `prep_recessive_impute` (`prep_qcatmix.cpp:196-208`): recessive coding
maps a MAF-0.05 variant to ≈0.25% ones, so in a small reference population the recoded genotype
is frequently constant and `SNP_Std_All_Pred_Recessive(i) == 0`.

**Fix:** return `NA_REAL` (or 0 with a counted warning) when the denominator is below a
tolerance, and filter such SNPs out before matrix assembly.

---

#### N-4 · `cal_af_norm_var()` divides by `mean*(1-mean)` and uses an unstable variance formula — **Medium** — difficulty **S**
**File:** [src/zmix.cpp:1206-1216](../src/zmix.cpp)

```cpp
double variance = sq_sum / n - mean * mean;      // catastrophic cancellation
double norm_var = variance / (mean * (1 - mean)); // 0/0 when the SNP is fixed across all pops
```
A SNP fixed in every reference population gives `mean ∈ {0,1}` and `norm_var = NaN`. That `NaN`
flows into `stats::quantile()` at `zmix.cpp:122`, which errors with *"missing values and NaN's not
allowed if 'na.rm' is FALSE"* — so `zmix()` aborts on a panel containing a single monomorphic
marker. `is.finite()` filtering in `R/zmix.R:52` happens too late to help.

The `E[X²] − E[X]²` form can also return a small negative variance for near-constant frequencies.
**Fix:** use a two-pass (or Welford) variance, skip SNPs with `mean` outside `(eps, 1-eps)`, and
pass `na.rm = TRUE` to `quantile`.

---

#### N-5 · `MakePosDef()` silently no-ops when the eigensolver fails — **Medium** — difficulty **S**
**File:** [src/util.cpp:299-315](../src/util.cpp)

```cpp
if(solver.info() != Eigen::Success){ return; }
```
The matrix is left unchanged — usually because it contained `NaN` (see N-3). The caller then
inverts a non-PD or `NaN` matrix via `fullPivLu().inverse()` and proceeds. `RmvPC` and `CountPC`
have the same silent early return.
**Fix:** propagate the failure (return a status, or `Rcpp::stop` with a diagnostic naming the window).

---

#### N-6 · `MakeSnpVec()`: division by zero and `ceil` used as "rounding" — **Low** — difficulty **S**
**File:** [src/gauss.cpp:585-591](../src/gauss.cpp)

`af1ref = allele_counter / (2 * num_subj)` is `NaN` when no population is flagged (`num_subj == 0`).
The next line, commented *"Round the allele frequency to 5 decimal places"*, is
`std::ceil(af1ref * 100000.0) / 100000.0` — a ceiling, not a round, biasing every reported
allele frequency upward by up to 1e-5 and (marginally) shifting which SNPs pass the MAF cutoff.
**Fix:** guard `num_subj > 0`; use `std::round`.

---

#### N-7 · `InvMat()` uses `fullPivLu().inverse()` on symmetric PD matrices — **Low** — difficulty **S**
**File:** [src/util.cpp:294-296](../src/util.cpp)

After `MakePosDef`, `B11`/`CovX`/`Cxx` are symmetric positive definite. `fullPivLu` is the most
expensive and least accurate option available. It is also used to invert the triangular `L` in
`run_qcat` (`qcat.cpp:207`), where a triangular solve is exact and `O(n²)`.
**Fix:** use `LLT` (or `LDLT`) and, wherever possible, `solve()` against the right-hand side
instead of forming an explicit inverse.

---

#### N-8 · `RmvPC` and `CountPC` disagree on the boundary case — **Low** — difficulty **S**
**File:** [src/util.cpp:317-380](../src/util.cpp)

`RmvPC` keeps components with `eig > cutoff`; `CountPC` drops components with `eig < cutoff`.
An eigenvalue exactly equal to `eig_cutoff` is counted by one and discarded by the other.
**Fix:** pick one comparison and share it.

---

### 2.5 Dead code and compiler warnings

#### D-1 · ~90 `-Wsign-compare` warnings and one uninitialised-variable warning — **Medium** — difficulty **M**
Verified with `clang++ -std=gnu++17 -Wall -Wextra`:

| File | count | File | count |
|---|---|---|---|
| `prep_qcatmix.cpp` | 13 | `gauss.cpp` | 9 |
| `zmix.cpp` | 8 | `distmix.cpp` | 8 |
| `qcatmix.cpp` | 8 | `gene.cpp` | 6 |
| `simulateLD.cpp` | 6 | `dist.cpp`, `qcat.cpp`, `prep_qcat.cpp` | 5 each |
| `computeLD.cpp` | 4 | `afmix.cpp`, `cpw2.cpp`, `jepegmix.cpp` | 1 each |

Plus `gauss.cpp:1336: variable 'categ_num' may be uninitialized` (see C-3) and
`bgzf.c:465: comparison of integers of different signs`.

The idiom `for(size_t i = 0; i < num_measured; i++)` where `num_measured` is `int` is pervasive.
These are benign today (all counts are small and non-negative) but they bury the one warning that
is a real bug, and CRAN's gcc/clang checks flag them.
**Fix:** make loop counters match the compared type (`std::size_t` throughout, or `int` with a
cast at the comparison), then build with `-Wall -Wextra -Werror` locally.

---

#### D-2 · Substantial dead and duplicated code — **Medium** — difficulty **M**

| Item | Location | Note |
|---|---|---|
| `prep_qcatmix()` | `prep_qcatmix.cpp:305-520` | ~215 lines commented out wholesale, including the `[[Rcpp::export]]` tag |
| `src/etc/` | 12 files | Older duplicates of `dist.cpp`, `jepeg.cpp`, `qcat.cpp`, `bgzf.c`, plus `gauss_main.cpp` and `RcppGSL.cpp` from the removed GSL era |
| `src/zmixrcpp.cpp_` | — | Stale renamed source |
| `read_input_cpw()` | `afmix.cpp:212-320` | Fully commented out |
| `read_input_jepeg`, `read_ref_index_jepeg` | `jepeg.cpp:157-280` | Fully commented out |
| `CalPopWgtVec()` | `util.cpp:497-560` | Commented-out GSL implementation |
| `CalMeanSumSq`, `CalCorMat`, `GetDiagMat`, `SubMatMat`, `AddNumMatDiag`, `SplitString` | `util.cpp` | Declared and defined, zero call sites |
| `RmvPC`, `CalBonfePval`, `CalSumUPval` | `util.cpp`, `gene.cpp` | Only reachable from commented-out lines |
| `Snp::flip_` + `FlipGenotype`/`FlipGenotypeVec` | `snp.h`, `util.cpp` | Flipping was disabled (`gauss.cpp:765-771`, `:849-857`); the field and helpers remain |
| `Arguments::mix_af1_cutoff`, `imp_info_cutoff` | `gauss.h:51, 68` | Set in the constructor, printed, never read |
| `Snp::PrintDistResult`, `PrintQcatResult`, `PrintSnpInfo(ofstream&)` | `snp.cpp` | Leftovers from the standalone CLI era |
| Triplicated `read_input_*` / `read_ref_index_*` | `afmix.cpp`, `cpw2.cpp`, `zmix.cpp` | Three near-identical copies of the `gauss.cpp` originals |
| `prep_zmix`, `prep_zmix2`, `prep_zmix3`, `prep_zmix4` | `zmix.cpp` | Superseded by `prep_zmix5`; all exported and documented |

**Consequence:** ~1,500 lines of code that must be read but never runs; three parser copies mean a
bug fix must be applied four times.
**Fix:** delete `src/etc/`, `zmixrcpp.cpp_` and all commented-out function bodies (git has them);
un-export the superseded `prep_zmix*` variants; collapse the three parser copies onto the
`gauss.cpp` implementations.

---

#### D-3 · Leftover debug output on user-facing paths — **Medium** — difficulty **S**

* [src/simulateLD.cpp](../src/simulateLD.cpp): `"Point 1"` … `"Point 11"` (lines 82–210), a loop
  printing **every** sampled genotype index (:154-157), and — worst — a
  `Rcpp::Rcout << k << " " << pop_num_sim_vec[k]` inside the `num_measured × num_pops` double loop
  (:168), i.e. tens of thousands of lines of console spam per call.
* [src/prep_qcat.cpp](../src/prep_qcat.cpp): `"push_vec"`, `"Make dataframe"`,
  `"release memory allocated for genotype"`, `"deletes snp map"`, `"return"`.
* `prep_qcatmix.cpp`, `qcat.cpp:212-214` (`VarZ1` / `VarLInvZ1` marked *"for testing"`).

**Fix:** remove, or move behind the existing `#ifdef *_Debug` convention.

---

#### D-4 · No `Rcpp::checkUserInterrupt()` anywhere — **Medium** — difficulty **S**
The `O(n²)`/`O(n³)` loops in `run_dist`, `run_distmix`, `run_qcat`, `computeLD`,
`prep_recessive_impute` and `prep_zmix5` can run for hours; Ctrl-C is ignored until the call
returns. **Fix:** call `Rcpp::checkUserInterrupt()` once per outer-loop iteration.

---

### 2.6 R package issues

#### R-1 · `data/` holds 63 MB of raw `.txt` files — **High** — difficulty **S**
`data/JEPEG_SNP_Annotation.v1.0.txt` (27 MB), `data/PGC3_SCZ_ilmn1M_Z.txt` (34 MB),
`PGC2_Chr22_ilmn1M_{Z,AF1}.txt`, `PGC2_3Mb.txt` — all tracked in git.

R treats `.txt` in `data/` as `data()`-loadable tables. Consequences:
* the source tarball is ~60 MB (CRAN's limit is 5 MB);
* `R CMD check` reports *"Undocumented data sets"* for all five (only `PGC2_SCZ_ANC_Prop` has an
  `.Rd`) and warns about the installed size;
* `data/.DS_Store` is also tracked.

**Fix:** move the example data to `inst/extdata/` (accessed via `system.file()`), keep only a
small subset in the repo, and host the rest alongside the reference panels. Add `data/.DS_Store`
and the other `.DS_Store` files to `.gitignore` and `git rm --cached` them.

---

#### R-2 · Missing `Imports` for base packages used at runtime — **High** — difficulty **S**
**Files:** `DESCRIPTION`, `NAMESPACE`

* `R/fiqt.R` calls `pnorm`, `qnorm`, `p.adjust` — **`stats`** is not imported.
* `R/zmix.R:60` calls `utils::read.table` — **`utils`** is not imported.
* `src/zmix.cpp:120` and `:273` call `Environment::namespace_env("stats")` — a runtime dependency
  on `stats` that neither `DESCRIPTION` nor `NAMESPACE` records.

`R CMD check` emits *"no visible global function definition for 'pnorm'"* and a note about the
undeclared dependency; the code only works today because `stats`/`utils` happen to be attached.
**Fix:** add `stats` and `utils` to `Imports:` and `importFrom()` the specific symbols
(or add `@importFrom stats pnorm qnorm p.adjust quantile` / `@importFrom utils read.table`).

---

#### R-3 · `gauss::dist()` masks `stats::dist()` — **High** — difficulty **M**
**File:** [src/dist.cpp:30](../src/dist.cpp)

`library(gauss)` shadows the base-R distance function with a nine-argument GWAS imputation
function. Any user (or dependent package) script calling `dist(x)` after attaching gauss gets a
confusing argument-matching error. `R CMD check` reports the conflict.
**Fix:** rename to something like `distimp()` / `gauss_dist()`, keeping `dist()` as a deprecated
alias for one release cycle.

---

#### R-4 · `exportPattern("^[[:alpha:]]+")` exports everything — **Medium** — difficulty **S**
**File:** `NAMESPACE` (generated from `R/gauss.R:2`)

All 20 functions are exported, including the experimental `prep_zmix`–`prep_zmix5_sup`, `cpw2`,
`prep_qcat` and `prep_recessive_impute` scaffolding. This freezes them as public API and forces
`.Rd` files for each.
**Fix:** replace the pattern with explicit `@export` tags on the ~8 user-facing functions.

---

#### R-5 · Vignettes cannot be built — **High** — difficulty **M**
**Files:** `vignettes/*.Rmd`, `DESCRIPTION`

Every vignette does `library(kableExtra)` / `library(data.table)` / `library(tidyverse)` /
`library(dplyr)` / `library(reshape2)`; `DESCRIPTION` lists only `rmarkdown` and `knitr` under
`Suggests`. Each vignette also reads `../data/...` and `../ref/33KG/...`, and `ref/` is
`.gitignore`d (multi-GB panel, not distributed).

`VignetteBuilder: knitr` means `R CMD build` **executes** these, so the build fails on a clean
checkout.
**Fix:** add the packages to `Suggests`, and either set `eval = FALSE` in
`knitr::opts_chunk$set()` with pre-rendered output, or gate every panel-dependent chunk on
`file.exists()`. Note the paths must be package-relative (`system.file`), not `../`.

---

#### R-6 · `PKG_CXXFLAGS +=` requires GNU make, which is not declared — **Medium** — difficulty **S**
**Files:** `src/Makevars.in`, `configure`, `configure.win`, `DESCRIPTION`

`+=` is a GNU make extension. Writing to Rtools/Solaris' make is the classic failure mode, and
`DESCRIPTION` has no `SystemRequirements: GNU make`.

Additionally the `configure` script does nothing but `cp src/Makevars.in src/Makevars` (and
`configure.win` copies to `Makevars.win`) — a whole configure step for a static copy. `Makevars`
is in `.Rbuildignore` but the generated `src/Makevars.win` is not, so a stale one can be shipped.
**Fix:** use `PKG_CXXFLAGS = -include rcpp_eigen_wrap.h`, ship `src/Makevars` and
`src/Makevars.win` directly, and delete `configure`, `configure.win` and `Makevars.in`.

---

#### R-7 · Stale and incomplete `DESCRIPTION` metadata — **Low** — difficulty **S**
* `Version: 1.0` / `Date: 2021-10-12`, though the package cites a 2024 *Bioinformatics* paper and
  has had substantive changes through 2026.
* No `Depends: R (>= 3.5.0)` despite `.RData` in `data/` and C++11+ requirements.
* No `NEWS.md`, so users have no changelog for the GSL removal or the `zmix` rewrite.
**Fix:** bump the version, refresh the date, add `Depends`, start a `NEWS.md`.

---

#### R-8 · `.DS_Store` files tracked in git — **Low** — difficulty **S**
`.DS_Store`, `R/.DS_Store`, `data/.DS_Store`, `docs/.DS_Store`, `man/.DS_Store`,
`src/.DS_Store`, `vignettes/.DS_Store` are all committed (and two show as modified in the current
working tree).
**Fix:** `git rm --cached` them and add `.DS_Store` to `.gitignore`.

---

### 2.7 Portability

#### P-1 · `simulateLD()` uses `std::random_device`/`std::mt19937` instead of R's RNG — **High** — difficulty **S**
**File:** [src/simulateLD.cpp:134-135](../src/simulateLD.cpp)

Results are not reproducible, ignore `set.seed()`, and are not covered by R's RNG stream. On some
platforms `std::random_device` is deterministic (a known libstdc++/MinGW issue), so on Windows the
"random" sample may be identical across runs. This violates *Writing R Extensions* §6.3.
**Fix:** use `GetRNGstate()` / `unif_rand()` / `PutRNGstate()` (or `Rcpp::sample`).

---

#### P-2 · `long long int bp` truncated to `int` in output — **Medium** — difficulty **S**
**Files:** [src/dist.cpp:89](../src/dist.cpp), `distmix.cpp:99`, and every `bp_vec.push_back()`
into an `IntegerVector` (`computeLD.cpp`, `simulateLD.cpp`, `qcat.cpp`, `qcatmix.cpp`,
`prep_qcat.cpp`, `prep_qcatmix.cpp`)

`Snp::bp_` is `long long`, but the filter variable is `int bp = (*it_sv)->GetBp();` and the R
column is an `IntegerVector`. Human chromosome 1 tops out around 2.5×10⁸ so this happens to fit,
but the truncation is silent and would break for any non-human or concatenated coordinate system.
The `Rcpp::export` signature also takes `long long int start_bp`, which R passes as a `double` —
positions above 2⁵³ would round.
**Fix:** keep `long long` internally and emit `NumericVector` (double) for `bp`.

---

#### P-3 · Vendored BGZF errors are not surfaced — **Medium** — difficulty **M**
**Files:** `src/bgzf.c`, all `bgzf_seek` / `bgzf_open` call sites

`report_error` writes into `fp->error` but nothing reads it (`bgzf.c:137` shows the `fprintf(stderr)`
was commented out to silence an R check warning, leaving no channel at all). No call site checks
`bgzf_seek`'s return value. A truncated or non-BGZF-compressed reference file therefore produces
plausible-looking but wrong numbers instead of an error.
**Fix:** check `bgzf_seek(...) < 0` and `fp->error` at every call site and `Rcpp::stop` with the
message. Longer term, consider depending on `Rhtslib` rather than a decade-old samtools fork.

---

#### P-4 · `::toupper` applied through `std::transform` — **Low** — difficulty **S**
**Files:** `distmix.cpp:52`, `computeLD.cpp:46`, `jepegmix.cpp:41`, `simulateLD.cpp:59`,
`prep_qcatmix.cpp:60`

`std::transform(s.begin(), s.end(), s.begin(), ::toupper)` passes a possibly-negative `char` to a
function whose argument must be representable as `unsigned char` — UB for non-ASCII input, and
locale-dependent (Turkish `i` in particular).
**Fix:** wrap in a lambda: `[](unsigned char c){ return std::toupper(c); }`.

---

### 2.8 Test coverage

#### T-1 · No automated tests at all — **High** — difficulty **L**
There is no `tests/` directory, no `testthat`, no CI workflow. `inst/dev-tests/` contains 27 manual
driver scripts that hard-code absolute paths and need the multi-GB 33KG panel.

**Consequence:** every finding above (especially C-1, C-3, C-5 and C-6, which change numeric
output without changing behaviour) is undetectable by any current check. Any refactor of the
duplicated parsers is unsafe.

**Recommended plan**
1. **Fixtures.** Build a tiny synthetic reference panel — say 3 populations × 20 samples × 200 SNPs
   — as `inst/extdata/`: a `bgzip`ped index, a `bgzip`ped genotype file, and a `pop_desc.txt`.
   Small enough to commit (<200 KB) and to run in seconds.
2. **Unit tests for the numeric kernels.** Expose (or test via a `testthat`-visible wrapper)
   `CalCor`, `CalWgtCov` and `CalCorSup`, and compare against R's `cor()`/`cov()` on the same
   genotype matrix. **This alone would have caught C-1.**
3. **Golden-output tests.** Run `dist()`, `distmix()`, `computeLD()`, `qcat()`, `jepeg()` and
   `zmix()` on the fixture panel and snapshot the results with `testthat::expect_snapshot_value`.
4. **Error-path tests.** Missing files, empty input, duplicate rows, `interval = 0`,
   unmatched population names, all-monomorphic SNPs — assert a clear `Rcpp::stop` message rather
   than a hang or a crash.
5. **Pure-R tests for `fiqt()`** against a hand-computed reference, including `NA` and length-0
   input.
6. **CI.** A GitHub Actions `R-CMD-check` matrix (ubuntu/macOS/windows, release + devel) plus one
   ASAN/UBSAN run (`rocker/r-devel-san`) — the latter would catch C-4 and C-5 directly.

---

### 2.9 Documentation

#### DOC-1 · `prep_zmix*` `@return` documentation is wrong — **Medium** — difficulty **S**
**Files:** [src/zmix.cpp:41, 199, 361, 509, 649, 938](../src/zmix.cpp) → `man/prep_zmix*.Rd`

Every one says *"R data frame containing population IDs and weights"*. They actually return a
**numeric matrix** whose first column is `z_i·z_j` for a SNP pair and whose remaining columns are
per-population (or per-superpopulation) LD correlations — the design matrix for the quadratic
program in `R/zmix.R`. A user following the documentation cannot use the result.

`prep_zmix5` is additionally **declared** as `NumericVector` (`zmix.cpp:44`) while returning a
`NumericMatrix`; `R/zmix.R:47` then relies on the `dim` attribute surviving the slice. It works
today, but the generated `.Rd` promises a vector and the code is one Rcpp change away from breaking.
**Fix:** change the declared return type to `NumericMatrix` and rewrite `@return` to describe the
actual matrix (dimensions, column meanings, column order).

---

#### DOC-2 · Percentile semantics contradict the code — **Medium** — difficulty **S**
**File:** [src/zmix.cpp:26-28, 127-135](../src/zmix.cpp), `R/zmix.R:22`

* The roxygen block says *"choose top 1% of SNPs"*; the C++ default is `pct = 0.99` but
  `R/zmix.R`'s default is `0.9` (top 10%) — the documented and actual defaults differ.
* The inline comment at `zmix.cpp:127` reads *"extract SNPs with normalized var of af1 **<** 99th
  percentile"* while the code keeps `snp_af_var_vec[i] > cutoff`.
* `interval` is documented as *"stepping distance … for selecting the first SNP of each pair"*
  (inherited from `prep_zmix4`), but in `prep_zmix5` it just decimates the SNP list; the
  all-pairs loop that follows ignores it.

**Fix:** align the defaults, fix the comparison comment, and rewrite `@param interval` for the
`prep_zmix5` algorithm.

---

#### DOC-3 · `README` typo and a workflow that does not compose — **Medium** — difficulty **S**
`README.Rmd:21` (and the generated `README.md:19`, `docs/index.md:13`) says **`zimx()`** —
no such function; it should be `zmix()`.

More substantively, the "Suggested Workflow" section presents `afmix()`/`zmix()` output as the
input to `distmix()`/`jepegmix()`/`computeLD()`, which does not work as written (C-7). Neither
the README nor the vignettes state that the first column must be dropped.
**Fix:** correct the typo in `README.Rmd` and re-knit; after fixing C-7, add a short worked
example showing the handoff.

---

#### DOC-4 · pkgdown site URL is the git repo, breaking generated links — **Low** — difficulty **S**
`_pkgdown.yml:1` sets `url: https://github.com/statsleelab/gauss`, so `docs/pkgdown.yml` records
`reference: https://github.com/statsleelab/gauss/reference` — a 404. The site itself is served
from `https://statsleelab.github.io/gauss`.
**Fix:** `url: https://statsleelab.github.io/gauss` and rebuild. Adding the same value to
`DESCRIPTION`'s `URL:` would also let roxygen generate correct cross-links.

---

#### DOC-5 · Undocumented behaviour and missing vignettes — **Low** — difficulty **M**
* No vignette or reference entry explains the reference-panel **file formats** (index columns,
  genotype record layout, `pop_desc` columns) — yet every function requires them and
  `read_ref_desc` will silently produce garbage on a mis-shaped `pop_desc.txt`.
* The `type` column (0 / 1 / 2) in `dist`/`distmix`/`qcat` output is never explained in any `.Rd`.
* `zmix()`, `qcat()`/`qcatmix()`, `simulateLD()` and `prep_recessive_impute()` have no vignette
  and are absent from the pkgdown navbar.
* `fiqt()`'s roxygen block has no `@export`, no `@examples` and no reference to the FIQT paper.

---

### 2.10 Performance

#### PERF-1 · `MakeSnpVec` / `MakeSnpVecMix` re-derive allele frequencies the index already contains — **High** — difficulty **M**
**Files:** [src/gauss.cpp:546-606](../src/gauss.cpp), `:634-692`; index parsing at `:331`, `:471`

`ReadReferenceIndex` parses the index's `af1ref` column and immediately discards it
(*"Reference allele frequency (not used in this function)"*). `MakeSnpVec` then seeks to each
SNP's `fpos` and reads its **entire** genotype record — for the 33KG panel that is ~33,000
characters per population × 29 populations ≈ 1 MB per SNP — purely to recompute an allele
frequency, then applies a MAF cutoff that usually discards the SNP anyway.

For genome-wide `jepeg()`/`jepegmix()` (which call `ReadInputZ(..., All=true)`) this is a full
scan of the multi-GB panel — one seek + one megabyte-scale decompression per SNP in the map,
including type-2 SNPs whose seek is invalid (C-5).

**Fix:** carry the index's `af1ref` onto the `Snp` in `ReadReferenceIndex*` and apply the cutoff
before touching the data file. For `MakeSnpVecMix`, add per-population `af1` columns to the index
(a one-off panel regeneration) or, failing that, store a separate small allele-frequency file.
Expected speedup for JEPEG-class workloads: one to two orders of magnitude.

---

#### PERF-2 · `Rcpp` vector `push_back` in output loops is O(n²) — **Medium** — difficulty **S**
**Files:** [src/dist.cpp:88-101](../src/dist.cpp), `distmix.cpp:98-111`, `computeLD.cpp:96-100` and
`:130-137`, `simulateLD.cpp:205-212`, `qcat.cpp`, `qcatmix.cpp`, `prep_qcat.cpp:145-156`,
`zmix.cpp:1183-1218` (`cal_af_norm_var`)

`Rcpp::Vector::push_back` reallocates and copies the whole vector on every call. `dist()` emits
one row per SNP in the window (10⁴–10⁵), across ten columns; `cal_af_norm_var` pushes once per
candidate SNP, which for `interval = 1` is the full measured set.
**Fix:** count first, pre-size (`NumericVector v(n)`) and assign by index — or accumulate in
`std::vector` and wrap once.

---

#### PERF-3 · `prep_zmix5` builds an `O(n²)` matrix with a 32-bit row count — **High** — difficulty **M**
**File:** [src/zmix.cpp:151-154](../src/zmix.cpp), `:314-315`, `:1017-1020`

```cpp
int total_rows = (snp_subvec_size * (snp_subvec_size - 1)) / 2;
NumericMatrix data_mat(total_rows, 1 + args.num_pops);
```
Two problems:
* **Overflow.** The multiplication is `int × int`. `snp_subvec_size = 65,536` already overflows
  signed 32-bit — undefined behaviour, and in practice a negative or wrapped row count passed to
  `NumericMatrix`.
* **Memory.** Even without overflow, 10,000 ancestry-informative SNPs → 5×10⁷ rows × 30 columns ×
  8 bytes ≈ **12 GB**, allocated as a single R object, before `R/zmix.R` even sees it. The
  `R/zmix.R` defaults (`percentile = 0.9`, `interval = 10`) make this reachable on a real
  genome-wide Z-score file.

The `crossprod(x)` in `R/zmix.R:88` only needs the `p×p` Gram matrix and the `p`-vector
`t(y) %*% x`, so the full design matrix never needs to exist.
**Fix:** accumulate `X'X` (30×30) and `X'y` (30×1) incrementally in C++ and return those; the
memory becomes `O(p²)` regardless of SNP count, and the overflow disappears. Interim mitigation:
use `R_xlen_t`/`double` for `total_rows` and `Rcpp::stop` above a size threshold.

---

#### PERF-4 · Genotype correlations recompute per-SNP sums `O(n²)` times — **Medium** — difficulty **M**
**Files:** [src/util.cpp:49-70](../src/util.cpp) (`CalCor`), `zmix.cpp:1221-1246` (`CalCorSup`);
callers `dist.cpp:166-175`, `distmix.cpp:187-200`, `qcat.cpp:183-199`, `computeLD.cpp:104-116`,
`gene.cpp:571-586`, `prep_qcatmix.cpp:167-240`

Each call re-walks both genotype strings accumulating `Σx`, `Σx²`, `Σy`, `Σy²`, `Σxy`. In an
`n × n` LD matrix, each SNP's `Σx` and `Σx²` are recomputed `n` times.
**Fix:** precompute per-SNP mean and standard deviation once (`O(n·m)`), then only `Σxy` in the
inner loop. Converting genotypes from ASCII strings to `std::vector<double>` (or `Eigen::VectorXd`)
once per window would additionally let Eigen vectorise the cross-products — the dominant cost in
`dist`/`distmix`/`computeLD`.

---

#### PERF-5 · `simulateLD` allocates two vectors per matrix entry — **Medium** — difficulty **S**
**File:** [src/simulateLD.cpp:190-199](../src/simulateLD.cpp)

```cpp
for(size_t j=i+1; j<num_measured; j++){
  NumericVector snp_i = geno_mat(i, _);   // re-allocated for every j
  NumericVector snp_j = geno_mat(j, _);
  ...
}
```
`snp_i` is invariant in the inner loop. **Fix:** hoist `snp_i`, and prefer an `Eigen::MatrixXd`
with column-major SNP layout so `CalCor` takes cheap `Eigen::Ref` views.

---

#### PERF-6 · Redundant `O(p³)` factorisations — **Low** — difficulty **S**
`run_dist`/`run_distmix` call `MakePosDef` (a full `SelfAdjointEigenSolver`) and then
`InvMat` (a `fullPivLu` inverse) on the same `n_measured × n_measured` matrix — two dense
`O(n³)` factorisations where a single `LDLT` with eigenvalue clipping, followed by `solve()`
per right-hand side, would suffice. See also N-7.

---

## 3. Suggested remediation order

| # | Finding | Why first |
|---|---|---|
| 1 | **C-2** `exit()` | One line; currently destroys user sessions |
| 2 | **C-3** uninitialised `categ_num` | One line; UB and silently wrong JEPEG results |
| 3 | **T-1** fixture panel + kernel tests | Prerequisite for safely changing any numeric code |
| 4 | **C-1** `CalWgtCov` | The single largest correctness defect; needs test 3 to validate |
| 5 | **C-5**, **C-4**, **V-1**, **V-2** | Invalid seeks, OOB reads, unvalidated input |
| 6 | **C-7**, **DOC-3** | The documented workflow does not run |
| 7 | **R-1**, **R-2**, **R-5**, **R-6** | Makes `R CMD check` pass — a gate on everything else |
| 8 | **PERF-3**, **PERF-1** | Overflow/OOM risk, then the dominant runtime cost |
| 9 | **M-1** RAII | Best done as one sweep, and it deletes ~15 cleanup loops |
| 10 | **D-1**, **D-2**, **D-3** | Warning-clean, dead-code-free build |
| 11 | **N-1**–**N-8**, **P-1**–**P-4** | Numerical and portability hardening |
| 12 | **DOC-1**, **DOC-2**, **DOC-4**, **DOC-5**, **R-3**, **R-4** | API and documentation cleanup for the next release |

---

## 4. Summary counts

| Severity | Count |
|---|---|
| Critical | 3 |
| High | 16 |
| Medium | 24 |
| Low | 13 |
| **Total** | **56** |

| Category | Findings | Count |
|---|---|---|
| Correctness | C-1 … C-11 | 11 |
| Memory / resources | M-1 … M-4 | 4 |
| Input validation | V-1 … V-5 | 5 |
| Numerical | N-1 … N-8 | 8 |
| Dead code / warnings | D-1 … D-4 | 4 |
| R packaging | R-1 … R-8 | 8 |
| Portability | P-1 … P-4 | 4 |
| Test coverage | T-1 | 1 |
| Documentation | DOC-1 … DOC-5 | 5 |
| Performance | PERF-1 … PERF-6 | 6 |
