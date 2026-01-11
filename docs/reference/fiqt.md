# FDR Inverse Quantile Transformation (FIQT)

FDR Inverse Quantile Transformation (FIQT)

## Usage

``` r
fiqt(z, min.p = 10^-300)
```

## Arguments

- z:

  association Z-score vector

- min.p:

  minimum p-value admitted (to avoid zero p-values/adjusted p-values
  which give troubles with inverse cdf)

## Value

R vector containing the Winner's Curse adjusted Z-scores
