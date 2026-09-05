//GAUSS : Genome Analysis Using Summary Statistics
//test_kernels.cpp
//
// Thin R-visible wrappers around the correlation/covariance kernels that every
// analysis in this package is built on.  They exist so that the kernels can be
// compared against R's own cor()/cov() on the same genotype matrix from the
// test suite; nothing in the package calls them.
//
// All of them are exported under a leading-dot name, so NAMESPACE's
// exportPattern("^[[:alpha:]]+") leaves them out of the public API.  Reach them
// from tests with gauss:::.cal_cor_pooled() and friends.
//
// Genotype arguments follow the same convention as the rest of the package: one
// string per population, one character per individual, each in {'0','1','2'}.

#include <string>
#include <vector>

#include <Rcpp.h>

#include "util.h"

// CalCorSup() lives in zmix.cpp and is not declared in any header.
double CalCorSup(const std::vector<std::string>& x,
                 const std::vector<std::string>& y,
                 const std::vector<int>& pop_idx);

// Correlation over the pooled samples of every supplied population.
// [[Rcpp::export(.cal_cor_pooled)]]
double cal_cor_pooled(std::vector<std::string> x, std::vector<std::string> y) {
  return CalCor(x, y);
}

// Correlation over a single genotype string.
// [[Rcpp::export(.cal_cor_one)]]
double cal_cor_one(std::string x, std::string y) {
  return CalCor(x, y);
}

// Population-weighted correlation (the mixture path).
// [[Rcpp::export(.cal_cor_wgt)]]
double cal_cor_wgt(std::vector<std::string> x, std::vector<std::string> y,
                   std::vector<double> pop_wgt) {
  return CalCor(x, y, pop_wgt);
}

// Population-weighted covariance (the mixture path).
// [[Rcpp::export(.cal_wgt_cov)]]
double cal_wgt_cov(std::vector<std::string> x, std::vector<std::string> y,
                   std::vector<double> pop_wgt) {
  return CalWgtCov(x, y, pop_wgt);
}

// Correlation over the pooled samples of the populations named by pop_idx.
// pop_idx is zero-based, matching the C++ caller in zmix.cpp.
// [[Rcpp::export(.cal_cor_sup)]]
double cal_cor_sup(std::vector<std::string> x, std::vector<std::string> y,
                   std::vector<int> pop_idx) {
  return CalCorSup(x, y, pop_idx);
}
