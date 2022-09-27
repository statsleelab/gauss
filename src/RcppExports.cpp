// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cal_pop_wgt
DataFrame cal_pop_wgt(std::string input_file, std::string reference_index_file, std::string reference_data_file, std::string reference_pop_desc_file, Rcpp::Nullable<int> interval);
RcppExport SEXP _gauss_cal_pop_wgt(SEXP input_fileSEXP, SEXP reference_index_fileSEXP, SEXP reference_data_fileSEXP, SEXP reference_pop_desc_fileSEXP, SEXP intervalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_index_file(reference_index_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_data_file(reference_data_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_pop_desc_file(reference_pop_desc_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type interval(intervalSEXP);
    rcpp_result_gen = Rcpp::wrap(cal_pop_wgt(input_file, reference_index_file, reference_data_file, reference_pop_desc_file, interval));
    return rcpp_result_gen;
END_RCPP
}
// dist
DataFrame dist(int chr, long long int start_bp, long long int end_bp, long long int wing_size, std::string study_pop, std::string input_file, std::string reference_index_file, std::string reference_data_file, std::string reference_pop_desc_file, Rcpp::Nullable<double> af1_cutoff);
RcppExport SEXP _gauss_dist(SEXP chrSEXP, SEXP start_bpSEXP, SEXP end_bpSEXP, SEXP wing_sizeSEXP, SEXP study_popSEXP, SEXP input_fileSEXP, SEXP reference_index_fileSEXP, SEXP reference_data_fileSEXP, SEXP reference_pop_desc_fileSEXP, SEXP af1_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< long long int >::type start_bp(start_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type end_bp(end_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type wing_size(wing_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type study_pop(study_popSEXP);
    Rcpp::traits::input_parameter< std::string >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_index_file(reference_index_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_data_file(reference_data_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_pop_desc_file(reference_pop_desc_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type af1_cutoff(af1_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(dist(chr, start_bp, end_bp, wing_size, study_pop, input_file, reference_index_file, reference_data_file, reference_pop_desc_file, af1_cutoff));
    return rcpp_result_gen;
END_RCPP
}
// distmix
DataFrame distmix(int chr, long long int start_bp, long long int end_bp, long long int wing_size, DataFrame pop_wgt_df, std::string input_file, std::string reference_index_file, std::string reference_data_file, std::string reference_pop_desc_file, Rcpp::Nullable<double> af1_cutoff);
RcppExport SEXP _gauss_distmix(SEXP chrSEXP, SEXP start_bpSEXP, SEXP end_bpSEXP, SEXP wing_sizeSEXP, SEXP pop_wgt_dfSEXP, SEXP input_fileSEXP, SEXP reference_index_fileSEXP, SEXP reference_data_fileSEXP, SEXP reference_pop_desc_fileSEXP, SEXP af1_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< long long int >::type start_bp(start_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type end_bp(end_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type wing_size(wing_sizeSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pop_wgt_df(pop_wgt_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_index_file(reference_index_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_data_file(reference_data_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_pop_desc_file(reference_pop_desc_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type af1_cutoff(af1_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(distmix(chr, start_bp, end_bp, wing_size, pop_wgt_df, input_file, reference_index_file, reference_data_file, reference_pop_desc_file, af1_cutoff));
    return rcpp_result_gen;
END_RCPP
}
// qcat
DataFrame qcat(int chr, long long int start_bp, long long int end_bp, long long int wing_size, std::string study_pop, std::string input_file, std::string reference_index_file, std::string reference_data_file, std::string reference_pop_desc_file, Rcpp::Nullable<double> af1_cutoff);
RcppExport SEXP _gauss_qcat(SEXP chrSEXP, SEXP start_bpSEXP, SEXP end_bpSEXP, SEXP wing_sizeSEXP, SEXP study_popSEXP, SEXP input_fileSEXP, SEXP reference_index_fileSEXP, SEXP reference_data_fileSEXP, SEXP reference_pop_desc_fileSEXP, SEXP af1_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< long long int >::type start_bp(start_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type end_bp(end_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type wing_size(wing_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type study_pop(study_popSEXP);
    Rcpp::traits::input_parameter< std::string >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_index_file(reference_index_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_data_file(reference_data_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_pop_desc_file(reference_pop_desc_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type af1_cutoff(af1_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(qcat(chr, start_bp, end_bp, wing_size, study_pop, input_file, reference_index_file, reference_data_file, reference_pop_desc_file, af1_cutoff));
    return rcpp_result_gen;
END_RCPP
}
// qcatmix
DataFrame qcatmix(int chr, long long int start_bp, long long int end_bp, long long int wing_size, DataFrame pop_wgt_df, std::string input_file, std::string reference_index_file, std::string reference_data_file, std::string reference_pop_desc_file, Rcpp::Nullable<double> af1_cutoff);
RcppExport SEXP _gauss_qcatmix(SEXP chrSEXP, SEXP start_bpSEXP, SEXP end_bpSEXP, SEXP wing_sizeSEXP, SEXP pop_wgt_dfSEXP, SEXP input_fileSEXP, SEXP reference_index_fileSEXP, SEXP reference_data_fileSEXP, SEXP reference_pop_desc_fileSEXP, SEXP af1_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type chr(chrSEXP);
    Rcpp::traits::input_parameter< long long int >::type start_bp(start_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type end_bp(end_bpSEXP);
    Rcpp::traits::input_parameter< long long int >::type wing_size(wing_sizeSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type pop_wgt_df(pop_wgt_dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type input_file(input_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_index_file(reference_index_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_data_file(reference_data_fileSEXP);
    Rcpp::traits::input_parameter< std::string >::type reference_pop_desc_file(reference_pop_desc_fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type af1_cutoff(af1_cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(qcatmix(chr, start_bp, end_bp, wing_size, pop_wgt_df, input_file, reference_index_file, reference_data_file, reference_pop_desc_file, af1_cutoff));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gauss_cal_pop_wgt", (DL_FUNC) &_gauss_cal_pop_wgt, 5},
    {"_gauss_dist", (DL_FUNC) &_gauss_dist, 10},
    {"_gauss_distmix", (DL_FUNC) &_gauss_distmix, 10},
    {"_gauss_qcat", (DL_FUNC) &_gauss_qcat, 10},
    {"_gauss_qcatmix", (DL_FUNC) &_gauss_qcatmix, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_gauss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
