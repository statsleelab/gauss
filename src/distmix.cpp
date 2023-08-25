#include <Rcpp.h>

using namespace Rcpp;

#include <string>
#include <vector>
#include <map>
#include <deque>
#include "snp.h"
#include "gauss.h"
#include "util.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>

void run_distmix(std::vector<Snp*>& snp_vec, Arguments& args);

//' Direct imputation of summary statistics for unmeasured SNPs from mixed ethnicity cohorts
//' 
//' @param chr chromosome number
//' @param start_bp start base pair position of prediction window
//' @param end_bp end base pair position of prediction window
//' @param wing_size the size of the area flanking the left and right of the prediction window
//' @param pop_wgt_df R data frame containing population IDs and weights
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2, af1, and z
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R data frame containing rsid, chr, bp, a1, a2, af1mix, z, pval, info, type 
// [[Rcpp::export]]
DataFrame distmix(int chr, 
                  long long int start_bp, 
                  long long int end_bp, 
                  long long int wing_size,
                  DataFrame pop_wgt_df,
                  std::string input_file, 
                  std::string reference_index_file, 
                  std::string reference_data_file,
                  std::string reference_pop_desc_file,
                  Rcpp::Nullable<double> af1_cutoff = R_NilValue){
  
  Arguments args;
  args.chr = chr; 
  args.start_bp = start_bp;
  args.end_bp = end_bp;
  args.wing_size = wing_size;
  
  // add pop_wgt_df info in args.pop_wgt_map
  std::vector<std::string> pop_vec_in = as<std::vector<std::string>>(pop_wgt_df[0]);
  std::vector<double> pop_wgt_vec_in = as<std::vector<double>>(pop_wgt_df[1]);
  for(int i=0; i<pop_vec_in.size(); i++){
    std::string pop = pop_vec_in[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper); //make capital
    args.pop_wgt_map[pop]=pop_wgt_vec_in[i];
  }  
  
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  
  if(af1_cutoff.isNotNull()){
    args.af1_cutoff = Rcpp::as<double>(af1_cutoff);
  } else {
    args.af1_cutoff = 0.01;
  }
  
  read_ref_desc(args);
  init_pop_flag_wgt_vec(args);
  //args.PrintArguments();
  
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  ReadInputZ(snp_map, args, false);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  ReadReferenceIndex(snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  
  // make a snp vector containing all SNPs
  std::vector<Snp*> snp_vec;
  MakeSnpVecMix(snp_vec, snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_vec.size() <<std::endl;
  ReadGenotype(snp_vec, args);

  // run DISTMIX imputation
  run_distmix(snp_vec, args);
  
  // release memory allocated for genotype
  FreeGenotype(snp_vec);
  
  StringVector rsid_vec;
  IntegerVector chr_vec;
  IntegerVector bp_vec;
  StringVector a1_vec;
  StringVector a2_vec;
  NumericVector af1mix_vec;
  NumericVector z_vec;
  NumericVector pval_vec;
  NumericVector info_vec;
  IntegerVector type_vec;
  
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int bp = (*it_sv)->GetBp();
    if(bp >= start_bp && bp <= end_bp){
      rsid_vec.push_back((*it_sv)->GetRsid());
      chr_vec.push_back((*it_sv)->GetChr());
      bp_vec.push_back((*it_sv)->GetBp());
      a1_vec.push_back((*it_sv)->GetA1());
      a2_vec.push_back((*it_sv)->GetA2());
      af1mix_vec.push_back((*it_sv)->GetAf1Mix());
      z_vec.push_back((*it_sv)->GetZ());
      pval_vec.push_back(2*gsl_cdf_ugaussian_Q(std::abs((*it_sv)->GetZ())));
      info_vec.push_back((*it_sv)->GetInfo());
      type_vec.push_back((*it_sv)->GetType());
    }
  }

  //deletes snp_map.
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end();){
    (it_sm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_sm->second;        // delete snp object
    snp_map.erase(it_sm++);      // delete map element
  }

  DataFrame df = DataFrame::create(Named("rsid")=rsid_vec,
                                   Named("chr")=chr_vec,
                                   Named("bp")=bp_vec,
                                   Named("a1")=a1_vec,
                                   Named("a2")=a2_vec,
                                   Named("af1mix")=af1mix_vec,
                                   Named("z")=z_vec,
                                   Named("pval")=pval_vec,
                                   Named("info")=info_vec,
                                   Named("type")=type_vec);
  return df;  
}


void run_distmix(std::vector<Snp*>& snp_vec, Arguments& args){
  std::deque<Snp*> sliding_window_measured;
  std::deque<Snp*> sliding_window_unmeasured;
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    long long int bp = (*it_sv)->GetBp();
    if(type == 0 && (bp >= args.start_bp && bp <= args.end_bp)){ // unmeasured and in the pred win.
      sliding_window_unmeasured.push_back(*it_sv);
    } else if(type == 1) { // if measured
      sliding_window_measured.push_back(*it_sv);
    } // if (type == 2) don't put the snp in the sliding window. Do nothing. 
  }
  
  int num_measured = sliding_window_measured.size();
  int num_unmeasured = sliding_window_unmeasured.size();
  
  if(sliding_window_measured.size() <= args.min_num_measured_snp || 
     sliding_window_unmeasured.size() <= args.min_num_unmeasured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured<<std::endl;
    Rcpp::Rcout<<"Number of unmeasured SNPs: "<<num_unmeasured<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - DISTMIX not performed");
  }
  
  /////////////////////////
  // run DIST imputation //
  /////////////////////////
  gsl_matrix* Z1 = gsl_matrix_calloc(num_measured, 1);
  gsl_vector* SNP_STD_VEC = gsl_vector_calloc(num_measured + num_unmeasured); //vector of SNP genotype standard deviations
  gsl_matrix* B11 = gsl_matrix_calloc(num_measured,num_measured); // correlation matrix B11 : correlation among measured SNPs 
  gsl_matrix* B11Inv = gsl_matrix_calloc(num_measured, num_measured);
  
  gsl_matrix* b21 = gsl_matrix_calloc(1, num_measured);
  gsl_matrix* b12 = gsl_matrix_calloc(num_measured, 1);
  gsl_matrix* b21B11Inv = gsl_matrix_calloc(1, num_measured);
  gsl_matrix* val = gsl_matrix_calloc(1, 1);
  
  // Init Z1 vector
  for(size_t i=0; i<num_measured; i++){   
    gsl_matrix_set(Z1, i, 0, (*sliding_window_measured[i]).GetZ());
  }
  // Init SNP_STD_VEC	
  for(size_t i=0; i < num_measured; i++){
    double v = CalWgtCov((*sliding_window_measured[i]).GetGenotypeVec(), (*sliding_window_measured[i]).GetGenotypeVec(), args.pop_wgt_vec);
    gsl_vector_set(SNP_STD_VEC, i, std::sqrt(v));
  }
  for(size_t i=0; i < num_unmeasured; i++){
    double v = CalWgtCov((*sliding_window_unmeasured[i]).GetGenotypeVec(), (*sliding_window_unmeasured[i]).GetGenotypeVec(), args.pop_wgt_vec);
    gsl_vector_set(SNP_STD_VEC, i + num_measured, std::sqrt(v));
  }
  // Init B11 matrix
  Rcpp::Rcout<<"Computing correlations between variants..."<<std::endl;
  for(size_t i=0; i<B11->size1; i++){
    gsl_matrix_set(B11, i, i, 1.0 + args.lambda); //add LAMBDA here (ridge regression trick)
    double stdi = gsl_vector_get(SNP_STD_VEC, i);
    for(size_t j=i+1; j<B11->size1; j++){
      double stdj = gsl_vector_get(SNP_STD_VEC, j);
      double cov = CalWgtCov((*sliding_window_measured[i]).GetGenotypeVec(), (*sliding_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov/(stdi*stdj);
      gsl_matrix_set(B11, i, j, cor);
      gsl_matrix_set(B11, j, i, cor);
    }
  }  
  // Init B11 Inverse matrix
  MakePosDef(B11, args.min_abs_eig);
  InvMat(B11Inv, B11);
  
  // Cal Z-score and Info
  Rcpp::Rcout<<"Imputing summary statistics of unmeasured SNPs..."<<std::endl;
  double prog_prev = 0;
  double prog = 0;
  for(size_t i=0; i<num_unmeasured; i++){
    double stdi = gsl_vector_get(SNP_STD_VEC, i + num_measured);
    for(size_t j=0; j<num_measured; j++){
      double stdj = gsl_vector_get(SNP_STD_VEC, j);
      double cov = CalWgtCov((*sliding_window_unmeasured[i]).GetGenotypeVec(),
                             (*sliding_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov/(stdi*stdj);
      gsl_matrix_set(b21, 0, j, cor);
    }
    gsl_matrix_transpose_memcpy(b12, b21);
    MpMatMat(b21B11Inv, b21, B11Inv);
    MpMatMat(val, b21B11Inv, Z1); //z2
    double z = gsl_matrix_get(val, 0, 0);
    //double pval = 2*gsl_cdf_ugaussian_Q(std::abs(z));
    MpMatMat(val, b21B11Inv, b12); // information of z2
    double info = std::abs(gsl_matrix_get(val, 0, 0));
    
    (*sliding_window_unmeasured[i]).SetZ(z/std::sqrt(info)); // use normalized z2 (imputed z-score)
    //(*sliding_window_unmeasured[i]).SetPval(pval); //pvalue
    (*sliding_window_unmeasured[i]).SetInfo(info); //info
    
    prog = 100*((i+1)/(double)(num_unmeasured));
    if(prog - prog_prev >= 1){
      int percent = static_cast<int>(prog);
      LoadProgressBar(percent);
      prog_prev++;
    }
  }
  
  gsl_matrix_free(Z1);
  gsl_vector_free(SNP_STD_VEC);
  gsl_matrix_free(B11);
  gsl_matrix_free(B11Inv);

  gsl_matrix_free(b21);
  gsl_matrix_free(b12);
  gsl_matrix_free(b21B11Inv);
  gsl_matrix_free(val);
  //////////////////////////
  // Imputation is done ! //
  //////////////////////////
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"Chromosome " <<args.chr<<" "<<args.start_bp<<"-"<<args.end_bp<<" locus successfully imputed!"<<std::endl; 		
  Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured<<std::endl;
  Rcpp::Rcout<<"Number of imputed SNPs: "<<num_unmeasured<<std::endl;

  //Delete sliding_window_measured.
  sliding_window_measured.clear();
  std::deque<Snp*>().swap(sliding_window_measured);

  //Delete sliding_window_unmeasured.
  sliding_window_unmeasured.clear();
  std::deque<Snp*>().swap(sliding_window_unmeasured);
}