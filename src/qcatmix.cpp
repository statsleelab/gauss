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

void run_qcatmix(std::vector<Snp*>& snp_vec, Arguments& args);

//' Testing causality of variants from mixed ethnicity cohorts
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
//' @return R dataframe containing rsid, chr, bp, a1, a2, af1mix, z, qcat_chisq, qcat_pval, type
// [[Rcpp::export]]
DataFrame qcatmix(int chr, 
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
    args.af1_cutoff = 0.05;
  }
  
  read_ref_desc(args);
  init_pop_flag_wgt_vec(args);
  //args.PrintArguments();
  
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  ReadInput(snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  ReadReferenceIndex(snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  
  // make a snp vector containing all SNPs
  std::vector<Snp*> snp_vec;
  MakeSnpVecMix(snp_vec, snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_vec.size() <<std::endl;
  ReadGenotype(snp_vec, args);

  // run QCATMIX 
  run_qcatmix(snp_vec, args);
  
  // release memory allocated for genotype
  FreeGenotype(snp_vec); 
  
  StringVector rsid_vec;
  IntegerVector chr_vec;
  IntegerVector bp_vec;
  StringVector a1_vec;
  StringVector a2_vec;
  NumericVector af1mix_vec;
  NumericVector z_vec;
  NumericVector qcat_chisq_vec;
  NumericVector qcat_pval_vec;
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
      qcat_chisq_vec.push_back((*it_sv)->GetQcatChisq());
      qcat_pval_vec.push_back(gsl_cdf_chisq_Q((*it_sv)->GetQcatChisq(), 1));
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
                                   Named("qcat_chisq")=qcat_chisq_vec,
                                   Named("qcat_pval")=qcat_pval_vec,
                                   Named("type")=type_vec);
  return df;  
}



void run_qcatmix(std::vector<Snp*>& snp_vec, Arguments& args){
  std::deque<Snp*> sliding_window_measured_ext;    //stores measured SNPs in the extended window
  std::deque<Snp*> sliding_window_unmeasured_pred; //stores unmeasured SNPs in the prediction window
  int num_measured_headwing = 0;
  int num_measured_pred = 0;
  
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    long long int bp = (*it_sv)->GetBp();
    if(type == 0 && (bp >= args.start_bp && bp <= args.end_bp)){ // unmeasured and in the pred win.
      sliding_window_unmeasured_pred.push_back(*it_sv);
    } else if(type == 1) { // measured
      sliding_window_measured_ext.push_back(*it_sv);
      if(bp < args.start_bp)
        num_measured_headwing++;
      else if(bp >= args.start_bp && bp <= args.end_bp)
        num_measured_pred++;
    } // if (type == 2) don't put the snp in the sliding window. Do nothing. 
  }
  
  int num_measured_ext = sliding_window_measured_ext.size();      // # of measured SNPs in ext win
  int num_unmeasured_pred = sliding_window_unmeasured_pred.size();// # of unmeasured SNPs in pred win
  
  if(num_measured_ext <= args.min_num_measured_snp || 
     num_unmeasured_pred <= args.min_num_unmeasured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured_ext<<std::endl;
    Rcpp::Rcout<<"Number of unmeasured SNPs: "<<num_unmeasured_pred<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - QCAT performed");
  }
  
  //////////////////////////
  // run QCATMIX analysis //
  //////////////////////////
  gsl_matrix* Z1 = gsl_matrix_calloc(num_measured_ext, 1);
  gsl_vector* SNP_STD_VEC = gsl_vector_calloc(num_measured_ext + num_unmeasured_pred); //vector of SNP genotype standard deviations
  gsl_matrix* B11 = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // correlation matrix B11 : correlation among measured SNPs 
  gsl_matrix* B21 = gsl_matrix_calloc(num_unmeasured_pred, num_measured_ext); // correlation matrix B21 : correlation btw measured and unmeasured SNPs
  
  gsl_matrix* L = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // L = square root (lower triangular matrix) of B11
  gsl_matrix* LInv = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // L Inverse
  gsl_matrix* LInvZ1 = gsl_matrix_calloc(num_measured_ext, 1); // L Inverse * Z1
  gsl_matrix* b11t = gsl_matrix_calloc(num_measured_ext, 1); // b11 (row vector of B11)
  gsl_matrix* LInvb11t = gsl_matrix_calloc(num_measured_ext, 1); // L Inverse * b11 (row vector of B11)
  gsl_matrix* b21t = gsl_matrix_calloc(num_measured_ext, 1); // b21 (row vector of B21)
  gsl_matrix* LInvb21t = gsl_matrix_calloc(num_measured_ext, 1); // L Inverse * b21 (row vector of B21)
  
  // Init Z1 vector
  for(size_t i=0; i<num_measured_ext; i++){   
    gsl_matrix_set(Z1, i, 0, (*sliding_window_measured_ext[i]).GetZ());
  }
  // Init SNP_STD_VEC	
  for(size_t i=0; i < num_measured_ext; i++){
    double v = CalWgtCov((*sliding_window_measured_ext[i]).GetGenotypeVec(), (*sliding_window_measured_ext[i]).GetGenotypeVec(), args.pop_wgt_vec);
    gsl_vector_set(SNP_STD_VEC, i, std::sqrt(v));
  }
  for(size_t i=0; i < num_unmeasured_pred; i++){
    double v = CalWgtCov((*sliding_window_unmeasured_pred[i]).GetGenotypeVec(), (*sliding_window_unmeasured_pred[i]).GetGenotypeVec(), args.pop_wgt_vec);
    gsl_vector_set(SNP_STD_VEC, i + num_measured_ext, std::sqrt(v));
  }
  // Init B11 matrix
  Rcpp::Rcout<<"Computing correlations between variants..."<<std::endl;
  for(size_t i=0; i<B11->size1; i++){
    gsl_matrix_set(B11, i, i, 1.0 + args.lambda); //add LAMBDA here (ridge regression trick)
    double stdi = gsl_vector_get(SNP_STD_VEC, i);
    for(size_t j=i+1; j<B11->size1; j++){
      double stdj = gsl_vector_get(SNP_STD_VEC, j);
      double cov = CalWgtCov((*sliding_window_measured_ext[i]).GetGenotypeVec(), (*sliding_window_measured_ext[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov/(stdi*stdj);
      gsl_matrix_set(B11, i, j, cor);
      gsl_matrix_set(B11, j, i, cor);
    }
  }  
  // Init B21 matrix
  for(size_t i=0; i<num_unmeasured_pred; i++){
    double stdi = gsl_vector_get(SNP_STD_VEC, i + num_measured_ext);
    for(size_t j=0; j<num_measured_ext; j++){
      double stdj = gsl_vector_get(SNP_STD_VEC, j);
      double cov = CalWgtCov((*sliding_window_unmeasured_pred[i]).GetGenotypeVec(),
                             (*sliding_window_measured_ext[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov/(stdi*stdj);
      gsl_matrix_set(B21, i, j, cor);
    }
  }
  
  //int num_eig = RmvPC(B11, args.eig_cutoff); //rmv PCs with eig value less than eig_cutoff
  int num_eig = CountPC(B11, args.eig_cutoff); //rmv PCs with eig value less than eig_cutoff
  CholeskyMat(L, B11); //Cholesky decomposion to get the square root (lower triangular matrix L) of B11.
  
  //MakePosDef(L, min_abs_eig_);
  InvMat(LInv, L);  
  MpMatMat(LInvZ1, LInv, Z1); // Cholesky transformation to uncorrelate Z1.
  
  //double varZ1 = CalVar(gsl_matrix_column(Z1,0));
  //double varLInvZ1 = CalVar(gsl_matrix_column(LInvZ1,0));  

  //Testing measured SNPs in the prediction window    
  Rcpp::Rcout<<"Testing measured variants..."<<std::endl;
  for(size_t i=0; i<num_measured_pred; i++){        
    gsl_matrix_view b11 = gsl_matrix_submatrix(B11, i+num_measured_headwing, 0, 1, num_measured_ext); //extract ith row vector b11 (as matrix image) from B11.
    gsl_matrix_transpose_memcpy(b11t, &b11.matrix);
    MpMatMat(LInvb11t, LInv, b11t); // Cholesky transformation to uncorrelate b11.
    double r = CalCor(gsl_matrix_column(LInvZ1,0), gsl_matrix_column(LInvb11t,0)); 
    //double r = CalCor(gsl_matrix_column(Z1,0), gsl_matrix_column(b11t,0));
    //double r = CalCor(gsl_matrix_column(LInvZ1,0), gsl_matrix_column(LInvb11t,0))*std::sqrt(varLInvZ1/(varLInvZ1+1));
    double chisq = (num_eig-3)*r*r;
    (*sliding_window_measured_ext[i+num_measured_headwing]).SetQcatChisq(chisq);
  }
  
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"Testing unmeasured variants..."<<std::endl;
  for(size_t i=0; i<num_unmeasured_pred; i++){  // for loop, unmeasured SNPs in the prediction window. 
    gsl_matrix_view b21 = gsl_matrix_submatrix(B21, i, 0, 1, num_measured_ext); //extract ith row vector b21 (as matrix image) from B21.
    gsl_matrix_transpose_memcpy(b21t, &b21.matrix);
    MpMatMat(LInvb21t, LInv, b21t); // Cholesky transformation to uncorrelate b21.
    double r = CalCor(gsl_matrix_column(LInvZ1,0), gsl_matrix_column(LInvb21t,0)); 
    //double r = CalCor(gsl_matrix_column(Z1,0), gsl_matrix_column(b21t,0));
    //double r = CalCor(gsl_matrix_column(LInvZ1,0), gsl_matrix_column(LInvb21t,0))*std::sqrt(varLInvZ1/(varLInvZ1+1));
    double chisq = (num_eig-3)*r*r;
    (*sliding_window_unmeasured_pred[i]).SetQcatChisq(chisq);
  }
  
  gsl_matrix_free(Z1);
  gsl_vector_free(SNP_STD_VEC);
  gsl_matrix_free(B11);
  gsl_matrix_free(B21);
  gsl_matrix_free(L);
  gsl_matrix_free(LInv);
  gsl_matrix_free(LInvZ1);
  gsl_matrix_free(b11t);
  gsl_matrix_free(LInvb11t);
  gsl_matrix_free(b21t);
  gsl_matrix_free(LInvb21t);
  
  ////////////////////
  // QCAT is done ! //
  ////////////////////
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"Chromosome " <<args.chr<<" "<<args.start_bp<<"-"<<args.end_bp<<" locus successfully tested!"<<std::endl; 		
  Rcpp::Rcout<<"Number of tested SNPs: "<<num_unmeasured_pred+num_measured_pred<<std::endl;

  //Delete sliding_window_measured_ext.
  sliding_window_measured_ext.clear();
  std::deque<Snp*>().swap(sliding_window_measured_ext);
  
  //Delete sliding_window_unmeasured_pred.
  sliding_window_unmeasured_pred.clear();
  std::deque<Snp*>().swap(sliding_window_unmeasured_pred);
}

