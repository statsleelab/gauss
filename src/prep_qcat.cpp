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

//void run_qcatR(std::vector<Snp*>& snp_vec, Arguments& args);

//' Prepare datasets for QCAT analysis
//' 
//' @param chr chromosome number
//' @param start_bp start base pair position of prediction window
//' @param end_bp end base pair position of prediction window
//' @param wing_size the size of the area flanking the left and right of the prediction window
//' @param study_pop study population group
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2, af1, and z
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R dataframe containing rsid, chr, bp, a1, a2, af1ref, z, qcat_chisq, qcat_pval, type
// [[Rcpp::export]]
List prep_qcat(int chr, 
               long long int start_bp, 
               long long int end_bp, 
               long long int wing_size,
               std::string study_pop,
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
  args.study_pop = study_pop;
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
  init_pop_flag_vec(args);
  //args.PrintArguments();
  
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  ReadInputZ(snp_map, args, false);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  ReadReferenceIndex(snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  
  // make a snp vector containing all SNPs
  std::vector<Snp*> snp_vec;
  MakeSnpVec(snp_vec, snp_map, args);
  ReadGenotype(snp_vec, args);
  
  // prep data for QCAT
  /*----------------------------------------------------*/
  std::deque<Snp*> sliding_window_measured_ext; //stores measured SNPs in the extended window
  std::deque<Snp*> sliding_window_all_pred;     //stores all SNPs in the prediction window
  
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    long long int bp = (*it_sv)->GetBp();
    if((type!=2)&(bp >= args.start_bp && bp <= args.end_bp)){ // all SNPs in the pred win. 
      sliding_window_all_pred.push_back(*it_sv);
    } else if(type == 1) { // measured
      sliding_window_measured_ext.push_back(*it_sv);
    } // if (type == 2) don't put the snp in the sliding window. type=2: measured SNP but not exist in rep. panel
  }
  
  int num_measured_ext = sliding_window_measured_ext.size();// # of measured SNPs in ext win
  int num_all_pred = sliding_window_all_pred.size();        // # of all SNPs in pred win

  Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured_ext<<std::endl;
  Rcpp::Rcout<<"Number of all SNPs in the prediction window: "<<num_all_pred<<std::endl;
    
  if(num_measured_ext <= args.min_num_measured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured_ext<<std::endl;
    Rcpp::Rcout<<"Number of all SNPs in the prediction window: "<<num_all_pred<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - QCAT not performed");
  }
  
  /*Prep data for QCAT*/
  //gsl_matrix* Z1 = gsl_matrix_calloc(num_measured_ext, 1);
  //gsl_matrix* B11 = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // LD among measured SNPs in ext win
  //gsl_matrix* B21 = gsl_matrix_calloc(num_all_pred, num_measured_ext); // LD btw all SNPs in pred win and measured SNPs in ext win
  NumericVector Z1(num_measured_ext);
  NumericMatrix B11(num_measured_ext,num_measured_ext);
  NumericMatrix B21(num_all_pred,num_measured_ext);
  
  // Init Z1 vector
  for(size_t i=0; i<num_measured_ext; i++){   
    Z1(i) = (*sliding_window_measured_ext[i]).GetZ();
  }
  
  // Init B11 matrix
  Rcpp::Rcout<<"Computing correlations between variants... B11"<<std::endl;
  for(size_t i=0; i<num_measured_ext; i++){
    B11(i,i) = 1.0; //diagonals
    for(size_t j=i+1; j<num_measured_ext; j++){
      double v = CalCor((*sliding_window_measured_ext[i]).GetGenotypeVec(), (*sliding_window_measured_ext[j]).GetGenotypeVec());
      B11(i,j) = v;
      B11(j,i) = v;
    }
  }
  
  // Init B21 matrix
  Rcpp::Rcout<<"Computing correlations between variants... B21"<<std::endl;
  for(size_t i=0; i<num_all_pred; i++){
    for(size_t j=0; j<num_measured_ext; j++){
      double v = CalCor((*sliding_window_all_pred[i]).GetGenotypeVec(),
                        (*sliding_window_measured_ext[j]).GetGenotypeVec());
      B21(i,j) = v;
    }
  }
  
  //Delete sliding_window_measured_ext.
  sliding_window_measured_ext.clear();
  std::deque<Snp*>().swap(sliding_window_measured_ext);
  
  //Delete sliding_window_all_pred.
  sliding_window_all_pred.clear();
  std::deque<Snp*>().swap(sliding_window_all_pred);
  
  
  /*----------------------------------------------------*/
  // release memory allocated for genotype
  Rcpp::Rcout<<"release memory allocated for genotype"<<std::endl;
  FreeGenotype(snp_vec); 
  
  StringVector rsid_vec;
  IntegerVector chr_vec;
  IntegerVector bp_vec;
  StringVector a1_vec;
  StringVector a2_vec;
  NumericVector af1ref_vec;
  NumericVector z_vec;
  IntegerVector type_vec;
  
  Rcpp::Rcout<<"push_vec"<<std::endl;
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int bp = (*it_sv)->GetBp();
    if(bp >= start_bp && bp <= end_bp){
      rsid_vec.push_back((*it_sv)->GetRsid());
      chr_vec.push_back((*it_sv)->GetChr());
      bp_vec.push_back((*it_sv)->GetBp());
      a1_vec.push_back((*it_sv)->GetA1());
      a2_vec.push_back((*it_sv)->GetA2());
      af1ref_vec.push_back((*it_sv)->GetAf1Ref());
      z_vec.push_back((*it_sv)->GetZ());
      type_vec.push_back((*it_sv)->GetType());
    }
  }
  
  Rcpp::Rcout<<"rsid:   "<<rsid_vec.length()<<std::endl;
  Rcpp::Rcout<<"chr :   "<<chr_vec.length()<<std::endl;
  Rcpp::Rcout<<"bp  :   "<<bp_vec.length()<<std::endl;
  Rcpp::Rcout<<"a1  :   "<<a1_vec.length()<<std::endl;
  Rcpp::Rcout<<"a2  :   "<<a2_vec.length()<<std::endl;
  Rcpp::Rcout<<"af1ref: "<<af1ref_vec.length()<<std::endl;
  Rcpp::Rcout<<"z   :   "<<z_vec.length()<<std::endl;
  Rcpp::Rcout<<"type:   "<<type_vec.length()<<std::endl;
  
  Rcpp::Rcout<<"Make dataframe"<<std::endl;
  DataFrame df = DataFrame::create(Named("rsid")=rsid_vec,
                                   Named("chr")=chr_vec,
                                   Named("bp")=bp_vec,
                                   Named("a1")=a1_vec,
                                   Named("a2")=a2_vec,
                                   Named("af1ref")=af1ref_vec,
                                   Named("z")=z_vec,
                                   Named("type")=type_vec);
  
  
  //deletes snp_map.
  Rcpp::Rcout<<"deletes snp map"<<std::endl;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end();){
    (it_sm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_sm->second;        // delete snp object
    snp_map.erase(it_sm++);      // delete map element
  }
  
  Rcpp::Rcout<<"return"<<std::endl;
  return List::create(Named("snplist")=df,
                      Named("z_vec") = Z1,
                      Named("cor_mat1") = B11,
                      Named("cor_mat2") = B21);
}



/*
void run_qcatR(std::vector<Snp*>& snp_vec, Arguments& args){
  
  std::deque<Snp*> sliding_window_measured_ext; //stores measured SNPs in the extended window
  std::deque<Snp*> sliding_window_all_pred;     //stores unmeasured SNPs in the prediction window
  
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    long long int bp = (*it_sv)->GetBp();
    if((type!=2)&(bp >= args.start_bp && bp <= args.end_bp)){ // all SNPs in the pred win. 
      sliding_window_all_pred.push_back(*it_sv);
    } else if(type == 1) { // measured
      sliding_window_measured_ext.push_back(*it_sv);
    } // if (type == 2) don't put the snp in the sliding window. type=2: measured SNP but not exist in rep. panel
  }
  
  int num_measured_ext = sliding_window_measured_ext.size();// # of measured SNPs in ext win
  int num_all_pred = sliding_window_all_pred.size();        // # of all SNPs in pred win
  
  if(num_measured_ext <= args.min_num_measured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured_ext<<std::endl;
    Rcpp::Rcout<<"Number of all SNPs in the prediction window: "<<num_all_pred<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - QCAT not performed");
  }
  
  //Prep data for QCAT
  gsl_matrix* Z1 = gsl_matrix_calloc(num_measured_ext, 1);
  gsl_matrix* B11 = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // LD among measured SNPs in ext win
  gsl_matrix* B21 = gsl_matrix_calloc(num_all_pred, num_measured_ext); // LD btw all SNPs in pred win and measured SNPs in ext win
  
  // Init Z1 vector
  for(size_t i=0; i<num_measured_ext; i++){   
    gsl_matrix_set(Z1, i, 0, (*sliding_window_measured_ext[i]).GetZ());
  }
  
  // Init B11 matrix
  Rcpp::Rcout<<"Computing correlations between variants..."<<std::endl;
  for(size_t i=0; i<B11->size1; i++){
    gsl_matrix_set(B11, i, i, 1.0 + args.lambda); //add LAMBDA here (ridge regression trick)
    for(size_t j=i+1; j<B11->size1; j++){
      double v = CalCor((*sliding_window_measured_ext[i]).GetGenotypeVec(), (*sliding_window_measured_ext[j]).GetGenotypeVec());
      gsl_matrix_set(B11, i, j, v);
      gsl_matrix_set(B11, j, i, v);
    }
  }
  
  // Init B21 matrix
  for(size_t i=0; i<num_all_pred; i++){
    for(size_t j=0; j<num_measured_ext; j++){
      double v = CalCor((*sliding_window_all_pred[i]).GetGenotypeVec(),
                             (*sliding_window_measured_ext[j]).GetGenotypeVec());
      gsl_matrix_set(B21, i, j, v);
    }
  }

  gsl_matrix_free(Z1);
  gsl_matrix_free(B11);
  gsl_matrix_free(B21);
  
  //Delete sliding_window_measured_ext.
  sliding_window_measured_ext.clear();
  std::deque<Snp*>().swap(sliding_window_measured_ext);
  
  //Delete sliding_window_all_pred.
  sliding_window_all_pred.clear();
  std::deque<Snp*>().swap(sliding_window_all_pred);
}
*/
