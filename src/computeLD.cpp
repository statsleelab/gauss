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

void cal_LD(std::vector<Snp*>& snp_vec, Arguments& args);

//' Compute LD for measured SNPs from mixed ethnicity cohorts
//' 
//' @param chr chromosome number
//' @param start_bp start base pair position of prediction window
//' @param end_bp end base pair position of prediction window
//' @param pop_wgt_df R data frame containing population IDs and weights
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2, af1, and z
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R list containing a data frame containing rsid, chr, bp, a1, a2 and af1mix and a correlation matrix 
// [[Rcpp::export]]
List computeLD(int chr, 
               long long int start_bp, 
               long long int end_bp, 
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
  args.wing_size = 0;
  
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

  /*------------------------------------------------------------*/
  // run calculate LD matrix
  std::vector<Snp*> snp_vec_measured;
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    if(type == 1) { // if measured
      snp_vec_measured.push_back(*it_sv);
    } // if (type == 2) don't put the snp in the sliding window. Do nothing. 
  }
  ReadGenotype(snp_vec_measured, args);
  
  int num_measured = snp_vec_measured.size();
  if(num_measured <= args.min_num_measured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - computeLD not performed");
  }
  
  // Calculate LD //
  NumericVector SNP_STD_VEC;
  NumericMatrix Cor_Mat(num_measured, num_measured);
  
  // Init SNP_STD_VEC	
  for(size_t i=0; i < num_measured; i++){
    double v = CalWgtCov((*snp_vec_measured[i]).GetGenotypeVec(), (*snp_vec_measured[i]).GetGenotypeVec(), args.pop_wgt_vec);
    SNP_STD_VEC.push_back(std::sqrt(v));
  }
  // Init B11 matrix
  Rcpp::Rcout<<"Computing correlations between variants..."<<std::endl;
  for(size_t i=0; i<num_measured; i++){
    Cor_Mat(i,i) = 1.0; //diagonals
    double stdi = SNP_STD_VEC(i);
    for(size_t j=i+1; j<num_measured; j++){
      double stdj = SNP_STD_VEC(j);
      double cov = CalWgtCov((*snp_vec_measured[i]).GetGenotypeVec(), (*snp_vec_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov/(stdi*stdj);
      Cor_Mat(i,j) = cor;
      Cor_Mat(j,i) = cor;
    }
  }  

  // LD computation is done ! //
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"Chromosome " <<args.chr<<" "<<args.start_bp<<"-"<<args.end_bp<<" locus LD computed!"<<std::endl; 		
  Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured<<std::endl;
  
  /*------------------------------------------------------------*/
  
  // release memory allocated for genotype
  FreeGenotype(snp_vec_measured);
  
  StringVector rsid_vec;
  IntegerVector chr_vec;
  IntegerVector bp_vec;
  StringVector a1_vec;
  StringVector a2_vec;
  NumericVector af1mix_vec;
  
  for(std::vector<Snp*>::iterator it_sv = snp_vec_measured.begin(); it_sv != snp_vec_measured.end(); ++it_sv){
    rsid_vec.push_back((*it_sv)->GetRsid());
    chr_vec.push_back((*it_sv)->GetChr());
    bp_vec.push_back((*it_sv)->GetBp());
    a1_vec.push_back((*it_sv)->GetA1());
    a2_vec.push_back((*it_sv)->GetA2());
    af1mix_vec.push_back((*it_sv)->GetAf1Mix());
  }
  DataFrame df = DataFrame::create(Named("rsid")=rsid_vec,
                                   Named("chr")=chr_vec,
                                   Named("bp")=bp_vec,
                                   Named("a1")=a1_vec,
                                   Named("a2")=a2_vec,
                                   Named("af1mix")=af1mix_vec);
  
  

  //Delete snp_vec_measured.
  snp_vec_measured.clear();
  std::vector<Snp*>().swap(snp_vec_measured);
  
  //deletes snp_map.
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end();){
    (it_sm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_sm->second;        // delete snp object
    snp_map.erase(it_sm++);      // delete map element
  }

  return List::create(Named("snplist") = df,
                      Named("cormat") = Cor_Mat);  
}
