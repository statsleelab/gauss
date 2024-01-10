#include <Rcpp.h>

using namespace Rcpp;

#include <string>
#include <vector>
#include <map>
#include <deque>
#include <random>

#include "snp.h"
#include "gauss.h"
#include "util.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>


double CalCor(NumericVector& x, NumericVector& y);

//' Simulate LD for user provided SNPs from mixed ethnicity cohorts
//' 
//' @param chr chromosome number
//' @param start_bp start base pair position of prediction window
//' @param end_bp end base pair position of prediction window
//' @param pop_wgt_df R data frame containing population IDs and weights
//' @param sim_size Number of simulated subjects
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2, af1, and z
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R list containing a data frame containing rsid, chr, bp, a1, a2 and af1mix and a correlation matrix 
// [[Rcpp::export]]
List simulateLD(int chr, 
                long long int start_bp, 
                long long int end_bp, 
                DataFrame pop_wgt_df,
                int sim_size,
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
  
  
  // Add pop_wgt_df info in args.pop_wgt_map
  // the code below is accessing the first column of pop_wgt_df, 
  // converting that column into a std::vector<std::string>, 
  // and storing the result in pop_vec_in. 
  std::vector<std::string> pop_vec_in = as<std::vector<std::string>>(pop_wgt_df[0]);
  std::vector<double> pop_wgt_vec_in = as<std::vector<double>>(pop_wgt_df[1]);
  std::vector<int> pop_num_sim_vec;
  for(int i=0; i<pop_vec_in.size(); i++){
    std::string pop = pop_vec_in[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper); //make capital
    args.pop_wgt_map[pop]=pop_wgt_vec_in[i];
    pop_num_sim_vec.push_back(static_cast<int>(pop_wgt_vec_in[i]*sim_size));
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
  // Simulate LD matrix
  /*------------------------------------------------------------*/
  
  // extract measured SNPs only
  std::vector<Snp*> snp_vec_measured;
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    if(type == 1) { // if measured
      snp_vec_measured.push_back(*it_sv);
    } // if (type == 2) don't put the snp in the sliding window. Do nothing. 
  }
  
  // read SNP genotypes
  ReadGenotype(snp_vec_measured, args);
  
  int num_measured = snp_vec_measured.size();
  if(num_measured <= args.min_num_measured_snp){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Number of measured SNPs: "<<num_measured<<std::endl;
    Rcpp::stop("Not enough number of SNPs loaded - computeLD not performed");
  }
  
  // To randomly draw with replacement a subject's genotypes from 
  // each combination of ethnicities, randomly select indexes 
  // of subjects in each ethnic group and store the information 
  // in geno_index_vec. The geno_index_vec will be used to simulate
  // genotypes
  
  std::random_device rd;
  std::mt19937 gen(rd());  
  std::vector<int> geno_index_vec;
  int num_pops = args.ref_pop_vec.size();
  int pop_counter=0;
  for(int k=0; k<num_pops; k++){
    if(args.pop_flag_vec[k]) {
      for(int j=0; j<pop_num_sim_vec[pop_counter]; j++){
        std::uniform_int_distribution<> dis(0, args.ref_pop_size_vec[k] - 1);
        int ran_index = dis(gen);
        geno_index_vec.push_back(ran_index);
      }
      pop_counter++;
    }
  }
  
  // Simulate genotypes using geno_index_vec
  NumericMatrix geno_mat(num_measured,sim_size);
  int subj_counter=0;
  
  for(size_t i=0; i<num_measured; i++){
    pop_counter=0;
    subj_counter=0;
    std::vector<std::string>& geno_vec = (*snp_vec_measured[i]).GetGenotypeVec();
    for(int k=0; k<num_pops; k++){
      if(args.pop_flag_vec[k]) {
        std::string geno_str = geno_vec[k];
        for(int j=0; j<pop_num_sim_vec[pop_counter]; j++){
          int ran_index = geno_index_vec[j+subj_counter];
          double geno = (double)(geno_str[ran_index] - '0');
          geno_mat(i,j+subj_counter)=geno;
        }
        subj_counter = subj_counter + pop_num_sim_vec[pop_counter];
        pop_counter++;
      }
    }
  }
  
  // Calculate LD
  NumericMatrix Cor_Mat(num_measured, num_measured);
  Rcpp::Rcout<<"Computing correlations between variants..."<<std::endl;
  for(size_t i=0; i<num_measured; i++){
    Cor_Mat(i,i) = 1.0; //diagonals
    for(size_t j=i+1; j<num_measured; j++){
      NumericVector snp_i = geno_mat(i, _);
      NumericVector snp_j = geno_mat(j, _);
      double cor = CalCor(snp_i, snp_j);
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


double CalCor(NumericVector& x, NumericVector& y){
  int n = x.size();
  double xi=0, yi=0, sumx=0, sumy=0, sumxsq=0, sumysq=0, sumxy=0;
  for(int i=0; i<n; i++){
    sumx += x[i];
    sumy += y[i];
    sumxsq += x[i]*x[i];
    sumysq += y[i]*y[i];
    sumxy += x[i]*y[i];
  }
  double numer = n*sumxy-sumx*sumy;
  double denor = std::sqrt((n)*sumxsq-sumx*sumx)*std::sqrt((n)*sumysq-sumy*sumy);
  double r = numer/denor;
  return r;
}
