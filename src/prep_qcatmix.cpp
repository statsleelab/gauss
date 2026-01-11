#include <Rcpp.h>

using namespace Rcpp;

#include <string>
#include <vector>
#include <map>
#include <deque>
#include "snp.h"
#include "gauss.h"
#include "util.h"
//#include <gsl/gsl_matrix.h> # delete gsl headers
//#include <gsl/gsl_cdf.h> # delete gsl headers


//' Prepare datasets for Recessive Imputation
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
//' @return A List containing:
//'         - snplist: A data frame of SNPs in the prediction window with columns rsid, chr, bp, a1, a2, af1mix, z, and type,
//'         - zvec: A numeric vector of Z-scores for measured SNPs in the extended window,
//'         - cormat: A correlation matrix (additive-coded) among measured SNPs in the extended window,
//'         - cormat_add: A correlation matrix between additive-coded predicted SNPs and additive-coded measured SNPs,
//'         - cormat_dom: A correlation matrix between dominant-coded predicted SNPs and additive-coded measured SNPs,
//'         - cormat_rec: A correlation matrix between recessive-coded predicted SNPs and additive-coded measured SNPs.
// [[Rcpp::export]]
List prep_recessive_impute(int chr, 
                           long long int start_bp, 
                           long long int end_bp, 
                           long long int wing_size,
                           DataFrame pop_wgt_df,
                           std::string input_file, 
                           std::string reference_index_file, 
                           std::string reference_data_file,
                           std::string reference_pop_desc_file,
                           Rcpp::Nullable<double> af1_cutoff = R_NilValue){
  
  // Initialize an Arguments object to hold the input parameters
  Arguments args;
  args.chr = chr; 
  args.start_bp = start_bp;
  args.end_bp = end_bp;
  args.wing_size = wing_size;
  
  // Convert population ID and weight data from the R DataFrame (pop_wgt_df) to C++ vectors
  std::vector<std::string> pop_vec_in = as<std::vector<std::string>>(pop_wgt_df[0]); // First column: Population IDs
  std::vector<double> pop_wgt_vec_in = as<std::vector<double>>(pop_wgt_df[1]); // Second column: Population Weights
  
  // Store population weights in args.pop_wgt_map, converting population IDs to uppercase
  for(int i=0; i<pop_vec_in.size(); i++){
    std::string pop = pop_vec_in[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper); // Convert ID to uppercase
    args.pop_wgt_map[pop] = pop_wgt_vec_in[i]; // Map population ID to its corresponding weight
  }
  
  // Set file paths for input and reference data in the Arguments object
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  
  // Set the allele frequency cutoff (af1_cutoff) if provided, otherwise use the default value (0.01)
  if(af1_cutoff.isNotNull()) {
    args.af1_cutoff = Rcpp::as<double>(af1_cutoff); // Convert R value to C++ double
  } else {
    args.af1_cutoff = 0.01; // Default value
  }
  
  // Read population description data and initialize population flag/weight vectors
  read_ref_desc(args);  // Load population description from the reference file
  init_pop_flag_wgt_vec(args); // Initialize population flags and weights
  
  // Initialize an empty map to store SNPs (Snp objects) indexed by a unique MapKey
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  
  // Read GWAS summary statistics (SNP Z-scores) into snp_map
  ReadInputZ(snp_map, args, false); 
  
  // Read the reference panel index file, updating snp_map with reference panel information
  ReadReferenceIndex(snp_map, args); 
  
  // Create an empty vector to store SNP pointers
  std::vector<Snp*> snp_vec;
  
  // Compute weighted reference allele frequencies (af1_mix) and store SNPs that meet the cutoff criteria in snp_vec
  MakeSnpVecMix(snp_vec, snp_map, args);
  
  // Read genotype data for the selected SNPs in snp_vec
  ReadGenotype(snp_vec, args);
  
  // Update each SNP in snp_vec so that the minor allele becomes the reference allele
  UpdateSnpToMinorAllele(snp_vec);
  
  
  /*----------------------------------------------------*/
  // Partition SNPs into two windows: 
  //  - ext_window_measured: measured SNPs (type==1) in the extended window
  //  - pred_window_all: all SNPs (excluding type==2) in the prediction window
  std::deque<Snp*> ext_window_measured;    // Measured SNPs in the extended window
  std::deque<Snp*> pred_window_all;        // All SNPs in the prediction window
  
  // The prediction window refers to the genomic region where the imputation will be conducted.
  // The extended window encompasses the prediction window, with additional buffer regions (wings) on both sides.
  // For instance, if the prediction window spans 1 MB and the wing size is 0.1 MB, the total size of the extended window
  // becomes 1.2 MB, with 0.1 MB of flanking regions added on either side of the prediction window. 
  // These extra regions help capture additional SNPs that may provide useful information for the imputation process.
  
  
  // Iterate through the SNP vector and categorize SNPs into measured or predicted windows
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    long long int bp = (*it_sv)->GetBp();
    
    // Add all SNPs within the prediction window to the prediction deque
    if((type != 2) && (bp >= args.start_bp && bp <= args.end_bp)){
      pred_window_all.push_back(*it_sv);
    }
    // Add measured SNPs to the extended window deque
    if(type == 1) { 
      ext_window_measured.push_back(*it_sv);
    }
  }
  
  // Count SNPs in each window
  int num_measured_ext = ext_window_measured.size();   // Number of measured SNPs in the extended window
  int num_all_pred = pred_window_all.size();           // Number of all SNPs in the prediction window
  
  // Print out the number of measured and predicted SNPs for reference
  Rcpp::Rcout<<"Number of measured SNPs in the extended window: "<<num_measured_ext<<std::endl;
  Rcpp::Rcout<<"Number of all SNPs in the prediction window: "<<num_all_pred<<std::endl;
  
  // If there are not enough measured SNPs, stop the process
  if(num_measured_ext <= args.min_num_measured_snp){
    Rcpp::stop("Not enough number of SNPs loaded - Recessive Imputation not performed");
  }
  
  // Z-scores for measured SNPs in the extended window
  NumericVector Z_vec(num_measured_ext);
  for(size_t i = 0; i < num_measured_ext; i++){   
    Z_vec(i) = (*ext_window_measured[i]).GetZ();
  }
  
  // Compute standard deviations for measured SNPs in the extended window
  NumericVector SNP_Std_Measured_Ext(num_measured_ext);
  for(size_t i = 0; i < num_measured_ext; i++){
    double v = CalWgtCov((*ext_window_measured[i]).GetGenotypeVec(), 
                         (*ext_window_measured[i]).GetGenotypeVec(), args.pop_wgt_vec);
    SNP_Std_Measured_Ext(i) = std::sqrt(v);
  }
  
  // Compute standard deviations for all additive-coded SNPs in the prediction window
  NumericVector SNP_Std_All_Pred(num_all_pred);
  for(size_t i = 0; i < num_all_pred; i++){
    double v = CalWgtCov((*pred_window_all[i]).GetGenotypeVec(), 
                         (*pred_window_all[i]).GetGenotypeVec(), args.pop_wgt_vec);
    SNP_Std_All_Pred(i) = std::sqrt(v);
  }  
  
  // ---- Compute cormat: Correlation matrix among measured SNPs in the extended window (additive-coded)
  NumericMatrix cormat(num_measured_ext, num_measured_ext);
  Rcpp::Rcout << "Computing correlations between additive-coded SNPs..." << std::endl;
  for(size_t i = 0; i < num_measured_ext; i++){
    cormat(i, i) = 1.0;
    double stdi = SNP_Std_Measured_Ext(i);
    for(size_t j = i + 1; j < num_measured_ext; j++){
      double stdj = SNP_Std_Measured_Ext(j);
      double cov = CalWgtCov((*ext_window_measured[i]).GetGenotypeVec(), 
                             (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov / (stdi * stdj);
      cormat(i, j) = cor;
      cormat(j, i) = cor;
    }
  }  
  
  
  // ---- Compute predicted SNP coding conversions and corresponding standard deviations ----
  // For additive coding, we use the genotype vector as is.
  // For dominant coding, convert using ConvertGenotypesToDominant.
  // For recessive coding, convert using ConvertGenotypesToRecessive.
  
  // Already computed standard deviations for additive-coded SNPs in pred win : SNP_Std_All_Pred
  
  // Compute standard deviations for all dominant-coded SNPs in the prediction window
  NumericVector SNP_Std_All_Pred_Dominant(num_all_pred);
  std::vector< std::vector<std::string> > pred_dominant(num_all_pred);
  for(size_t i = 0; i < num_all_pred; i++){
    pred_dominant[i] = ConvertGenotypesToDominant((*pred_window_all[i]).GetGenotypeVec());
    double v = CalWgtCov(pred_dominant[i], pred_dominant[i], args.pop_wgt_vec);
    SNP_Std_All_Pred_Dominant(i) = std::sqrt(v);
  }
  
  // Compute standard deviations for all recessive-coded SNPs in the prediction window
  NumericVector SNP_Std_All_Pred_Recessive(num_all_pred);
  std::vector< std::vector<std::string> > pred_recessive(num_all_pred);
  for(size_t i = 0; i < num_all_pred; i++){
    pred_recessive[i] = ConvertGenotypesToRecessive((*pred_window_all[i]).GetGenotypeVec());
    double v = CalWgtCov(pred_recessive[i], pred_recessive[i], args.pop_wgt_vec);
    SNP_Std_All_Pred_Recessive(i) = std::sqrt(v);
  }
  
  // ---- Compute cormat.add: Correlation matrix between additive-coded SNPs (in pred) vs. additive-coded measured SNPs (in ext) ----
  NumericMatrix cormat_add(num_all_pred, num_measured_ext);
  Rcpp::Rcout << "Computing correlations between additive-coded SNPs (pred) and additive-coded measured SNPs (ext)..." << std::endl;
  for (size_t i = 0; i < num_all_pred; i++){
    double stdi = SNP_Std_All_Pred(i);
    for (size_t j = 0; j < num_measured_ext; j++){
      double stdj = SNP_Std_Measured_Ext(j);
      double cov = CalWgtCov((*pred_window_all[i]).GetGenotypeVec(), 
                             (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov / (stdi * stdj);
      cormat_add(i, j) = cor;
    }
  }
  
  // ---- Compute cormat.dom: Correlation matrix between dominant-coded SNPs (in pred) vs. additive-coded measured SNPs (in ext) ----
  NumericMatrix cormat_dom(num_all_pred, num_measured_ext);
  Rcpp::Rcout << "Computing correlations between dominant-coded SNPs (pred) and additive-coded measured SNPs (ext)..." << std::endl;
  for (size_t i = 0; i < num_all_pred; i++){
    double stdi = SNP_Std_All_Pred_Dominant(i);
    for (size_t j = 0; j < num_measured_ext; j++){
      double stdj = SNP_Std_Measured_Ext(j);
      double cov = CalWgtCov(pred_dominant[i], (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov / (stdi * stdj);
      cormat_dom(i, j) = cor;
    }
  }
  
  // ---- Compute cormat.rec: Correlation matrix between recessive-coded SNPs (in pred) vs. additive-coded measured SNPs (in ext) ----
  NumericMatrix cormat_rec(num_all_pred, num_measured_ext);
  Rcpp::Rcout << "Computing correlations between recessive-coded SNPs (pred) and additive-coded measured SNPs (ext)..." << std::endl;
  for (size_t i = 0; i < num_all_pred; i++){
    double stdi = SNP_Std_All_Pred_Recessive(i);
    for (size_t j = 0; j < num_measured_ext; j++){
      double stdj = SNP_Std_Measured_Ext(j);
      double cov = CalWgtCov(pred_recessive[i], (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
      double cor = cov / (stdi * stdj);
      cormat_rec(i, j) = cor;
    }
  }
  
  
  // Build a data frame from SNPs in the prediction window.
  StringVector rsid_vec, a1_vec, a2_vec;
  IntegerVector chr_vec, bp_vec, type_vec;
  NumericVector af1mix_vec, z_vec_pred;
  Rcpp::Rcout << "Building SNP data frame..." << std::endl;
  for(std::deque<Snp*>::iterator it = pred_window_all.begin(); it != pred_window_all.end(); ++it){
    rsid_vec.push_back((*it)->GetRsid());
    chr_vec.push_back((*it)->GetChr());
    bp_vec.push_back((*it)->GetBp());
    a1_vec.push_back((*it)->GetA1());
    a2_vec.push_back((*it)->GetA2());
    af1mix_vec.push_back((*it)->GetAf1Mix());
    z_vec_pred.push_back((*it)->GetZ());
    type_vec.push_back((*it)->GetType());
  }
  
  DataFrame df = DataFrame::create(Named("rsid") = rsid_vec,
                                   Named("chr") = chr_vec,
                                   Named("bp") = bp_vec,
                                   Named("a1") = a1_vec,
                                   Named("a2") = a2_vec,
                                   Named("af1mix") = af1mix_vec,
                                   Named("z") = z_vec_pred,
                                   Named("type") = type_vec);
  
  // Clear the sliding windows
  ext_window_measured.clear();
  std::deque<Snp*>().swap(ext_window_measured);
  pred_window_all.clear();
  std::deque<Snp*>().swap(pred_window_all);
  
  // Release memory allocated for genotype vectors
  Rcpp::Rcout << "Releasing genotype memory..." << std::endl;
  FreeGenotype(snp_vec);
  
  // Delete SNP objects and clear the SNP map
  Rcpp::Rcout << "Deleting SNP map..." << std::endl;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end(); ){
    (it_sm->second)->ClearSnp();
    delete it_sm->second;
    snp_map.erase(it_sm++);
  }
  
  Rcpp::Rcout << "Returning results..." << std::endl;
  return List::create(Named("snplist") = df,         // Data frame of SNPs in the prediction window
                      Named("zvec") = Z_vec,           // Z-score vector of measured SNPs in the extended window
                      Named("cormat") = cormat,        // Correlation matrix among measured SNPs (additive) in the extended window
                      Named("cormat_add") = cormat_add,// Correlation matrix: additive-coded SNPs (in prediction window) vs. additive-coded measured SNPs (in extended window)
                      Named("cormat_dom") = cormat_dom,// Correlation matrix: dominant-coded SNPs (in prediction window) vs. additive-coded measured SNPs (in extended window)
                      Named("cormat_rec") = cormat_rec // Correlation matrix: recessive-coded SNPs (in prediction window) vs. additive-coded measured SNPs (in extended window)
  );
}


/*
 //' Prepare datasets for QCATMIX analysis
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
 List prep_qcatmix(int chr, 
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
 
 
 // prep data for QCAT
 
 //----------------------------------------------------//
 std::deque<Snp*> ext_window_measured;    //stores measured SNPs in the extended window
 std::deque<Snp*> pred_window_all;        //stores all SNPs in the prediction window
 
 for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
 int type = (*it_sv)->GetType();
 long long int bp = (*it_sv)->GetBp();
 if((type!=2)&(bp >= args.start_bp && bp <= args.end_bp)){ // all SNPs in the pred win. 
 pred_window_all.push_back(*it_sv);
 }
 if(type == 1) { // measured
 ext_window_measured.push_back(*it_sv);
 } // if (type == 2) don't put the snp in the sliding window. type=2: measured SNP but not exist in rep. panel
 }
 
 int num_measured_ext = ext_window_measured.size();      // # of measured SNPs in ext win
 int num_all_pred = pred_window_all.size();       // # of all SNPs in pred win
 
 Rcpp::Rcout<<"Number of measured SNPs in the ext window: "<<num_measured_ext<<std::endl;
 Rcpp::Rcout<<"Number of all SNPs in the pred window: "<<num_all_pred<<std::endl;
 
 if(num_measured_ext <= args.min_num_measured_snp){
 Rcpp::Rcout<<std::endl;
 Rcpp::Rcout<<"Number of measured SNPs in the ext window: "<<num_measured_ext<<std::endl;
 Rcpp::Rcout<<"Number of all SNPs in the pred window: "<<num_all_pred<<std::endl;
 Rcpp::stop("Not enough number of SNPs loaded - QCAT not performed");
 }
 
 //////////////////////////
 // run QCATMIX analysis //
 //////////////////////////
 //gsl_matrix* Z1 = gsl_matrix_calloc(num_measured_ext, 1);
 //gsl_vector* SNP_STD_VEC = gsl_vector_calloc(num_measured_ext + num_unmeasured_pred); //vector of SNP genotype standard deviations
 //gsl_matrix* B11 = gsl_matrix_calloc(num_measured_ext,num_measured_ext); // correlation matrix B11 : correlation among measured SNPs 
 //gsl_matrix* B21 = gsl_matrix_calloc(num_unmeasured_pred, num_measured_ext); // correlation matrix B21 : correlation btw measured and unmeasured SNPs
 
 NumericVector Z1(num_measured_ext);
 NumericVector SNP_Std_Measured_Ext(num_measured_ext);
 NumericVector SNP_Std_All_Pred(num_all_pred);
 NumericMatrix B11(num_measured_ext,num_measured_ext);
 NumericMatrix B21(num_all_pred,num_measured_ext);
 
 // Init Z1 vector
 for(size_t i=0; i<num_measured_ext; i++){   
 Z1(i) = (*ext_window_measured[i]).GetZ();
 }
 // Init SNP_Std_Measured_Ext	
 for(size_t i=0; i < num_measured_ext; i++){
 double v = CalWgtCov((*ext_window_measured[i]).GetGenotypeVec(), (*ext_window_measured[i]).GetGenotypeVec(), args.pop_wgt_vec);
 SNP_Std_Measured_Ext(i) = std::sqrt(v);
 }
 // Init SNP_Std_All_Pred
 for(size_t i=0; i < num_all_pred; i++){
 double v = CalWgtCov((*pred_window_all[i]).GetGenotypeVec(), (*pred_window_all[i]).GetGenotypeVec(), args.pop_wgt_vec);
 SNP_Std_All_Pred(i) = std::sqrt(v);
 }
 // Init B11 matrix
 Rcpp::Rcout<<"Computing correlations between variants... B11"<<std::endl;
 for(size_t i=0; i<num_measured_ext; i++){
 B11(i,i) = 1.0; //diagonals
 double stdi = SNP_Std_Measured_Ext(i);
 for(size_t j=i+1; j<num_measured_ext; j++){
 double stdj = SNP_Std_Measured_Ext(j);
 double cov = CalWgtCov((*ext_window_measured[i]).GetGenotypeVec(), 
 (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
 double cor = cov/(stdi*stdj);
 B11(i,j) = cor;
 B11(j,i) = cor;
 }
 }  
 // Init B21 matrix
 Rcpp::Rcout<<"Computing correlations between variants... B21"<<std::endl;
 for(size_t i=0; i<num_all_pred; i++){
 double stdi = SNP_Std_All_Pred(i);
 for(size_t j=0; j<num_measured_ext; j++){
 double stdj = SNP_Std_Measured_Ext(j);
 double cov = CalWgtCov((*pred_window_all[i]).GetGenotypeVec(),
 (*ext_window_measured[j]).GetGenotypeVec(), args.pop_wgt_vec);
 double cor = cov/(stdi*stdj);
 B21(i,j) = cor;
 }
 }
 
 StringVector rsid_vec;
 IntegerVector chr_vec;
 IntegerVector bp_vec;
 StringVector a1_vec;
 StringVector a2_vec;
 NumericVector af1mix_vec;
 NumericVector z_vec;
 IntegerVector type_vec;
 
 Rcpp::Rcout<<"push_vec"<<std::endl;
 for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
 rsid_vec.push_back((*it_sv)->GetRsid());
 chr_vec.push_back((*it_sv)->GetChr());
 bp_vec.push_back((*it_sv)->GetBp());
 a1_vec.push_back((*it_sv)->GetA1());
 a2_vec.push_back((*it_sv)->GetA2());
 af1mix_vec.push_back((*it_sv)->GetAf1Mix());
 z_vec.push_back((*it_sv)->GetZ());
 type_vec.push_back((*it_sv)->GetType());
 }
 
 Rcpp::Rcout<<"rsid:   "<<rsid_vec.length()<<std::endl;
 Rcpp::Rcout<<"chr :   "<<chr_vec.length()<<std::endl;
 Rcpp::Rcout<<"bp  :   "<<bp_vec.length()<<std::endl;
 Rcpp::Rcout<<"a1  :   "<<a1_vec.length()<<std::endl;
 Rcpp::Rcout<<"a2  :   "<<a2_vec.length()<<std::endl;
 Rcpp::Rcout<<"af1mix: "<<af1mix_vec.length()<<std::endl;
 Rcpp::Rcout<<"z   :   "<<z_vec.length()<<std::endl;
 Rcpp::Rcout<<"type:   "<<type_vec.length()<<std::endl;
 
 Rcpp::Rcout<<"Make dataframe"<<std::endl;
 DataFrame df = DataFrame::create(Named("rsid")=rsid_vec,
 Named("chr")=chr_vec,
 Named("bp")=bp_vec,
 Named("a1")=a1_vec,
 Named("a2")=a2_vec,
 Named("af1mix")=af1mix_vec,
 Named("z")=z_vec,
 Named("type")=type_vec);  
 
 //Delete ext_window_measured.
 ext_window_measured.clear();
 std::deque<Snp*>().swap(ext_window_measured);
 
 //Delete pred_window_all.
 pred_window_all.clear();
 std::deque<Snp*>().swap(pred_window_all);
 
 //----------------------------------------------------//
 // release memory allocated for genotype
 Rcpp::Rcout<<"release memory allocated for genotype"<<std::endl;
 FreeGenotype(snp_vec); 
 
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
 */
