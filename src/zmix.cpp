#include <Rcpp.h>

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include "util.h"
#include "bgzf.h"
#include "snp.h"
#include "gauss.h"

using namespace Rcpp;

//#define CPW2_Debug  

//Forward declaration
void read_input_zmix(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
void read_ref_index_zmix(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);


//' Calculate population weights using association Z-scores
//' 
//' @param input_file file name of input data containing rsid, chr, bp, a1, a2, and z 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
////' @param interval number of non-overlapping SNP sets used in calculating population weights 
//' @return R data frame containing population IDs and weights 
// [[Rcpp::export]]
NumericMatrix zmix4(std::string input_file,
                    std::string reference_index_file,
                    std::string reference_data_file,
                    std::string reference_pop_desc_file,
                    Rcpp::Nullable<int> interval = R_NilValue,
                    Rcpp::Nullable<int> offset = R_NilValue){
  
  Arguments args;
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  //args.interval = 1000;
  
  if(interval.isNotNull()){
    args.interval = Rcpp::as<int>(interval);
  } else {
    args.interval = 1000;
  }
  
  int offset_value = 0;
  if(offset.isNotNull()){
    offset_value = Rcpp::as<int>(offset);
  } else {
    offset_value = 3;
  }
  
  read_ref_desc(args);
  
#ifdef ZMIX_Debug  
  args.PrintArguments();
#endif
  
  std::map<MapKey, Snp*, LessThanMapKey> measured_snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_msm;
  std::vector<Snp*> snp_vec;
  
  read_input_zmix(measured_snp_map, args);
  read_ref_index_zmix(measured_snp_map, args);
  
#ifdef ZMIX_Debug    
  Rcpp::Rcout<<"Measured snp map size: "<< measured_snp_map.size() <<std::endl;
#endif
  
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end(); ++it_msm){
    int type = (it_msm->second)->GetType();
    if(type == 1)
      snp_vec.push_back(it_msm->second);
    //(it_msm->second)->PrintSnpInfo();
  }
  
#ifdef ZMIX_Debug  
  Rcpp::Rcout<<"Num of SNPs used for calculations: "<< snp_vec.size() <<std::endl;
#endif
  
  
  //////////////////////////
  Rcpp::Rcout<<"Num of SNPs: "<<snp_vec.size()<<std::endl;
  Rcpp::Rcout<<"Interval length: "<<args.interval<<std::endl;
  Rcpp::Rcout<<"Num of populations: "<<args.num_pops<<std::endl;
  
  Rcpp::Rcout<<"Calculating population weights..."<<std::endl;
  
  // std::vector<Snp*> snp_subvec;
  // for(int i=0;;i++){
  //   int index = i*args.interval;
  //   if(index < measured_snp_vec.size()){
  //     snp_subvec.push_back(measured_snp_vec[index]);
  //   }else {
  //     break;
  //   }
  // }
  
  int snp_vec_size = snp_vec.size();
  //Rcpp::Rcout<<"snp_subvec_size: "<<snp_subvec_size<<std::endl;
  
  for(int i=0; i<args.num_pops; i++){
    args.pop_flag_vec.push_back(1);
  }
  ReadGenotype(snp_vec, args);
  
  // Determine the total number of rows that will be needed in the matrix.
  int total_rows = 0;
  for (int h = 0; h < args.interval; h++) {  // Loop over 
    for (int i = h; i < snp_vec_size; i += args.interval) {
      if (i + offset_value < snp_vec_size) {
        total_rows++;
      }
    }
  }
  // Create the output matrix with the appropriate number of rows and columns
  NumericMatrix data_mat(total_rows, 2 + args.num_pops);
  int row_index = 0;
  for (int h = 0; h < args.interval; h++) {  // Loop over
    for (int i = h; i < snp_vec_size; i += args.interval) {
      int index_i = i;
      int index_j = i + offset_value;
      if (index_j < snp_vec_size) {
        std::vector<std::string>& snpi_geno_vec = snp_vec[index_i]->GetGenotypeVec();
        double snpi_z = snp_vec[index_i]->GetZ();
        std::vector<std::string>& snpj_geno_vec = snp_vec[index_j]->GetGenotypeVec();
        double snpj_z = snp_vec[index_j]->GetZ();
        // The first column contains h index
        data_mat(row_index, 0) = h;
        // The second column contains the product of snpi_z and snpj_z
        data_mat(row_index, 1) = snpi_z * snpj_z;
        // The rest of the columns contain the correlations
        for (int k = 0; k < args.num_pops; k++) {
          double cor = CalCor(snpi_geno_vec[k], snpj_geno_vec[k]);
          data_mat(row_index, k + 2) = cor;
        }
        row_index++;
      } else {
        break;
      }
    }
  }
  
  // release memory allocated for genotype
  FreeGenotype(snp_vec);
  
  //deletes measured_snp_map.
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end();){
    (it_msm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_msm->second;        // delete snp object
    measured_snp_map.erase(it_msm++);      // delete map element
  }
  
  return data_mat;
}

//' Calculate population weights using association Z-scores
//' 
//' @param input_file file name of input data containing rsid, chr, bp, a1, a2, and z 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
////' @param interval number of non-overlapping SNP sets used in calculating population weights 
//' @return R data frame containing population IDs and weights 
// [[Rcpp::export]]
NumericMatrix zmix3(std::string input_file,
                   std::string reference_index_file,
                   std::string reference_data_file,
                   std::string reference_pop_desc_file,
                   Rcpp::Nullable<int> interval = R_NilValue,
                   Rcpp::Nullable<int> steps = R_NilValue){
  
  Arguments args;
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  //args.interval = 1000;
  
  if(interval.isNotNull()){
    args.interval = Rcpp::as<int>(interval);
  } else {
    args.interval = 1000;
  }
  
  int steps_value = 0;
  if(steps.isNotNull()){
    steps_value = Rcpp::as<int>(steps);
  } else {
    steps_value = 5;
  }
  
  
  read_ref_desc(args);
  
#ifdef ZMIX_Debug  
  args.PrintArguments();
#endif
  
  std::map<MapKey, Snp*, LessThanMapKey> measured_snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_msm;
  std::vector<Snp*> measured_snp_vec;
  
  read_input_zmix(measured_snp_map, args);
  read_ref_index_zmix(measured_snp_map, args);
  
#ifdef ZMIX_Debug    
  Rcpp::Rcout<<"Measured snp map size: "<< measured_snp_map.size() <<std::endl;
#endif
  
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end(); ++it_msm){
    int type = (it_msm->second)->GetType();
    if(type == 1)
      measured_snp_vec.push_back(it_msm->second);
    //(it_msm->second)->PrintSnpInfo();
  }
  
#ifdef ZMIX_Debug  
  Rcpp::Rcout<<"Num of measured SNPs used for calculations: "<< measured_snp_vec.size() <<std::endl;
#endif
  
  
  //////////////////////////
  Rcpp::Rcout<<"Num of SNPs: "<<measured_snp_vec.size()<<std::endl;
  Rcpp::Rcout<<"Interval length: "<<args.interval<<std::endl;
  Rcpp::Rcout<<"Num of populations: "<<args.num_pops<<std::endl;
  
  
  Rcpp::Rcout<<"Calculating population weights..."<<std::endl;
  std::vector<Snp*> snp_subvec;
  for(int i=0;;i++){
    int index = i*args.interval;
    if(index < measured_snp_vec.size()){
      snp_subvec.push_back(measured_snp_vec[index]);
    }else {
      break;
    }
  }
  int snp_subvec_size = snp_subvec.size();
  Rcpp::Rcout<<"snp_subvec_size: "<<snp_subvec_size<<std::endl;
  
  for(int i=0; i<args.num_pops; i++){
    args.pop_flag_vec.push_back(1);
  }
  ReadGenotype(snp_subvec, args);
  
  // Determine the total number of rows that will be needed in the matrix.
  int total_rows = 0;
  for (int i = 0; i < snp_subvec_size; i++) {
    total_rows += std::min(steps_value, snp_subvec_size - (i + 1));
  }
  
  // Create the output matrix with the appropriate number of rows and columns
  NumericMatrix data_mat(total_rows, 1 + args.num_pops); 
  int row_index = 0;
  
  for(int i = 0; i < snp_subvec_size; i++) {
    std::vector<std::string>& snpi_geno_vec = snp_subvec[i]->GetGenotypeVec();
    double snpi_z = snp_subvec[i]->GetZ();
    for(int j = i + 1; j < std::min(i + 1 + steps_value, snp_subvec_size); j++) {
      std::vector<std::string>& snpj_geno_vec = snp_subvec[j]->GetGenotypeVec();
      double snpj_z = snp_subvec[j]->GetZ();
      
      // The first column contains the product of snpi_z and snpj_z
      data_mat(row_index, 0) = snpi_z * snpj_z;
      
      // The rest of the columns contain the correlations
      for(int k = 0; k < args.num_pops; k++) {
        double cor = CalCor(snpi_geno_vec[k], snpj_geno_vec[k]);
        data_mat(row_index, k + 1) = cor;
      }
      row_index++;
    }
  }

  // release memory allocated for genotype
  FreeGenotype(snp_subvec);
  
  //deletes measured_snp_map.
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end();){
    (it_msm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_msm->second;        // delete snp object
    measured_snp_map.erase(it_msm++);      // delete map element
  }
  
  return data_mat;
}

//' Calculate population weights using association Z-scores
//' 
//' @param input_file file name of input data containing rsid, chr, bp, a1, a2, and z 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
////' @param interval number of non-overlapping SNP sets used in calculating population weights 
//' @return R data frame containing population IDs and weights 
// [[Rcpp::export]]
NumericMatrix zmix2(std::string input_file,
                   std::string reference_index_file,
                   std::string reference_data_file,
                   std::string reference_pop_desc_file,
                   Rcpp::Nullable<int> interval = R_NilValue,
                   Rcpp::Nullable<int> offset = R_NilValue){
  
  Arguments args;
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  //args.interval = 1000;
  
  if(interval.isNotNull()){
    args.interval = Rcpp::as<int>(interval);
  } else {
    args.interval = 1000;
  }
  
  int offset_value = 0;
  if(offset.isNotNull()){
    offset_value = Rcpp::as<int>(offset);
  } else {
    offset_value = 3;
  }
  
  read_ref_desc(args);
  
#ifdef ZMIX_Debug  
  args.PrintArguments();
#endif
  
  std::map<MapKey, Snp*, LessThanMapKey> measured_snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_msm;
  std::vector<Snp*> snp_vec;
  
  read_input_zmix(measured_snp_map, args);
  read_ref_index_zmix(measured_snp_map, args);
  
#ifdef ZMIX_Debug    
  Rcpp::Rcout<<"Measured snp map size: "<< measured_snp_map.size() <<std::endl;
#endif
  
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end(); ++it_msm){
    int type = (it_msm->second)->GetType();
    if(type == 1)
      snp_vec.push_back(it_msm->second);
    //(it_msm->second)->PrintSnpInfo();
  }
  
#ifdef ZMIX_Debug  
  Rcpp::Rcout<<"Num of SNPs used for calculations: "<< snp_vec.size() <<std::endl;
#endif
  
  
  //////////////////////////
  Rcpp::Rcout<<"Num of SNPs: "<<snp_vec.size()<<std::endl;
  Rcpp::Rcout<<"Interval length: "<<args.interval<<std::endl;
  Rcpp::Rcout<<"Num of populations: "<<args.num_pops<<std::endl;
  
  Rcpp::Rcout<<"Calculating population weights..."<<std::endl;

  // std::vector<Snp*> snp_subvec;
  // for(int i=0;;i++){
  //   int index = i*args.interval;
  //   if(index < measured_snp_vec.size()){
  //     snp_subvec.push_back(measured_snp_vec[index]);
  //   }else {
  //     break;
  //   }
  // }

  int snp_vec_size = snp_vec.size();
  //Rcpp::Rcout<<"snp_subvec_size: "<<snp_subvec_size<<std::endl;
  
  for(int i=0; i<args.num_pops; i++){
    args.pop_flag_vec.push_back(1);
  }
  ReadGenotype(snp_vec, args);
  
  // Determine the total number of rows that will be needed in the matrix.
  int total_rows = 0;
  // int offset = 3; // to get moderate LD  
  for(int i = 0; i < snp_vec_size; i += args.interval) {
    if(i + offset_value < snp_vec_size){total_rows++;}
  }

  // Create the output matrix with the appropriate number of rows and columns
  NumericMatrix data_mat(total_rows, 1 + args.num_pops); 
  int row_index = 0;
  for(int i = 0; i < snp_vec_size; i += args.interval) {
    int index_i = i;
    int index_j = i + offset_value;
    if(index_j < snp_vec_size){
      std::vector<std::string>& snpi_geno_vec = snp_vec[index_i]->GetGenotypeVec();
      double snpi_z = snp_vec[index_i]->GetZ();
      std::vector<std::string>& snpj_geno_vec = snp_vec[index_j]->GetGenotypeVec();
      double snpj_z = snp_vec[index_j]->GetZ();    
      // The first column contains the product of snpi_z and snpj_z
      data_mat(row_index, 0) = snpi_z * snpj_z;
      // The rest of the columns contain the correlations
      for(int k = 0; k < args.num_pops; k++) {
        double cor = CalCor(snpi_geno_vec[k], snpj_geno_vec[k]);
        data_mat(row_index, k + 1) = cor;
      }
      row_index++;
    } else{
      break;
    }
  }
 
  // release memory allocated for genotype
  FreeGenotype(snp_vec);
  
  //deletes measured_snp_map.
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end();){
    (it_msm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_msm->second;        // delete snp object
    measured_snp_map.erase(it_msm++);      // delete map element
  }
  
  return data_mat;
}

//' Calculate population weights using association Z-scores
//' 
//' @param input_file file name of input data containing rsid, chr, bp, a1, a2, and z 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
////' @param interval number of non-overlapping SNP sets used in calculating population weights 
//' @return R data frame containing population IDs and weights 
// [[Rcpp::export]]
NumericMatrix zmix(std::string input_file,
                   std::string reference_index_file,
                   std::string reference_data_file,
                   std::string reference_pop_desc_file,
                   Rcpp::Nullable<int> interval = R_NilValue){
  
  Arguments args;
  //args.chr = chr;
  //args.start_bp = start_bp;
  //args.end_bp = end_bp;
  args.input_file = input_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  //args.interval = 1000;
  
  if(interval.isNotNull()){
    args.interval = Rcpp::as<int>(interval);
  } else {
    args.interval = 1;
  }
  
  read_ref_desc(args);
  
#ifdef ZMIX_Debug  
  args.PrintArguments();
#endif
  
  std::map<MapKey, Snp*, LessThanMapKey> measured_snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_msm;
  std::vector<Snp*> measured_snp_vec;

  read_input_zmix(measured_snp_map, args);
  read_ref_index_zmix(measured_snp_map, args);
  
#ifdef ZMIX_Debug    
  Rcpp::Rcout<<"Measured snp map size: "<< measured_snp_map.size() <<std::endl;
#endif
  
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end(); ++it_msm){
    int type = (it_msm->second)->GetType();
    if(type == 1)
      measured_snp_vec.push_back(it_msm->second);
      //(it_msm->second)->PrintSnpInfo();
  }
  
#ifdef ZMIX_Debug  
  Rcpp::Rcout<<"Num of measured SNPs used for calculations: "<< measured_snp_vec.size() <<std::endl;
#endif
  

  //////////////////////////
  Rcpp::Rcout<<"Num of SNPs: "<<measured_snp_vec.size()<<std::endl;
  Rcpp::Rcout<<"Interval length: "<<args.interval<<std::endl;
  Rcpp::Rcout<<"Num of populations: "<<args.num_pops<<std::endl;

    
  Rcpp::Rcout<<"Calculating population weights..."<<std::endl;
  std::vector<Snp*> snp_subvec;
  for(int i=0;;i++){
    int index = i*args.interval;
    if(index < measured_snp_vec.size()){
      snp_subvec.push_back(measured_snp_vec[index]);
    }else {
      break;
    }
  }
  int snp_subvec_size = snp_subvec.size();
  Rcpp::Rcout<<"snp_subvec_size: "<<snp_subvec_size<<std::endl;

  for(int i=0; i<args.num_pops; i++){
    args.pop_flag_vec.push_back(1);
  }
  ReadGenotype(snp_subvec, args);

  // Determine the total number of rows that will be needed in the matrix.
  int total_rows = (snp_subvec_size * (snp_subvec_size - 1)) / 2;
  
  // Create the output matrix with the appropriate number of rows and columns
  NumericMatrix data_mat(total_rows, 1 + args.num_pops); 
  int row_index = 0;
  
  for(int i = 0; i < snp_subvec_size; i++) {
    std::vector<std::string>& snpi_geno_vec = snp_subvec[i]->GetGenotypeVec();
    double snpi_z = snp_subvec[i]->GetZ();
    for(int j = i + 1; j < snp_subvec_size; j++) {
      std::vector<std::string>& snpj_geno_vec = snp_subvec[j]->GetGenotypeVec();
      double snpj_z = snp_subvec[j]->GetZ();
      
      // The first column contains the product of snpi_z and snpj_z
      data_mat(row_index, 0) = snpi_z * snpj_z;
      
      // The rest of the columns contain the correlations
      for(int k = 0; k < args.num_pops; k++) {
        double cor = CalCor(snpi_geno_vec[k], snpj_geno_vec[k]);
        data_mat(row_index, k + 1) = cor;
      }
      row_index++;
    }
  }
  
  
/*    
  std::ofstream data_out;
  data_out.open(output_file.c_str());

  for(int i=0; i<snp_subvec_size; i++){
    std::vector<std::string>& snpi_geno_vec = (*snp_subvec[i]).GetGenotypeVec();
    double snpi_z = (*snp_subvec[i]).GetZ(); 
    for(int j=i+1; j<snp_subvec_size; j++){
      std::vector<std::string>& snpj_geno_vec = (*snp_subvec[j]).GetGenotypeVec();
      double snpj_z = (*snp_subvec[j]).GetZ();
      data_out << snpi_z*snpj_z <<" ";  
      for(int k=0; k<args.num_pops; k++){
        double cor = CalCor(snpi_geno_vec[k],snpj_geno_vec[k]);
        data_out << cor <<" ";
      }
      data_out<<std::endl;
    }
  }
  data_out.close(); //close output filestream
*/ 
 
  // release memory allocated for genotype
  FreeGenotype(snp_subvec);

  //deletes measured_snp_map.
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end();){
    (it_msm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_msm->second;        // delete snp object
    measured_snp_map.erase(it_msm++);      // delete map element
  }
  
  return data_mat;
}


void read_input_zmix(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  Rcpp::Rcout<<"Reading input...";
  Rcpp::Rcout.flush();
  
  std::string input_file = args.input_file;
  std::ifstream in_input(input_file.c_str());
  if(!in_input){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open input file '"+input_file+"'");
  }
  
  std::string line;
  std::string rsid, a1, a2;
  int chr;
  long long int bp;
  double zscore;
  Snp* snp;
  
  std::getline(in_input, line); //read header of input file.  
  while(std::getline(in_input, line)){
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> zscore;
    
    snp = new Snp();
    snp->SetRsid(rsid);
    snp->SetChr(chr);
    snp->SetBp(bp);
    snp->SetA1(a1);
    snp->SetA2(a2);
    snp->SetZ(zscore); 
    MapKey mkey(chr, bp, a1, a2);
    snp_map[mkey]=snp;
  }//while
  in_input.close();
  Rcpp::Rcout<<std::endl;
}

void read_ref_index_zmix(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  Rcpp::Rcout<<"Reading reference index...";
  Rcpp::Rcout.flush();
  
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;
  std::string reference_index_file = args.reference_index_file;
  BGZF* fp = bgzf_open(reference_index_file.c_str(), "r");
  if(!fp){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference index file '"+reference_index_file+"'");
  }
  
  int last_char;
  std::string line;
  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  
  while(true){
    
    last_char = BgzfGetLine(fp, line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    MapKey mkey1(chr, bp, a1, a2);
    MapKey mkey2(chr, bp, a2, a1);
    
    it1 = snp_map.find(mkey1);
    it2 = snp_map.find(mkey2);
    
    if((it1 != snp_map.end()) && (it2 == snp_map.end())){ // snp exists in input and a1=a1 & a2=a2.
      
      (it1->second)->SetRsid(rsid);
      (it1->second)->SetType(1); // change to "measured and exists in reference panel"
      //(it1->second)->SetAf1ref(af1ref);
      (it1->second)->SetFpos(fpos); // index of snp genotypes
      
    } else if((it1 == snp_map.end()) && (it2 != snp_map.end())){ // snp exists in input but a1=a2 & a2=a1.
      
      (it2->second)->SetRsid(rsid);
      (it2->second)->SetA1(a1);
      (it2->second)->SetA2(a2);
      (it2->second)->SetZ( (it2->second)->GetZ()*(-1) ); //change the sign of z-score
      (it2->second)->SetType(1); // change to "measured and exists in reference panel"
      (it2->second)->SetAf1Study( 1 - (it2->second)->GetAf1Study() );
      //(it2->second)->SetAf1ref(af1ref);
      (it2->second)->SetFpos(fpos); // index of snp genotypes
      
      MapKey new_key(chr, bp, a1, a2); //makes a new key for the modified snp object.
      snp_map[new_key] = it2->second;  //makes a map element with the new key. 
      snp_map.erase(it2);              //delete snp with the old key. 
      
    } else if(it1 != snp_map.end() && it2 != snp_map.end()){ //Throw error message! This should not happen.
      Rcpp::Rcout<<std::endl;
      Rcpp::stop("ERROR: input file contains duplicates");
    }//if
    
  }//while
  bgzf_close(fp);
  Rcpp::Rcout<<std::endl;
}
