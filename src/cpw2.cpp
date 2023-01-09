#include <Rcpp.h>

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "util.h"
#include "bgzf.h"
#include "snp.h"
#include "gauss.h"

using namespace Rcpp;

//#define CPW2_Debug  

//Forward declaration
void cpw2_vec(std::vector<Snp*>& snp_vec, Arguments& args);
void read_input_cpw2(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
void read_ref_index_cpw2(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);


//' Calculate population weights using RAF
//' 
//' @param input_file file name of input data containing rsid, chr, bp, a1, a2, and af1 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param interval number of non-overlapping SNP sets used in calculating population weights 
//' @return R data frame containing population IDs and weights 
// [[Rcpp::export]]
DataFrame cpw2(std::string input_file,
              std::string reference_index_file,
              std::string reference_data_file,
              std::string reference_pop_desc_file,
              Rcpp::Nullable<int> interval = R_NilValue){
  
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
  
  
  read_ref_desc(args);
  
#ifdef CPW2_Debug  
  args.PrintArguments();
#endif
  
  std::map<MapKey, Snp*, LessThanMapKey> measured_snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_msm;
  std::vector<Snp*> measured_snp_vec;

  read_input_cpw2(measured_snp_map, args);
  read_ref_index_cpw2(measured_snp_map, args);
  
#ifdef CPW2_Debug    
  Rcpp::Rcout<<"Measured snp map size: "<< measured_snp_map.size() <<std::endl;
#endif
  
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end(); ++it_msm){
    int type = (it_msm->second)->GetType();
    if(type == 1)
      measured_snp_vec.push_back(it_msm->second);
      //(it_msm->second)->PrintSnpInfo();
  }
  
#ifdef CPW2_Debug  
  Rcpp::Rcout<<"Num of measured SNPs used for calculations: "<< measured_snp_vec.size() <<std::endl;
#endif
  
  cpw2_vec(measured_snp_vec, args);
  
  //deletes measured_snp_map.
  for(it_msm = measured_snp_map.begin(); it_msm != measured_snp_map.end();){
    (it_msm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_msm->second;        // delete snp object
    measured_snp_map.erase(it_msm++);      // delete map element
  }
  
  std::vector<std::string> pop_vec;
  std::vector<double> wgt_vec;
  
  for(int i=0; i<args.ref_pop_vec.size(); i++){
    if(args.pop_wgt_vec[i]>0){
      pop_vec.push_back(args.ref_pop_vec[i]);
      wgt_vec.push_back(args.pop_wgt_vec[i]);
    }
  }
  
  StringVector pop(pop_vec.size());
  NumericVector wgt(pop_vec.size());
  
  pop = pop_vec;
  wgt = wgt_vec;
  
  DataFrame df = DataFrame::create(Named("pop")=pop,
                                   Named("wgt")=wgt);
  return df;  
}


void cpw2_vec(std::vector<Snp*>& snp_vec, Arguments& args){
  
  int size = snp_vec.size(); 
  int interval = args.interval;
  
  Rcpp::Rcout<<"Num of SNPs: "<<size<<std::endl;
  Rcpp::Rcout<<"Interval length: "<<interval<<std::endl;
  Rcpp::Rcout<<"Num of populations: "<<args.num_pops<<std::endl;

  gsl_matrix* af1_cor_mat = gsl_matrix_calloc(args.num_pops+1, args.num_pops+1);
  gsl_matrix* W_mat_i = gsl_matrix_calloc(args.num_pops, 1);
  gsl_matrix* W_mat = gsl_matrix_calloc(args.num_pops, 1);
  gsl_matrix_view Cxy;
  gsl_matrix_view Cxx;
  gsl_matrix* CxxInv = gsl_matrix_calloc(args.num_pops, args.num_pops);
  
  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+args.reference_data_file+"'");
  }
  
  Rcpp::Rcout<<"Calculating population weights..."<<std::endl;
  for(int i=0; i<interval; i++){
    std::vector<Snp*> snp_subvec;
    for(int j=0;;j++){
      int index = i+j*interval;
      if(index < size){
        snp_subvec.push_back(snp_vec[index]);
      }else {
        break;
      }
    }
    int snp_subvec_size = snp_subvec.size();
    //Rcpp::Rcout<<"snp_subvec_size: "<<snp_subvec_size<<std::endl;
    
    gsl_matrix* af1_mat = gsl_matrix_calloc(snp_subvec_size, args.num_pops+1);
    for(int j=0; j<snp_subvec_size; j++){
      gsl_matrix_set(af1_mat, j, 0, std::asin(std::sqrt(snp_subvec[j]->GetAf1Study()))); //use arcsine square root transformation to stabilize the variance    
      //Rcpp::Rcout<<"af1study :"<<snp_subvec[j]->GetAf1Study();
      std::string line;
      bgzf_seek(fp, (*snp_subvec[j]).GetFpos(), SEEK_SET);
      BgzfGetLine(fp, line);
      //Rcpp::Rcout<<line<<std::endl;
      std::istringstream buffer(line);
      //skipping genotypes
      for(int k=0; k<args.num_pops; k++){
        std::string tmp;
        buffer >> tmp;
        //Rcpp::Rcout<<tmp<<std::endl;
      }
      //reading allele frequencies
      int kk=1;
      for(int k=0; k<args.num_pops; k++){
        double af1;
        buffer >> af1;
        //Rcpp::Rcout<<af1<<std::endl;
        gsl_matrix_set(af1_mat, j, kk, std::asin(std::sqrt(af1))); // use arcsine square root transformation to stabilize the variance
        kk++;
      }
    }
    
    CalCovMat(af1_cor_mat, af1_mat);
    
    //Calculate W
    Cxy = gsl_matrix_submatrix(af1_cor_mat, 1, 0, args.num_pops, 1);
    Cxx = gsl_matrix_submatrix(af1_cor_mat, 1, 1, args.num_pops, args.num_pops);
    MakePosDef(&Cxx.matrix, args.min_abs_eig);
    InvMat(CxxInv, &Cxx.matrix);
    MpMatMat(W_mat_i, CxxInv, &Cxy.matrix);
    
    for(int j=0; j<args.num_pops; j++){
      double wval = gsl_matrix_get(W_mat, j, 0) + gsl_matrix_get(W_mat_i, j, 0)/interval;
      gsl_matrix_set(W_mat, j, 0, wval);
    }
    
    gsl_matrix_free(af1_mat);
    
    if((i%10)==9){
      int percent = static_cast<int>(0.1+100*(i/(double)(interval)));
      LoadProgressBar(percent);
    }
    
  }
  // if w is less than zero then make it zero.
  for(int i=0; i<args.num_pops; i++){
    double wval = gsl_matrix_get(W_mat, i, 0);
    if(wval < 0)
      gsl_matrix_set(W_mat, i, 0, 0);
    else
      gsl_matrix_set(W_mat, i, 0, floor(wval*1000+0.5)/1000);
  }

  double sum_pop_wgt = 0.0;
  for(int i=0; i<args.num_pops; i++){
    double wval = gsl_matrix_get(W_mat, i, 0);
    args.pop_wgt_vec.push_back(wval);
    sum_pop_wgt += wval;
  }
  args.sum_pop_wgt = sum_pop_wgt;

  gsl_matrix_free(af1_cor_mat);
  gsl_matrix_free(W_mat_i);
  gsl_matrix_free(W_mat);
  gsl_matrix_free(CxxInv);
  
  //closes file connections
  bgzf_close(fp); //closes BGZF file connnection.
  Rcpp::Rcout<<std::endl;
}


void read_input_cpw2(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
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
  double af1study;
  Snp* snp;
  
  std::getline(in_input, line); //read header of input file.  
  while(std::getline(in_input, line)){
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1study;
    
    snp = new Snp();
    snp->SetRsid(rsid);
    snp->SetChr(chr);
    snp->SetBp(bp);
    snp->SetA1(a1);
    snp->SetA2(a2);
    snp->SetAf1Study(af1study); 
    MapKey mkey(chr, bp, a1, a2);
    snp_map[mkey]=snp;
  }//while
  in_input.close();
  Rcpp::Rcout<<std::endl;
}

void read_ref_index_cpw2(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
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
