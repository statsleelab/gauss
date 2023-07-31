#include <Rcpp.h>

using namespace Rcpp;

#include <string>
#include <vector>
#include <map>
#include <deque>
#include "snp.h"
#include "gene.h"
#include "gauss.h"
#include "util.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cdf.h>

void read_input_jepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
void read_ref_index_jepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);

//' Joint effect on phenotype of eQTLs/functional SNPs associated with a gene
//' 
//' @param study_pop study population group
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2 and z
//' @param annotation file name of the SNP annotation data set 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R dataframe containing geneid, chisq, df, jepeg_pval, num_snp, top_categ, top_categ_pval, top_snp, and top_snp_pval   
// [[Rcpp::export]]
DataFrame jepeg(std::string study_pop,
                std::string input_file,
                std::string annotation_file,
                std::string reference_index_file,
                std::string reference_data_file,
                std::string reference_pop_desc_file,
                Rcpp::Nullable<double> af1_cutoff = R_NilValue){
  
  Arguments args;
  args.study_pop = study_pop;
  args.input_file = input_file;
  args.annotation_file = annotation_file;
  args.reference_index_file = reference_index_file;
  args.reference_data_file = reference_data_file;
  args.reference_pop_desc_file = reference_pop_desc_file;
  
  if(af1_cutoff.isNotNull()){
    args.af1_cutoff = Rcpp::as<double>(af1_cutoff);
  } else {
    args.af1_cutoff = 0.01;
  }

  read_ref_desc(args);
  init_pop_flag_vec(args);
  //args.PrintArguments();
  
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  
  read_input_jepeg(snp_map, args);
  Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  read_ref_index_jepeg(snp_map, args);
  Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  
  ReadAnnotation(snp_map, args);
  
  // make a snp vector containing SNPs with af1ref > af1_cutoff & af1ref < 1-af1_cutoff
  std::vector<Snp*> snp_vec;
  std::vector<Snp*> gene_snp_vec;
  std::vector<Snp*>::iterator it_sv;
  MakeSnpVec(snp_vec, snp_map, args);
  
  Rcpp::Rcout<<"snp_map size : "<<snp_map.size()<<std::endl;
  Rcpp::Rcout<<"snp_vec size : "<<snp_vec.size()<<std::endl;
  // use only measured SNPs with geneid
  for(it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    if( (*it_sv)->GetGeneid() != "." && 
        (*it_sv)->GetType() == 1){
      gene_snp_vec.push_back(*it_sv);
      //(it_sm->second)->PrintSnpInfo();
    }
  }
  Rcpp::Rcout<<"gene_snp_vec size : "<<gene_snp_vec.size()<<std::endl;
  
  // read genotypes of SNPs in gene_snp_vec
  ReadGenotype(gene_snp_vec, args);
  
  //Sort gene_snp_vec by geneid
  //TODO: make LessThanGeneid() to sort snps by (geneid, chr, bp)
  std::sort(gene_snp_vec.begin(), gene_snp_vec.end(), LessThanGeneid());
 
  
  StringVector geneid_vec;
  NumericVector chisq_vec;
  IntegerVector df_vec;
  NumericVector jepeg_pval_vec;
  IntegerVector num_snp_vec;
  StringVector top_categ_vec;
  NumericVector top_categ_pval_vec;
  StringVector top_snp_vec;
  NumericVector top_snp_pval_vec;
  
  //Make a vector containing start iterator and end iterator of each gene block.
  std::vector<StartEnd> gene_start_end_vec;
  std::vector<StartEnd>::iterator it_gsev;
  Rcpp::Rcout<<"gene_start_end_vec size  (before adding): "<<gene_start_end_vec.size()<<std::endl;
  MakeGeneStartEndVec(gene_start_end_vec, gene_snp_vec);
  Rcpp::Rcout<<"gene_start_end_vec size  (after adding): "<<gene_start_end_vec.size()<<std::endl;
  for(it_gsev = gene_start_end_vec.begin(); it_gsev != gene_start_end_vec.end(); ++it_gsev){
    std::vector<Snp*> gvec(it_gsev->start_it, it_gsev->end_it);
    Gene gene(args);
    gene.RunJepeg(gvec);
    geneid_vec.push_back(gene.GetGeneid());
    chisq_vec.push_back(gene.GetChisq());
    df_vec.push_back(gene.GetDf());
    jepeg_pval_vec.push_back(gene.GetJepegPval());
    num_snp_vec.push_back(gene.GetNumSnp());
    top_categ_vec.push_back(gene.GetTopCategName());
    top_categ_pval_vec.push_back(gene.GetTopCategPval());
    top_snp_vec.push_back(gene.GetTopSnpId());
    top_snp_pval_vec.push_back(gene.GetTopSnpPval());
  }
  
  // release memory allocated for genotype
  FreeGenotype(gene_snp_vec);
  
  //deletes snp_map.
  for(it_sm = snp_map.begin(); it_sm != snp_map.end();){
    (it_sm->second)->ClearSnp(); // clear categ map in each snp object
    delete it_sm->second;        // delete snp object
    snp_map.erase(it_sm++);      // delete map element
  }
  
  DataFrame df = DataFrame::create(Named("geneid")=geneid_vec,
                                   Named("chisq")=chisq_vec,
                                   Named("df")=df_vec,
                                   Named("jepeg_pval")=jepeg_pval_vec,
                                   Named("num_snp")=num_snp_vec,
                                   Named("top_categ")=top_categ_vec,
                                   Named("top_categ_pval")=top_categ_pval_vec,
                                   Named("top_snp_vec")=top_snp_vec,
                                   Named("top_snp_pval")=top_snp_pval_vec);
  return df;  
}


void read_input_jepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  
  Rcpp::Rcout<<"Reading input...";
  Rcpp::Rcout.flush();
  
  std::string input_file = args.input_file;
  std::ifstream in_input(input_file.c_str());
  
  if(!in_input){
    Rcpp::stop("ERROR: can't open input file '"+input_file+"'");
  }
  
  std::string line;
  std::string rsid, a1, a2;
  int chr;
  long long int bp;
  double z;
  double info = 1.0;
  Snp* snp;
  
  std::getline(in_input, line); //read header of input file.  
  while(std::getline(in_input, line)){
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> z;
    
    snp = new Snp();
    snp->SetRsid(rsid);
    snp->SetChr(chr);
    snp->SetBp(bp);
    snp->SetA1(a1);
    snp->SetA2(a2);
    snp->SetZ(z);
    snp->SetInfo(info);
    snp->SetType(2); // 2: measured SNP but does not exist in reference data.
    
    MapKey mkey(chr, bp, a1, a2);
    snp_map[mkey]=snp;
  }//while
  in_input.close();
  Rcpp::Rcout<<std::endl;
}

void read_ref_index_jepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
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
      //(it2->second)->SetAf1Study( 1 - (it2->second)->GetAf1Study() );
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
