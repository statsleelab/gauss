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


//' JEPEGMIX: gene-level joint analysis of functional SNPs in cosmopolitan cohorts
//' 
//' @param pop_wgt_df R data frame containing population IDs and weights
//' @param input_file file name of GWAS summary statistics data containing rsid, chr, bp, a1, a2 and z
//' @param annotation file name of the SNP annotation data set 
//' @param reference_index_file file name of reference panel index data
//' @param reference_data_file  file name of reference panel data
//' @param reference_pop_desc_file file name of reference panel population description data
//' @param af1_cutoff cutoff of reference allele, a1, frequency
//' @return R dataframe containing geneid, chisq, df, jepeg_pval, num_snp, top_categ, top_categ_pval, top_snp, and top_snp_pval   
// [[Rcpp::export]]
DataFrame jepegmix(DataFrame pop_wgt_df,
                   std::string input_file,
                   std::string annotation_file,
                   std::string reference_index_file,
                   std::string reference_data_file,
                   std::string reference_pop_desc_file,
                   Rcpp::Nullable<double> af1_cutoff = R_NilValue){
  
  Arguments args;

  // add pop_wgt_df info in args.pop_wgt_map
  std::vector<std::string> pop_vec_in = as<std::vector<std::string>>(pop_wgt_df[0]);
  std::vector<double> pop_wgt_vec_in = as<std::vector<double>>(pop_wgt_df[1]);
  for(int i=0; i<pop_vec_in.size(); i++){
    std::string pop = pop_vec_in[i];
    std::transform(pop.begin(), pop.end(), pop.begin(), ::toupper); //make capital
    args.pop_wgt_map[pop]=pop_wgt_vec_in[i];
    //Rcpp::Rcout<<pop<<" "<<pop_wgt_vec_in[i]<<std::endl;
  }  
  
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
  init_pop_flag_wgt_vec(args);
  //args.PrintArguments();
  
  std::map<MapKey, Snp*, LessThanMapKey> snp_map;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  
  ReadInputZ(snp_map, args, true);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  ReadReferenceIndexAll(snp_map, args);
  //Rcpp::Rcout<<"size: "<< snp_map.size() <<std::endl;
  
  ReadAnnotation(snp_map, args);
  
  // make a snp vector containing SNPs with af1ref > af1_cutoff & af1ref < 1-af1_cutoff
  std::vector<Snp*> snp_vec;
  std::vector<Snp*> gene_snp_vec;
  std::vector<Snp*>::iterator it_sv;
  MakeSnpVecMix(snp_vec, snp_map, args);
  
  //Rcpp::Rcout<<"snp_map size : "<<snp_map.size()<<std::endl;
  //Rcpp::Rcout<<"snp_vec size : "<<snp_vec.size()<<std::endl;
  // use only measured SNPs with geneid
  for(it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    if( (*it_sv)->GetGeneid() != "." && 
        (*it_sv)->GetType() == 1){
      gene_snp_vec.push_back(*it_sv);
      //(it_sm->second)->PrintSnpInfo();
    }
  }
  //Rcpp::Rcout<<"gene_snp_vec size : "<<gene_snp_vec.size()<<std::endl;
  
  // read genotypes of SNPs in gene_snp_vec
  ReadGenotype(gene_snp_vec, args);
  
  //Sort gene_snp_vec by geneid
  //TODO: make LessThanGeneid() to sort snps by (geneid, chr, bp)
  std::sort(gene_snp_vec.begin(), gene_snp_vec.end(), LessThanGeneid());
  
  /*
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############### Print SNPs in gene_snp_vec " <<std::endl;
  for(std::vector<Snp*>::iterator it_gsv=gene_snp_vec.begin(); it_gsv != gene_snp_vec.begin()+40; ++it_gsv){
    (*it_gsv)->PrintSnpInfo();
  } 
  */
  
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
  //Rcpp::Rcout<<"gene_start_end_vec size  (before adding): "<<gene_start_end_vec.size()<<std::endl;
  MakeGeneStartEndVec(gene_start_end_vec, gene_snp_vec);
  //Rcpp::Rcout<<"gene_start_end_vec size  (after adding): "<<gene_start_end_vec.size()<<std::endl;
  //int counter=0;
  for(it_gsev = gene_start_end_vec.begin(); it_gsev != gene_start_end_vec.end(); ++it_gsev){
    std::vector<Snp*> gvec(it_gsev->start_it, it_gsev->end_it);
    //Rcpp::Rcout<<gvec.size()<<std::endl;
    //counter = counter + gvec.size();
    
    Gene gene(args);
    gene.RunJepegmix(gvec);
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
  //Rcpp::Rcout<<"Total: "<< counter <<std::endl;
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
                                   Named("top_snp")=top_snp_vec,
                                   Named("top_snp_pval")=top_snp_pval_vec);
  return df;  
}
