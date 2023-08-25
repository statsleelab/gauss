//GAUSS: Genome Analysis Using Summary Statistics
//gene.h

#ifndef GENE_H
#define GENE_H

struct Arguments;
class Snp;

#include <fstream>
#include <string>
#include <vector>
#include <gsl/gsl_matrix.h>


struct Categ{
public:
  Categ(int, double, bool);
  int GetNum() {return num_;}
  double GetPval() {return pval_;}
  bool GetRmv() {return rmv_;}
  std::string GetName();
  
  void SetPval(double pval) { pval_ = pval; }
  void SetRmv(bool rmv) { rmv_ = rmv; }
  
  void PrintInfo();
  
private:
  int num_;
  double pval_;    //p-value of the category
  bool rmv_;    //true means a categ not being used in JEPEG pval calculation
};


class Gene{
  
public:
  Gene(Arguments& args);
  void RunJepeg(std::vector<Snp*>& gene_snp_vec);
  void RunJepegmix(std::vector<Snp*>& gene_snp_vec);
  
  //get
  std::string GetGeneid() {return geneid_;}
  double GetChisq() {return chisq_;}
  int GetDf() {return df_;}
  double GetJepegPval() {return jepeg_pval_;}
  int GetNumSnp() {return num_snp_;}
  std::string GetTopCategName() {return top_categ_name_;}
  double GetTopCategPval() {return top_categ_pval_;}
  std::string GetTopSnpId() {return top_snp_id_;}
  double GetTopSnpPval() {return top_snp_pval_;}
  
  //set
  void SetGeneid(std::string geneid) {geneid_=geneid;}
  void SetChisq(double chisq) {chisq_=chisq;}
  void SetDf(int df) {df_=df;}
  void SetJepegPval(double jepeg_pval) {jepeg_pval_=jepeg_pval;}
  void SetNumSnp(int num_snp) {num_snp_=num_snp;}
  void SetTopCategName(std::string top_categ_name) {top_categ_name_=top_categ_name;}
  void SetTopCategPval(double top_categ_pval) {top_categ_pval_=top_categ_pval;}
  void SetTopSnpId(std::string top_snp_id) {top_snp_id_=top_snp_id;}
  void SetTopSnpPval(double top_snp_pval) {top_snp_pval_=top_snp_pval;}
  
private:  
  
  //JEPEG/MIX output
  std::string geneid_;
  double chisq_;
  int df_;
  double jepeg_pval_;
  int num_snp_;
  std::string top_categ_name_;
  double top_categ_pval_;
  std::string top_snp_id_;
  double top_snp_pval_;
  
  int num_avail_categ_;
  int num_rmv_categ_;
  
  int num_protein_; //Protein Effect stopcodon
  int num_tfbs_;    //TFBS
  int num_wthhair_; //within hairpin(micro RNA)
  int num_wthtar_;  //within target (micro RNA)
  int num_cis_;     //cis eqtls
  int num_trans_;   //trans eqtls
  
  std::vector<int> categ_count_vec_;
  std::vector<Categ> categ_vec_;  
  
  double bonfe_pval_bf_;
  double bonfe_pval_af_;
  double sumU_pval_;
  
  double lambda_;
  double min_abs_eig_;
  int total_num_categ_;
  double categ_cor_cutoff_;
  int denorm_norm_w_;
  
  std::vector<bool> rmv_categ_list_debug_; //for Debugging
  std::vector<double> U_list_debug_;  //for Debugging
  
  std::vector<int> pop_flag_vec_;
  std::vector<double> pop_wgt_vec_;
  
  void CalJepegPval(std::vector<Snp*>& gene_snp_vec);
  void CalJepegmixPval(std::vector<Snp*>& gene_snp_vec);
  
  void GetW(gsl_matrix*, std::vector<Snp*>& gene_snp_vec);
  void GetZ(gsl_matrix*, std::vector<Snp*>& gene_snp_vec);
  double CalBonfePval(gsl_matrix*);
  double CalSumUPval(gsl_matrix*, gsl_matrix*);
  
  Categ GetTopCateg(std::vector<Categ>& categ_vec);
  Snp* GetTopSNP(std::vector<Snp*>& gene_snp_vec);
  
};

#endif
