//GAUSS: Genome Analysis Using Summary Statistics
//gene.cpp

#include "gene.h"

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <Rmath.h>


#include "snp.h"
#include "gauss.h"
#include "util.h"
//#include "bgzf.h"

//#define RunJepeg_Debug  //should move to gauss.h
//#define CalJepegPval_Debug
//#define GetW_Debug

Categ::Categ(int num, double pval, bool rmv){
  num_ = num;
  pval_ = pval;
  rmv_ = rmv;
}

std::string Categ::GetName(){
  std::string name;
  if(num_ == 0){
    name = "PFS";
  } else if(num_ == 1){
    name = "TFB";
  } else if(num_ == 2){
    name = "STR";
  } else if(num_ == 3){
    name = "TAR";
  } else if(num_ == 4){
    name = "CIS";
  } else if(num_ == 5){
    name = "TRN";
  }
  return name;
}

void Categ::PrintInfo(){
  Rcpp::Rcout<<num_<<" "<<pval_<<" "<<rmv_<<std::endl;
}

Gene::Gene(Arguments& args){
  
  //JEPEG/MIX output
  geneid_ = ".";
  chisq_ = -1.0;
  df_ = 0;
  jepeg_pval_ = -1.0;
  num_snp_ = 0;
  top_categ_name_ = ".";
  top_categ_pval_ = -1.0;
  top_snp_id_ = ".";
  top_snp_pval_ = -1.0;
  
  num_avail_categ_ = 0;
  num_rmv_categ_ = 0;
  
  num_protein_ = 0;
  num_tfbs_ = 0;
  num_wthhair_ = 0;
  num_wthtar_ = 0;
  num_cis_ = 0;
  num_trans_ = 0;  
  
  bonfe_pval_bf_ = -1.0;
  bonfe_pval_af_ = -1.0;
  sumU_pval_ = -1.0;
  
  lambda_ = args.lambda;
  min_abs_eig_ = args.min_abs_eig;
  total_num_categ_ = args.total_num_categ;
  categ_cor_cutoff_ = args.categ_cor_cutoff;
  denorm_norm_w_ = args.denorm_norm_w;
  
  pop_flag_vec_= args.pop_flag_vec;
  pop_wgt_vec_= args.pop_wgt_vec;
}


void Gene::RunJepeg(std::vector<Snp*>& gene_snp_vec){
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"################################ RunJepeg() Start " <<std::endl;
  if(gene_snp_vec.size()){
    Rcpp::Rcout<<"### Gene Name: "<<(*gene_snp_vec[0]).GetGeneid()<<std::endl;
  }
  Rcpp::Rcout<<std::endl;
#endif
  
  // count the number of SNPs belong to each categ.
  // TODO: make a function for the following
  std::vector<Snp*>::iterator it_gsv;
  for(it_gsv = gene_snp_vec.begin(); it_gsv != gene_snp_vec.end(); ++it_gsv){
    std::map< int, double >& categ_map = (*it_gsv)->GetCategMap();
    std::map< int, double >::iterator it_cm;
    for(it_cm = categ_map.begin(); it_cm != categ_map.end(); ++it_cm){
      if( (it_cm->first) == 0 ){       // Protein
        if(num_protein_ == 0)
          num_avail_categ_++;
        num_protein_++;
      } else if( (it_cm->first) == 1){ // TFBS
        if(num_tfbs_ == 0)
          num_avail_categ_++;
        num_tfbs_++;
      } else if( (it_cm->first) == 2){ // Within Hairpin
        if(num_wthhair_ == 0)
          num_avail_categ_++;
        num_wthhair_++;
      } else if( (it_cm->first) == 3){ // Within Target
        if(num_wthtar_ == 0)
          num_avail_categ_++;
        num_wthtar_++;
      } else if( (it_cm->first) == 4){ // cis eqtl
        if(num_cis_ == 0)
          num_avail_categ_++;
        num_cis_++;
      } else if( (it_cm->first) == 5){ // trans eqtl
        if(num_trans_ == 0)
          num_avail_categ_++;
        num_trans_++;
      }
    }	
  }//for


  num_snp_ = gene_snp_vec.size();
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# num_snp_: " <<num_snp_<<std::endl;
  Rcpp::Rcout<<"############# num_avail_categ_: " <<num_avail_categ_<<std::endl;
#endif
  
  
  categ_count_vec_.push_back(num_protein_);
  categ_count_vec_.push_back(num_tfbs_);
  categ_count_vec_.push_back(num_wthhair_);
  categ_count_vec_.push_back(num_wthtar_);
  categ_count_vec_.push_back(num_cis_);
  categ_count_vec_.push_back(num_trans_);
  
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# categ_count_vec_" << std::endl;
  PrintVector(categ_count_vec_);
#endif
  
  for(int i=0; i<categ_count_vec_.size(); i++){
    if(categ_count_vec_[i]){
      Categ categ(i, 0.0, false);
      categ_vec_.push_back(categ);
    }
  }
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# categ_vec_" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif
  
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############### Print SNPs in gene_snp_vec " <<std::endl;
  for(std::vector<Snp*>::iterator it_gsv=gene_snp_vec.begin(); it_gsv != gene_snp_vec.end(); ++it_gsv){
    (*it_gsv)->PrintSnpInfo();
  }
#endif
  
  //Calculate Jepeg p-value.
  CalJepegPval(gene_snp_vec);
  
}

void Gene::RunJepegmix(std::vector<Snp*>& gene_snp_vec){

#ifdef RunJepegmix_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"################################ RunJepegmix() Start " <<std::endl;
  if(gene_snp_vec.size()){
    Rcpp::Rcout<<"### Gene Name: "<<(*gene_snp_vec[0]).GetGeneid()<<std::endl;
  }
  Rcpp::Rcout<<std::endl;
#endif
  
  // count the number of SNPs belong to each categ.
  // TODO: make a function for the following
  std::vector<Snp*>::iterator it_gsv;
  for(it_gsv = gene_snp_vec.begin(); it_gsv != gene_snp_vec.end(); ++it_gsv){
    std::map< int, double >& categ_map = (*it_gsv)->GetCategMap();
    std::map< int, double >::iterator it_cm;
    for(it_cm = categ_map.begin(); it_cm != categ_map.end(); ++it_cm){
      if( (it_cm->first) == 0 ){       // Protein
        if(num_protein_ == 0)
          num_avail_categ_++;
        num_protein_++;
      } else if( (it_cm->first) == 1){ // TFBS
        if(num_tfbs_ == 0)
          num_avail_categ_++;
        num_tfbs_++;
      } else if( (it_cm->first) == 2){ // Within Hairpin
        if(num_wthhair_ == 0)
          num_avail_categ_++;
        num_wthhair_++;
      } else if( (it_cm->first) == 3){ // Within Target
        if(num_wthtar_ == 0)
          num_avail_categ_++;
        num_wthtar_++;
      } else if( (it_cm->first) == 4){ // cis eqtl
        if(num_cis_ == 0)
          num_avail_categ_++;
        num_cis_++;
      } else if( (it_cm->first) == 5){ // trans eqtl
        if(num_trans_ == 0)
          num_avail_categ_++;
        num_trans_++;
      }
    }	
  }//for
  
  
  num_snp_ = gene_snp_vec.size();
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# num_snp_: " <<num_snp_<<std::endl;
  Rcpp::Rcout<<"############# num_avail_categ_: " <<num_avail_categ_<<std::endl;
#endif
  
  
  categ_count_vec_.push_back(num_protein_);
  categ_count_vec_.push_back(num_tfbs_);
  categ_count_vec_.push_back(num_wthhair_);
  categ_count_vec_.push_back(num_wthtar_);
  categ_count_vec_.push_back(num_cis_);
  categ_count_vec_.push_back(num_trans_);
  
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# categ_count_vec_" << std::endl;
  PrintVector(categ_count_vec_);
#endif
  
  for(int i=0; i<categ_count_vec_.size(); i++){
    if(categ_count_vec_[i]){
      Categ categ(i, 0.0, false);
      categ_vec_.push_back(categ);
    }
  }
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############# categ_vec_" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif
  
  
#ifdef RunJepeg_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"############### Print SNPs in gene_snp_vec " <<std::endl;
  for(std::vector<Snp*>::iterator it_gsv=gene_snp_vec.begin(); it_gsv != gene_snp_vec.end(); ++it_gsv){
    (*it_gsv)->PrintSnpInfo();
  }
#endif
  
  //Calculate Jepeg p-value.
  
  CalJepegmixPval(gene_snp_vec);
  
}


void Gene::CalJepegPval(std::vector<Snp*>& gene_snp_vec){
  
  Eigen::MatrixXd CorG = Eigen::MatrixXd::Zero(gene_snp_vec.size(), gene_snp_vec.size()); // =CorZ
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(num_avail_categ_, gene_snp_vec.size()); // weight matrix
  Eigen::MatrixXd Wt = Eigen::MatrixXd::Zero(gene_snp_vec.size(), num_avail_categ_);
  Eigen::MatrixXd WWt = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_);
  Eigen::MatrixXd W_CorG = Eigen::MatrixXd::Zero(num_avail_categ_, gene_snp_vec.size());
  Eigen::MatrixXd CovU = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_); // CovU (covariance matrix of U)= W * CorG * Wt
  //To remove collinear categs
  Eigen::MatrixXd CorU = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_); // CorU (correlation matrix of U)
  
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(gene_snp_vec.size(), 1);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(num_avail_categ_, 1); 
  Eigen::MatrixXd Ut = Eigen::MatrixXd::Zero(1, num_avail_categ_); 
  Eigen::MatrixXd normU = Eigen::MatrixXd::Zero(num_avail_categ_, 1); //normalized U
  
  
  //Get CorG (correlation matrix[n by n] btw snp genotypes)  
  for(size_t i=0; i<static_cast<size_t>(CorG.rows()); i++){
    CorG(i, i) = 1.0 + lambda_; //add LAMBDA here (ridge regression trick)
    for(size_t j=i+1; j<static_cast<size_t>(CorG.rows()); j++){
      double v = CalCor((*gene_snp_vec[i]).GetGenotypeVec(), 
                        (*gene_snp_vec[j]).GetGenotypeVec());
      //v = floor(v*DECIMAL+0.5)/DECIMAL;
      CorG(i, j) = v;
      CorG(j, i) = v;
    }
  }
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# CorG" << std::endl;
  PrintMatrix(CorG);
#endif
  
  //Get W (weight score matrix[k by n])
  GetW(W, gene_snp_vec);
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# W" << std::endl;
  PrintMatrix(W);
#endif
  
  //Get Wt
  Wt = W.transpose();
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# Wt" << std::endl;
  PrintMatrix(Wt);
#endif 
  
  //Get WWt
  MpMatMat(WWt, W, Wt);
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# WWt" << std::endl;
  PrintMatrix(WWt);
#endif 
  
  //Get Cov U  
  MpMatMat(W_CorG, W, CorG);
  MpMatMat(CovU, W_CorG, Wt); //covariance matrix of U
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# CovU (W*CorG*Wt)" << std::endl;
  PrintMatrix(CovU);
#endif 
  
  //Get Cor U 
  CnvrtCovToCor(CorU, CovU); // covariance (correlation) matrix of U = correlation matrix of V
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# CorU" << std::endl;
  PrintMatrix(CorU);
#endif 
  
  //Get Z
  GetZ(Z, gene_snp_vec);
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# Z" << std::endl;
  PrintMatrix(Z);
#endif 
  
  //Get U
  MpMatMat(U, W, Z);
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# U (=W*Z)" << std::endl;
  PrintMatrix(U);
#endif 
  
  //Get normU (normalized U)
  for(int i=0; i< normU.rows(); i++){
    double var = CovU(i, i);
    double u = U(i, 0) / std::sqrt(var);
    normU(i, 0) = u;
    categ_vec_[i].SetPval(2 * R::pnorm5(std::abs(u), 0.0, 1.0, 0, 0));
  }
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# normU (nomalized U) = U/sqrt(varU)" << std::endl;
  PrintMatrix(normU);
#endif 
  
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector before removing" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif   
  
  //TODO: find collinear categs (to be removed).
  for(int j=num_avail_categ_-1; j>0; j--){
    for(int i=0; i<j; i++){
      double cor = CorU(i, j);
      if(std::abs(cor) > categ_cor_cutoff_){
        categ_vec_[j].SetRmv(true);
        break;
      }
    }
  }
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector after removing collinear categs" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif   
  
  //TODO: find categs with low variance (to be removed).  
  for(int i=0; i<num_avail_categ_; i++){
    double varU = CovU(i, i); 
    double normW = WWt(i, i)/denorm_norm_w_;
    if(varU < normW){
      categ_vec_[i].SetRmv(true);
    }
  }
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector after removing categs with low variance" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
  //PrintVector(rmv_categ_list);
#endif   
  
  //The number of available categs after removing collinear or low-variance categs
  for(int i=0; i<num_avail_categ_; i++){
    if(categ_vec_[i].GetRmv())
      num_rmv_categ_++;
  }
  df_ = num_avail_categ_ - num_rmv_categ_;
  
#ifdef CalJepegPval_Debug
  Rcpp::Rcout<<"############# The number of available/rmv categs after removing collinear or low-variance categs" << std::endl;
  Rcpp::Rcout<<"num_avail_categ_: "<<num_avail_categ_<<std::endl;
  Rcpp::Rcout<<"num_rmv_categ_: "<<num_rmv_categ_<<std::endl;
  Rcpp::Rcout<<"df_: "<<df_<<std::endl;
#endif   
  
  //////////////////////////////////////////////////////////////////////////////////// 
  //START
  //1. make new U (X) and CovU (CovX) after removing collinear or low variance categs
  //2. calculate JEPEG pvalue 
  if(df_){
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(df_, 1);
    Eigen::MatrixXd Xt = Eigen::MatrixXd::Zero(1, df_);
    Eigen::MatrixXd CovX = Eigen::MatrixXd::Zero(df_, df_);
    Eigen::MatrixXd CovXInv = Eigen::MatrixXd::Zero(df_, df_);
    Eigen::MatrixXd Xt_CovXInv = Eigen::MatrixXd::Zero(1, df_);
    Eigen::MatrixXd Xt_CovXInv_X = Eigen::MatrixXd::Zero(1, 1); //chisq
    
    //Remove rmv categs from U.
    //Make new vector X.
    int ii=0;
    for(int i=0; i<num_avail_categ_; i++){
      if(!(categ_vec_[i].GetRmv())){
        X(ii, 0) = U(i, 0);
        ii++;
      }
    }
    
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# X" << std::endl;
    PrintMatrix(X);
#endif
    
    //Remove rmv categs from CovU.
    //Make new matrix CovX.
    int n = 0;
    int m = 0;
    for(int i=0; i<num_avail_categ_; i++){
      if(!(categ_vec_[i].GetRmv())){
        for(int j=0; j<num_avail_categ_; j++){
          if(!(categ_vec_[j].GetRmv())){
            CovX(n, m) = CovU(i, j);
            m++;
          }
        }
        m=0;
        n++;
      }
    }
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# CovX" << std::endl;
    PrintMatrix(CovX);
#endif
    
    //Get Xt
    Xt = X.transpose();
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# Xt" << std::endl;
    PrintMatrix(Xt);
#endif
    
    //Get inverce matrix of CovX (CovXInv)
    MakePosDef(CovX, min_abs_eig_);
    InvMat(CovXInv, CovX);
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# CovXInv" << std::endl;
    PrintMatrix(CovXInv);
#endif 
    
    //Get JEPEG Chisq value (Xt * CovXInv * X)
    MpMatMat(Xt_CovXInv, Xt, CovXInv);
    MpMatMat(Xt_CovXInv_X, Xt_CovXInv, X);
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# JEPEG Chisq value = Xt * CovXInv * X" << std::endl;
    PrintMatrix(Xt_CovXInv_X);
#endif 
    
    chisq_ = Xt_CovXInv_X(0, 0);
    jepeg_pval_ = R::pchisq(chisq_, df_, 0, 0);
    
    //bonfe_pval_bf_ = CalBonfePval(U);
    //bonfe_pval_af_ = CalBonfePval(X);
    //sumU_pval_ = CalSumUPval(X, CovX);
    
    //////////////////////////////
    Categ top_categ = GetTopCateg(categ_vec_);
    top_categ_name_ = top_categ.GetName();
    top_categ_pval_ = top_categ.GetPval();
    
    Snp* top_snp = GetTopSNP(gene_snp_vec);
    top_snp_id_ = (*top_snp).GetRsid();
    top_snp_pval_ = 2 * R::pnorm5(std::abs((*top_snp).GetZ()), 0.0, 1.0, 0, 0);
    
    geneid_ = (*gene_snp_vec[0]).GetGeneid();
    
#ifdef CalJepegPval_Debug
    Rcpp::Rcout<<"############# JEPEG results" << std::endl;
    Rcpp::Rcout<<"geneid_: "<<geneid_<<std::endl;
    Rcpp::Rcout<<"chisq_: "<<chisq_<<std::endl;
    Rcpp::Rcout<<"df_: "<< df_ <<std::endl;
    Rcpp::Rcout<<"jepeg_pval_: "<<jepeg_pval_<<std::endl;
    
    Rcpp::Rcout<<"top_categ_name_: "<<top_categ_name_<<std::endl;
    Rcpp::Rcout<<"top_categ_pval_: "<<top_categ_pval_<<std::endl; // Z-score
    
    Rcpp::Rcout<<"top_snp_id_: "<<top_snp_id_<<std::endl;
    Rcpp::Rcout<<"top_snp_pval_: "<<top_snp_pval_<<std::endl;
    
    //Rcpp::Rcout<<"bonfe_pval_bf_: "<<bonfe_pval_bf_<<std::endl;
    //Rcpp::Rcout<<"bonfe_pval_af_: "<<bonfe_pval_af_<<std::endl;
    //Rcpp::Rcout<<"sumU_pval_: "<<sumU_pval_<<std::endl<<std::endl;
#endif   
    
    //END 
    //1. make new U (X) and CovU (CovX) after removing collinear or low variance categs
    //2. calculate JEPEG pvalue 
    ////////////////////////////////////////////////////////////////////////////////////
  }
  
}


void Gene::CalJepegmixPval(std::vector<Snp*>& gene_snp_vec){
  
  Eigen::MatrixXd CorG = Eigen::MatrixXd::Zero(gene_snp_vec.size(), gene_snp_vec.size()); // =CorZ
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(num_avail_categ_, gene_snp_vec.size()); // weight matrix
  Eigen::MatrixXd Wt = Eigen::MatrixXd::Zero(gene_snp_vec.size(), num_avail_categ_);
  Eigen::MatrixXd WWt = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_);
  Eigen::MatrixXd W_CorG = Eigen::MatrixXd::Zero(num_avail_categ_, gene_snp_vec.size());
  Eigen::MatrixXd CovU = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_); // CovU (covariance matrix of U)= W * CorG * Wt
  //To remove collinear categs
  Eigen::MatrixXd CorU = Eigen::MatrixXd::Zero(num_avail_categ_, num_avail_categ_); // CorU (correlation matrix of U)
  
  Eigen::MatrixXd Z = Eigen::MatrixXd::Zero(gene_snp_vec.size(), 1);
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(num_avail_categ_, 1); 
  Eigen::MatrixXd Ut = Eigen::MatrixXd::Zero(1, num_avail_categ_); 
  Eigen::MatrixXd normU = Eigen::MatrixXd::Zero(num_avail_categ_, 1); //normalized U
  
  Eigen::VectorXd SNP_STD_VEC = Eigen::VectorXd::Zero(gene_snp_vec.size()); //vector of standard deviations for SNPs in a gene
  
  for(size_t i=0; i < gene_snp_vec.size(); i++){
    double v = CalWgtCov((*gene_snp_vec[i]).GetGenotypeVec(), (*gene_snp_vec[i]).GetGenotypeVec(), pop_wgt_vec_);
    SNP_STD_VEC(i) = std::sqrt(v);
  }
  //Cal CorG (correlation matrix[n by n] btw snp genotypes)  
  for(size_t i=0; i<static_cast<size_t>(CorG.rows()); i++){
    CorG(i, i) = 1.0 + lambda_; //add LAMBDA here (ridge regression trick)
    double stdi = SNP_STD_VEC(i);
    for(size_t j=i+1; j<static_cast<size_t>(CorG.rows()); j++){
      double stdj = SNP_STD_VEC(j);
      double cov = CalWgtCov((*gene_snp_vec[i]).GetGenotypeVec(), (*gene_snp_vec[j]).GetGenotypeVec(), pop_wgt_vec_);
      double cor = cov/(stdi*stdj);
      CorG(i, j) = cor;
      CorG(j, i) = cor;
    }
  }
  
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# CorG" << std::endl;
  PrintMatrix(CorG);
#endif
  
  //Get W (weight score matrix[k by n])
  GetW(W, gene_snp_vec);
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# W" << std::endl;
  PrintMatrix(W);
#endif
  
  //Get Wt
  Wt = W.transpose();
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# Wt" << std::endl;
  PrintMatrix(Wt);
#endif 
  
  //Get WWt
  MpMatMat(WWt, W, Wt);
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# WWt" << std::endl;
  PrintMatrix(WWt);
#endif 
  
  //Get Cov U  
  MpMatMat(W_CorG, W, CorG);
  MpMatMat(CovU, W_CorG, Wt); //covariance matrix of U
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# CovU (W*CorG*Wt)" << std::endl;
  PrintMatrix(CovU);
#endif 
  
  //Get Cor U 
  CnvrtCovToCor(CorU, CovU); // covariance (correlation) matrix of U = correlation matrix of V
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# CorU" << std::endl;
  PrintMatrix(CorU);
#endif 
  
  //Get Z
  GetZ(Z, gene_snp_vec);
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# Z" << std::endl;
  PrintMatrix(Z);
#endif 
  
  //Get U
  MpMatMat(U, W, Z);
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# U (=W*Z)" << std::endl;
  PrintMatrix(U);
#endif 
  
  //Get normU (normalized U)
  for(int i=0; i< normU.rows(); i++){
    double var = CovU(i, i);
    double u = U(i, 0) / std::sqrt(var);
    normU(i, 0) = u;
    categ_vec_[i].SetPval(2 * R::pnorm5(std::abs(u), 0.0, 1.0, 0, 0));
  }
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# normU (nomalized U) = U/sqrt(varU)" << std::endl;
  PrintMatrix(normU);
#endif 
  
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector before removing" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif   
  
  //TODO: find collinear categs (to be removed).
  for(int j=num_avail_categ_-1; j>0; j--){
    for(int i=0; i<j; i++){
      double cor = CorU(i, j);
      if(std::abs(cor) > categ_cor_cutoff_){
        categ_vec_[j].SetRmv(true);
        break;
      }
    }
  }
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector after removing collinear categs" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
#endif   
  
  //TODO: find categs with low variance (to be removed).  
  for(int i=0; i<num_avail_categ_; i++){
    double varU = CovU(i, i); 
    double normW = WWt(i, i)/denorm_norm_w_;
    if(varU < normW){
      categ_vec_[i].SetRmv(true);
    }
  }
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# Print Categ Vector after removing categs with low variance" << std::endl;
  for(int i=0; i<num_avail_categ_; i++){
    categ_vec_[i].PrintInfo();
  }
  //PrintVector(rmv_categ_list);
#endif   
  
  //The number of available categs after removing collinear or low-variance categs
  for(int i=0; i<num_avail_categ_; i++){
    if(categ_vec_[i].GetRmv())
      num_rmv_categ_++;
  }
  df_ = num_avail_categ_ - num_rmv_categ_;
  
#ifdef CalJepegmixPval_Debug
  Rcpp::Rcout<<"############# The number of available/rmv categs after removing collinear or low-variance categs" << std::endl;
  Rcpp::Rcout<<"num_avail_categ_: "<<num_avail_categ_<<std::endl;
  Rcpp::Rcout<<"num_rmv_categ_: "<<num_rmv_categ_<<std::endl;
  Rcpp::Rcout<<"df_: "<<df_<<std::endl;
#endif   
  
  //////////////////////////////////////////////////////////////////////////////////// 
  //START
  //1. make new U (X) and CovU (CovX) after removing collinear or low variance categs
  //2. calculate JEPEG pvalue 
  if(df_){
    Eigen::MatrixXd X = Eigen::MatrixXd::Zero(df_, 1);
    Eigen::MatrixXd Xt = Eigen::MatrixXd::Zero(1, df_);
    Eigen::MatrixXd CovX = Eigen::MatrixXd::Zero(df_, df_);
    Eigen::MatrixXd CovXInv = Eigen::MatrixXd::Zero(df_, df_);
    Eigen::MatrixXd Xt_CovXInv = Eigen::MatrixXd::Zero(1, df_);
    Eigen::MatrixXd Xt_CovXInv_X = Eigen::MatrixXd::Zero(1, 1); //chisq
    
    //Remove rmv categs from U.
    //Make new vector X.
    int ii=0;
    for(int i=0; i<num_avail_categ_; i++){
      if(!(categ_vec_[i].GetRmv())){
        X(ii, 0) = U(i, 0);
        ii++;
      }
    }
    
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# X" << std::endl;
    PrintMatrix(X);
#endif
    
    //Remove rmv categs from CovU.
    //Make new matrix CovX.
    int n = 0;
    int m = 0;
    for(int i=0; i<num_avail_categ_; i++){
      if(!(categ_vec_[i].GetRmv())){
        for(int j=0; j<num_avail_categ_; j++){
          if(!(categ_vec_[j].GetRmv())){
            CovX(n, m) = CovU(i, j);
            m++;
          }
        }
        m=0;
        n++;
      }
    }
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# CovX" << std::endl;
    PrintMatrix(CovX);
#endif
    
    //Get Xt
    Xt = X.transpose();
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# Xt" << std::endl;
    PrintMatrix(Xt);
#endif
    
    //Get inverce matrix of CovX (CovXInv)
    MakePosDef(CovX, min_abs_eig_);
    InvMat(CovXInv, CovX);
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# CovXInv" << std::endl;
    PrintMatrix(CovXInv);
#endif 
    
    //Get JEPEG Chisq value (Xt * CovXInv * X)
    MpMatMat(Xt_CovXInv, Xt, CovXInv);
    MpMatMat(Xt_CovXInv_X, Xt_CovXInv, X);
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# JEPEG Chisq value = Xt * CovXInv * X" << std::endl;
    PrintMatrix(Xt_CovXInv_X);
#endif 
    
    chisq_ = Xt_CovXInv_X(0, 0);
    jepeg_pval_ = R::pchisq(chisq_, df_, 0, 0);
    
    //bonfe_pval_bf_ = CalBonfePval(U);
    //bonfe_pval_af_ = CalBonfePval(X);
    //sumU_pval_ = CalSumUPval(X, CovX);
    
    //////////////////////////////
    Categ top_categ = GetTopCateg(categ_vec_);
    top_categ_name_ = top_categ.GetName();
    top_categ_pval_ = top_categ.GetPval();
    
    Snp* top_snp = GetTopSNP(gene_snp_vec);
    top_snp_id_ = (*top_snp).GetRsid();
    top_snp_pval_ = 2 * R::pnorm5(std::abs((*top_snp).GetZ()), 0.0, 1.0, 0, 0);
    
    geneid_ = (*gene_snp_vec[0]).GetGeneid();
    
#ifdef CalJepegmixPval_Debug
    Rcpp::Rcout<<"############# JEPEG results" << std::endl;
    Rcpp::Rcout<<"geneid_: "<<geneid_<<std::endl;
    Rcpp::Rcout<<"chisq_: "<<chisq_<<std::endl;
    Rcpp::Rcout<<"df_: "<< df_ <<std::endl;
    Rcpp::Rcout<<"jepeg_pval_: "<<jepeg_pval_<<std::endl;
    
    Rcpp::Rcout<<"top_categ_name_: "<<top_categ_name_<<std::endl;
    Rcpp::Rcout<<"top_categ_pval_: "<<top_categ_pval_<<std::endl; // Z-score
    
    Rcpp::Rcout<<"top_snp_id_: "<<top_snp_id_<<std::endl;
    Rcpp::Rcout<<"top_snp_pval_: "<<top_snp_pval_<<std::endl;
    
    //Rcpp::Rcout<<"bonfe_pval_bf_: "<<bonfe_pval_bf_<<std::endl;
    //Rcpp::Rcout<<"bonfe_pval_af_: "<<bonfe_pval_af_<<std::endl;
    //Rcpp::Rcout<<"sumU_pval_: "<<sumU_pval_<<std::endl<<std::endl;
#endif   
    
    //END 
    //1. make new U (X) and CovU (CovX) after removing collinear or low variance categs
    //2. calculate JEPEG pvalue 
    ////////////////////////////////////////////////////////////////////////////////////
  }
  
}

double Gene::CalBonfePval(const Eigen::MatrixXd& U){
  Eigen::VectorXd pval_vec = Eigen::VectorXd::Zero(U.rows());
  double u = 0.0;
  for(int i=0; i<U.rows(); i++){
    u = U(i, 0);
    pval_vec(i) = 2 * R::pnorm5(std::abs(u), 0.0, 1.0, 0, 0) * U.rows();
  }
  double p = pval_vec.minCoeff();
  if(p > 1)
    p = 1.0;
  return p;
}

double Gene::CalSumUPval(const Eigen::MatrixXd& U, const Eigen::MatrixXd& CovU){
  double sumU=0;
  double var_sumU=0;
  for(int i=0; i<CovU.rows(); i++){
    sumU = sumU + U(i, 0);
    for(int j=0; j<CovU.rows(); j++){
      var_sumU = var_sumU + CovU(i ,j);
    }
  }
  double p = 2 * R::pnorm5(std::abs(sumU), 0.0, std::sqrt(var_sumU), 0, 0);
  return p;
}



void Gene::GetZ(Eigen::MatrixXd& Z, std::vector<Snp*>& gene_snp_vec){
  for(int i=0; i<gene_snp_vec.size(); i++){
    Z(i, 0) = (*gene_snp_vec[i]).GetZ();
  }
}


void Gene::GetW(Eigen::MatrixXd& W, std::vector<Snp*>& gene_snp_vec){
  
#ifdef GetW_Debug
  Rcpp::Rcout<<"############# categ_count_vec_" << std::endl;
  Rcpp::Rcout<<"size: "<<categ_count_vec_.size()<< std::endl;
  PrintVector(categ_count_vec_);
  Rcpp::Rcout<< std::endl;
#endif  
  for(int i=0; i<gene_snp_vec.size(); i++){
    int k=0;
    for(int j=0; j<categ_count_vec_.size(); j++){
      if(categ_count_vec_[j]!=0){
        W(k, i) = ((*gene_snp_vec[i]).GetCategWgt(j)) * std::sqrt(((*gene_snp_vec[i]).GetInfo()));
        //gsl_matrix_set(W, k, i, ((*gene_snp_vec[i]).GetCategWgt(j)));
        k++;
      }
    }
  }
}


Categ Gene::GetTopCateg(std::vector<Categ>& categ_vec){
  int top_index = 0;
  for(int i=0; i<num_avail_categ_; i++){
    double top_p = categ_vec[top_index].GetPval();
    double p = categ_vec[i].GetPval();
    bool rmv = categ_vec[i].GetRmv();
    if( (top_p > p) & !rmv ){
      top_index = i;
    }
  }
  return categ_vec[top_index];
}

//TODO: need to select Top SNP from only non-collinear categs
Snp* Gene::GetTopSNP(std::vector<Snp*>& gene_snp_vec){
  int top_index = 0;
  for(int i=0; i<gene_snp_vec.size(); i++){
    double top_absz = std::abs((*gene_snp_vec[top_index]).GetZ());
    double absz = std::abs((*gene_snp_vec[i]).GetZ());
    if(top_absz < absz){
      top_index = i;
    }
  }
  return gene_snp_vec[top_index];
}
