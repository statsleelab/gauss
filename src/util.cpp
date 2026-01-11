//GAUSS : Genome Analysis Using Summary Statistics
//util.cpp


#include "util.h"

#include <cstdlib> //exit(), EXIT_FAILURE ..
#include <cmath>   //std::sqrt()
#include <iostream>//endl
#include <sstream> //istringstream
#include <iomanip> //setw()
//#include <map>
//#include <sstream>
#include <string>
#include <vector>

#include <RcppEigen.h>

#include "gauss.h"
//#include "bgzf.h"
#include "snp.h"


/*
double CalCor(std::string& x, std::string& y){
 
 int n = x.length();
 double ex=0,ey=0,xt=0,yt=0,sxx=0,syy=0,sxy=0;
 
 for(int i=0; i<n; i++){
 ex += (int)(x[i]-'0');
 ey += (int)(y[i]-'0');
 }
 ex /= n;
 ey /= n;
 
 for(int i=0; i<n; i++){
 xt = (int)(x[i]-'0') - ex;
 yt = (int)(y[i]-'0') - ey;
 sxx += xt * xt;
 syy += yt * yt;
 sxy += xt * yt;
 }
 return sxy/std::sqrt(sxx*syy);
}
*/

//This function is used when --mix flag is not specified. 
double CalCor(std::vector<std::string>& x, std::vector<std::string>& y){
  int num_samples = 0;
  int n = x.size();
  double xij=0, yij=0, sumx=0, sumy=0, sumxsq=0, sumysq=0, sumxy=0;
  for(int i=0; i<n; i++){
    int m = x[i].length();
    for(int j=0; j<m; j++){
      xij = (double)(x[i][j]-'0');
      yij = (double)(y[i][j]-'0');
      sumx += xij;
      sumy += yij;
      sumxsq += xij*xij;
      sumysq += yij*yij;
      sumxy += xij*yij;
    }
    num_samples += m;
  }
  double numer = num_samples*sumxy-sumx*sumy;
  double denor = std::sqrt((num_samples)*sumxsq-sumx*sumx)*std::sqrt((num_samples)*sumysq-sumy*sumy);
  double r = numer/denor;
  return r;
}

/**
 * @brief Calculates the weighted covariance between two genotype data vectors.
 * 
 * This function computes a weighted sum of covariances between the genotype strings 
 * provided in vectors x and y across multiple populations. Each element in x and y is 
 * a string representing the genotype data for a particular population, where each 
 * character (typically '0', '1', or '2') represents the count of the reference allele 
 * for an individual.
 * 
 * For each population, the function:
 *  - Computes the sum of genotype counts from x and y.
 *  - Computes the sum of cross-products of the genotype counts.
 *  - Uses the formula with an unbiased correction factor: (m/(m-1)) * (m * sumxy - sumx * sumy),
 *    where m is the number of individuals (length of the genotype string).
 *  - Multiplies the result by the corresponding population weight from pop_wgt_vec.
 * 
 * The function then adjusts for the weighted means of the populations by combining the 
 * accumulated values and returns the final weighted covariance.
 * 
 * @param x A vector of strings, each containing genotype calls (as digits) for a SNP in a specific population.
 * @param y A vector of strings, each containing genotype calls (as digits) for a paired SNP corresponding to the populations in x.
 * @param pop_wgt_vec A vector of population weights corresponding to each population in x and y.
 * 
 * @return double The computed weighted covariance between the genotype data vectors.
 * 
 * @note The function assumes that each genotype string has a length greater than 1. 
 * 
 * @warning Ensure that the genotype strings only contain valid digits ('0', '1', or '2').
 */
//This function is used when --mix flag is specified.
//DISTMIX (ImputeMix) uses this function.
double CalWgtCov(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec){
  int n = x.size();
  double wsumcov=0, wsum_mi_mj=0, wsum_mi=0, wsum_mj=0;
  for(int i=0; i<n; i++){
    int m = x[i].length();
    double xij=0, yij=0, sumx=0, sumy=0, sumxy=0;
    double wgt_val = pop_wgt_vec[i];
    for(int j=0; j<m; j++){
      xij = (double)(x[i][j]-'0');
      yij = (double)(y[i][j]-'0');
      sumx += xij;
      sumy += yij;
      sumxy += xij*yij;
    }
    double factor = ((double)m) / (m - 1); // ensure floating-point division
    wsumcov += wgt_val*factor*(m*sumxy-sumx*sumy);
    wsum_mi_mj += wgt_val*(sumx/m)*(sumy/m);
    wsum_mi += wgt_val*(sumx/m);
    wsum_mj += wgt_val*(sumy/m);
  }
  return (wsumcov + wsum_mi_mj - wsum_mi*wsum_mj);
}

//This function is used when --mix flag is specified.
//DISTMIX (ImputeMix) uses this function.
double CalCor(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec){
  //int num_samples = 0;
  int n = x.size();
  double xij=0, yij=0, sumwx=0, sumwy=0, sumwxsq=0, sumwysq=0, sumwxy=0;
  for(int i=0; i<n; i++){
    int m = x[i].length();
    double wgt_val = pop_wgt_vec[i]/m;
    for(int j=0; j<m; j++){
      xij = (double)(x[i][j]-'0');
      yij = (double)(y[i][j]-'0');
      sumwx += wgt_val*xij;
      sumwy += wgt_val*yij;
      sumwxsq += wgt_val*xij*xij;
      sumwysq += wgt_val*yij*yij;
      sumwxy += wgt_val*xij*yij;
    }
    //num_samples += m;
  }
  double numer = sumwxy-sumwx*sumwy;
  double denor = std::sqrt(sumwxsq-sumwx*sumwx)*std::sqrt(sumwysq-sumwy*sumwy);
  double r = numer/denor;
  return r;
}


double CalCor(std::string& x, std::string& y){
  int n = x.length();
  double xi=0, yi=0, sumx=0, sumy=0, sumxsq=0, sumysq=0, sumxy=0;
  for(int i=0; i<n; i++){
    xi = (double)(x[i]-'0');
    yi = (double)(y[i]-'0');
    sumx += xi;
    sumy += yi;
    sumxsq += xi*xi;
    sumysq += yi*yi;
    sumxy += xi*yi;
  }
  double numer = n*sumxy-sumx*sumy;
  double denor = std::sqrt((n)*sumxsq-sumx*sumx)*std::sqrt((n)*sumysq-sumy*sumy);
  double r = numer/denor;
  return r;
}

double CalCor(std::vector<unsigned char>& x, std::vector<unsigned char>& y){
  
  int n = x.size();
  double ex=0,ey=0,xt=0,yt=0,sxx=0,syy=0,sxy=0;
  
  for(int i=0; i<n; i++){
    ex += (int)(x[i]-'0');
    ey += (int)(y[i]-'0');
  }
  ex /= n;
  ey /= n;
  
  for(int i=0; i<n; i++){
    xt = (int)(x[i]-'0') - ex;
    yt = (int)(y[i]-'0') - ey;
    sxx += xt * xt;
    syy += yt * yt;
    sxy += xt * yt;
  }
  return sxy/std::sqrt(sxx*syy);
}


double CalCor(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y){
  double mean_x = x.mean();
  double mean_y = y.mean();
  Eigen::ArrayXd dx = x.array() - mean_x;
  Eigen::ArrayXd dy = y.array() - mean_y;
  double sxx = (dx * dx).sum();
  double syy = (dy * dy).sum();
  double sxy = (dx * dy).sum();
  return sxy / std::sqrt(sxx * syy);
}

double CalCov(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y){
  double mean_x = x.mean();
  double mean_y = y.mean();
  Eigen::ArrayXd dx = x.array() - mean_x;
  Eigen::ArrayXd dy = y.array() - mean_y;
  double sxy = (dx * dy).sum();
  double denom = static_cast<double>(x.size() - 1);
  return sxy / denom;
}

double CalVar(const Eigen::Ref<const Eigen::VectorXd>& x){
  double mean_x = x.mean();
  Eigen::ArrayXd dx = x.array() - mean_x;
  double denom = static_cast<double>(x.size() - 1);
  return (dx * dx).sum() / denom;
}

double CalMeanSumSq(const Eigen::Ref<const Eigen::MatrixXd>& m){
  int n = static_cast<int>(m.rows());
  double ss = 0;
  for(int i=0; i < n; i++){
    ss += m(i, 0) * m(i, 0);
  }
  return ss / n;
}

void CalCorMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m){ //calculate correlation between columns of m
  double cor;
  for(int i=0; i < m.cols(); i++){
    for(int j=i; j < m.cols(); j++){
      cor = CalCor(m.col(i), m.col(j));
      result(i, j) = cor;
      if(i != j)
        result(j, i) = cor;
    }
  }
}

void CalCovMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m){ //calculate covariance between columns of m
  double cov;
  for(int i=0; i < m.cols(); i++){
    for(int j=i; j < m.cols(); j++){
      cov = CalCov(m.col(i), m.col(j));
      result(i, j) = cov;
      if(i != j)
        result(j, i) = cov;
    }
  }
}

void GetDiagMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m){
  result.setZero();
  for(int i=0; i < m.rows(); i++){
    result(i, i) = m(i, i);
  }
}

void MpMatMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m1, const Eigen::Ref<const Eigen::MatrixXd>& m2){
  result.noalias() = m1 * m2;
}

void MpNumMat(Eigen::MatrixXd& result, double num, const Eigen::Ref<const Eigen::MatrixXd>& m){
  result = m * num;
}

//Cholesky decomposition (M=L*Lt) return L 
void CholeskyMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m){
  Eigen::LLT<Eigen::MatrixXd> llt(m);
  result = llt.matrixL();
}

void SubMatMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& x1, const Eigen::Ref<const Eigen::MatrixXd>& x2){
  result = x1 - x2;
}

void AddNumMatDiag(Eigen::MatrixXd& x1, double num){
  x1.diagonal().array() += num;
}

void CnvrtCovToCor(Eigen::MatrixXd& cor_mat, const Eigen::Ref<const Eigen::MatrixXd>& cov_mat){
  double std1, std2, cor;
  for(int i=0; i<cov_mat.rows(); i++){
    for(int j=i; j<cov_mat.cols(); j++){
      std1 = std::sqrt(cov_mat(i, i));
      std2 = std::sqrt(cov_mat(j, j));
      cor = cov_mat(i, j) / (std1 * std2);
      cor_mat(i, j) = cor;
      if(i!=j)
        cor_mat(j, i) = cor;
    }
  }
}

void InvMat(Eigen::MatrixXd& inverse, const Eigen::Ref<const Eigen::MatrixXd>& m1){
  inverse = m1.fullPivLu().inverse();
}

void MakePosDef(Eigen::MatrixXd& m1, double min_abs_eig){
  int size = static_cast<int>(m1.rows());
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m1);
  if(solver.info() != Eigen::Success){
    return;
  }
  Eigen::VectorXd eig_vals = solver.eigenvalues();
  Eigen::MatrixXd eig_vecs = solver.eigenvectors();
  if(eig_vals.minCoeff() < min_abs_eig){
    for(int i=0; i<size; i++){
      if(eig_vals(i) < min_abs_eig){
        eig_vals(i) = min_abs_eig;
      }
    }
    m1 = eig_vecs * eig_vals.asDiagonal() * eig_vecs.transpose();
  }
}

int RmvPC(Eigen::MatrixXd& m1, double eig_cutoff){
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"RmvPC start!"<<std::endl;
#endif
  
  int size = static_cast<int>(m1.rows());
  int num_eig = size;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m1);
  if(solver.info() != Eigen::Success){
    return num_eig;
  }
  Eigen::VectorXd eig_vals = solver.eigenvalues();
  Eigen::MatrixXd eig_vecs = solver.eigenvectors();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(size, size);
  
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"min_eig_val: "<<eig_vals(0)<<std::endl;
  Rcpp::Rcout<<"eig_cutoff: "<<eig_cutoff<<std::endl;
#endif
  
  if(eig_vals(0) < eig_cutoff){  // if the smallest eig value is less than eig_cutoff, investigate all eig values.
    for(int i=0; i<size; i++){
      if(eig_vals(i) > eig_cutoff){	        
        res.noalias() += eig_vals(i) * eig_vecs.col(i) * eig_vecs.col(i).transpose();
      } else {
        num_eig--;
      }
    }
    m1 = res;
  }
  return num_eig;
}

int CountPC(const Eigen::Ref<const Eigen::MatrixXd>& m1, double eig_cutoff){
#ifdef CountPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"CountPC start!"<<std::endl;
#endif
  
  int size = static_cast<int>(m1.rows());
  int num_eig = size;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m1);
  if(solver.info() != Eigen::Success){
    return num_eig;
  }
  Eigen::VectorXd eig_vals = solver.eigenvalues();
  
  //gsl_matrix* res = gsl_matrix_calloc(size, size);
  //gsl_matrix* eig_vec_sq = gsl_matrix_calloc(size, size);
  //gsl_matrix_view eig_vec;
  //gsl_matrix* eig_vec_trans = gsl_matrix_calloc(1, size);
  
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"min_eig_val: "<<eig_vals(0)<<std::endl;
  Rcpp::Rcout<<"eig_cutoff: "<<eig_cutoff<<std::endl;
#endif
  
  if(eig_vals(0) < eig_cutoff){  // if the smallest eig value is less than eig_cutoff, investigate all eig values.
    for(int i=0; i<size; i++){
      if(eig_vals(i) < eig_cutoff){	        
        num_eig--;
      }
    }
  }
  return num_eig;
}


void PrintVector(const Eigen::VectorXd& x){
  for(int i=0 ; i<x.size() ; i++){
    Rcpp::Rcout<<" "<<x(i);
  }
  Rcpp::Rcout<<std::endl;
}

void PrintVector(const Eigen::VectorXd& x, std::ofstream& outFile){
  for(int i=0 ; i<x.size() ; i++){
    outFile<<" "<<x(i);
  }
  outFile<<std::endl;
}

void PrintVector(std::vector<bool>& x){
  for(size_t i=0 ; i<x.size() ; i++){
    Rcpp::Rcout<<" "<<x[i];
  }
  Rcpp::Rcout<<std::endl;
}

void PrintVector(std::vector<int>& x){
  for(size_t i=0 ; i<x.size() ; i++){
    Rcpp::Rcout<<" "<<x[i];
  }
  Rcpp::Rcout<<std::endl;
}

void PrintMatrix(const Eigen::MatrixXd& x){
  for(int i=0 ; i< x.rows() ; i++){
    for(int j=0 ; j<x.cols() ; j++){
      Rcpp::Rcout<<" "<<x(i,j);
    }
    Rcpp::Rcout<<std::endl;  
  }
  Rcpp::Rcout<<std::endl;
}

void PrintMatrix(const Eigen::MatrixXd& x, std::ofstream& outFile){
  for(int i=0 ; i< x.rows() ; i++){
    for(int j=0 ; j<x.cols() ; j++){
      outFile<<" "<<x(i,j);
    }
    outFile<<std::endl;  
  }
  outFile<<std::endl;
}

void SplitString(std::vector<std::string>& vec, 
                 std::string& str, char delim){
  std::stringstream ss(str);
  std::string token;
  while(std::getline(ss, token, delim)){
    vec.push_back(token);
  }
}


void LoadProgressBar( int percent ){
  std::string bar;
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  Rcpp::Rcout<< "\r" <<percent << "%  " <<"[" << bar << "] " <<std::flush;
}

void FlipGenotype(std::string& genotype){
  
  int n = genotype.length();
  for(int i=0; i<n; i++){
    if(genotype[i] == '0')
      genotype[i] = '2';
    else if (genotype[i] == '2')
      genotype[i] = '0';
  }
}

void FlipGenotypeVec(std::vector<std::string>& genotype_vec){
  int n = genotype_vec.size();
  for(int i=0; i<n; i++){
    int m = genotype_vec[i].length();
    for(int j=0; j<m; j++){
      if(genotype_vec[i][j] == '0')
        genotype_vec[i][j] = '2';
      else if (genotype_vec[i][j] == '2')
        genotype_vec[i][j] = '0';
    }
  }
}


int BgzfGetLine(BGZF* fp, std::string& line){
  line.erase();
  int i=0;
  int c;
  while(true){
    i++;
    c = bgzf_getc(fp);
    if(c == -2){
      Rcpp::stop("Error: can't read "+ std::to_string(i) +"-th character");
    }
    if(c == -1){ // end of file                                                                                 
      break;
    }
    if(c == 10){ // end of line                                                                                 
      break;
    }
    line += static_cast<char>(c);
  }
  return c;
}

/*
void CalPopWgtVec(std::vector<Snp*>& snp_vec, BGZF* fp, Options& opts){
  
  
  int size = snp_vec.size(); 
  int interval = opts.interval;
  
#ifdef CalPopWgt_Debug
  //Rcpp::Rcout<<"CalPopWgtVec 1"<<std::endl;
  Rcpp::Rcout<<"snp_vec size: "<<size<<std::endl;
  Rcpp::Rcout<<"interval: "<<interval<<std::endl;
  Rcpp::Rcout<<"num_pops: "<<opts.num_pops<<std::endl;  
#endif
  
  gsl_matrix* af1_cor_mat = gsl_matrix_calloc(opts.num_pops+1, opts.num_pops+1);
  gsl_matrix* W_mat_i = gsl_matrix_calloc(opts.num_pops, 1);
  gsl_matrix* W_mat = gsl_matrix_calloc(opts.num_pops, 1);
  gsl_matrix_view Cxy;
  gsl_matrix_view Cxx;
  gsl_matrix* CxxInv = gsl_matrix_calloc(opts.num_pops, opts.num_pops);
  
  //Rcpp::Rcout<<"CalPopWgtVec 2"<<std::endl;
  
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
    gsl_matrix* af1_mat = gsl_matrix_calloc(snp_subvec_size, opts.num_pops+1);
    for(int j=0; j<snp_subvec_size; j++){
      gsl_matrix_set(af1_mat, j, 0, snp_subvec[j]->GetAf1study());    
      
      //Rcpp::Rcout<<"af1study :"<<snp_subvec[j]->GetAf1study();
      
      std::string line;
      bgzf_seek(fp, (*snp_subvec[j]).GetFpos(), SEEK_SET);
      BgzfGetLine(fp, line);
      std::istringstream buffer(line);
      //skipping genotypes
      for(int k=0; k<opts.pop_flag_vec.size(); k++){
        std::string tmp;
        buffer >> tmp;
      }
      //reading allele frequencies
      int kk=1;
      for(int k=0; k<opts.pop_flag_vec.size(); k++){
        double af1;
        buffer >> af1;
        if(opts.pop_flag_vec[k]){
          gsl_matrix_set(af1_mat, j, kk, af1);
          kk++;
          //Rcpp::Rcout<<" "<<af1<<" ";
        }
      }
      //Rcpp::Rcout<<std::endl;
    }
    
#ifdef CalPopWgt_Debug
    Rcpp::Rcout<<"af1_mat:"<<std::endl;
    PrintMatrix(af1_mat);
    Rcpp::Rcout<<std::endl;
#endif
    
    //Calculate corr matrix of afs.
    //CalCorMat(af1_cor_mat, af1_mat);
    CalCovMat(af1_cor_mat, af1_mat);
    
#ifdef CalPopWgt_Debug
    Rcpp::Rcout<<"af1_cor_mat:"<<std::endl;
    PrintMatrix(af1_cor_mat);
    Rcpp::Rcout<<std::endl;
#endif
    
    //Calculate W
    Cxy = gsl_matrix_submatrix(af1_cor_mat, 1, 0, opts.num_pops, 1);
    Cxx = gsl_matrix_submatrix(af1_cor_mat, 1, 1, opts.num_pops, opts.num_pops);
    MakePosDef(&Cxx.matrix, opts.min_abs_eig);
    InvMat(CxxInv, &Cxx.matrix);
    MpMatMat(W_mat_i, CxxInv, &Cxy.matrix);
    
#ifdef CalPopWgt_Debug
    Rcpp::Rcout<<"Cxy:"<<std::endl;
    PrintMatrix(&Cxy.matrix);
    Rcpp::Rcout<<std::endl;
    
    Rcpp::Rcout<<"Cxx:"<<std::endl;
    PrintMatrix(&Cxx.matrix);
    Rcpp::Rcout<<std::endl;
    
    Rcpp::Rcout<<"CxxInv:"<<std::endl;
    PrintMatrix(CxxInv);
    Rcpp::Rcout<<std::endl;
    
    Rcpp::Rcout<<"W_mat_i:"<<std::endl;
    PrintMatrix(W_mat_i);
    Rcpp::Rcout<<std::endl;
#endif
    
    for(int j=0; j<opts.num_pops; j++){
      double wval = gsl_matrix_get(W_mat, j, 0) + gsl_matrix_get(W_mat_i, j, 0)/interval;
      gsl_matrix_set(W_mat, j, 0, wval);
    }
    
#ifdef CalPopWgt_Debug
    Rcpp::Rcout<<"W_mat:"<<std::endl;
    PrintMatrix(W_mat);
    Rcpp::Rcout<<std::endl;
#endif
    
    gsl_matrix_free(af1_mat);
  }
  // if w is less than zero then make it zero.
  for(int i=0; i<opts.num_pops; i++){
    double wval = gsl_matrix_get(W_mat, i, 0);
    if(wval < 0)
      gsl_matrix_set(W_mat, i, 0, 0);
    else
      gsl_matrix_set(W_mat, i, 0, floor(wval*1000+0.5)/1000);
  }
  
  double sum_pop_wgt = 0.0;
  for(int i=0; i<opts.pop_flag_vec.size(); i++){
    double wval = gsl_matrix_get(W_mat, i, 0);
    if(wval > 0){
      opts.pop_wgt_vec.push_back(wval);
      sum_pop_wgt += wval;
    } else {
      opts.pop_flag_vec[i]=0;
      opts.num_pops--;
    }
    //if(i==0) opts.pop_wgt_map["ASW"]=wval;
    //else if (i==1) opts.pop_wgt_map["CEU"]=wval;
    //else if (i==2) opts.pop_wgt_map["CHB"]=wval;
    //else if (i==3) opts.pop_wgt_map["CHS"]=wval;
    //else if (i==4) opts.pop_wgt_map["CLM"]=wval;
    //else if (i==5) opts.pop_wgt_map["FIN"]=wval;
    //else if (i==6) opts.pop_wgt_map["GBR"]=wval;
    //else if (i==7) opts.pop_wgt_map["IBS"]=wval;
    //else if (i==8) opts.pop_wgt_map["JPT"]=wval;
    //else if (i==9) opts.pop_wgt_map["LWK"]=wval;
    //else if (i==10) opts.pop_wgt_map["MXL"]=wval;
    //else if (i==11) opts.pop_wgt_map["PUR"]=wval;
    //else if (i==12) opts.pop_wgt_map["TSI"]=wval;
    //else if (i==13) opts.pop_wgt_map["YRI"]=wval;
  }
  opts.sum_pop_wgt = sum_pop_wgt;
  

//   //calculate sum of wgt
//   double wgt_sum = 0.0;
//   for(int i=0; i<opts.num_pops; i++){
//   wgt_sum += gsl_matrix_get(W_mat, i, 0);
//   }
//   //renormalize them
//   double wgt_sum_after = 0.0;
//   for(int i=0; i<opts.num_pops; i++){
//   double wval = gsl_matrix_get(W_mat, i, 0)/wgt_sum;
//   opts.pop_wgt_vec.push_back(wval);
//   wgt_sum_after += wval;
//   }

#ifdef CalPopWgt_Debug
  Rcpp::Rcout<<"sum_pop_wgt: "<< sum_pop_wgt << std::endl;
  //Rcpp::Rcout<<"wgt_sum before: "<< wgt_sum << std::endl;
  //Rcpp::Rcout<<"wgt_sum after : "<< wgt_sum_after << std::endl;
#endif
  
///   double wgt_sum = 0.0;
///   for(int i=0; i<opts.num_pops; i++)
///   wgt_sum += gsl_matrix_get(W_mat, i, 0);
   
///   double wgt_sum_after = 0.0;
///   for(int i=0; i<opts.num_pops; i++){
///   double wval = gsl_matrix_get(W_mat, i, 0)/wgt_sum;
///   if(wval < 0)
///   wval=0;
///   pop_wgt_vec.push_back(wval);
///   wgt_sum_after += wval;
///   }
   
///   Rcpp::Rcout<<"wgt_sum before: "<< wgt_sum << std::endl;
///   Rcpp::Rcout<<"wgt_sum after : "<< wgt_sum_after << std::endl;

  
  gsl_matrix_free(af1_cor_mat);
  gsl_matrix_free(W_mat_i);
  gsl_matrix_free(W_mat);
  gsl_matrix_free(CxxInv);
}
*/
