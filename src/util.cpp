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

//#include <RcppGSL.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics_double.h>

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
//This function is used when --mix flag is specified.
//DISTMIX (ImputeMix) uses this function.
//This function calculates weighted sum of Covariances
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
    wsumcov += wgt_val*(m/(m-1))*(m*sumxy-sumx*sumy);
    wsum_mi_mj += wgt_val*(sumx/m)*(sumy/m);
    wsum_mi += wgt_val*(sumx/m);
    wsum_mj += wgt_val*(sumy/m);
  }
  return (wsumcov + wsum_mi_mj - wsum_mi*wsum_mj);
  //return wsumcov;
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


double CalCor(gsl_vector_view x, gsl_vector_view y){
  
  size_t n = (&x.vector)->size;
  double cor = gsl_stats_correlation( x.vector.data, x.vector.stride, 
                                      y.vector.data, y.vector.stride,
                                      n);
  return cor;
  
}


double CalCov(gsl_vector_view x, gsl_vector_view y){
  
  size_t n = (&x.vector)->size;
  double cov =  gsl_stats_covariance( x.vector.data, x.vector.stride, 
                                      y.vector.data, y.vector.stride,
                                      n);
  return cov;
  
}

double CalVar(gsl_vector_view x){
  
  size_t n = (&x.vector)->size;
  double var =  gsl_stats_variance(x.vector.data, x.vector.stride, n);
  return var;
}
double CalMeanSumSq(gsl_matrix* m){
  size_t n = m->size1;
  double ss = 0;
  for(int i=0; i < n; i++){
    ss += gsl_matrix_get(m, i, 0)*gsl_matrix_get(m, i, 0);
  }
  return ss/n;
}

void CalCorMat(gsl_matrix* result, gsl_matrix* m){ //calculate correlation between columns of m
  
  double cor;
  for(int i=0; i < m->size2; i++){
    for(int j=i; j < m->size2; j++){
      cor = CalCor(gsl_matrix_column(m,i), gsl_matrix_column(m,j));
      gsl_matrix_set(result, i, j, cor);
      if(i != j)
        gsl_matrix_set(result, j, i, cor);
    }
  } 
}

void CalCovMat(gsl_matrix* result, gsl_matrix* m){ //calculate covariance between columns of m
  
  double cov;
  for(int i=0; i < m->size2; i++){
    for(int j=i; j < m->size2; j++){
      cov = CalCov(gsl_matrix_column(m,i), gsl_matrix_column(m,j));
      gsl_matrix_set(result, i, j, cov);
      if(i != j)
        gsl_matrix_set(result, j, i, cov);
    }
  } 
}

void GetDiagMat(gsl_matrix* result, gsl_matrix* m){
  for(int i=0; i < m->size1; i++){
    gsl_matrix_set(result, i, i, gsl_matrix_get(m, i, i));
  } 
}

void MpMatMat(gsl_matrix* result, gsl_matrix* m1, gsl_matrix* m2){
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m1, m2, 0.0, result);
}

void MpNumMat(gsl_matrix* result, double num, gsl_matrix* m){
  for(int i=0; i < m->size1; i++){
    for(int j=0; j < m->size2; j++){
      gsl_matrix_set(result, i, j, gsl_matrix_get(m, i, j)*num);
    }
  } 
}

//Cholesky decomposition (M=L*Lt) return L 
void CholeskyMat(gsl_matrix* result, gsl_matrix* m){
  gsl_matrix_memcpy(result, m);
  gsl_linalg_cholesky_decomp(result); 
  //make upper triangular part zero to make result as a lower triangular matrix.  
  for(int i=0; i < m->size1; i++){
    for(int j=i+1; j < m->size2; j++){
      gsl_matrix_set(result, i, j, 0);
    }
  }
}
void SubMatMat(gsl_matrix* result, gsl_matrix* x1, gsl_matrix* x2){
  gsl_matrix_memcpy(result, x1);
  gsl_matrix_sub(result, x2);
}

void AddNumMatDiag(gsl_matrix* x1, double num){
  for(int i=0; i < x1->size1; i++){
    gsl_matrix_set(x1, i, i, gsl_matrix_get(x1, i, i) + num);
  }
}

void CnvrtCovToCor(gsl_matrix* cor_mat, gsl_matrix* cov_mat){
  
  double std1, std2, cor;
  for(int i=0; i<cov_mat->size1; i++){
    for(int j=i; j<cov_mat->size2; j++){
      std1 = std::sqrt(gsl_matrix_get(cov_mat, i, i));
      std2 = std::sqrt(gsl_matrix_get(cov_mat, j, j));
      cor = gsl_matrix_get(cov_mat, i, j)/(std1*std2);
      gsl_matrix_set(cor_mat, i, j, cor);
      if(i!=j)
        gsl_matrix_set(cor_mat, j, i, cor);
    }
  }
}

void InvMat(gsl_matrix* inverse, gsl_matrix* m1){
  int s;
  int size = m1->size1;
  gsl_matrix* tmp = gsl_matrix_calloc(m1->size1, m1->size2);
  gsl_matrix_memcpy(tmp, m1);
  //makePosDef(tmp, 1e-3); //
  gsl_permutation* perm = gsl_permutation_alloc(size);
  gsl_linalg_LU_decomp(tmp, perm, &s);
  gsl_linalg_LU_invert(tmp, perm, inverse);
  gsl_permutation_free(perm);
  gsl_matrix_free(tmp);
}


void MakePosDef(gsl_matrix* m1, double min_abs_eig){
  int size = m1->size1;
  gsl_matrix* tmp = gsl_matrix_calloc(m1->size1, m1->size2);
  gsl_matrix_memcpy(tmp, m1);
  
  gsl_vector* eig_vals = gsl_vector_calloc(size);
  gsl_matrix* eig_vecs = gsl_matrix_calloc(size, size);
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(tmp, eig_vals, eig_vecs, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eig_vals, eig_vecs, GSL_EIGEN_SORT_VAL_ASC); // ascending order in numerical value.
  
  gsl_matrix* res = gsl_matrix_calloc(size, size);
  gsl_matrix* eig_vec_sq = gsl_matrix_calloc(size, size);
  gsl_matrix_view eig_vec;
  gsl_matrix* eig_vec_trans = gsl_matrix_calloc(1, size);
  
  if(gsl_vector_get(eig_vals,0) < min_abs_eig){
    
    for(int i=0; i<size; i++){
      if(gsl_vector_get(eig_vals,i) < min_abs_eig){	      
        gsl_vector_set(eig_vals, i, min_abs_eig);    
      }  
      eig_vec = gsl_matrix_submatrix(eig_vecs, 0, i, size, 1);
      gsl_matrix_transpose_memcpy(eig_vec_trans, &eig_vec.matrix);
      MpMatMat(eig_vec_sq, &eig_vec.matrix, eig_vec_trans);
      gsl_matrix_scale(eig_vec_sq, gsl_vector_get(eig_vals, i));
      gsl_matrix_add(res, eig_vec_sq);
    }
    gsl_matrix_memcpy(m1, res);
  }
  
  gsl_matrix_free(tmp);
  gsl_vector_free(eig_vals);
  gsl_matrix_free(eig_vecs);
  gsl_matrix_free(res);
  gsl_matrix_free(eig_vec_sq);
  gsl_matrix_free(eig_vec_trans);
}

int RmvPC(gsl_matrix* m1, double eig_cutoff){
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"RmvPC start!"<<std::endl;
#endif
  
  int size = m1->size1;
  int num_eig = size;
  gsl_matrix* tmp = gsl_matrix_calloc(m1->size1, m1->size2);
  gsl_matrix_memcpy(tmp, m1);
  
  gsl_vector* eig_vals = gsl_vector_calloc(size);
  gsl_matrix* eig_vecs = gsl_matrix_calloc(size, size);
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(tmp, eig_vals, eig_vecs, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eig_vals, eig_vecs, GSL_EIGEN_SORT_VAL_ASC); // ascending order in numerical value.
  
  gsl_matrix* res = gsl_matrix_calloc(size, size);
  gsl_matrix* eig_vec_sq = gsl_matrix_calloc(size, size);
  gsl_matrix_view eig_vec;
  gsl_matrix* eig_vec_trans = gsl_matrix_calloc(1, size);
  
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"min_eig_val: "<<gsl_vector_get(eig_vals,0)<<std::endl;
  Rcpp::Rcout<<"eig_cutoff: "<<eig_cutoff<<std::endl;
#endif
  
  if(gsl_vector_get(eig_vals,0) < eig_cutoff){  // if the smallest eig value is less than eig_cutoff, investigate all eig values.
    for(int i=0; i<size; i++){
      if(gsl_vector_get(eig_vals,i) > eig_cutoff){	        
        eig_vec = gsl_matrix_submatrix(eig_vecs, 0, i, size, 1);
        gsl_matrix_transpose_memcpy(eig_vec_trans, &eig_vec.matrix);
        MpMatMat(eig_vec_sq, &eig_vec.matrix, eig_vec_trans);
        gsl_matrix_scale(eig_vec_sq, gsl_vector_get(eig_vals, i));
        gsl_matrix_add(res, eig_vec_sq);
      } else {
        num_eig--;
      }
    }
    gsl_matrix_memcpy(m1, res);
  }
  gsl_matrix_free(tmp);
  gsl_vector_free(eig_vals);
  gsl_matrix_free(eig_vecs);
  gsl_matrix_free(res);
  gsl_matrix_free(eig_vec_sq);
  gsl_matrix_free(eig_vec_trans);
  return num_eig;
}

int CountPC(gsl_matrix* m1, double eig_cutoff){
#ifdef CountPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"CountPC start!"<<std::endl;
#endif
  
  int size = m1->size1;
  int num_eig = size;
  gsl_matrix* tmp = gsl_matrix_calloc(m1->size1, m1->size2);
  gsl_matrix_memcpy(tmp, m1);
  
  gsl_vector* eig_vals = gsl_vector_calloc(size);
  gsl_matrix* eig_vecs = gsl_matrix_calloc(size, size);
  gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(size);
  gsl_eigen_symmv(tmp, eig_vals, eig_vecs, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eig_vals, eig_vecs, GSL_EIGEN_SORT_VAL_ASC); // ascending order in numerical value.
  
  //gsl_matrix* res = gsl_matrix_calloc(size, size);
  //gsl_matrix* eig_vec_sq = gsl_matrix_calloc(size, size);
  //gsl_matrix_view eig_vec;
  //gsl_matrix* eig_vec_trans = gsl_matrix_calloc(1, size);
  
#ifdef RmvPC_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"min_eig_val: "<<gsl_vector_get(eig_vals,0)<<std::endl;
  Rcpp::Rcout<<"eig_cutoff: "<<eig_cutoff<<std::endl;
#endif
  
  if(gsl_vector_get(eig_vals,0) < eig_cutoff){  // if the smallest eig value is less than eig_cutoff, investigate all eig values.
    for(int i=0; i<size; i++){
      if(gsl_vector_get(eig_vals,i) < eig_cutoff){	        
        num_eig--;
      }
    }
  }
  gsl_matrix_free(tmp);
  gsl_vector_free(eig_vals);
  gsl_matrix_free(eig_vecs);
  //gsl_matrix_free(res);
  //gsl_matrix_free(eig_vec_sq);
  //gsl_matrix_free(eig_vec_trans);
  return num_eig;
}


void PrintVector(gsl_vector* x){
  for(size_t i=0 ; i<x->size ; i++){
    Rcpp::Rcout<<" "<<gsl_vector_get(x,i);
  }
  Rcpp::Rcout<<std::endl;
}

void PrintVector(gsl_vector* x, std::ofstream& outFile){
  for(size_t i=0 ; i<x->size ; i++){
    outFile<<" "<<gsl_vector_get(x,i);
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

void PrintMatrix(gsl_matrix* x){
  for(size_t i=0 ; i< x->size1 ; i++){
    for(size_t j=0 ; j<x->size2 ; j++){
      Rcpp::Rcout<<" "<<gsl_matrix_get(x,i,j);
    }
    Rcpp::Rcout<<std::endl;  
  }
  Rcpp::Rcout<<std::endl;
}

void PrintMatrix(gsl_matrix* x, std::ofstream& outFile){
  for(size_t i=0 ; i< x->size1 ; i++){
    for(size_t j=0 ; j<x->size2 ; j++){
      outFile<<" "<<gsl_matrix_get(x,i,j);
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
