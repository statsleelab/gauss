//GAUSS : Genome Analysis Using Summary Statistics
//util.h

#ifndef UTIL_H
#define UTIL_H

//#include <RcppGSL.h>
#include <Rcpp.h>

// forward declaration
//struct BGZF;

//struct Options;
//struct LessThanMapKey;
//class Snp;
//class MapKey;
 
#include <string>
#include <vector>
//#include <map>
#include <fstream> //ofstream

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "bgzf.h"

double CalCor(std::vector<std::string>& x, std::vector<std::string>& y); 
double CalWgtCov(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec); 
double CalCor(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec); 
double CalCor(std::string& x, std::string& y);
double CalCor(std::vector<unsigned char>& x, std::vector<unsigned char>& y);
double CalCor(gsl_vector_view x, gsl_vector_view y);


double CalCov(gsl_vector_view x, gsl_vector_view y);
double CalVar(gsl_vector_view x);
double CalMeanSumSq(gsl_matrix* m);
void CalCorMat(gsl_matrix* result, gsl_matrix* m);
void CalCovMat(gsl_matrix* result, gsl_matrix* m);
void GetDiagMat(gsl_matrix* result, gsl_matrix* m);
void MpMatMat(gsl_matrix* result, gsl_matrix* m1, gsl_matrix* m2);
void MpNumMat(gsl_matrix* result, double num, gsl_matrix* m);
void CholeskyMat(gsl_matrix* result, gsl_matrix* m);
void SubMatMat(gsl_matrix* result, gsl_matrix* x1, gsl_matrix* x2);
void AddNumMatDiag(gsl_matrix* x1, double num);
void CnvrtCovToCor(gsl_matrix* corMat, gsl_matrix* covMat);
void InvMat(gsl_matrix* inverse, gsl_matrix* m1);
void MakePosDef(gsl_matrix* m1, double minAbsEig);
int RmvPC(gsl_matrix* m1, double eigCutoff);
int CountPC(gsl_matrix* m1, double eig_cutoff);
void PrintVector(gsl_vector* x);
void PrintVector(std::vector<bool>& x);
void PrintVector(std::vector<int>& x);
void PrintVector(gsl_vector* x, std::ofstream& outFile);
void PrintMatrix(gsl_matrix* x);
void PrintMatrix(gsl_matrix* x, std::ofstream& outFile);

void SplitString(std::vector<std::string>& vec, std::string& str, char delim);  

void LoadProgressBar( int percent );

void FlipGenotype(std::string& genotype);
void FlipGenotypeVec(std::vector<std::string>& genotype_vec);

int BgzfGetLine(BGZF* fp, std::string& line);

//void CalPopWgtVec(std::vector<Snp*>& snp_vec, BGZF* fp, Options& opts);


#endif
