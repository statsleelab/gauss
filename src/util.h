//GAUSS : Genome Analysis Using Summary Statistics
//util.h

#ifndef UTIL_H
#define UTIL_H

#include <RcppEigen.h>

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

#include <Eigen/Dense>

#include "bgzf.h"

double CalCor(std::vector<std::string>& x, std::vector<std::string>& y); 
double CalWgtCov(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec); 
double CalCor(std::vector<std::string>& x, std::vector<std::string>& y, std::vector<double>& pop_wgt_vec); 
double CalCor(std::string& x, std::string& y);
double CalCor(std::vector<unsigned char>& x, std::vector<unsigned char>& y);
double CalCor(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y);


double CalCov(const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y);
double CalVar(const Eigen::Ref<const Eigen::VectorXd>& x);
double CalMeanSumSq(const Eigen::Ref<const Eigen::MatrixXd>& m);
void CalCorMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m);
void CalCovMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m);
void GetDiagMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m);
void MpMatMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m1, const Eigen::Ref<const Eigen::MatrixXd>& m2);
void MpNumMat(Eigen::MatrixXd& result, double num, const Eigen::Ref<const Eigen::MatrixXd>& m);
void CholeskyMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& m);
void SubMatMat(Eigen::MatrixXd& result, const Eigen::Ref<const Eigen::MatrixXd>& x1, const Eigen::Ref<const Eigen::MatrixXd>& x2);
void AddNumMatDiag(Eigen::MatrixXd& x1, double num);
void CnvrtCovToCor(Eigen::MatrixXd& corMat, const Eigen::Ref<const Eigen::MatrixXd>& covMat);
void InvMat(Eigen::MatrixXd& inverse, const Eigen::Ref<const Eigen::MatrixXd>& m1);
void MakePosDef(Eigen::MatrixXd& m1, double minAbsEig);
int RmvPC(Eigen::MatrixXd& m1, double eigCutoff);
int CountPC(const Eigen::Ref<const Eigen::MatrixXd>& m1, double eig_cutoff);
void PrintVector(const Eigen::VectorXd& x);
void PrintVector(std::vector<bool>& x);
void PrintVector(std::vector<int>& x);
void PrintVector(const Eigen::VectorXd& x, std::ofstream& outFile);
void PrintMatrix(const Eigen::MatrixXd& x);
void PrintMatrix(const Eigen::MatrixXd& x, std::ofstream& outFile);

void SplitString(std::vector<std::string>& vec, std::string& str, char delim);  

void LoadProgressBar( int percent );

void FlipGenotype(std::string& genotype);
void FlipGenotypeVec(std::vector<std::string>& genotype_vec);

int BgzfGetLine(BGZF* fp, std::string& line);

//void CalPopWgtVec(std::vector<Snp*>& snp_vec, BGZF* fp, Options& opts);


#endif
