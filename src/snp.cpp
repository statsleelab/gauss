//GAUSS : Genome Analysis Using Summary Statistics
//snp.cpp

#include "snp.h"

#include <string>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <iostream>

Snp::Snp()
{
  chr_ = -1;
	bp_ = -1;
	rsid_ = ".";
	a1_ = ".";
	a2_ = ".";
	af1mix_ = -1.0;
  af1study_ = -1.0;
	af1ref_ = -1.0;
	z_ = 0.0;
	info_ = -1.0;

	qcat_chisq_ = 0.0;

	type_ = -1;
	fpos_ = -1;
  flip_ = false;

	geneid_ = ".";
}

double Snp::GetCategWgt(int categ_num){
  if(categ_map_.find(categ_num)!=categ_map_.end())
    return categ_map_.find(categ_num)->second;
    //return 1;
  else
    return 0;
}

void Snp::RmvCateg(int categ_num){
  std::map < int, double >::iterator it = categ_map_.find(categ_num);
  categ_map_.erase(it);
}

void Snp::PrintSnpInfo(){
  std::cout<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1study_<<" "<<af1ref_<<" ";
  std::cout<<z_<<" "<<info_<<" "<<type_<<" "<<fpos_<<" "<<flip_<<" "<<geneid_<<std::endl;
  //print categ_map_
  for(std::map< int, double >::iterator it = categ_map_.begin(); it != categ_map_.end(); ++it){
    std::cout<<it->first<<" "<<it->second<<" ";
  }
  std::cout<<std::endl;  
  //print genotype strings
  for(std::vector<std::string>::iterator it = genotype_vec_.begin(); it != genotype_vec_.end(); ++it){
    std::cout<<*it<<std::endl;
  }
}

void Snp::PrintSnpInfo(std::ofstream& outfile){
  outfile<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1study_<<" "<<af1ref_<<" ";
  outfile<<z_<<" "<<info_<<" "<<type_<<" "<<fpos_<<" "<<flip_<<" "<<geneid_<<std::endl;
  //print categ_map_
  for(std::map< int, double >::iterator it = categ_map_.begin(); it != categ_map_.end(); ++it){
    outfile<<it->first<<" "<<it->second<<" ";
  }
  outfile<<std::endl;  
  //print genotype strings
  for(std::vector<std::string>::iterator it = genotype_vec_.begin(); it != genotype_vec_.end(); ++it){
    outfile<<*it<<std::endl;
  }
}

void Snp::PrintDistResult(){
  std::cout<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1ref_<<" ";
  std::cout<<z_<<" "<<info_<<" "<<type_<<std::endl;
}

void Snp::PrintDistResult(std::ofstream& outfile){
  outfile<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1ref_<<" ";
  outfile<<z_<<" "<<info_<<" "<<type_<<std::endl;
}

void Snp::PrintQcatResult(){
  std::cout<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1ref_<<" ";
  std::cout<<z_<<" "<<qcat_chisq_<<" "<<qcat_pval_<<" "<<type_<<std::endl;
}


void Snp::PrintQcatResult(std::ofstream& outfile){
  outfile<<rsid_<<" "<<chr_<<" "<<bp_<<" "<<a1_<<" "<<a2_<<" "<<af1ref_<<" ";
  outfile<<z_<<" "<<qcat_chisq_<<" "<<qcat_pval_<<" "<<type_<<std::endl;
}


void Snp::ClearSnp(){
  categ_map_.clear();
  //genotype_vec_.clear();
  //std::vector<std::string>().swap(genotype_vec_);
}
