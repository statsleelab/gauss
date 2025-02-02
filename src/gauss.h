//GAUSS : Genome Analysis Using Summary Statistics
//gauss.h

#ifndef GAUSS_H
#define GAUSS_H

#include <Rcpp.h>

//forward declaration
class Snp;

#include <string>
#include <vector>
#include <map>
#include "bgzf.h"


struct Arguments{
  
  Arguments();
  void PrintArguments();

  //JEPEGMIX
  //bool imputation_flag;

  int chr;
  long long int start_bp;
  long long int end_bp;
  long long int wing_size;
  std::string study_pop;
  
  std::string input_file;
  std::string reference_index_file; //contains snpid, chr, bp, a1, a2, fp 
  std::string reference_data_file;  //contains genotype strings and af1s
  std::string reference_pop_desc_file; //contains population descriptions  
  std::string annotation_file;

  std::vector<std::string> ref_pop_vec;
  std::vector<int> ref_pop_size_vec;
  std::vector<std::string> ref_sup_pop_vec;
  
  /////##//////////////
  // Hidden Arguments//
  ///////##////////////
  
  double lambda; // ridge regression estimate
  double min_abs_eig; // minimum absolute eigen value used in spectral decomposition.
  double eig_cutoff; // for rmvPC

  // for Mix
  double mix_af1_cutoff;
  int interval;
  std::vector<int> pop_flag_vec;
  std::vector<double> pop_wgt_vec;
  std::map<std::string, double> pop_wgt_map;
  double sum_pop_wgt;
  int num_pops;
  int num_samples;
  
  double af1_cutoff; // allele1 frequency cutoff
  int min_num_measured_snp; // minimum number of meanusred SNPs
  int min_num_unmeasured_snp; // minimum number of unmeasured SNPs
  
  //JEPEG/MIX
  int total_num_categ;
  double categ_cor_cutoff;
  int denorm_norm_w; // denorminator of norm of W
  double imp_info_cutoff; // imputation information cutoff for JEPEG p-value calculation.
};


class MapKey{
 public:
  MapKey(int chr, long long int bp, std::string a1, std::string a2)
    : chr_(chr), bp_(bp), a1_(a1), a2_(a2) {}

  bool operator<(const MapKey &right) const{
    if( chr_ == right.chr_ ){
        if( bp_ == right.bp_ ){
	        if(a1_ == right.a1_){
	          return a2_ < right.a2_;
	        } else {
	          return a1_ < right.a1_;
	        }
        } else {
	        return bp_ < right.bp_;
        }
    } else {
      return chr_ < right.chr_;
    }
  }//operator<

 private:
  int chr_;
  long long int bp_;//base pair position
  std::string a1_;  //reference allele
  std::string a2_;  //alternative allele

};//class MapKey


struct LessThanMapKey{
  inline bool operator() (const MapKey& key1, const MapKey& key2) const{
    return key1.operator<(key2);
  }
};

struct StartEnd{
  std::vector<Snp*>::iterator start_it;
  std::vector<Snp*>::iterator end_it;
};

// used in dist&distmix with All=false, jepeg&jepegmix with All=true
void ReadInputZ(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args, bool All);
// used in afmix
void ReadInputAf(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
// used in dist, distmix, qcat, qcatmix Read SNPs in a genomic region
void ReadReferenceIndex(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
// used in afmix, jepeg, jepegmix Read all SNPs
void ReadReferenceIndexAll(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);

void MakeSnpVec(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
void MakeSnpVecMix(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args);
void ReadGenotype(std::vector<Snp*>& snp_vec, Arguments& args);
void ReadGenotypeOne(Snp* snp, Arguments& args);
void FreeGenotype(std::vector<Snp*>& snp_vec);
void FreeGenotypeOne(Snp* snp);
void read_ref_desc(Arguments& args);
void init_pop_flag_vec(Arguments& args);
void init_pop_flag_wgt_vec(Arguments& args);

void UpdateSnpToMinorAllele(std::vector<Snp*>& snp_vec);
std::vector<std::string> ConvertGenotypesToRecessive(const std::vector<std::string>& genoVec);

//JEPEG & JEPEGMIX  
void ReadAnnotation(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args); 
void MakeGeneStartEndVec(std::vector<StartEnd>& gene_start_end_vec, std::vector<Snp*>& snp_vec);

 #endif
