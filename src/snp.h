//GAUSS : Genome Analysis Using Summary Statistics
//snp.h

#ifndef SNP_H
#define SNP_H

#include <string>
#include <map>
#include <vector>
#include <fstream>


class Snp{

public:
  Snp();

  //get
  std::string GetRsid() {return rsid_; }
  int GetChr() { return chr_; }
  long long int GetBp() const { return bp_; }
  std::string GetA1() { return a1_; }
  std::string GetA2() { return a2_; }
  double GetAf1Mix() { return af1mix_; }    // DISTMIX, QCATMIX, JEPEGMIX output
  double GetAf1Study() {return af1study_; } // only used in calculating pop wgt, cal_pop_wgt.cpp 
  double GetAf1Ref() { return af1ref_; }    // DIST, QCAT, JEPEG output
  double GetZ() {return z_; }
  double GetInfo() {return info_; }

  double GetQcatChisq() { return qcat_chisq_; } //QCAT chisq;
  
  int GetType() const {return type_;}
  long long int GetFpos() { return fpos_; }
  bool GetFlip() { return flip_; }

  std::string GetGeneid() const { return geneid_; }
  double GetCategWgt(int categ_num);
  std::map<int, double>& GetCategMap() { return categ_map_; }
  std::vector<std::string>& GetGenotypeVec() { return genotype_vec_; }
  

  //set
  void SetRsid(std::string rsid) { rsid_ = rsid; }
  void SetChr(int chr) { chr_ = chr; }
  void SetBp(long long int bp) { bp_ = bp; }
  void SetA1(std::string a1) { a1_ = a1; }
  void SetA2(std::string a2) { a2_ = a2; }
  void SetAf1Mix(double af1mix) { af1mix_ = af1mix; } // DISTMIX, QCATMIX, JEPEGMIX output
  void SetAf1Study(double af1study) { af1study_ = af1study; }// only used in calculating pop wgt, cal_pop_wgt.cpp
  void SetAf1Ref(double af1ref) { af1ref_ = af1ref; }// DIST, QCAT, JEPEG output
  void SetZ(double z) { z_ = z; }
  void SetInfo(double info) { info_ = info; }

  void SetQcatChisq(double qcat_chisq) { qcat_chisq_ = qcat_chisq; } //QCAT

  void SetType(int type) { type_ = type; } //0: unmeasured/ref, 1: measured/ref, 2: measured/no ref
  void SetFpos(long long int fpos) { fpos_ = fpos; }  
  void SetFlip(bool flip) { flip_ = flip; }  

  void SetGeneid(std::string geneid) { geneid_ = geneid; }
  void SetCateg(int categ_num, double wgt) { categ_map_[categ_num]=wgt; }

  void SetGenotypeVec(std::vector<std::string>& genotype_vec) { genotype_vec_ = genotype_vec; }

  //remove
  void RmvCateg(int categ_num);

  //print
  void PrintSnpInfo();
  void PrintSnpInfo(std::ofstream&);
  void PrintDistResult();
  void PrintDistResult(std::ofstream&);

  void PrintQcatResult();
  void PrintQcatResult(std::ofstream&);


  //clear
  void ClearSnp();


private:

  std::string rsid_;
  int chr_;
  long long int bp_;
  std::string a1_,a2_;
  double af1mix_; // DISTMIX, QCATMIX, JEPEGMIX output
  double af1study_;// only used in calculating pop wgt, cal_pop_wgt.cpp
  double af1ref_;// DIST, QCAT, JEPEG output
  double z_;
  double info_;

  double qcat_chisq_;

  int type_; // 0:unmeasured snp/1KG, 1:measured snp/1KG, 2:measured snp/No 1KG.
  long long int fpos_;
  bool flip_;

  std::string geneid_;
  std::map< int, double > categ_map_;
  std::vector< std::string > genotype_vec_;
};


struct LessThanBp{
  inline bool operator() (const Snp& obj1, const Snp& obj2){
    return (obj1.GetBp() < obj2.GetBp());
  }	
};

struct LessThanType{	
  inline bool operator() (const Snp* obj1, const Snp* obj2){
    return ( (*obj1).GetType() < (*obj2).GetType() );
  }
};

struct GreaterThanType{	
  inline bool operator() (const Snp* obj1, const Snp* obj2){
    return ( (*obj1).GetType() > (*obj2).GetType() );
  }
};

struct LessThanGeneid{
  inline bool operator() (const Snp* obj1, const Snp* obj2){
    return ( (*obj1).GetGeneid() < (*obj2).GetGeneid() );
  }
};

#endif
