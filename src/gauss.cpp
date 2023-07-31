//GAUSS : Genome Analysis Using Summary Statistics
//gauss.cpp

#include "gauss.h"

#include <iostream>
#include <fstream> // ifstream
#include <string>
#include <sstream>
#include <algorithm>
#include "util.h"
#include "bgzf.h"
#include "snp.h"

//#define ReadInput_Debug
//#define ReadReferenceIndex_Debug

Arguments::Arguments(){

  lambda = 0.1;
  min_abs_eig = 1e-5;
  eig_cutoff = 0.01;
  mix_af1_cutoff = 0.05;
  interval = 1000;
  sum_pop_wgt = 1;

  min_num_measured_snp = 10;
  min_num_unmeasured_snp = 10;
  
  // JEPEG/MIX
  total_num_categ = 6;
  categ_cor_cutoff = 0.8;
  denorm_norm_w = 3;
  imp_info_cutoff = 0.3;
}

void Arguments::PrintArguments(){

  Rcpp::Rcout << "chromosome: "<< chr << std::endl;
  Rcpp::Rcout << "start_bp: " << start_bp <<std::endl;
  Rcpp::Rcout << "end_bp: " << end_bp << std::endl;
  Rcpp::Rcout << "wing_size: "<< wing_size << std::endl;

  Rcpp::Rcout << "input_file: "<< input_file << std::endl;
  Rcpp::Rcout << "reference_index_file: "<< reference_index_file << std::endl;
  Rcpp::Rcout << "reference_data_file: "<< reference_data_file << std::endl;
  Rcpp::Rcout << "reference_pop_desc_file: "<< reference_pop_desc_file << std::endl;

  //JEPEGMIX
  //Rcpp::Rcout << "annotation: "<<annotation_file << std::endl;

  //Hidden
  Rcpp::Rcout << "lambda: " << lambda << std::endl;
  Rcpp::Rcout << "min_abs_eig: " << min_abs_eig << std::endl;
  Rcpp::Rcout << "eig_cutoff: " << eig_cutoff << std::endl;

  // for mix
  Rcpp::Rcout << "mix_af1_cutoff: " << mix_af1_cutoff << std::endl;
  Rcpp::Rcout << "interval: " << interval << std::endl; 

  Rcpp::Rcout << "ref_pop_vec: ";
  for(int i=0; i<ref_pop_vec.size(); i++){
    Rcpp::Rcout << ref_pop_vec[i] << " ";
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "ref_pop_size_vec: ";
  for(int i=0; i<ref_pop_size_vec.size(); i++){
    Rcpp::Rcout << ref_pop_size_vec[i] << " ";
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "ref_sup_pop_vec: ";
  for(int i=0; i<ref_sup_pop_vec.size(); i++){
    Rcpp::Rcout << ref_sup_pop_vec[i] << " ";
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "pop_flag_vec: ";
  for(int i=0; i<pop_flag_vec.size(); i++){
    Rcpp::Rcout << pop_flag_vec[i] << " ";
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "pop_wgt_vec: ";
  for(int i=0; i<pop_wgt_vec.size(); i++){
    Rcpp::Rcout << pop_wgt_vec[i] << " ";
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "pop_wgt_map: ";
  for(std::map<std::string, double>::iterator it = pop_wgt_map.begin(); it != pop_wgt_map.end(); ++it){
    Rcpp::Rcout<<it->first<<" "<<it->second<<std::endl;
  }
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout << "sum_pop_wgt: " << sum_pop_wgt << std::endl;
  Rcpp::Rcout << "num_pops: " << num_pops << std::endl; 
  
  Rcpp::Rcout << "af1_cutoff: " << af1_cutoff << std::endl;

  // for JEPEG
  Rcpp::Rcout << "total_num_categ: " << total_num_categ << std::endl;
  Rcpp::Rcout << "categ_cor_cutoff: " << categ_cor_cutoff << std::endl;
  Rcpp::Rcout << "denorm_norm_w: " << denorm_norm_w << std::endl;
  Rcpp::Rcout << "imp_info_cutoff: " << imp_info_cutoff << std::endl;
  
}


void ReadInput(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  
  Rcpp::Rcout<<"Reading input...";
  Rcpp::Rcout.flush();

  std::string input_file = args.input_file;
  std::ifstream in_input(input_file.c_str());
  
  if(!in_input){
    Rcpp::stop("ERROR: can't open input file '"+input_file+"'");
  }

  std::string line;
  std::string rsid, a1, a2;
  int chr;
  long long int bp;
  double z;
  double info = 1.0;
  Snp* snp;

  std::getline(in_input, line); //read header of input file.  
  while(std::getline(in_input, line)){
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> z;

    if( (args.chr > 0) && (args.chr != chr) )  // if only one chromosome is specified with -c option and
      continue;                                   // the user specified chr is not equal to chr, skip it.

    if((args.start_bp - args.wing_size) > bp || (args.end_bp + args.wing_size) < bp)
      continue;
    
    snp = new Snp();
    snp->SetRsid(rsid);
    snp->SetChr(chr);
    snp->SetBp(bp);
    snp->SetA1(a1);
    snp->SetA2(a2);
    snp->SetZ(z);
    snp->SetInfo(info);
    snp->SetType(2); // 2: measured SNP but does not exist in reference data.

    MapKey mkey(chr, bp, a1, a2);
    snp_map[mkey]=snp;
  }//while
  in_input.close();
  Rcpp::Rcout<<std::endl;
}


void ReadReferenceIndex(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  Rcpp::Rcout<<"Reading reference index...";
  Rcpp::Rcout.flush();
  
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;
  
  std::string reference_index_file = args.reference_index_file;
  BGZF* fp = bgzf_open(reference_index_file.c_str(), "r");
  
  if(!fp){
    Rcpp::stop("ERROR: can't open reference index file '"+reference_index_file+"'");
  }
  
  int last_char;
  std::string line;

  std::string rsid, a1, a2;
  int chr;
  double af1ref;
  long long int bp, fpos;
  Snp* snp;
  
  while(true){
    
    last_char = BgzfGetLine(fp, line);
    if(last_char == -1) //EOF
      break;
    
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    if( (args.chr > 0) && (args.chr != chr) )  // if only one chromosome is specified with -c option and
      continue;                                   // the user specified chr is not equal to chr, skip it.
    
    if((args.start_bp - args.wing_size) > bp || (args.end_bp + args.wing_size) < bp)
      continue;
      
      MapKey mkey1(chr, bp, a1, a2);
      MapKey mkey2(chr, bp, a2, a1);
      
      it1 = snp_map.find(mkey1);
      it2 = snp_map.find(mkey2);
      
      if((it1 != snp_map.end()) && (it2 == snp_map.end())){ // snp exists in input and a1=a1 & a2=a2.
        
        (it1->second)->SetRsid(rsid);
        (it1->second)->SetType(1); // change Type to 1:"measured and exists in reference panel"
        //(it1->second)->SetAf1Ref(af1ref);
        (it1->second)->SetFpos(fpos); // position of snp genotype string in reference panel genotype data
        
      } else if((it1 == snp_map.end()) && (it2 != snp_map.end())){ // if snp exists in input but a1=a2 & a2=a1.
        
        (it2->second)->SetRsid(rsid);
        (it2->second)->SetA1(a1);
        (it2->second)->SetA2(a2);
        (it2->second)->SetZ( (it2->second)->GetZ()*(-1) ); //change the sign of z-score
        (it2->second)->SetType(1); // change Type to 1:"measured and exists in reference panel"
        //(it2->second)->SetAf1study( 1 - (it2->second)->GetAf1study() );
        //(it2->second)->SetAf1Ref(af1ref);
        (it2->second)->SetFpos(fpos); // position of snp genotype string in reference panel genotype data
        
        MapKey new_key(chr, bp, a1, a2); //makes a new key for the modified snp object.
        snp_map[new_key] = it2->second;  //makes a map element with the new key. 
        snp_map.erase(it2);              //delete snp with the old key. 
        
      } else if(it1 == snp_map.end() && it2 == snp_map.end()){ // if the snp does not exist in input
        
        //make new snp
        snp = new Snp();
        snp->SetRsid(rsid);
        snp->SetChr(chr);
        snp->SetBp(bp);
        snp->SetA1(a1);
        snp->SetA2(a2);
        snp->SetType(0); // change Type to 0:"unmeasured and exists in reference panel"
        //snp->SetAf1Ref(af1ref);
        snp->SetFpos(fpos);
        snp_map[mkey1]=snp; // add the new snp in snp_map.
        
      } else { //Throw error message! This should not happen.
        Rcpp::stop("ERROR: input file contains duplicates");
      }//if
  }//while
  
  bgzf_close(fp);
  Rcpp::Rcout<<std::endl;
}


// Used in DIST, QCAT, jepeg
void MakeSnpVec(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+args.reference_data_file+"'");
  }
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  double af1ref = 0;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end(); ++it_sm){
    bgzf_seek(fp, (it_sm->second)->GetFpos(), SEEK_SET);
    std::string line;
    BgzfGetLine(fp, line);
    std::istringstream buffer(line);
    double allele_counter=0, num_subj=0;
    for(int k=0; k<args.num_pops; k++){
      std::string geno_str;
      buffer >> geno_str;
      if(args.pop_flag_vec[k]){
        num_subj += geno_str.length();
        for(int i=0; i<geno_str.length(); i++){
          allele_counter += (double)(geno_str[i]-'0');
        }
      }
    }
    af1ref = allele_counter/(2*num_subj);
    af1ref = std::ceil(af1ref*100000.0)/100000.0;  //round up to 5 decimal places 
    (it_sm->second)->SetAf1Ref(af1ref);
    if( (af1ref > args.af1_cutoff) && (af1ref < (1-args.af1_cutoff)) ){
      snp_vec.push_back(it_sm->second);
    }      
  }
  bgzf_close(fp); //closes BGZF file connnection.
}


// Used in DISTMIX, QCATMIX, jepegMIX
// This function computes weighted sum of the reference allele frequencies (af1_mix) of each SNP and 
// store SNPs with af1_mix > af1_cutoff or af1_mix < 1-af1_cutoff and store them in SNP vector.   

void MakeSnpVecMix(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+args.reference_data_file+"'");
  }
  
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end(); ++it_sm){
    std::vector<double> af1_vec;
    double af1_mix = 0;
    std::string line;
    bgzf_seek(fp, (it_sm->second)->GetFpos(), SEEK_SET);
    BgzfGetLine(fp, line);
    std::istringstream buffer(line);
    //skipping genotypes
    for(int k=0; k<args.pop_flag_vec.size(); k++){
      std::string tmp;
      buffer >> tmp;
    }
    //reading reference allele frequencies
    for(int k=0; k<args.pop_flag_vec.size(); k++){
      double af1;
      buffer >> af1;
      if(args.pop_flag_vec[k])
        af1_vec.push_back(af1);
    }
    for(int k=0; k<af1_vec.size(); k++){
      af1_mix += af1_vec[k]*args.pop_wgt_vec[k];
    }
    if( (af1_mix > args.af1_cutoff) && (af1_mix < (1-args.af1_cutoff)) ){
      //Rcpp::Rcout<<"af1_mix: "<<af1_mix<<", af1_cutoff: "<<args.af1_cutoff<<std::endl; 
      (it_sm->second)->SetAf1Mix(af1_mix);
      snp_vec.push_back(it_sm->second);
    }
  }
  bgzf_close(fp); //closes BGZF file connnection.
}

void ReadGenotype(std::vector<Snp*>& snp_vec, Arguments& args){
  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+args.reference_data_file+"'");
  }
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    int type = (*it_sv)->GetType();
    if(type == 0 || type == 1){ // if unmeasured & exists in ref or measured & exists in ref
      std::string line;
      bgzf_seek(fp, (*it_sv)->GetFpos(), SEEK_SET);
      BgzfGetLine(fp, line);
      std::istringstream buffer(line);
      std::vector<std::string> geno_vec;
      for(int i=0; i<args.num_pops; i++){
        std::string geno_str;
        buffer >> geno_str;
        if(args.pop_flag_vec[i])
          geno_vec.push_back(geno_str); 
      }
      //check flip variable
      if((*it_sv)->GetFlip()){ //if flip is on, flip the genotypes.
        FlipGenotypeVec(geno_vec); //std::string& genotype
      }
      (*it_sv)->SetGenotypeVec(geno_vec);
    } // if (type == 2) don't put the snp in the sliding window. Do nothing. 
  }
  bgzf_close(fp); //closes BGZF file connnection.
}

void FreeGenotype(std::vector<Snp*>& snp_vec){
  //Release memory allocated for genotypes of all SNPs in the window, 
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    std::vector<std::string>& geno_vec = (*it_sv)->GetGenotypeVec();
    std::vector<std::string>().swap(geno_vec); //release memory allocated for genotype
  }
}


void read_ref_desc(Arguments& args){
  
  std::string ref_desc_file = args.reference_pop_desc_file;
  std::ifstream in_ref_desc(ref_desc_file.c_str());
  
  if(!in_ref_desc){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference population description file '"+ref_desc_file+"'");
  }
  
  std::string line;
  std::string pop_abb, sup_pop_abb;
  int pop_num_subj;
  
  std::getline(in_ref_desc, line); //read header of input file.  
  while(std::getline(in_ref_desc, line)){
    std::istringstream buffer(line);
    buffer >> pop_abb >> pop_num_subj >> sup_pop_abb;
    args.ref_pop_vec.push_back(pop_abb);
    args.ref_pop_size_vec.push_back(pop_num_subj);
    args.ref_sup_pop_vec.push_back(sup_pop_abb);
  }//while
  
  args.num_pops=args.ref_pop_vec.size();
  
  in_ref_desc.close();
  Rcpp::Rcout<<std::endl;
}

// used in DIST 
void init_pop_flag_vec(Arguments& args){
  
  int in_pop = std::count(args.ref_pop_vec.begin(), args.ref_pop_vec.end(), args.study_pop);
  int in_sup_pop = std::count(args.ref_sup_pop_vec.begin(), args.ref_sup_pop_vec.end(), args.study_pop);
  
#ifdef init_pop_flag_vec_test
  for(int i=0; i<args.num_pops; i++){
    Rcpp::Rcout<<args.ref_pop_vec[i]<<" "<<args.ref_pop_size_vec[i]<<" "<<args.ref_sup_pop_vec[i]<<std::endl;
  }
  Rcpp::Rcout<<"in_pop: "<<in_pop<<std::endl;
  Rcpp::Rcout<<"in_sup_pop: "<<in_sup_pop<<std::endl;
#endif
  
  std::vector<std::string> pop_vec;
  if(in_pop!=0 && in_sup_pop==0)
    pop_vec = args.ref_pop_vec;
  if(in_pop==0 && in_sup_pop!=0)
    pop_vec = args.ref_sup_pop_vec;
  
  if(in_pop==0 && in_sup_pop==0){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: invalid population name '"+args.study_pop+"'");
  }
  
  int sample_counter = 0;
  for(int i=0; i<args.num_pops; i++){
    if(pop_vec[i]==args.study_pop){
      args.pop_flag_vec.push_back(1);
      sample_counter += args.ref_pop_size_vec[i];
    } else
      args.pop_flag_vec.push_back(0);
  }
  args.num_samples = sample_counter;
}

// used in DISTMIX
void init_pop_flag_wgt_vec(Arguments& args){
  std::string pop;
  for(int i=0; i<args.num_pops; i++){
    pop = args.ref_pop_vec[i];
    if(args.pop_wgt_map.find(pop) != args.pop_wgt_map.end()){ // if pop is found in pop_wgt_map
      args.pop_flag_vec.push_back(1);
      args.pop_wgt_vec.push_back(args.pop_wgt_map[pop]);
    } else {
      args.pop_flag_vec.push_back(0);
      //args.pop_wgt_vec.push_back(0);
    }
  }
}
  

  /*
   void ReadPopulationWeight(Arguments& args){
   Rcpp::Rcout<<"Reading population weight data...";
   Rcpp::Rcout.flush();  
   
   Rcpp::Rcout<<std::endl;
   
   std::string pop_wgt_file = args.pop_wgt_file;
   std::ifstream in_pop_wgt(pop_wgt_file.c_str());
   if(!in_pop_wgt){
   Rcpp::Rcout<<"ERROR: can't open population weight file '"<<pop_wgt_file<<"'"<<std::endl;
   exit(EXIT_FAILURE);
   }
   
   std::map<std::string, double> pop_wgt_map;
   std::string line;
   std::string pop;
   double wgt;
   std::getline(in_pop_wgt, line); //read header of pop wgt file.
   while(std::getline(in_pop_wgt, line)){
   std::istringstream buffer(line);
   buffer >> pop >> wgt;
#ifdef ReadPopulationWeight_Debug    
   Rcpp::Rcout<<"pop: "<<pop<<"  wgt: "<<wgt<<std::endl;
#endif
   //validate population name
   boost::regex pattern ("asw|ceu|chb|chs|clm|fin|gbr|ibs|jpt|lwk|mxl|pur|tsi|yri");
   //boost::regex pattern ("acb|asw|beb|cdx|ceu|chb|chs|clm|esn|fin|gbr|gih|gwd|ibs|itu|jpt|khv|lwk|msl|mxl|pel|pjl|pur|stu|tsi|yri");
   boost::smatch matches;
   std::string pop_low = pop;
   std::transform(pop_low.begin(), pop_low.end(), pop_low.begin(), ::tolower); //make capital 
   if(!boost::regex_match(pop_low, matches, pattern)){
   Rcpp::Rcout << "WARNING: invalid population name '" << pop << "'" << std::endl;
   } else {
   if(wgt < 0 || wgt > 1.5){
   Rcpp::Rcout << "WARNING: invalid population wgt '" << wgt << "'" << std::endl;
   } else {
   pop_wgt_map[pop_low]=wgt; //add pop and wgt to pop_wgt_map
   
   if(pop_low=="asw") args.pop_flag_vec[0]=1;
   else if(pop_low=="ceu") args.pop_flag_vec[1]=1;
   else if(pop_low=="chb") args.pop_flag_vec[2]=1;
   else if(pop_low=="chs") args.pop_flag_vec[3]=1;
   else if(pop_low=="clm") args.pop_flag_vec[4]=1;
   else if(pop_low=="fin") args.pop_flag_vec[5]=1;
   else if(pop_low=="gbr") args.pop_flag_vec[6]=1;
   else if(pop_low=="ibs") args.pop_flag_vec[7]=1;
   else if(pop_low=="jpt") args.pop_flag_vec[8]=1;
   else if(pop_low=="lwk") args.pop_flag_vec[9]=1;
   else if(pop_low=="mxl") args.pop_flag_vec[10]=1;
   else if(pop_low=="pur") args.pop_flag_vec[11]=1;
   else if(pop_low=="tsi") args.pop_flag_vec[12]=1;
   else if(pop_low=="yri") args.pop_flag_vec[13]=1;
   
   //args.pop_wgt_vec(wgt);
   //sum_pop_wgt += wgt; 
   }
   }
   }//while
   in_pop_wgt.close();
   
#ifdef ReadPopulationWeight_Debug    
   //print pop_wgt_map
   Rcpp::Rcout << "pop_wgt_map: "<< std::endl;
   for(std::map< std::string, double >::iterator it = pop_wgt_map.begin(); it != pop_wgt_map.end(); ++it){
   Rcpp::Rcout<<it->first<<" "<<it->second<<std::endl;
   }
   Rcpp::Rcout<<std::endl;  
#endif
   
   //init pop_wgt_vec
   double sum_pop_wgt=0.0;
   for(std::map< std::string, double>::iterator it = pop_wgt_map.begin(); it != pop_wgt_map.end(); ++it){
   args.pop_wgt_vec.push_back(it->second);
   sum_pop_wgt += it->second; 
   }
   
#ifdef VERBOSE
   if(sum_pop_wgt > 1.0){
   Rcpp::Rcout << "WARNING: sum of all population weights exceeds one" << std::endl;
   }
#endif
   if(args.pop_wgt_vec.size() < 1){
   Rcpp::Rcout << "in option 'populationWeight': at least one population needed" << std::endl;
   exit(EXIT_FAILURE); 
   }
   args.sum_pop_wgt = sum_pop_wgt;
   args.num_pops = args.pop_wgt_vec.size();
   
   Rcpp::Rcout<<std::endl;
   }
*/  




/*
void ReadAnnotation(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  Rcpp::Rcout<<"Reading SNP Annotation data...";
  Rcpp::Rcout.flush();

  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;

  std::string annotation_file = args.annotation_file;
  std::ifstream in_annotation(annotation_file.c_str());
  if(!in_annotation){
    Rcpp::Rcout << "ERROR: can't open snp annotation data file '"<<annotation_file<<"'"<<std::endl;
    exit(EXIT_FAILURE);
  }
  std::string line;
  std::string rsid, a1, a2, geneid, categ;
  int chr;
  int categ_num;
  long long int bp;
  double wgt;
 
  std::getline(in_annotation, line); // read header
  while(std::getline(in_annotation, line)){
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> geneid >> categ >> wgt;

    if( (args.chr > 0) && (args.chr != chr) )  // if only one chromosome is specified with -c option and
      continue;                                   // the user specified chr is not equal to chr, skip it.

    MapKey mkey1(chr, bp, a1, a2);
    MapKey mkey2(chr, bp, a2, a1);
    it1 = snp_map.find(mkey1);
    it2 = snp_map.find(mkey2);

    if(categ == "PROTEIN")
      categ_num = 0;
    else if (categ == "TFBS")
      categ_num = 1;
    else if (categ == "WTH_HAIR")
      categ_num = 2;
    else if (categ == "WTH_TARGET")
      categ_num = 3;
    else if (categ == "CIS_EQTL")
      categ_num = 4;
    else if (categ == "TRANS_EQTL")
      categ_num = 5;

    if((it1 != snp_map.end()) && (it2 == snp_map.end())){ // snp exists in snp_map and a1=a1 & a2=a2.

      (it1->second)->SetGeneid(geneid);
      (it1->second)->SetCateg(categ_num, wgt);

    }else if(it1 == snp_map.end() && it2 != snp_map.end() ){ // snp exists in snp_map but a1=a2 & a2=a1.

      // flip genotypes and the sign of z-score of this snp in snp_map based on the annotation info.
      // we don't touch the weight value.

      (it2->second)->SetA1(a1);
      (it2->second)->SetA2(a2);
      (it2->second)->SetAf1study( 1 - (it2->second)->GetAf1study() );
      (it2->second)->SetAf1ref( 1 - (it2->second)->GetAf1ref() );
      (it2->second)->SetZ( (it2->second)->GetZ()*(-1) );
      (it2->second)->SetFlip(true); // true means genotypes of this snp needs to be flipped.

      (it2->second)->SetGeneid(geneid);
      (it2->second)->SetCateg(categ_num, wgt);

      //TODO: when deleting snp using it2, check it2 changes or not. If this changes,
      //then it makes a problem.
      MapKey new_key(chr, bp, a1, a2); //makes a new key for the modified snp object.
      snp_map[new_key] = it2->second;  //makes a map element with the new key. 
      snp_map.erase(it2);              //delete snp with the old key. 

    }//if
  }//while

  in_annotation.close();
  Rcpp::Rcout<<std::endl;
} 

void DeleteUnusedSnp(std::map<MapKey, Snp*, LessThanMapKey>& snp_map){
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end();){
    int type = (it_sm->second)->GetType();
    std::string geneid = (it_sm->second)->GetGeneid();
    if((type == 0 && geneid == ".") || (type == 2)){
      (it_sm->second)->ClearSnp(); // clear categ map in each snp object
      delete it_sm->second;        // delete snp object
      snp_map.erase(it_sm++);      // delete map element	  
    } else {
      ++it_sm;
    }
  }
  //To release allocated memory of deleted Snps immediately, 
  //use swap trick here.
  std::map<MapKey, Snp*, LessThanMapKey> tmp_map;
  tmp_map.insert(snp_map.begin(), snp_map.end());
  tmp_map.swap(snp_map);
}




void MakeGeneStartEndVec(std::vector<StartEnd>& gene_start_end_vec, std::vector<Snp*>& snp_vec){
  std::vector<Snp*>::iterator gene_start = snp_vec.begin();
  std::vector<Snp*>::iterator gene_end = snp_vec.begin();

  std::vector<Snp*>::iterator it_sv;
  std::string current_gene;
  bool first_gene_snp = true;

  for(it_sv=snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    std::string geneid = (*it_sv)->GetGeneid();
    if(first_gene_snp){
      first_gene_snp = false;
      gene_start = gene_end = it_sv;
      current_gene = geneid;
    } else {
      gene_end = it_sv;
      if(current_gene != geneid){
	StartEnd gene_start_end;
	gene_start_end.start_it = gene_start;
	gene_start_end.end_it = gene_end;
	gene_start_end_vec.push_back(gene_start_end);
	--it_sv;
	first_gene_snp = true;
      }//if(current_gene...
    }//if(first_gene...

    if( it_sv == snp_vec.end()-1 ){
      	StartEnd gene_start_end;
	gene_start_end.start_it = gene_start;
	gene_start_end.end_it = snp_vec.end();
	gene_start_end_vec.push_back(gene_start_end);
    }//if(it_sv...
  }//for
}
*/

/*
void Executejepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){

  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  std::vector<Snp*> snp_vec;
  std::vector<Snp*>::iterator it_sv;

  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Error: can't open reference data file '" << args.reference_file << "'" <<std::endl;
    exit(EXIT_FAILURE);
  }
  //opens output data file.
  std::ofstream out_result;
  out_result.open(args.output_file.c_str());
  if(!out_result){
    Rcpp::Rcout<<"Error: can't open output data file '" << args.output_file << "'" << std::endl;
    exit(EXIT_FAILURE);
  }

  ReadAnnotation(snp_map, args);
  Rcpp::Rcout<<"snp_map size : "<<snp_map.size()<<std::endl;
  Rcpp::Rcout<<"snp_vec size (before adding): "<<snp_vec.size()<<std::endl;
  for(it_sm = snp_map.begin(); it_sm != snp_map.end(); ++it_sm){
    if( (it_sm->second)->GetGeneid() != "." && 
	(it_sm->second)->GetType() == 1 &&
	(it_sm->second)->GetInfo() > args.imp_info_cutoff){
      snp_vec.push_back(it_sm->second);
      //(it_sm->second)->PrintSnpInfo();
    }
  }
  Rcpp::Rcout<<"snp_vec size  (after adding): "<<snp_vec.size()<<std::endl;
  
  //Sort snp_vec by geneid
  //TODO: fix LessThanGeneid() to sort snps by (geneid, chr, bp)
  std::sort(snp_vec.begin(), snp_vec.end(), LessThanGeneid());
  
  //Make a table (vector) containing start iterator and end iterator of each gene block.
  std::vector<StartEnd> gene_start_end_vec;
  std::vector<StartEnd>::iterator it_gsev;
  MakeGeneStartEndVec(gene_start_end_vec, snp_vec);
  out_result << "geneid chisq df jepeg_pval num_categ num_protein num_tfbs num_wthhair num_wthtar num_cis num_trans rmvc0 rmvc1 rmvc2 rmvc3 rmvc4 rmvc5 u0 u1 u2 u3 u4 u5"<<std::endl; //print output header.
  for(it_gsev = gene_start_end_vec.begin(); it_gsev != gene_start_end_vec.end(); ++it_gsev){
    std::vector<Snp*> gene_snp_vec(it_gsev->start_it, it_gsev->end_it);
    jepeg jepeg(args, fp);
    jepeg.Runjepeg(gene_snp_vec);
    jepeg.PrintjepegResult(out_result, gene_snp_vec);
  }

  //closes file connections
  bgzf_close(fp); //closes BGZF file connnection.
  out_result.close(); // closes output file connection.
}
*/

/*
void ExecuteQcat(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
   
  std::vector<Snp*> snp_vec;
  std::vector<Snp*>::iterator it_sv;
  std::vector<StartEnd> chr_arm_start_end_vec;
  std::vector<StartEnd>::iterator it_casev;

  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Error: can't open reference data file '" << args.reference_file << "'" <<std::endl;
    exit(EXIT_FAILURE);
  }
  //opens output data file.
  std::ofstream out_result;
  out_result.open(args.output_file.c_str());
  if(!out_result){
    Rcpp::Rcout<<"Error: can't open output data file '" << args.output_file << "'" << std::endl;
    exit(EXIT_FAILURE);
  }

  //Rcpp::Rcout<<"ExecuteQcat 1"<<std::endl;

  MakeSnpVecMix(snp_vec, snp_map, args, fp);

  //Rcpp::Rcout<<"Num of SNPs for Imputation: "<< snp_vec.size() <<std::endl;  

#ifdef main_Debug
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    (*it_sv)->PrintQcatResult();
  }
#endif
  //Rcpp::Rcout<<"ExecuteQcat 2"<<std::endl;
  MakeChrArmStartEndVec(chr_arm_start_end_vec, snp_vec);
  //Rcpp::Rcout<<"ExecuteQcat 3"<<std::endl;

  //out_result << "rsid chr bp a1 a2 waf1 af1 z info pval type"<<std::endl; //print output header.
  out_result << "rsid chr bp a1 a2 waf1 af1 assoc_pval qcat_chisq qcat_pval type"<<std::endl; //print output header.
  for(it_casev = chr_arm_start_end_vec.begin(); it_casev != chr_arm_start_end_vec.end(); ++it_casev){
    std::vector<Snp*> chr_arm_snp_vec(it_casev->start_it, it_casev->end_it);
    Qcat qcat(args, fp);
    qcat.RunQcat(chr_arm_snp_vec);
    qcat.PrintQcatResult(out_result, chr_arm_snp_vec);
    //qcat.PrintQcatResult(chr_arm_snp_vec);
  }
  //Rcpp::Rcout<<"ExecuteQcat 4"<<std::endl; 

  //closes file connections
  bgzf_close(fp); //closes BGZF file connnection.
  out_result.close(); // closes output file connection.

}


void ExecuteQcatLocal(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){

#ifdef ExecuteQcatLocal_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"========================"<<std::endl;
  Rcpp::Rcout<<"ExecuteQcatLocal() Start"<<std::endl;
#endif
   
  std::vector<Snp*> snp_vec;
  std::vector<Snp*>::iterator it_sv;
  
  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Error: can't open reference data file '" << args.reference_file << "'" <<std::endl;
    exit(EXIT_FAILURE);
  }
  //opens output data file.
  std::ofstream out_result;
  out_result.open(args.output_file.c_str());
  if(!out_result){
    Rcpp::Rcout<<"Error: can't open output data file '" << args.output_file << "'" << std::endl;
    exit(EXIT_FAILURE);
  }

#ifdef ExecuteQcatLocal_Debug
  Rcpp::Rcout<<"Before MakeSnpVecMix()"<<std::endl;
#endif
  MakeSnpVecMix(snp_vec, snp_map, args, fp);

#ifdef ExecuteQcatLocal_Debug
  Rcpp::Rcout<<"After MakeSnpVecMix()"<<std::endl;
#endif

#ifdef ExecuteQcatLocal_Debug
  Rcpp::Rcout<<"============="<<std::endl;
  Rcpp::Rcout<<"snp_vec print"<<std::endl;
  for(std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv){
    (*it_sv)->PrintSnpInfo();
  }
#endif
  out_result << "rsid chr bp a1 a2 waf1 af1 assoc_pval qcat_chisq qcat_pval type"<<std::endl; //print output header.

  Qcat qcat(args, fp);
  qcat.RunQcatLocal(snp_vec);
  qcat.PrintQcatLocalResult(out_result, snp_vec);
  //qcat.PrintQcatResult(out_result, snp_vec);
  //qcat.PrintQcatResult(snp_vec);
  
  //closes file connections
  bgzf_close(fp); //closes BGZF file connnection.
  out_result.close(); // closes output file connection.

#ifdef ExecuteQcatLocal_Debug
  Rcpp::Rcout<<std::endl;
  Rcpp::Rcout<<"ExecuteQcatLocal() End"<<std::endl;
  Rcpp::Rcout<<"======================"<<std::endl;
#endif
}

void ExecuteDistjepeg(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){

  std::vector<Snp*> snp_vec;
  std::vector<Snp*>::iterator it_sv;
  std::vector<StartEnd> chr_arm_start_end_vec;
  std::vector<StartEnd>::iterator it_casev;

  //opens reference genotype matrix file (BGZF)
  BGZF* fp = bgzf_open(args.reference_file.c_str(), "r");
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::Rcout<<"Error: can't open reference data file '" << args.reference_file << "'" <<std::endl;
    exit(EXIT_FAILURE);
  }
  //opens output data file.
  std::ofstream out_result;
  out_result.open(args.output_file.c_str());
  if(!out_result){
    Rcpp::Rcout<<"Error: can't open output data file '" << args.output_file << "'" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  ReadAnnotation(snp_map, args);
  //Deletes all unmeasured SNPs without geneid from snp_map.
  Rcpp::Rcout<<"snp_map size (before deleting): "<<snp_map.size()<<std::endl;
  DeleteUnusedSnp(snp_map);
  Rcpp::Rcout<<"snp_map size  (after deleting): "<<snp_map.size()<<std::endl;
  
  MakeSnpVecMix(snp_vec, snp_map, args, fp);
  
  //if(!args.mix_flag)
  //  MakeSnpVec(snp_vec, snp_map, args, fp);
  //else
  //  MakeSnpVecMix(snp_vec, snp_map, args, fp);
  

  MakeChrArmStartEndVec(chr_arm_start_end_vec, snp_vec);
  
  //Runs DIST. 
  //Here DIST imputes only unmeasured snps with geneid.
  for(it_casev = chr_arm_start_end_vec.begin(); it_casev != chr_arm_start_end_vec.end(); ++it_casev){
    std::vector<Snp*> chr_arm_snp_vec(it_casev->start_it, it_casev->end_it);
    Dist dist(args, fp);
    dist.RunDist(chr_arm_snp_vec);
    //dist.PrintDistResult(out_result, chr_arm_snp_vec);
    //dist.PrintDistResult(chr_arm_snp_vec);
  }
  
  //Remove from snp_vec SNPs 1)that do not have geneid, 2)whose type is 2, 3)cannot be imputed (info=-1) 
  //or 4) whose imputation information is less than args->imp_info_cutoff.
  for(it_sv = snp_vec.begin(); it_sv != snp_vec.end();){
    if((*it_sv)->GetGeneid() == "." || (*it_sv)->GetType() == 2 || (*it_sv)->GetInfo() < args.imp_info_cutoff){
      it_sv = snp_vec.erase(it_sv);
    } else {
      ++it_sv;
    }
  }
  //Sort snp_vec by geneid
  //TODO: fix LessThanGeneid() to sort snps by (geneid, chr, bp)
  std::sort(snp_vec.begin(), snp_vec.end(), LessThanGeneid());
  
  //Make a table (vector) containing start iterator and end iterator of each gene block.
  std::vector<StartEnd> gene_start_end_vec;
  std::vector<StartEnd>::iterator it_gsev;
  MakeGeneStartEndVec(gene_start_end_vec, snp_vec);
  out_result << "geneid chisq df jepeg_pval num_categ num_protein num_tfbs num_wthhair num_wthtar num_cis num_trans rmvc0 rmvc1 rmvc2 rmvc3 rmvc4 rmvc5 u0 u1 u2 u3 u4 u5"<<std::endl; //print output header.
  for(it_gsev = gene_start_end_vec.begin(); it_gsev != gene_start_end_vec.end(); ++it_gsev){
    std::vector<Snp*> gene_snp_vec(it_gsev->start_it, it_gsev->end_it);
    Jepeg jepeg(args, fp);
    jepeg.RunJepeg(gene_snp_vec);
    jepeg.PrintJepegResult(out_result, gene_snp_vec);
  } 

  //closes file connections
  bgzf_close(fp); //closes BGZF file connnection.
  out_result.close(); // closes output file connection.
}

void ExecuteDistSnpx(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){

}
*/