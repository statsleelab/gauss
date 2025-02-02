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


// Function: ReadInputZ
// Purpose:
//   This function reads GWAS summary statistics from an input file
//   and stores the SNP information in a map (`snp_map`). Each SNP is represented by a unique key 
//   (a combination of chromosome, base pair position, and alleles) and the SNP object contains details 
//   like the rsid, chromosome, base pair, alleles, and Z-score.
//
// Parameters:
//   - snp_map: A reference to a map where SNP data will be stored. The key is a MapKey object 
//     (which uniquely identifies the SNP based on chromosome, base pair, and alleles), and the 
//     value is a pointer to the corresponding Snp object.
//   - args: A reference to an Arguments object, which contains parameters such as file paths, 
//     chromosome number, start and end base pair positions, and optional filtering criteria.
//   - All: A boolean flag that determines whether to filter the SNPs by chromosome and position. 
//     If `All` is false, only SNPs within the specified chromosome and base pair window are processed.
//.    IF `All` is true, all SNPs in the summary statistics data are processed.
void ReadInputZ(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args, bool All) {
  
  Rcpp::Rcout << "Reading input...";
  Rcpp::Rcout.flush(); // Ensures the message is displayed immediately
  
  // Get the input file path from the args structure
  std::string input_file = args.input_file;
  
  // Open the input file as an input stream
  std::ifstream in_input(input_file.c_str());
  
  // Check if the file was successfully opened, otherwise throw an error and stop
  if (!in_input) {
    Rcpp::stop("ERROR: can't open input file '" + input_file + "'");
  }
  
  std::string line;       // Store each line of the file
  std::string rsid, a1, a2; // SNP ID (rsid), allele 1 (a1), allele 2 (a2)
  int chr;                 // Chromosome number
  long long int bp;        // Base pair position
  double z;                // Z-score from GWAS summary statistics
  double info = 1.0;       // Default imputation information value 
  Snp* snp;                // Pointer to a SNP object
  
  // Read the header line of the input file (this is skipped as it's not used)
  std::getline(in_input, line); // Read the first line (header)
  
  // Loop over each subsequent line in the input file
  while (std::getline(in_input, line)) {
    // Use a string stream to parse the line into SNP variables
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> z; // Extract rsid, chromosome, base pair, alleles, and Z-score
    
    // If not processing all SNPs (All is false), apply filtering based on chromosome and position
    if (!All) {
      // If a specific chromosome is specified and the current SNP's chromosome doesn't match, skip it
      if ((args.chr > 0) && (args.chr != chr)) 
        continue; // Skip this SNP and go to the next one
      
      // If the SNP's base pair position is outside the specified window (start_bp - wing_size to end_bp + wing_size), skip it
      if ((args.start_bp - args.wing_size) > bp || (args.end_bp + args.wing_size) < bp)
        continue; // Skip this SNP and go to the next one
    }
    
    // Create a new SNP object dynamically
    snp = new Snp();
    
    // Set the values of the SNP object using the extracted data
    snp->SetRsid(rsid);  // Set the rsid
    snp->SetChr(chr);     // Set the chromosome number
    snp->SetBp(bp);       // Set the base pair position
    snp->SetA1(a1);       // Set allele 1
    snp->SetA2(a2);       // Set allele 2
    snp->SetZ(z);         // Set the Z-score
    snp->SetInfo(info);   // Set the info (defaulted to 1.0)
    snp->SetType(2);      // Set the SNP type to 2, meaning it is a measured SNP but does not exist in the reference data
    
    // Create a MapKey object that uniquely identifies the SNP based on chr, bp, a1, and a2
    MapKey mkey(chr, bp, a1, a2);
    
    // Insert the SNP object into the map with the MapKey as the key
    snp_map[mkey] = snp;
  } // End of while loop
  
  // Close the input file after processing all lines
  in_input.close();
  
  // Print a newline to the R console to indicate that the process is finished
  Rcpp::Rcout << std::endl;
}


// Function: ReadInputAf
// Purpose:
//   This function reads an input file containing study allele frequency data (`af1study`) for SNPs and stores the SNPs
//   in a map (`snp_map`). The map uses `MapKey` (a combination of chromosome, base pair position, and alleles) 
//   as the key and `Snp` objects as the values.
//
//   The function is used in `afmix()` to populate the `snp_map` with SNP data, including allele frequencies, 
//   from the input file.
//
// Parameters:
//   - snp_map: A map where each SNP is stored using a `MapKey`. The key is based on the chromosome, base pair position, and alleles.
//   - args: An `Arguments` object that holds relevant information such as the input file path.
//
// Function Details:
//   1. Opens the input file that contains SNP data (including allele frequencies).
//   2. Iterates through the file line by line, reading SNP details (rsid, chromosome, base pair, alleles, and allele frequency).
//   3. Creates an `Snp` object for each SNP and stores it in `snp_map` using `MapKey`.
//   4. Closes the input file after reading all the lines.
void ReadInputAf(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  
  Rcpp::Rcout<<"Reading input...";
  Rcpp::Rcout.flush();
  
  // Get the input file path from the arguments object
  std::string input_file = args.input_file;
  
  // Open the input file as an input stream
  std::ifstream in_input(input_file.c_str());
  
  // Check if the input file was opened successfully, if not, print an error message and stop
  if(!in_input){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open input file '"+input_file+"'");
  }
  
  std::string line;
  std::string rsid, a1, a2;  // Variables to store SNP ID and alleles
  int chr;                   // Variable to store the chromosome number
  long long int bp;           // Variable to store the base pair position
  double af1study;            // Variable to store the study allele frequency
  Snp* snp;                   // Pointer to an Snp object
  
  // Read and discard the first line of the file (header)
  std::getline(in_input, line);  // Read header of input file.  
  
  // Loop through the rest of the file, reading each SNP's data
  while(std::getline(in_input, line)){
    
    // Use a string stream to parse the line into individual variables
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1study;
    
    // Create a new SNP object and set its properties (rsid, chr, bp, a1, a2, af1study)
    snp = new Snp();
    snp->SetRsid(rsid);
    snp->SetChr(chr);
    snp->SetBp(bp);
    snp->SetA1(a1);
    snp->SetA2(a2);
    snp->SetAf1Study(af1study);
    
    // Create a MapKey for the SNP (using chr, bp, a1, and a2) and store the SNP in snp_map
    MapKey mkey(chr, bp, a1, a2);
    snp_map[mkey] = snp;  // Add the SNP to the map with the generated key
  } // End of while loop
  
  // Close the input file after reading all the data
  in_input.close();
  Rcpp::Rcout<<std::endl;  // Print a newline to indicate completion
}



// Function: ReadReferenceIndex
// Purpose:
//   This function reads the reference index file (compressed BGZF format) and compares each SNP from the 
//   reference file to the SNPs in `snp_map` (which holds SNPs from the input GWAS summary statistics). 
//   It updates SNPs in `snp_map` with information from the reference file, such as rsid, allele order, and 
//   the file position (`fpos`) of the SNP's genotype string in the reference panel data.
//   If a SNP from the reference index is not found in `snp_map`, it is added as a new SNP entry.
//
// Parameters:
//   - snp_map: A map containing SNPs from the input GWAS file, keyed by MapKey (a combination of chr, bp, a1, a2).
//     The values are pointers to Snp objects containing detailed SNP information.
//   - args: An Arguments object containing parameters such as file paths, chromosome number, base pair positions,
//     and filtering options (e.g., `chr`, `start_bp`, `end_bp`, and `wing_size`).
//
// Function Details:
//   1. Opens the reference index file using the BGZF format.
//   2. Iterates through each line of the reference index file, extracting SNP data such as chromosome, position,
//      alleles, and `fpos` (which stores the file position in the reference panel data for the SNP).
//   3. Filters SNPs by chromosome and base pair range based on the options provided in `args`.
//   4. For each SNP, checks if it exists in `snp_map` (from the input GWAS):
//      - If found, updates the SNP with reference information (rsid, allele info, file position in reference data, etc.).
//      - If not found, adds it to `snp_map` as a new unmeasured SNP with the `fpos` recorded.
//   5. If the SNP has the reverse alleles (a1 and a2 swapped), updates the Z-score, swaps alleles, and modifies the key.
//
// Error Handling:
//   If duplicates are found in the input or reference files, the function stops with an error message.

void ReadReferenceIndex(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args) {
  
  Rcpp::Rcout << "Reading reference index...";
  Rcpp::Rcout.flush();  // Ensure immediate printing
  
  // Declare iterators for finding SNPs in the map
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;
  
  // Open the reference index file using BGZF format (compressed)
  std::string reference_index_file = args.reference_index_file;
  BGZF* fp = bgzf_open(reference_index_file.c_str(), "r"); // Open for reading
  
  // Check if the file was successfully opened, otherwise stop with an error
  if (!fp) {
    Rcpp::stop("ERROR: can't open reference index file '" + reference_index_file + "'");
  }
  
  int last_char;
  std::string line;
  
  // SNP data variables
  std::string rsid, a1, a2;    // SNP ID and alleles
  int chr;                     // Chromosome number
  double af1ref;               // Reference allele frequency (not used in this function)
  long long int bp, fpos;      // Base pair position and `fpos` stores the file position in the reference panel data
  Snp* snp;                    // Pointer to an Snp object
  
  // Loop to read each line from the reference index file
  while (true) {
    
    last_char = BgzfGetLine(fp, line); // Read a line from the file
    if (last_char == -1)  // Check if end of file (EOF)
      break;              // Exit loop if EOF is reached
    
    // Parse the line into SNP variables (rsid, chr, bp, a1, a2, af1ref, fpos)
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    // Filter SNPs by chromosome if specified (args.chr > 0)
    if ((args.chr > 0) && (args.chr != chr)) 
      continue; // Skip SNPs that are not on the specified chromosome
    
    // Filter SNPs by base pair position (start_bp - wing_size to end_bp + wing_size)
    if ((args.start_bp - args.wing_size) > bp || (args.end_bp + args.wing_size) < bp)
      continue; // Skip SNPs outside the desired base pair range
    
    // Create keys for SNP (both possible allele orders: a1/a2 and a2/a1)
    MapKey mkey1(chr, bp, a1, a2);
    MapKey mkey2(chr, bp, a2, a1);
    
    // Find SNP in snp_map using both possible keys
    it1 = snp_map.find(mkey1);
    it2 = snp_map.find(mkey2);
    
    // Case 1: SNP found with key mkey1 (alleles in correct order, i.e., a1=a1 & a2=a2)
    if ((it1 != snp_map.end()) && (it2 == snp_map.end())) {
      // Update existing SNP object with reference data
      (it1->second)->SetRsid(rsid);        // Set the reference rsid
      (it1->second)->SetType(1);           // Type 1: Measured SNP that exists in the reference panel
      (it1->second)->SetFpos(fpos);        // Set file position (fpos) in reference panel data
      
      // Case 2: SNP found with key mkey2 (alleles in reverse order, i.e., a1=a2 & a2=a1)
      // Here, I am updating the alleles in the snp_map to match the reference panel alleles (i.e., a1 and a2),
      // Therefore, I do not need to flip the genotypes in the reference panel. 
    } else if ((it1 == snp_map.end()) && (it2 != snp_map.end())) {
      // Update SNP object with swapped alleles
      (it2->second)->SetRsid(rsid);        // Set the reference rsid
      (it2->second)->SetA1(a1);            // Set allele 1
      (it2->second)->SetA2(a2);            // Set allele 2
      (it2->second)->SetZ((it2->second)->GetZ() * (-1));  // Reverse Z-score
      (it2->second)->SetType(1);           // Type 1: Measured SNP that exists in the reference panel
      (it2->second)->SetFpos(fpos);        // Set file position (fpos) in reference panel data
      
      // Modify the key in the map
      MapKey new_key(chr, bp, a1, a2);     // Create new key with updated allele order
      snp_map[new_key] = it2->second;      // Insert SNP with new key
      snp_map.erase(it2);                  // Remove the old key (reverse alleles)
      
      // Case 3: SNP not found in snp_map (new SNP from reference panel)
    } else if (it1 == snp_map.end() && it2 == snp_map.end()) {
      // Create a new SNP object for the unmeasured SNP in the reference panel
      snp = new Snp();
      snp->SetRsid(rsid);                  // Set the reference rsid
      snp->SetChr(chr);                    // Set the chromosome number
      snp->SetBp(bp);                      // Set the base pair position
      snp->SetA1(a1);                      // Set allele 1
      snp->SetA2(a2);                      // Set allele 2
      snp->SetType(0);                     // Type 0: Unmeasured SNP that exists in the reference panel
      snp->SetFpos(fpos);                  // Set file position (fpos) in reference panel data
      
      // Add the new SNP to the map
      snp_map[mkey1] = snp;
      
      // Case 4: Duplicates found (this should not happen, error case)
    } else {
      // Throw error message if duplicates are found
      Rcpp::stop("ERROR: input file contains duplicates");
    }
  } // End of while loop
  
  // Close the BGZF file after reading
  bgzf_close(fp);
  
  // Print newline to indicate completion
  Rcpp::Rcout << std::endl;
}



// Function: ReadReferenceIndexAll
// Purpose:
//   This function reads the reference index file (compressed in BGZF format) and updates SNPs
//   from the reference file into `snp_map`, which stores SNPs from the input GWAS summary statistics.
//   It focuses on updating existing/measured SNPs in `snp_map` without adding new SNPs that are not present
//   in the input GWAS data.
//
// Difference from ReadReferenceIndex:
//   Unlike `ReadReferenceIndex`, this function does not add new SNPs from the reference panel that
//   are not already present in the input GWAS file. `ReadReferenceIndexAll` is focused on updating 
//   SNPs that already exist in the `snp_map`. This distinction makes it more efficient when you only need 
//   to update existing SNPs without introducing unmeasured ones from the reference panel.
//
// Parameters:
//   - snp_map: A map containing SNPs from the input GWAS file, keyed by MapKey (a combination of chr, bp, a1, a2).
//     The values are pointers to Snp objects that hold detailed SNP information.
//   - args: An Arguments object containing parameters such as file paths and options for filtering.
//
// Function Details:
//   1. Opens the reference index file using the BGZF format.
//   2. Iterates through each line of the reference index file to extract SNP information.
//   3. Checks if each SNP exists in `snp_map`:
//      - If found with matching alleles, it updates SNP details (rsid, fpos).
//      - If found with swapped alleles, it modifies Z-score, alleles, and allele frequencies, and updates the map.
//   4. Throws an error if the same SNP appears more than once (duplicate).
//
// Error Handling:
//   If duplicates are found in the input or reference files, the function stops with an error message.
void ReadReferenceIndexAll(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args) {
  
  Rcpp::Rcout << "Reading reference index...";
  Rcpp::Rcout.flush();  // Ensure immediate printing
  
  // Declare iterators for finding SNPs in the map
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;
  
  // Open the reference index file using BGZF format
  std::string reference_index_file = args.reference_index_file;
  BGZF* fp = bgzf_open(reference_index_file.c_str(), "r"); // Open for reading
  
  // Error handling if the file cannot be opened
  if (!fp) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: can't open reference index file '" + reference_index_file + "'");
  }
  
  int last_char;
  std::string line;
  
  // SNP data variables
  std::string rsid, a1, a2;    // SNP ID and alleles
  int chr;                     // Chromosome number
  double af1ref;               // Reference allele frequency (not used in this function)
  long long int bp, fpos;      // Base pair position and file position (fpos) in reference data
  
  // Loop to read each line from the reference index file
  while (true) {
    
    // Read a line from the reference index file
    last_char = BgzfGetLine(fp, line);
    if (last_char == -1)  // Check if end of file (EOF)
      break;              // Exit loop if EOF is reached
    
    // Parse the line into SNP variables (rsid, chr, bp, a1, a2, af1ref, fpos)
    std::istringstream buffer(line);
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> af1ref >> fpos;
    
    // Create MapKey for both possible allele orders (a1/a2 and a2/a1)
    MapKey mkey1(chr, bp, a1, a2);
    MapKey mkey2(chr, bp, a2, a1);
    
    // Find SNP in snp_map using both possible keys (mkey1 and mkey2)
    it1 = snp_map.find(mkey1);
    it2 = snp_map.find(mkey2);
    
    // Case 1: SNP found with key mkey1 (alleles in correct order)
    if ((it1 != snp_map.end()) && (it2 == snp_map.end())) {
      // Update existing SNP object with reference data
      it1->second->SetRsid(rsid);        // Set the reference rsid
      it1->second->SetType(1);           // Type 1: SNP measured and exists in reference panel
      it1->second->SetFpos(fpos);        // Set file position in reference panel data
      
    // Case 2: SNP found with key mkey2 (alleles in reverse order, i.e., a1=a2 & a2=a1)
    // Here, I am updating the alleles in the snp_map to match the reference panel alleles (i.e., a1 and a2),
    // Therefore, I do not need to flip the genotypes in the reference panel.   
    } else if ((it1 == snp_map.end()) && (it2 != snp_map.end())) {
      // Update SNP object with swapped alleles and reverse Z-score
      it2->second->SetRsid(rsid);        // Set the reference rsid
      it2->second->SetA1(a1);            // Set allele 1
      it2->second->SetA2(a2);            // Set allele 2
      it2->second->SetZ(it2->second->GetZ() * (-1));  // Reverse Z-score
      it2->second->SetType(1);           // Type 1: SNP measured and exists in reference panel
      it2->second->SetAf1Study(1 - it2->second->GetAf1Study());  // Adjust allele frequency
      it2->second->SetFpos(fpos);        // Set file position in reference panel data
      
      // Reinsert the SNP with the corrected key and remove the old one
      MapKey new_key(chr, bp, a1, a2);     // Create new key with updated allele order
      snp_map[new_key] = it2->second;      // Insert SNP with new key
      snp_map.erase(it2);                  // Remove the old key (reverse alleles)
      
      // Case 3: Duplicate SNP found (this should not happen)
    } else if (it1 != snp_map.end() && it2 != snp_map.end()) {
      // Print an error message and stop if duplicates are found
      Rcpp::Rcout << std::endl;
      Rcpp::stop("ERROR: input file contains duplicates");
    }
    
  } // End of while loop
  
  // Close the BGZF file after reading
  bgzf_close(fp);
  
  // Print newline to indicate completion
  Rcpp::Rcout << std::endl;
}



// Function: MakeSnpVec
// Purpose:
//   This function reads the reference genotype data for each SNP in the `snp_map` and calculates the reference 
//   allele frequency (`af1ref`) for each SNP. If the allele frequency meets the specified cutoff criteria, 
//   the SNP is added to the `snp_vec`.
//
//   The function processes genotype data from multiple populations, summing the alleles and subject counts to calculate 
//   the allele frequency. The result is a list of SNPs (`snp_vec`) that satisfy the allele frequency cutoff condition.
//
// Parameters:
//   - snp_vec: A vector to store pointers to `Snp` objects that meet the allele frequency cutoff criteria.
//   - snp_map: A map where each SNP is stored using a `MapKey`. The key is based on the chromosome, base pair position, and alleles.
//   - args: An `Arguments` object that holds information such as the reference data file path, population flags, and allele frequency cutoff.
//
// Function Details:
//   1. Opens the reference genotype data file (BGZF format).
//   2. Iterates through the SNPs in `snp_map` and reads their genotype data from the file.
//   3. Calculates the reference allele frequency (`af1ref`) for each SNP.
//   4. Adds SNPs that meet the allele frequency cutoff criteria to `snp_vec`.
//   5. Closes the reference data file after processing.
// Used in DIST, QCAT, JEPEG
void MakeSnpVec(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  
  // Open the reference genotype data file (BGZF format)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  
  // If the file cannot be opened, print an error message and stop
  if(fp == NULL){
    Rcpp::Rcout<<std::endl;
    Rcpp::stop("ERROR: can't open reference data file '"+args.reference_data_file+"'");
  }
  
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;  // Iterator for snp_map
  double af1ref = 0;  // Variable to store the reference allele frequency
  
  // Loop through each SNP in snp_map
  for(it_sm = snp_map.begin(); it_sm != snp_map.end(); ++it_sm){
    
    // Seek to the SNP's position in the reference file using its stored file position (fpos)
    bgzf_seek(fp, (it_sm->second)->GetFpos(), SEEK_SET);
    
    std::string line;
    BgzfGetLine(fp, line);  // Read the line containing the SNP's genotype data
    
    std::istringstream buffer(line);  // Create a string stream to parse the line
    
    double allele_counter = 0;  // Counter for the number of alleles
    double num_subj = 0;        // Counter for the number of subjects (individuals)
    
    // Loop through each population to read and process genotype data
    for(int k = 0; k < args.num_pops; k++){
      std::string geno_str;
      buffer >> geno_str;  // Read the genotype string for the population
      
      // Process genotype data only for populations that are flagged (pop_flag_vec[k] == 1)
      if(args.pop_flag_vec[k]){
        num_subj += geno_str.length();  // Add the number of subjects (length of the genotype string)
        
        // Sum the alleles for each individual in the population
        for(int i = 0; i < geno_str.length(); i++){
          allele_counter += (double)(geno_str[i] - '0');  // Convert character to number and sum
        }
      }
    }
    
    // Calculate the reference allele frequency (af1ref) by dividing allele count by total alleles (2 * num_subj)
    af1ref = allele_counter / (2 * num_subj);
    
    // Round the allele frequency to 5 decimal places
    af1ref = std::ceil(af1ref * 100000.0) / 100000.0;
    
    // Set the reference allele frequency for the SNP
    (it_sm->second)->SetAf1Ref(af1ref);
    
    // If the allele frequency meets the cutoff criteria, add the SNP to snp_vec
    if( (af1ref > args.af1_cutoff) && (af1ref < (1 - args.af1_cutoff)) ){
      snp_vec.push_back(it_sm->second);
    }      
  }
  
  // Close the BGZF file after processing
  bgzf_close(fp);
}



// Function: MakeSnpVecMix
// Purpose:
//   This function reads allele frequencies from the reference panel data, computes the weighted 
//   sum of reference allele frequencies (af1_mix) for each SNP, and stores the SNPs in a vector if 
//   they meet the specified allele frequency cutoff criteria. 
//
//   The function is used in DISTMIX, QCATMIX, and JEPEGMIX
//
// Parameters:
//   - snp_vec: A vector that will store SNPs that meet the allele frequency cutoff criteria.
//   - snp_map: A map of SNPs, where the key is a MapKey (chr, bp, alleles) and the value is a pointer 
//     to an Snp object. These SNPs come from the input GWAS summary statistics.
//   - args: An Arguments object that contains relevant information like file paths, population weights, 
//     allele frequency cutoff, and population flags.
//
// Function Details:
//   1. Opens the reference genotype matrix file using the BGZF format.
//   2. For each SNP in `snp_map`, retrieves the SNP's file position (fpos) and reads the reference 
//      allele frequencies from the reference panel.
//   3. Computes the weighted sum of reference allele frequencies (af1_mix) using population-specific weights.
//   4. Adds the SNP to `snp_vec` if its allele frequency (af1_mix) is greater than `af1_cutoff` 
//      and less than `1 - af1_cutoff`.
//   5. Closes the BGZF file after processing.
void MakeSnpVecMix(std::vector<Snp*>& snp_vec, std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args) {
  
  // Opens the reference genotype matrix file (BGZF format)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if (fp == NULL) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: can't open reference data file '" + args.reference_data_file + "'");
  }
  
  // Iterator to traverse through the SNP map (snp_map)
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it_sm;
  
  // Loop through each SNP in the snp_map
  for (it_sm = snp_map.begin(); it_sm != snp_map.end(); ++it_sm) {
    
    std::vector<double> af1_vec; // Vector to store allele frequencies for each population
    double af1_mix = 0;          // Weighted allele frequency (af1_mix) initialization
    std::string line;            // String to hold a line from the reference data file
    
    // Use the stored file position (fpos) to seek to the SNP's location in the BGZF file
    bgzf_seek(fp, (it_sm->second)->GetFpos(), SEEK_SET);
    
    // Read the line containing genotype data and allele frequencies
    BgzfGetLine(fp, line);
    
    // Create a stringstream buffer to parse the line
    std::istringstream buffer(line);
    
    // Skipping genotype data (the first part of the line corresponds to genotype strings)
    for (int k = 0; k < args.pop_flag_vec.size(); k++) {
      std::string tmp;
      buffer >> tmp; // Skip the genotype data
    }
    
    // Reading reference allele frequencies from the file for each population
    for (int k = 0; k < args.pop_flag_vec.size(); k++) {
      double af1;
      buffer >> af1; // Read allele frequency
      
      // Only consider allele frequencies from populations that are flagged (pop_flag_vec[k] == 1)
      if (args.pop_flag_vec[k]) {
        af1_vec.push_back(af1); // Store the allele frequency
      }
    }
    
    // Compute the weighted sum of allele frequencies (af1_mix) using population-specific weights (pop_wgt_vec)
    for (int k = 0; k < af1_vec.size(); k++) {
      af1_mix += af1_vec[k] * args.pop_wgt_vec[k]; // Weighted allele frequency
    }
    
    // Check if the weighted allele frequency falls within the specified cutoff range
    if ((af1_mix > args.af1_cutoff) && (af1_mix < (1 - args.af1_cutoff))) {
      // Set the computed af1_mix value in the SNP object
      (it_sm->second)->SetAf1Mix(af1_mix);
      
      // Add the SNP to the snp_vec vector if it passes the cutoff criteria
      snp_vec.push_back(it_sm->second);
    }
  }
  
  // Closes the BGZF file after processing all SNPs
  bgzf_close(fp);
}



// Function: ReadGenotype
// Purpose:
//   This function reads genotype data from the reference panel data file for each SNP in the SNP vector (`snp_vec`).
//   It extracts the genotype strings for the populations specified by the population flags (`pop_flag_vec`)
//   and stores them in the corresponding SNP objects.
//
//   The function is used to load the genotype data into the SNP objects for SNPs that are either unmeasured 
//   but exist in the reference panel (type 0) or are measured and exist in the reference panel (type 1).
//   It skips SNPs that are measured but do not exist in the reference panel (type 2).
//
// Parameters:
//   - snp_vec: A vector of pointers to `Snp` objects.
//   - args: An Arguments object that contains relevant information like file paths, population flags, and 
//     the number of populations.
//
// Function Details:
//   1. Opens the reference panel data file using the BGZF format.
//   2. For each SNP in `snp_vec`, if its type is 0 (unmeasured but exists in reference) or 1 (measured and exists in reference):
//      - Seeks to the SNP's file position (`fpos`) in the reference file and reads the genotype strings.
//      - Extracts genotype data for populations where `pop_flag_vec` is true (i.e., relevant populations).
//      - If the SNP requires flipping, it flips the genotype strings using the `FlipGenotypeVec` function.
//      - Stores the genotype data in the SNP object.
//   3. Closes the BGZF file after processing all SNPs.
void ReadGenotype(std::vector<Snp*>& snp_vec, Arguments& args) {
  
  // Opens the reference genotype matrix file (BGZF format)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if (fp == NULL) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: can't open reference data file '" + args.reference_data_file + "'");
  }
  
  // Loop through each SNP in the SNP vector (snp_vec)
  for (std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv) {
    
    // Get the type of the SNP: 
    // 0 - unmeasured but exists in reference panel, 
    // 1 - measured and exists in reference panel, 
    // 2 - measured but does not exist in reference panel.
    int type = (*it_sv)->GetType();
    
    // Process SNPs of type 0 (unmeasured but exists in reference) and type 1 (measured and exists in reference)
    if (type == 0 || type == 1) {
      
      std::string line;  // String to store a line from the reference data file
      
      // Use the stored file position (fpos) to seek to the SNP's location in the reference file
      bgzf_seek(fp, (*it_sv)->GetFpos(), SEEK_SET);
      
      // Read the line corresponding to the SNP from the reference file
      BgzfGetLine(fp, line);
      
      // Create a stringstream buffer to parse the line
      std::istringstream buffer(line);
      
      std::vector<std::string> geno_vec;  // Vector to store genotype strings for each population
      
      // Loop through each population and extract the genotype strings
      for (int i = 0; i < args.num_pops; i++) {
        std::string geno_str;
        buffer >> geno_str;  // Read the genotype string
        
        // Only add genotype data for populations that are flagged (pop_flag_vec[i] == 1)
        if (args.pop_flag_vec[i]) {
          geno_vec.push_back(geno_str);
        }
      }
      
      // If the SNP has the "flip" flag set to true, flip the genotypes
      // In ReadReferenceIndex or ReadReferenceIndexAll function, I am updating the alleles 
      // in the snp_map to match the reference panel alleles (i.e., a1 and a2),
      // Therefore, I do not need to flip the genotypes in the reference panel. 
      /*
      if ((*it_sv)->GetFlip()) {
        // Flip the genotype strings if needed
        FlipGenotypeVec(geno_vec);  //defined in util.cpp
      }
      */
      
      // Store the genotype data in the SNP object
      (*it_sv)->SetGenotypeVec(geno_vec);
      
    } // If type == 2 (measured but does not exist in reference), do nothing. 
    
  }
  
  // Close the BGZF file after processing all SNPs
  bgzf_close(fp);
}


// Function: ReadGenotypeOne
// Purpose:
//   This function reads the genotype data for a single SNP from the reference genotype matrix file
//   (in BGZF format) and stores the genotype data in the SNP object.
//   
//   It processes SNPs that either:
//   - Type 0: Unmeasured but exist in the reference panel.
//   - Type 1: Measured and exist in the reference panel.
//   
//   SNPs with Type 2 (measured but do not exist in the reference panel) are ignored. 
//
// Parameters:
//   - snp: A pointer to a `Snp` object, which contains the SNP's data (type, file position, etc.).
//   - args: An `Arguments` object that contains relevant information like file paths, population flags, 
//     and the number of populations.
//
// Function Details:
//   1. Opens the reference genotype matrix file in BGZF format.
//   2. For SNPs of Type 0 or Type 1:
//      - Seeks to the SNP's file position (fpos) in the reference file and reads the genotype strings.
//      - Extracts the genotype data for populations where `pop_flag_vec` is true (i.e., relevant populations).
//      - Stores the genotype data in the SNP object.
//   3. Closes the BGZF file after processing.
void ReadGenotypeOne(Snp* snp, Arguments& args) {
  
  // Opens the reference genotype matrix file (BGZF format)
  BGZF* fp = bgzf_open(args.reference_data_file.c_str(), "r");
  if (fp == NULL) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: can't open reference data file '" + args.reference_data_file + "'");
  }
  
  // Get the type of the SNP: 
  // 0 - unmeasured but exists in reference panel, 
  // 1 - measured and exists in reference panel, 
  // 2 - measured but does not exist in reference panel.
  int type = snp->GetType();
  
  // Process SNPs of type 0 (unmeasured but exists in reference) and type 1 (measured and exists in reference)
  if (type == 0 || type == 1) {
    
    std::string line;  // String to store a line from the reference data file
    
    // Use the stored file position (fpos) to seek to the SNP's location in the reference file
    bgzf_seek(fp, snp->GetFpos(), SEEK_SET);
    
    // Read the line corresponding to the SNP from the reference file
    BgzfGetLine(fp, line);
    
    // Create a stringstream buffer to parse the line
    std::istringstream buffer(line);
    
    std::vector<std::string> geno_vec;  // Vector to store genotype strings for each population
    
    // Loop through each population and extract the genotype strings
    for (int i = 0; i < args.num_pops; i++) {
      std::string geno_str;
      buffer >> geno_str;  // Read the genotype string
      
      // Only add genotype data for populations that are flagged (pop_flag_vec[i] == 1)
      if (args.pop_flag_vec[i]) {
        geno_vec.push_back(geno_str);
      }
    }

    // In the ReadReferenceIndex or ReadReferenceIndexAll function, the alleles in snp_map
    // are updated to match the reference panel alleles (i.e., a1 and a2),
    // Therefore, there is no need to flip the genotypes in the reference panel.
    // This code can be safely removed after verification through testing.
    // If the SNP has the "flip" flag set to true, flip the genotypes
    /*
    if (snp->GetFlip()) {
      // Flip the genotype strings if needed
      FlipGenotypeVec(geno_vec);  // Flip the genotype vector in-place
    }
    */
    
    // Store the genotype data in the SNP object
    snp->SetGenotypeVec(geno_vec);
    
  } // If type == 2 (measured but does not exist in reference), do nothing. 
  
  // Close the BGZF file after processing the SNP
  bgzf_close(fp);
}

// Function: FreeGenotype
// Purpose:
//   This function releases memory that was allocated for the genotype vectors (`geno_vec`) 
//   of all SNPs in the provided SNP vector (`snp_vec`). 
//   
//   It ensures that memory used to store genotype data is properly freed, which is particularly 
//   important when working with large datasets to prevent memory leaks and reduce memory usage.
//
// Parameters:
//   - snp_vec: A vector of pointers to `Snp` objects. Each `Snp` object contains a genotype vector 
//     (`geno_vec`), which holds genotype strings for the SNP.
//
// Function Details:
//   1. Iterates through each SNP in the `snp_vec`.
//   2. For each SNP, retrieves its genotype vector (`geno_vec`).
//   3. Swaps the `geno_vec` with an empty vector to release the memory allocated for genotypes.
//      This effectively frees the memory that was used to store genotype strings.
void FreeGenotype(std::vector<Snp*>& snp_vec) {
  
  // Loop through each SNP in the SNP vector (snp_vec)
  for (std::vector<Snp*>::iterator it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv) {
    
    // Retrieve the genotype vector (geno_vec) for the current SNP
    std::vector<std::string>& geno_vec = (*it_sv)->GetGenotypeVec();
    
    // Release memory allocated for genotype data by swapping the vector with an empty one
    // Swapping with an empty vector ensures that all memory is fully released back to the system.
    // This method is better than using clear() function.
    std::vector<std::string>().swap(geno_vec); // Fully deallocates memory used by the genotype vector
  }
}


// Function: FreeGenotypeOne
// Purpose:
//   This function releases memory that was allocated for the genotype vector (`geno_vec`) 
//   of a single SNP object (`snp`). 
//
//   It is used when memory allocated for a specific SNP's genotype data needs to be freed, ensuring 
//   that the memory is deallocated once the genotype data is no longer needed.
//
// Parameters:
//   - snp: A pointer to an `Snp` object. The `Snp` object contains a genotype vector (`geno_vec`), 
//     which holds genotype strings for that SNP.
//
// Function Details:
//   1. Retrieves the genotype vector (`geno_vec`) for the SNP.
//   2. Swaps the `geno_vec` with an empty vector to fully deallocate the memory used for storing genotypes.
void FreeGenotypeOne(Snp* snp) {
  
  // Retrieve the genotype vector (geno_vec) for the SNP
  std::vector<std::string>& geno_vec = snp->GetGenotypeVec();
  
  // Release memory allocated for genotype data by swapping the vector with an empty one
  std::vector<std::string>().swap(geno_vec); // Fully deallocates memory used by the genotype vector
}



// Function: read_ref_desc
// Purpose:
//   This function reads the reference population description file, which contains information about 
//   different populations in the reference panel. It extracts the population abbreviation, number of subjects, 
//   and superpopulation abbreviation from the file and stores them in the respective vectors in the `args` object.
//
//   The function is used to initialize population-related data, including the total number of populations, 
//   and store this information for use in downstream analyses.
//
// Parameters:
//   - args: An `Arguments` object that holds relevant vectors (`ref_pop_vec`, `ref_pop_size_vec`, `ref_sup_pop_vec`)
//     and the path to the reference population description file.
//
// Function Details:
//   1. Opens the reference population description file.
//   2. Reads the file line by line, skipping the header.
//   3. Extracts population data (abbreviation, size, superpopulation) from each line and stores it in the respective vectors in `args`.
//   4. Stores the total number of populations and closes the file.
void read_ref_desc(Arguments& args) {
  
  // Get the file path for the reference population description file
  std::string ref_desc_file = args.reference_pop_desc_file;
  
  // Open the file as an input stream
  std::ifstream in_ref_desc(ref_desc_file.c_str());
  
  // Check if the file was opened successfully, if not, stop the program and show an error
  if (!in_ref_desc) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: can't open reference population description file '" + ref_desc_file + "'");
  }
  
  std::string line;
  std::string pop_abb, sup_pop_abb;  // Population abbreviation and superpopulation abbreviation
  int pop_num_subj;                  // Number of subjects in the population
  
  // Read the first line (header) and discard it
  std::getline(in_ref_desc, line); // Read header of input file.
  
  // Loop through each line in the file (after the header)
  while (std::getline(in_ref_desc, line)) {
    
    // Use a string stream to parse the line
    std::istringstream buffer(line);
    
    // Read the population abbreviation, number of subjects, and superpopulation abbreviation
    buffer >> pop_abb >> pop_num_subj >> sup_pop_abb;
    
    // Store the parsed data into the respective vectors in args
    args.ref_pop_vec.push_back(pop_abb);           // Add population abbreviation
    args.ref_pop_size_vec.push_back(pop_num_subj); // Add population size
    args.ref_sup_pop_vec.push_back(sup_pop_abb);   // Add superpopulation abbreviation
  } // while
  
  // After reading all lines, store the total number of populations
  args.num_pops = args.ref_pop_vec.size();
  
  // Close the file after processing
  in_ref_desc.close();
  Rcpp::Rcout << std::endl;
}





// Function: init_pop_flag_vec
// Purpose:
//   This function initializes the population flag vector (`pop_flag_vec`) and determines the total number of samples 
//   (`num_samples`) for a specific study population (`study_pop`). It checks if the population exists in the reference
//   population list (`ref_pop_vec`) or the superpopulation list (`ref_sup_pop_vec`).
//
//   The function is used in DIST to set up flags that indicate whether each population in the reference panel is part 
//   of the study population, and it calculates the total number of samples.
//
// Parameters:
//   - args: An `Arguments` object that contains the reference population vectors, study population name, 
//     and other related data.
//
// Function Details:
//   1. Counts how many times the `study_pop` is found in the reference population list (`ref_pop_vec`) and 
//      superpopulation list (`ref_sup_pop_vec`).
//   2. Determines which population vector to use for comparison (either `ref_pop_vec` or `ref_sup_pop_vec`).
//   3. Initializes the `pop_flag_vec` by setting a flag (1) for populations that match the study population and 
//      (0) for those that do not.
//   4. Calculates the total number of samples (`num_samples`).
void init_pop_flag_vec(Arguments& args) {
  
  // Count how many times the study population is found in the reference population vector (ref_pop_vec)
  int in_pop = std::count(args.ref_pop_vec.begin(), args.ref_pop_vec.end(), args.study_pop);
  
  // Count how many times the study population is found in the superpopulation vector (ref_sup_pop_vec)
  int in_sup_pop = std::count(args.ref_sup_pop_vec.begin(), args.ref_sup_pop_vec.end(), args.study_pop);
  
#ifdef init_pop_flag_vec_test
  // Debugging: Print each population's abbreviation, size, and superpopulation abbreviation
  for (int i = 0; i < args.num_pops; i++) {
    Rcpp::Rcout << args.ref_pop_vec[i] << " " << args.ref_pop_size_vec[i] << " " << args.ref_sup_pop_vec[i] << std::endl;
  }
  Rcpp::Rcout << "in_pop: " << in_pop << std::endl;
  Rcpp::Rcout << "in_sup_pop: " << in_sup_pop << std::endl;
#endif
  
  std::vector<std::string> pop_vec; // Vector to hold the population comparison vector
  
  // If study population is found in ref_pop_vec and not in ref_sup_pop_vec, use ref_pop_vec
  if (in_pop != 0 && in_sup_pop == 0)
    pop_vec = args.ref_pop_vec;
  
  // If study population is found in ref_sup_pop_vec and not in ref_pop_vec, use ref_sup_pop_vec
  if (in_pop == 0 && in_sup_pop != 0)
    pop_vec = args.ref_sup_pop_vec;
  
  // If study population is not found in either list, show an error and stop
  if (in_pop == 0 && in_sup_pop == 0) {
    Rcpp::Rcout << std::endl;
    Rcpp::stop("ERROR: invalid population name '" + args.study_pop + "'");
  }
  
  int sample_counter = 0; // Initialize a counter for the total number of samples
  
  // Loop through each population and set the population flag vector (pop_flag_vec)
  for (int i = 0; i < args.num_pops; i++) {
    if (pop_vec[i] == args.study_pop) {
      args.pop_flag_vec.push_back(1);  // Set flag to 1 if the population matches the study population
      sample_counter += args.ref_pop_size_vec[i];  // Add the number of samples for the matching population
    } else {
      args.pop_flag_vec.push_back(0);  // Set flag to 0 if the population does not match
    }
  }
  
  // Store the total number of samples
  args.num_samples = sample_counter;
}



// Function: init_pop_flag_wgt_vec
// Purpose:
//   This function iterates through each population in `args.ref_pop_vec` (the list of all population IDs from the reference panel) 
//   and checks whether each population exists in `args.pop_wgt_map` (which holds population weights from the input data).
//
//   For populations found in `args.pop_wgt_map`, it appends a flag (1) to `pop_flag_vec` and the corresponding population 
//   weight to `pop_wgt_vec`. For populations not found in `pop_wgt_map`, it appends a flag (0) to `pop_flag_vec`.
//
//   The length of `pop_flag_vec` is the same as the number of populations in the reference panel (`ref_pop_vec`), 
//   while the length of `pop_wgt_vec` matches the number of populations present in `pop_wgt_map` (i.e., the populations 
//   present in the input data).
//
//   This function is used in DISTMIX, QCATMIX, and JEPEGMIX to initialize population flags and weights for the analysis.
//
// Parameters:
//   - args: An `Arguments` object that contains the reference population list (`ref_pop_vec`), 
//           population weight map (`pop_wgt_map`), and other related data.
//
// Function Details:
//   1. Iterates through all populations in the reference panel (`ref_pop_vec`).
//   2. For each population, it checks if the population exists in the population weight map (`pop_wgt_map`).
//   3. If the population is found in `pop_wgt_map`, a flag value of 1 is appended to `pop_flag_vec` and the corresponding weight is added to `pop_wgt_vec`.
//   4. If the population is not found, a flag value of 0 is appended to `pop_flag_vec`.
void init_pop_flag_wgt_vec(Arguments& args) {
  
  std::string pop;
  
  // Iterate through each population in the reference panel (ref_pop_vec)
  for (int i = 0; i < args.num_pops; i++) {
    pop = args.ref_pop_vec[i];  // Get the population ID from the reference panel
    
    // Check if the population exists in the population weight map (pop_wgt_map)
    if (args.pop_wgt_map.find(pop) != args.pop_wgt_map.end()) {  // If population is found in pop_wgt_map
      
      // Add a flag value of 1 to pop_flag_vec and add the population's weight to pop_wgt_vec
      args.pop_flag_vec.push_back(1);  // Population is part of the study, set flag to 1
      args.pop_wgt_vec.push_back(args.pop_wgt_map[pop]);  // Add the corresponding population weight to pop_wgt_vec
      
    } else {  // If population is not found in pop_wgt_map
      
      // Add a flag value of 0 to pop_flag_vec (population not part of the study)
      args.pop_flag_vec.push_back(0);
      
      // No weight is added to pop_wgt_vec for populations that are not found
      // args.pop_wgt_vec.push_back(0);  // Optional: This line is commented out, as it isn't needed
    }
  }
}


/**
 * @brief Updates each SNP in snp_vec so that the minor allele becomes the reference allele.
 *
 * This function iterates over the SNP vector and for any SNP where the weighted allele frequency
 * (af1_mix) is greater than 0.5 (indicating that the current reference allele is actually the major allele),
 * it performs the following updates:
 *   - Sets the new allele frequency to 1 - af1_mix.
 *   - Multiplies the Z-score by -1.
 *   - Swaps the allele labels (i.e., exchanges a1 and a2).
 *   - Updates the genotype vector: each genotype string is processed character by character,
 *     converting each digit 'd' (assumed to be '0', '1', or '2') to '2 - d'.
 *
 * This ensures that downstream computations (such as LD calculations or imputation)
 * consistently use the minor allele as the reference allele.
 *
 * @param snp_vec A vector of pointers to Snp objects to update.
 */
void UpdateSnpToMinorAllele(std::vector<Snp*>& snp_vec) {
  // Loop through each SNP in the vector.
  for (size_t i = 0; i < snp_vec.size(); ++i) {
    Snp* snp = snp_vec[i];
    double af = snp->GetAf1Mix();
    
    // If the weighted allele frequency is greater than 0.5,
    // then the current reference allele is the major allele.
    if (af > 0.5) {
      // Update the allele frequency so that the minor allele becomes the reference allele.
      snp->SetAf1Mix(1 - af);
      
      // Flip the Z-score.
      snp->SetZ(-snp->GetZ());
      
      // Swap allele labels.
      std::string origA1 = snp->GetA1();
      std::string origA2 = snp->GetA2();
      snp->SetA1(origA2);
      snp->SetA2(origA1);
      
      // Update the genotype vector.
      // We assume that each genotype is stored as a string of digits representing the count of the reference allele.
      // For example, a genotype "0", "1", or "2" in additive coding.
      // To flip the coding (e.g., from additive to its complement), replace each digit 'd' with '2 - d'.
      std::vector<std::string>& genoVec = snp->GetGenotypeVec();
      for (size_t j = 0; j < genoVec.size(); ++j) {
        std::string oldGeno = genoVec[j];
        std::string newGeno;
        newGeno.reserve(oldGeno.size());
        for (size_t k = 0; k < oldGeno.size(); ++k) {
          char c = oldGeno[k];
          // Check if the character is a digit between '0' and '2'.
          if (c >= '0' && c <= '2') {
            int val = c - '0';
            int flipped = 2 - val; // Flip the value.
            newGeno.push_back('0' + flipped);
          } else {
            // If the character is not a digit (e.g., a missing value marker), copy it unchanged.
            newGeno.push_back(c);
          }
        }
        // Replace the old genotype string with the updated (flipped) version.
        genoVec[j] = newGeno;
      }
    }
  }
}


/**
 * @brief Converts an additive-coded genotype vector to recessive-coded.
 *
 * For each genotype string in the vector, each character (assumed to be '0', '1', or '2')
 * is converted such that if the value is '2' (i.e. homozygous minor), it becomes '1'; otherwise, it becomes '0'.
 *
 * @param genoVec The original genotype vector (additive coding).
 * @return std::vector<std::string> A new genotype vector with recessive coding.
 */
std::vector<std::string> ConvertGenotypesToRecessive(const std::vector<std::string>& genoVec) {
  std::vector<std::string> recessiveVec;
  recessiveVec.reserve(genoVec.size());
  for (size_t i = 0; i < genoVec.size(); i++) {
    std::string oldGeno = genoVec[i];
    std::string newGeno;
    newGeno.reserve(oldGeno.size());
    for (size_t j = 0; j < oldGeno.size(); j++) {
      char c = oldGeno[j];
      if (c >= '0' && c <= '2') {
        int val = c - '0';
        // In recessive coding, only genotype '2' becomes 1; others become 0.
        int recVal = (val == 2) ? 1 : 0;
        newGeno.push_back('0' + recVal);
      } else {
        // For any non-digit (e.g., missing data), copy unchanged.
        newGeno.push_back(c);
      }
    }
    recessiveVec.push_back(newGeno);
  }
  return recessiveVec;
}



// Function: ReadAnnotation
// Purpose:
//   This function reads SNP annotation data from a file and updates the SNPs in `snp_map` with 
//   the corresponding gene IDs, categories, and weights. It also handles cases where alleles are 
//   reversed (`a1` and `a2` swapped) and ensures that the `snp_map` is updated accordingly.
//
//   The function supports different SNP categories (e.g., "PROTEIN", "TFBS", etc.) and assigns a numerical 
//   category value (`categ_num`) to each SNP, along with a weight (`wgt`).
//
// Parameters:
//   - snp_map: A map where each SNP is stored using a `MapKey`. The key is based on the chromosome, base pair position, and alleles.
//   - args: An `Arguments` object that contains relevant file paths, including the annotation file path (`annotation_file`).
//
// Function Details:
//   1. Opens the SNP annotation file.
//   2. Iterates through the file line by line, reading SNP annotation details.
//   3. For each SNP in the annotation file, checks if it exists in `snp_map` (both normal and reversed allele order).
//   4. Updates the SNP with gene ID, category, and weight.
//   5. If alleles are reversed, updates the SNP, flips the reference allele frequency, and reverses the Z-score.
//   6. Ensures `snp_map` is updated with the correct allele order after modifying the SNP.
// Used in JEPEG and JEPEGMIX
void ReadAnnotation(std::map<MapKey, Snp*, LessThanMapKey>& snp_map, Arguments& args){
  
  Rcpp::Rcout<<"Reading SNP Annotation data...";
  Rcpp::Rcout.flush();
  
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it1;  // Iterator for snp_map (normal allele order)
  std::map<MapKey, Snp*, LessThanMapKey>::iterator it2;  // Iterator for snp_map (reversed allele order)
  
  // Get the file path for the SNP annotation data
  std::string annotation_file = args.annotation_file;
  
  // Open the annotation file as an input stream
  std::ifstream in_annotation(annotation_file.c_str());
  
  // Check if the annotation file was opened successfully, if not, stop the program and print an error
  if(!in_annotation){
    Rcpp::Rcout << "ERROR: can't open snp annotation data file '"<<annotation_file<<"'"<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::string line;  // Variable to store each line from the file
  std::string rsid, a1, a2, geneid, categ;  // Variables to store SNP details
  int chr;                                  // Chromosome number
  int categ_num;                            // Category number
  long long int bp;                         // Base pair position
  double wgt;                               // Weight of the category
  
  std::getline(in_annotation, line);  // Read and discard the header line
  
  // Loop through the rest of the file, reading SNP annotation data
  while(std::getline(in_annotation, line)){
    
    std::istringstream buffer(line);  // Use a string stream to parse the line
    buffer >> rsid >> chr >> bp >> a1 >> a2 >> geneid >> categ >> wgt;  // Parse SNP details
    
    // Create map keys for both normal and reversed allele order
    MapKey mkey1(chr, bp, a1, a2);  // Normal allele order (a1, a2)
    MapKey mkey2(chr, bp, a2, a1);  // Reversed allele order (a2, a1)
    
    // Find SNP in the map using both possible allele orders
    it1 = snp_map.find(mkey1);  
    it2 = snp_map.find(mkey2);
    
    // Assign numerical category values based on the annotation category
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
    
    // If SNP exists in snp_map with normal allele order
    if((it1 != snp_map.end()) && (it2 == snp_map.end())){
      // Update the SNP with gene ID, category number, and weight
      (it1->second)->SetGeneid(geneid);
      (it1->second)->SetCateg(categ_num, wgt);
      
      // If SNP exists in snp_map but with reversed allele order
    }else if(it1 == snp_map.end() && it2 != snp_map.end()){
      
      // Update the alleles and flip reference allele frequency and Z-score
      (it2->second)->SetA1(a1);
      (it2->second)->SetA2(a2);
      (it2->second)->SetAf1Ref(1 - (it2->second)->GetAf1Ref());  // Flip allele frequency
      (it2->second)->SetZ((it2->second)->GetZ() * (-1));  // Flip the Z-score
      
      // Update the SNP with gene ID, category number, and weight
      (it2->second)->SetGeneid(geneid);
      (it2->second)->SetCateg(categ_num, wgt);
      
      // Create a new map key for the modified SNP and update the snp_map
      MapKey new_key(chr, bp, a1, a2);  // Create a new key with the updated allele order
      snp_map[new_key] = it2->second;   // Insert the updated SNP with the new key
      snp_map.erase(it2);               // Erase the old entry with the reversed alleles
    }
  }
  
  // Close the annotation file after processing all lines
  in_annotation.close();
  Rcpp::Rcout<<std::endl;  // Print a newline to indicate completion
}


// Function: MakeGeneStartEndVec
// Purpose:
//   This function processes a vector of SNPs (`snp_vec`) and groups the SNPs by their associated gene IDs.
//   For each gene, it records the start and end iterators of the SNPs that belong to that gene, and stores
//   this information in `gene_start_end_vec`. Each gene is represented by a `StartEnd` structure that contains
//   the iterators pointing to the first and last SNPs for the gene.
//
//   This is used for downstream analysis where gene-specific SNP ranges need to be identified.
//
// Parameters:
//   - gene_start_end_vec: A vector of `StartEnd` structures that will store the start and end iterators for each gene.
//   - snp_vec: A vector of SNP pointers, where each SNP is associated with a gene via its `geneid`.
//
// Function Details:
//   1. Iterates through the SNPs in `snp_vec`, grouping them by their associated gene IDs (`geneid`).
//   2. Records the start and end iterators for each group of SNPs that belong to the same gene.
//   3. For each gene, a `StartEnd` structure is created, which stores the start and end iterators.
//   4. Handles edge cases such as the first SNP of a gene and the last SNP in the vector.

void MakeGeneStartEndVec(std::vector<StartEnd>& gene_start_end_vec, std::vector<Snp*>& snp_vec) {
  
  // Iterator to store the start and end positions of SNPs for each gene
  std::vector<Snp*>::iterator gene_start = snp_vec.begin();
  std::vector<Snp*>::iterator gene_end = snp_vec.begin();
  
  std::vector<Snp*>::iterator it_sv;  // Iterator to loop through the SNP vector
  std::string current_gene;  // String to store the current gene ID being processed
  bool first_gene_snp = true;  // Flag to indicate the first SNP of a gene
  
  // Loop through all SNPs in the snp_vec
  for(it_sv = snp_vec.begin(); it_sv != snp_vec.end(); ++it_sv) {
    
    // Get the gene ID for the current SNP
    std::string geneid = (*it_sv)->GetGeneid();
    
    // If processing the first SNP for a gene
    if(first_gene_snp) {
      first_gene_snp = false;  // Reset the flag
      gene_start = gene_end = it_sv;  // Set both start and end iterators to the current SNP
      current_gene = geneid;  // Store the current gene ID
      
    } else {  // If it's not the first SNP for the gene
      
      gene_end = it_sv;  // Update the end iterator to the current SNP
      
      // If the gene ID changes, we've reached the end of the current gene's SNPs
      if(current_gene != geneid) {
        
        // Create a new StartEnd structure to store the start and end iterators for the current gene
        StartEnd gene_start_end;
        gene_start_end.start_it = gene_start;
        gene_start_end.end_it = gene_end;
        
        // Add the gene's SNP range to the gene_start_end_vec
        gene_start_end_vec.push_back(gene_start_end);
        
        --it_sv;  // Move the iterator back one position to process the new gene correctly
        first_gene_snp = true;  // Set flag to true to start processing the new gene
      }  // End of if(current_gene != geneid)
      
    }  // End of if(first_gene_snp)
    
    // Handle the last SNP in snp_vec, which may not trigger the "gene change" logic
    if(it_sv == snp_vec.end() - 1) {
      // Create a StartEnd structure for the last gene in the SNP vector
      StartEnd gene_start_end;
      gene_start_end.start_it = gene_start;
      gene_start_end.end_it = gene_end + 1;  // End iterator is one past the last SNP for the gene
      
      // Add the last gene's SNP range to the gene_start_end_vec
      gene_start_end_vec.push_back(gene_start_end);
    }  // End of if(it_sv == snp_vec.end() - 1)
    
  }  // End of for loop
  
}


