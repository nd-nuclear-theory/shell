/****************************************************************
  read_wavefunctions.cpp

  adapted read functions from group_read.cpp by Patrick J. Fasano,  University of Notre Dame

  Shwetha Vittal
  University of Notre Dame

  + 01/07/25 (slv): First version of code to read many body states from mfdn_MBgroupsxxx
  + 01/22/25 (slv): Create ReadCoefficients function to read the amplitudes, of the 
    many body states, from mfdn_smwf001 file 
****************************************************************/

#ifndef PARTITIONUTILS_GROUP_READ_H_
#define PARTITIONUTILS_GROUP_READ_H_

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <execution>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <map>
#include <valarray>

#include "basis/nlj_orbital.h"
#include "mcutils/fortran_io.h"
#include "mcutils/parsing.h"
#include "mcutils/profiling.h"


// namespace
// {

#define SYMMETRIC_PARTITIONS

struct MBGroupsMetadata
{
  int32_t version_number, num_classes, num_particles, reserved3;
  int32_t parity, twoM, max_Nmax, num_diag;
  int32_t num_groupids, num_states, reserved10, num_blocks;
  int32_t reserved12, reserved13, reserved14, reserved15;
};

struct MFDnPartitioning
{
  std::vector<int> proton_partitions, neutron_partitions;
};

struct MFDnSMWFInfo
{
  int version;
  int Z, N, twoM;
  std::string uniqueID;
  basis::OrbitalSpacePN orbital_space;
  int num_proton_states, num_neutron_states;
  MFDnPartitioning partitioning;
  int parity;
  float weight_max;
  std::size_t dimension;
  int num_diag;
  int num_eigenvectors;
  std::vector<float> energy;
  std::vector<float> J;
  std::vector<float> T;
};

MFDnSMWFInfo ReadMFDnSMWFInfo(std::string filename)
{
  MFDnSMWFInfo info;
  std::string line;
  int line_count = 0;
  mcutils::FileExistCheck(
      filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
    );
  auto stream = std::ifstream(filename, std::ios_base::in);

  // line 1: version
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.version;
    mcutils::ParsingCheck(line_stream, line_count, line);
    assert(info.version == 15200);
  }

  // line 2: Z, N, 2Mj
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.Z >> info.N >> info.twoM;
    mcutils::ParsingCheck(line_stream, line_count, line);
    assert((info.Z + info.N) % 2 == info.twoM % 2);
  }

  // line 3: unique identifier
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.uniqueID;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  // line 4+: number of proton, neutron orbitals, orbital listing body
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    int num_orbitals_p, num_orbitals_n;
    line_stream >> num_orbitals_p >> num_orbitals_n;
    mcutils::ParsingCheck(line_stream, line_count, line);

    std::size_t num_orbitals = num_orbitals_p + num_orbitals_n;
    std::string orbital_info_str;
    for (std::size_t orbital_line_count = 0; orbital_line_count < num_orbitals;
         ++orbital_line_count)
    {
      mcutils::GetLine(stream, line, line_count);
      orbital_info_str.append(line);
      orbital_info_str.append("\n");  // need to restore newline to input line
    }
    std::istringstream orbital_info_stream(orbital_info_str);
    auto orbital_list = basis::ParseOrbitalPNStream(
        orbital_info_stream,
        /*standalone=*/false,
        basis::MFDnOrbitalFormat::kVersion15200
      );
    info.orbital_space = basis::OrbitalSpacePN(orbital_list);
  }

  // derive number of proton and neutron states
  info.num_proton_states = 0;
  info.num_neutron_states = 0;
  const auto& proton_subspace = info.orbital_space.GetSubspace(0);
  for (int index = 0; index < proton_subspace.size(); ++index)
  {
    info.num_proton_states += TwiceValue(proton_subspace.GetState(index).j()) + 1;
  }
  const auto& neutron_subspace = info.orbital_space.GetSubspace(0);
  for (int index = 0; index < neutron_subspace.size(); ++index)
  {
    info.num_neutron_states +=
        TwiceValue(neutron_subspace.GetState(index).j()) + 1;
  }


  // partitioning
  {
    int num_proton_partitions, num_neutron_partitions;
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> num_proton_partitions >> num_neutron_partitions;
    mcutils::ParsingCheck(line_stream, line_count, line);

    auto& partitioning = info.partitioning;
    partitioning.proton_partitions.reserve(num_proton_partitions);
    for (int i = 0; i < num_proton_partitions; ++i)
    {
      int partition;
      stream >> partition;
      partitioning.proton_partitions.push_back(partition);
    }
    partitioning.neutron_partitions.reserve(num_neutron_partitions);
    for (int i = 0; i < num_neutron_partitions; ++i)
    {
      int partition;
      stream >> partition;
      partitioning.neutron_partitions.push_back(info.num_proton_states + partition);
    }
  }

  // basis information
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.parity >> info.weight_max >> info.dimension
        >> info.num_diag;
    mcutils::ParsingCheck(line_stream, line_count, line);
  }

  //Number of eigen vectors in the smwf file. Read eigen energies, J and T
  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> info.num_eigenvectors;
    mcutils::ParsingCheck(line_stream, line_count, line);

    for (int i =0; i < info.num_eigenvectors ; i++){
      float Eb, two_j, t, res; 
      int index, nJ;
      stream >> index >> two_j >> nJ >> t >> Eb >> res;
      info.energy.push_back(Eb);
      info.J.push_back(two_j/2);
      info.T.push_back(t);
    }
    //mcutils::ParsingCheck(line_stream, line_count, line);
  }
  return info;
}

std::vector<int> ReadMBGroups(
    std::string filename_pattern,
    const MFDnSMWFInfo& smwf_info,
    std::vector<std::vector<uint16_t> >& groupid_list,
    bool verbose = false
  )
  {
    std::vector<int> numStatesPerFile;
    // convenience mode variable
    static constexpr std::ios_base::openmode mode_argument =
        std::ios_base::in | std::ios_base::binary;

    // read metadata
    {
      // open stream
      std::string filename = fmt::format(filename_pattern, 1);
      mcutils::FileExistCheck(
          filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
        );
      auto stream = std::ifstream(filename, mode_argument);
      mcutils::StreamCheck(
          bool(stream), filename, "Failure opening groups file for input"
        );

      const auto metadata_vector =
          mcutils::ReadFortranRecord<MBGroupsMetadata>(stream);
      const auto metadata = metadata_vector[0];

      assert(metadata.num_particles == smwf_info.Z + smwf_info.N);
      assert(metadata.twoM == smwf_info.twoM);
  
    }

    for (std::size_t i = 1; i <= smwf_info.num_diag; ++i)
    {

      const std::string filename = fmt::format(filename_pattern, i);
      if (verbose)
        fmt::print("reading file {} ({}/{})\n", filename, i, smwf_info.num_diag);
      auto stream = std::ifstream(filename, mode_argument);
      mcutils::StreamCheck(
          bool(stream), filename, "Failure opening groups file for input"
        );

      const auto metadata_vector =
          mcutils::ReadFortranRecord<MBGroupsMetadata>(stream);
      const auto& metadata = metadata_vector[0];
      if (verbose)
        fmt::print("  num_groups: {}\n", metadata.num_groupids);
      numStatesPerFile.push_back(metadata.num_states);
      //fmt::print("  num_states: {}\n", metadata.num_states);
      mcutils::SkipFortranRecord(stream);  // nblksNm
      const auto Mstateptr = mcutils::ReadFortranRecord<int32_t>(stream);
      const auto groupIDs = mcutils::ReadFortranRecord<int16_t>(stream);

      assert(metadata.num_groupids * metadata.num_particles == groupIDs.size());
      assert(Mstateptr.size() == metadata.num_groupids + 1);

      for (std::size_t j = 0; j < metadata.num_groupids; ++j)
      {
        std::vector<uint16_t> groupid(metadata.num_particles, 0);

        for (int p = 0; p < metadata.num_particles; ++p)
          groupid[p] = groupIDs[metadata.num_particles * j + p] - 1;
        assert(Mstateptr[j + 1] - Mstateptr[j] > 0);
        groupid_list.push_back(groupid);
      }
    }
  return numStatesPerFile;
  }

int setLastMj(const uint16_t numParticles, int num_sp_states, 
              std::vector<int> &mj2_sp, int two_mj, 
              std::vector<int> &next_bin, 
              std::vector<uint16_t> &tempState, int flag){
                //int flag = 1;
                int delta_mj = two_mj;
                for (uint16_t i = 0; i < numParticles; i++)
                  { //std::cout << mj2_sp[tempState[i] ] << std::endl; 
                    delta_mj -= mj2_sp[tempState[i] ]; // because mj2 indices begin from 0 in c++ whereas they begin from 1 in Fortran
                    }
                
                if (delta_mj< 0) {
                  //std::cout<< "Flag from setLastMj : " << flag << std::endl;
                  flag = 1;
                  return flag;}
                
                int iLast = tempState[numParticles-1] + delta_mj/2;
                  //std::cout << "iLast = " << iLast << std::endl;
                
                if (iLast< next_bin[tempState[numParticles -1]]){
                  //std::cout << "next_bin[tempState[numParticles -1] ] " << next_bin[tempState[numParticles -1] ] << std::endl;
                  tempState[numParticles-1]= iLast;
                  //fmt::print("New tempState  {:>4d}\n", fmt::join(tempState," "));
                  flag = 0;
                }
                else flag = -1; 
                //std::cout<< "Flag from setLastMj : " << flag << std::endl;
                return flag;
              }

int incrementMj(const uint16_t numParticles, int num_sp_states, 
                  std::vector<int> &mj2_sp, int two_mj, 
                  std::vector<int> &next_bin, 
                  std::vector<uint16_t> &mbGroup,
                  std::vector<uint16_t> &mbState, int flag){
                  
                  for(int i = numParticles-2; i >= 0; i--){ // index i goes from 4 -> 0 when there are 6 particles
                    if (mbState[i] < (next_bin[mbGroup[i]] -1)){
                      mbState[i] += 1;
                      for(int j = i+1; j < numParticles; j++){
                        //std::cout<< "mbGroup[j] " << mbGroup[j] << std::endl;
                        //std::cout<< "mbState[j-1] " << mbState[j-1] << " i "<< i << " j "<< j <<std::endl;

                        if ((mbState[j-1] + 1) >= (next_bin[mbGroup[j]])){
                        
                          //std::cout<< "Flag is 2 "<<std::endl;
                          //fmt::print(" mbGroup {:>4d}\n", fmt::join(mbGroup," "));
                          //fmt::print(" mbState {:>4d}\n", fmt::join(mbState," "));
                          flag = 2;
                          break; 
                        }
                        else{
                          mbState[j] = mbState[j-1] +1;
                          mbState[j] = std::max(mbState[j], mbGroup[j]);
                        }
                      }
                      if (flag==2){
                        flag = 0;
                        //std::cout<< "flag is 2 here----------------------------" << " i " << i << std::endl;
                        continue;
                      }
                      flag = setLastMj(numParticles, num_sp_states, mj2_sp, two_mj, next_bin, mbState, flag);

                      if (flag ==1) continue;
                        else {
                          //std::cout<< "Flag from incrementMj : " << flag << std::endl;
                          return flag;}
                    }
                  }
                  flag = 1;
                  //std::cout<< "Flag from incrementMj : " << flag << std::endl;
                  return flag;
                }

int MjStatesGen(const uint16_t numParticles, int num_sp_states, 
                  std::vector<int> &mj2_sp, int two_mj, 
                  std::vector<int> &next_bin, 
                  std::vector<uint16_t> &tempVar, int num_states, 
                  std::vector<std::vector<uint16_t> > &mb_state_list,
                  int currentState){
  /****************************************************************
    Functions just like the subroutine mfdn_transitions/src/module_MjStates/MjStatesGen

    numParticles : Total number of particles
    num_sp_states : number of single particle states (it is assumed that there are same number of
                    sp states in both species)
    mj2_sp : contains a list of 2*mj values for each single particle state
    two_mj : the selected 2*mJ
    next_bin : contains the index of the next partition a list of size 1 less than sizeof(mj2_sp)
    tempVar : single GroupID from a list of GroupIDs, the state that marks the beginning of the 
              group(the set of many body states in the partition). It is called ID but it is a 
              vector (of size the numParticles) with sp state indices.
    num_states : Total number of states 
    mb_state_list : List to be updated with the states matching the criteria of having same 2*mj
    currentState : count of number of states matching the criteria of having same 2*mj

    *****************************************************************/

          std::vector<uint16_t> mbstate = tempVar;
          std::vector<uint16_t> mbgroup = mbstate;
          //std::cout<< "MjStatesGen is running .. " << std::endl;
          int flag = 0; 
          flag = setLastMj(numParticles, num_sp_states, mj2_sp, two_mj, next_bin, mbstate, flag);

          if (flag==0){
            //std::cout<< "Adding a state to mb_state_list---------------------------" << currentState <<std::endl;

            // Sanity check
            int total2mj =0;
            for(int i =0; i< numParticles; i++){
              total2mj +=mj2_sp[mbstate[i]];
            }
            if (total2mj != two_mj) fmt::print(" Wrong MJ {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            
            //fmt::print(" mbState {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            mb_state_list[currentState] = mbstate;
            currentState += 1;
            }
          flag = 0;
          while(flag==0){
            flag = -1;
            while(flag == -1){
              //std::cout << "flag is -1" << std::endl;
              flag = incrementMj(numParticles, num_sp_states, mj2_sp, two_mj, next_bin, mbgroup, mbstate, flag); // flag = 1 is the end of loop condition
              //std::cout << "flag is " << flag << std::endl;
            }
            if(flag ==0){
              //std::cout<< "Adding a state to mb_state_list---------------------------" << currentState <<std::endl;
              //std::cout << "flag is 0" << std::endl;
              // Sanity check
            int total2mj =0;
            for(int i =0; i< numParticles; i++){
              total2mj +=mj2_sp[mbstate[i]];
            }
            if (total2mj != two_mj) fmt::print(" Wrong MJ {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            
            //fmt::print(" mbState {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            mb_state_list[currentState] = mbstate;
            currentState += 1;
            }
          }
          return currentState;
        }

std::vector<double> ReadCoefficients(std::string filename_pattern,
                                    const MFDnSMWFInfo& smwf_info,
                                    int state, 
                                    std::vector<int>& numStatesPerFile,
                                    bool verbose = true){
  /****************************************************************
  Reads the coefficients of a wavefunction
  
  filename_pattern : "mfdn_smwf{:03d}"
  smwf_info        : consists of wavefunction information from mfdn_smwf.info
  state            : 0 for ground state, 1 for next eigen state and so on.
  num_statesPerFile: list of number of many body basis states in each file
  
  *****************************************************************/
  std::vector<double> coeffs;

  for (std::size_t i = 1; i <= smwf_info.num_diag; ++i)
  {
    const std::string filename = fmt::format(filename_pattern, i);
    if (verbose)
      fmt::print("reading file {} ({}/{})\n", filename, i, smwf_info.num_diag);
      mcutils::FileExistCheck(
        filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
      );
    
    auto stream = std::ifstream(filename, std::ios_base::binary);
    float buffer;
    int count =0;
      
    // quite unlikely but here is a sanity check to ensure number of bytes in the file 
    // matches the expected number 
    int file_size = stream.tellg();
    stream.seekg(0, std::ios_base::end);
    file_size = int(stream.tellg()) - file_size;
    fmt::print("file_size {:d}   numStates  {:d}    numBytes  {:d}\n", file_size, numStatesPerFile[i-1], (numStatesPerFile[i-1] + 2) * sizeof(float));
    if (file_size % ((numStatesPerFile[i-1] + 2) * sizeof(float)) != 0){ // |1 byte|wf coeffs|1 byte|
      fmt::print("Corrupted file : Unexpected size \n ");
      std::exit(1); // Perhaps a different kind of error must be thrown
    }

    // TO DO (slv): Need to document this 
    int offset = sizeof(float) * (2* state +1);
    stream.seekg(std::ios_base::beg + numStatesPerFile[i-1] * sizeof(float) * state + offset); // set position back to beginning of the state in the stream
    
    while(stream.read(reinterpret_cast<char*>(&buffer), sizeof(float))){
      if(count < numStatesPerFile[i-1]){
            coeffs.push_back(buffer);
          }
      else break;
      count++;
    }
    fmt::print("Count after reading {:d}th file   : {:d} \n", i, count );
  }  
  return coeffs;
}

// Discontinuing work on this idea.
int generateSPStates(//const MFDnSMWFInfo& smwf_info,
                    int numShells,
                    std::vector<std::vector<int16_t> > &spstates_list){
  /*
  Generates single particle states in the BIGSTICK order 
  T--> wt -->2mj --> 2j --> n,l 

  */
  int16_t lmax = numShells - 1; //(int)smwf_info.weight_max;
  int16_t maxWt = numShells - 1; // for clarity
  int16_t twoJMax = 2 * lmax + 1;
  int16_t twoMjMax = twoJMax;
  std::vector<int16_t> twoMjList;
  for(int16_t i = 1 ; i <= twoJMax ; i+=2)
    twoMjList.push_back(i);
  std::cout << "Num elements in twoMjList : " <<twoMjList.size() << std::endl;
  
  int16_t index = 1;
  for(int16_t t = 1 ; t >= -1; t-=2){
    for(int i = -1; i<= 1; i+=2)
    for(int16_t w = 0; w <= maxWt; w++){
      for(int16_t mj = 0; mj < twoMjList.size(); mj++){
                     
          for(int16_t n = 0; n <=w/2; n++){
              for(int16_t l = w; l>=0; l--){        
                if (2*n + l == w) {//smwf_info.weight_max){
                  for(int16_t j = twoJMax; j > 0; j-=2){
                    if(std::abs(twoMjList[mj])> j) continue;   
                  if (l>0)
                    if(j< (2*l - 1)) continue;  

                  if(j> (2*l + 1)) continue;
                  spstates_list.push_back({index, n, l, j, t, (int16_t)(2*n + l), (int16_t)(i * twoMjList[mj])});
                  index++;
                }
              }
            }
          }
      }
    }
  }
  
  return spstates_list.size();
}


int readTrwfn(std::string filename,
                std::vector<std::vector<int16_t> > &sp_state_list_temp,
                std::vector<std::vector<uint16_t> > &mb_state_list_temp){
  /* Reads list of SP states  from a model trwfn file
    
    filename            : Trwfn filename including the path
    sp_state_list_temp  : List of single particle states in trwfn file with columns {n, l, twice_j, twice_mj, species_code}
    mb_state_list_temp  : List of many body states in trwfn file with proton and neutron SP indices in ascending order

  */

  fmt::print(" Reading trwfn file .. \n");
  int count =0;

  std::string line;
  int line_count = 0;
  mcutils::FileExistCheck(
      filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
    );
  auto stream = std::ifstream(filename, std::ios_base::in);
  int Z, N;
  for(int i =0; i<5; i++)
  {
    mcutils::GetLine(stream, line, line_count);
    if(i ==0){
      std::istringstream line_stream(line);
      line_stream >> Z;
    }
    if(i ==1){
      std::istringstream line_stream(line);
      line_stream >> N;
    }
  }
  //std::cout << "First 5 lines "<<line_count<< std::endl;

  std::size_t numSpstates;
  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> numSpstates;
  }
  //std::cout <<"read numSPstates " <<line_count<< std::endl;
  for(int i =0; i<4; i++)
  {
    mcutils::GetLine(stream, line, line_count);
    }
    //std::cout <<"Ignored 4 lines " <<line_count<< std::endl;
  std::size_t num_eigenvectors;
  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> num_eigenvectors;
  }
  //std::cout <<"read numEigenvectors " << num_eigenvectors << "  " << line_count<< std::endl;
  for(int i = 0; i < num_eigenvectors; i++)
  {
    mcutils::GetLine(stream, line, line_count);
  }
  //std::cout <<"Read eigenvectors " <<line_count<< std::endl;

  {
    std::string sps_info_str;
    for (std::size_t sps_line_count = 0; sps_line_count < numSpstates;
         ++sps_line_count)
    {
      mcutils::GetLine(stream, line, line_count);
      //std::cout<< line << std::endl;
      sps_info_str.append(line);
      sps_info_str.append("\n");  // need to restore newline to input line
    }
    std::cout <<"read sps" <<line_count<< std::endl;
    std::istringstream spstates_str(sps_info_str);
    int sp_line_count = 0;
    
    while (mcutils::GetLine(spstates_str,line,sp_line_count))
      {
        // set up for parsing
        std::istringstream line_stream(line);

        int index;
        int16_t n, l, twice_j, twice_mj;
            
        int16_t species_code;
        line_stream >> index >> n >> l >> twice_j >> twice_mj >>species_code;
    
        mcutils::ParsingCheck(line_stream, sp_line_count, line);
        std::vector<int16_t> state({n, l, twice_j, twice_mj, species_code});
        sp_state_list_temp.push_back(state);
        //spstate_list[count][1]= index;
        count++;
      }
  }
  {
    while(mcutils::GetLine(stream, line, line_count)){
      // set up for parsing
      std::istringstream line_stream(line);
      std::vector<uint16_t> spIndices(Z+N, 0);
      for(int i = 0; i< Z+N; i++ ){
        line_stream>>spIndices[i];
      }
      mb_state_list_temp.push_back(spIndices);
      mcutils::GetLine(stream, line, line_count); // Ignore next line
    }
  }
  return count;
}

void findAndReplaceSPIndices(std::vector<std::vector<int16_t> > &sp_state_list, 
                            std::vector<std::vector<int16_t> > &sp_state_list_temp,
                            std::vector<std::vector<uint16_t> > &mb_state_list,
                            int num_sp_states){
  /*
  Looks up indices of states in sp_states_list from sp_state_list_temp (temp for template)
  Replace sp states in mb_state_list with indices found by lookup.
  
  sp_state_list       : List of single particle states in the order that MFDn operates.
                        (n, l, j ,mj, t) 5 columns in the same order as in trwfn
  sp_state_list_temp  : List of single particle states in the order of trwfn (Look at docstring for readTrwfn)
  mb_state_list       : List of many body states from MFDn
  num_sp_states       : Number of single particle states

  */
  fmt::print(" Replacing SP indices .. \n ");
  std::vector<int> index_list(num_sp_states, 0);
  for(int i =0; i < num_sp_states; i++){
    auto it =  std::find(sp_state_list_temp.begin(), sp_state_list_temp.end(), sp_state_list[i]);
    if(it != sp_state_list_temp.end()){
      index_list[i] = it - sp_state_list_temp.begin() + 1; // finding index of the element
    } 
  }
  /*
  // Test
  for(int i =0; i< num_sp_states; i++)
    fmt::print("{:>4d}  \n", index_list[i]);
  */
  for(int i = 0 ; i<mb_state_list.size(); i++){
    for(int j = 0; j<(mb_state_list[i]).size(); j++){
      mb_state_list[i][j] = index_list[mb_state_list[i][j] ];
    }
  }
}

void sortMBBasisStates(std::vector<std::vector<uint16_t> > &mb_state_list,
                      std::map<std::vector<uint16_t> , std::vector<double> > &mb_states_B,
                      std::vector<std::vector<double> > &coefficients_list,
                      std::vector<int16_t> &phaseFactor,// Not needed anymore
                      uint16_t Z, uint16_t N,
                      int num_eigenvectors){
  /*
  Sorts the single particle state(SPS) indices in the mbstates 
  
  mb_state_list     : List of many body states with SPS indices of trwfn
  mb_states_B       : Map with sorted SPS indices as *key* and coefficients of 
                      the corresponding basis state multiplied by the phase factor as *value*
  coefficients_list : Array of coefficients with row -> eigenvector and column -> basis state
  phaseFactor       : phase factor generated to preserve antisymmetry after swapping of states for sorting
  Z                 : Number of protons
  N                 : Number of neutrons
  num_eigenvectors  : Number of eigen vectors in mfdn_smwf*** file

  */

  fmt::print("Sorting SP Indices .. \n ");
  for(std::vector<std::vector<uint16_t> >::iterator it= mb_state_list.begin(); it !=mb_state_list.end(); it++ ) {
    std::vector<uint16_t> mbstate = *it; 
    std::vector<double> Coeffs(num_eigenvectors, 0);
    int phase = 1;
    for(int i =0; i < Z; i++){  
      for(int j = i+1 ; j < Z; j++){
        if( mbstate[i]>mbstate[j]){
          phase *= -1;
          std::swap(mbstate[i], mbstate[j]);
        }
      }
    }
    for(int i =Z; i < Z+N; i++){  
      for(int j = i+1 ; j < Z+N; j++){
        if( mbstate[i]>mbstate[j]){
          phase *= -1;
          std::swap(mbstate[i], mbstate[j]);
        }
      }
    }
    
    for(int i =0; i< num_eigenvectors; i++){
      Coeffs[i] = (coefficients_list[i])[it - mb_state_list.begin()] * phase;
    }
    
    phaseFactor.push_back(phase);
    mb_states_B[mbstate] = Coeffs;
  }
}

int main(int argc, char* argv[])
{
  // header
  std::cout << std::endl;
  std::cout << "read MFDn wavefunctions " << std::endl;
  std::cout << std::endl;

  // usage message
  if (argc-1 < 4)
    {
      std::cout << "Syntax: read_wavefunctions state runmode template_filename output_filename" << std::endl;
      std::exit(EXIT_SUCCESS);
    }
  
  int state; // 0 for lowest eigen wavefunction aka ground state in mfdn_smwf001
  std::istringstream parameter_1(argv[1]);
  parameter_1 >> state;
  if (!parameter_1)
    {
      std::cerr << "Expecting numeric value for Nmax argument" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  
  int runmode;
  std::istringstream parameter_2(argv[2]); // If runmode is 0 -> generate MB states list for low Nmax ;
  // 1 -> generate MBO file
  // 2 -> generate trwfn (..WIP)
  parameter_2 >> runmode;

  // trwfn filename
  std::string template_filename = argv[3];

  // output filename
  std::string out_filename = argv[4];
  
  MBGroupsMetadata metadata{};
  const auto smwf_info = ReadMFDnSMWFInfo("mfdn_smwf.info");
  const uint16_t N = smwf_info.N;
  const uint16_t Z = smwf_info.Z;
  const uint16_t numParticles = N + Z;

  int num_sp_states = smwf_info.num_proton_states + smwf_info.num_neutron_states;
  int two_mj = smwf_info.twoM;
  std::size_t num_states = smwf_info.dimension;
    
  fmt::print(
      "partitions_p: {:>4d}\n",
      fmt::join(smwf_info.partitioning.proton_partitions, " ")
    );
  fmt::print(
      "partitions_n: {:>4d}\n",
      fmt::join(smwf_info.partitioning.neutron_partitions, " ")
    );
  
  fmt::print("dimension: {:d}\n", num_states);
  fflush(stdout);

  std::vector<std::vector<uint16_t> > groupid_list; 
  
  //ReadMBGroups("mfdn_MBgroups{:03d}", smwf_info, groupid_list, false);
  std::vector<int> numStatesPerFile = ReadMBGroups("mfdn_MBgroups{:03d}", smwf_info, groupid_list, false);

  fmt::print("number of groups: {:d}\n", groupid_list.size());
  fflush(stdout);
  
  std::vector<std::vector<uint16_t> > mb_state_list(num_states, std::vector<uint16_t>( numParticles,0)); 
  std::vector<std::vector<int16_t> > sp_state_list; 
  // Create the mj2_sp vector that contains the 2M values of all the single particle states
  std::vector<int> mj2_sp;
  // Create the n_sp, l_sp and twoJ_sp vectors which needs to be printed for runmode=2 i.e., smwf --> trwfn
  std::vector<int> n_sp;
  std::vector<int> l_sp;
  std::vector<int> twoJ_sp;
  const auto& proton_subspace = smwf_info.orbital_space.GetSubspace(0);
  for (int index = 0; index < proton_subspace.size(); ++index)
  {
    auto orb_2j = TwiceValue(proton_subspace.GetState(index).j());
    auto orb_n = proton_subspace.GetState(index).n();
    auto orb_l = proton_subspace.GetState(index).l();
    for(int mj = -1 * orb_2j ; mj <= orb_2j; mj+=2){
      mj2_sp.push_back(mj);
      n_sp.push_back(orb_n);
      l_sp.push_back(orb_l);
      twoJ_sp.push_back(orb_2j);
      std::vector<int16_t> state({(int16_t)orb_n, (int16_t)orb_l, (int16_t)orb_2j, (int16_t)mj, 1});
      // (n, l, j ,mj, t) 5 columns in the same order
      sp_state_list.push_back(state);
    }
  }
  int num_proton_sp_states = mj2_sp.size();
  const auto& neutron_subspace = smwf_info.orbital_space.GetSubspace(0);

  for (int index = 0; index < neutron_subspace.size(); ++index)
  {
    auto orb_2j = TwiceValue(neutron_subspace.GetState(index).j());
    auto orb_n = neutron_subspace.GetState(index).n();
    auto orb_l = neutron_subspace.GetState(index).l();    
    for(int mj = -1 * orb_2j ; mj <= orb_2j; mj+=2){
      mj2_sp.push_back(mj);
      n_sp.push_back(orb_n);
      l_sp.push_back(orb_l);
      twoJ_sp.push_back(orb_2j);
      std::vector<int16_t> state({(int16_t)orb_n, (int16_t)orb_l, (int16_t)orb_2j, (int16_t)mj, -1});
       // (n, l, j ,mj, t) 5 columns in the same order
      sp_state_list.push_back(state);
    }
  }

  // Create the next_bin vector that contains the next partition bin corresponding to each single particle state
  std::vector<int> next_bin;
  std::vector<int> partition;
  partition.reserve(smwf_info.partitioning.proton_partitions.size() + smwf_info.partitioning.neutron_partitions.size());
  partition.insert(partition.end(), smwf_info.partitioning.proton_partitions.begin(), smwf_info.partitioning.proton_partitions.end());
  partition.insert(partition.end(), smwf_info.partitioning.neutron_partitions.begin(), smwf_info.partitioning.neutron_partitions.end());
  partition.push_back(mj2_sp.size()+1);
  
  int next = 1;
  for(int i = 0; i< mj2_sp.size(); i++)
  {
    if((i+1)< partition[next]){
    next_bin.push_back(partition[next]-1);
    }
    else{
      i--;
      next++;}
  
  }
  // Check for proper creation of next_bin vector
  //int count = 0;
  //for(std::vector<int>::iterator it = next_bin.begin(); it !=next_bin.end(); it++)
    //{std::cout<<count++ <<"  "<< *it << std::endl;
    //}

  int currentnumstates = 0;
  for(std::vector<std::vector<uint16_t> >::iterator it = groupid_list.begin(); it != groupid_list.end(); it++ )
  {
    std::vector<uint16_t> groupID = *it;
    // Imitating Fortran subroutine MjStatesGen
    //tempVar = {0, 1, 4, 40, 44, 54}; // test case
    currentnumstates = MjStatesGen(numParticles, num_sp_states, mj2_sp, two_mj, next_bin, 
                       groupID, num_states, mb_state_list, currentnumstates); 

  }  
  

  // write list of MB states to a text file
  if (runmode == 0)  {
    auto stream1 = std::ofstream(out_filename, std::ios_base::out);
    
    stream1 << fmt::format("  {:>4d}   {:>4d}  \n", numParticles, num_states);
    for (int i = 0; i< num_states; i++){
      stream1 << fmt::format("  {:>4d}   \n ", fmt::join(mb_state_list[i],"  "));
    }
    std::cout<< "Writing list of MB states to output file .. " << num_states << " states" << std::endl;
  }


  if (runmode == 1){
    std::vector<double> coefficients = ReadCoefficients("mfdn_smwf{:03d}", smwf_info, state, numStatesPerFile); 
    // for(std::vector<double>::iterator it = coefficients.begin(); it !=coefficients.end(); it++)
    //   std::cout<< *it<<std::endl;

    fmt::print("number of states: {:d}\n", coefficients.size());

    auto stream2 = std::ofstream(out_filename, std::ios_base::binary);
    std::vector<uint16_t> buffer(numParticles , 0);
    
    mcutils::WriteBinary(stream2, &numParticles, 1);
    mcutils::WriteBinary(stream2, &num_states, 1);

    // Here I need to implement the change of indices
    for(int i = 0; i< num_states; i++){
      for(int j =0; j<numParticles; j++){
        
        if(j < Z)
          buffer[j+N] = (mb_state_list[i])[j] + smwf_info.num_proton_states; // change the indices of proton states and move them to positions after the neutron indices
        else
          buffer[j-Z] = (mb_state_list[i])[j] - smwf_info.num_neutron_states; // change the indices of neutron states and move them to positions before the proton indices
        
      }
  /*
      0 1 2 40 41 47               0 1 2 40 41 47
        \_  \__                         __/   __/
          \___ \___                  __/ ____/
              \    \                /   /
        0 1 7 40 41 42             0 1 7 40 41 42
  */
      mcutils::WriteBinary(stream2, buffer.data(), numParticles);
      mcutils::WriteBinary(stream2, &(coefficients[i]), 1);
      //stream2.write(reinterpret_cast<const char*>(&(coefficients[i])),sizeof(double));

    }
  }

  if(runmode==2){

    /*
    //Generate SP states in Calvin's method
    int numShells = 2;
    std::vector<std::vector<int16_t> > spstates_list;
    //int numSPStates = generateSPStates(smwf_info, numShells,spstates_list);
    int numSPStates = generateSPStates(numShells,spstates_list);
    
    for(int i =0; i<numSPStates; i++)
      fmt::print(" {:>4d} \n",fmt::join(spstates_list[i],"  "));
    std::exit(0);
    */

    //Reading a model trwfn file
    std::vector<std::vector<int16_t> > sp_state_list_temp; // temp for template 
    std::vector<std::vector<uint16_t> > mb_state_list_temp; // temp for template
    int numSPStates = readTrwfn(template_filename, sp_state_list_temp, mb_state_list_temp);

    /*
    //Test
    for(std::map<std::vector<int16_t> , int>::iterator it = sp_state_list_temp.begin(); it!=sp_state_list_temp.end(); it++){
      fmt::print(" {:>4d}  ",fmt::join(it->first,"  "));
      std::cout << it->second << std::endl;}
    
    for(std::vector<std::vector<uint16_t> >::iterator it = mb_state_list_temp.begin(); it != mb_state_list_temp.end(); it++ )
      fmt::print(" {:>4d}  \n",fmt::join(*it,"  "));
    std::exit(0);
    */
    
    //Sanity check
    if (numSPStates != num_sp_states){
      std::cout<< "numSPStates BIGSTICK : " << numSPStates << " num_sp_states : " << num_sp_states << std::endl;
      fmt::print(" Incorrect template file . \n");
      std::exit(1);}


    std::vector<std::vector<double> > coefficients_list(state, std::vector<double>(1, 0)); // 1-> dimension
    if (state <= smwf_info.num_eigenvectors){
      for(int i =0; i < state; i++){
        coefficients_list[i] = ReadCoefficients("mfdn_smwf{:03d}", smwf_info, i, numStatesPerFile);
      }
    }
    std::map<std::vector<uint16_t> , std::vector<double> > mb_states_B; // B for BIGSTICK
    findAndReplaceSPIndices(sp_state_list, sp_state_list_temp, mb_state_list, num_sp_states);
    std::vector<int16_t> phaseFactor;
    sortMBBasisStates(mb_state_list, mb_states_B, coefficients_list, phaseFactor, Z, N, state);
    
    /*
      //print mb_states
      auto stream1 = std::ofstream(out_filename, std::ios_base::out);
      
      stream1 << fmt::format("  {:>4d}   {:>4d}  \n", numParticles, num_states);
      for (int i = 0; i< num_states; i++){
        stream1 << fmt::format("  {:>4d}   || ", fmt::join(mb_state_list[i],"  "));
        stream1 << fmt::format("  {:>4d}   \n", fmt::join(mb_state_list_B[i],"  "));
      }
      std::cout<< "Writing list of MB states to output file .. " << num_states << " states" << std::endl;

    */

    auto stream3 = std::ofstream(out_filename, std::ios_base::out);
    std::vector<std::vector<double> > Coeffs;
    stream3<< fmt::format("  {:>4d} ! Z  \n", Z);
    stream3<< fmt::format("  {:>4d} ! N  \n", N);
    stream3<< fmt::format("  ! interaction file  \n");
    stream3<< fmt::format("  ! hw  \n");
    stream3<< fmt::format("  ! # of shell  \n");
    stream3<< fmt::format("  {:>4d} ! total number of p,n s.p. states  \n", num_sp_states);
    stream3<< fmt::format("  ! Nmax  \n");
    stream3<< fmt::format("  ! # of many-body configurations  \n");
    stream3<< fmt::format("  {:>4d} ! parity  \n", smwf_info.parity);
    stream3<< fmt::format("  {:>4d} ! 2 x Jz  \n", smwf_info.twoM);
    stream3<< fmt::format("  {:>4d} ! # of eigenstates  \n", state); 

    for(int i =0; i < state; i++){
          stream3 << fmt::format("  {:>4f}   {:>4f}  {:>4f}  \n", smwf_info.energy[i], smwf_info.J[i], smwf_info.T[i]);
    }
    // Print the state labels and n , l, 2J, Jz
    for(int nsp = 0 ; nsp < num_sp_states; nsp++){
      stream3 << fmt::format("  {:>4d}  ", nsp + 1);
      stream3 << fmt::format(" {:>4d}   \n", fmt::join(sp_state_list_temp[nsp] , "  "));
    }
    
    int mbstateCount = 0;
    for(int index =0; index< num_states; index++){
      auto it = mb_states_B.find(mb_state_list_temp[index]);
      if(it !=mb_states_B.end()){
        stream3 << fmt::format("  {:>4d}   \n", fmt::join(it->first ,"  "));
        stream3 << fmt::format("{:>4e}   \n", fmt::join(it->second, "   "));
        mbstateCount++;
      }
    }
      
    if (mbstateCount == num_states)
      std::cout<<"Validation successful .. " << std::endl;
  
  }

  return 0;
}

#endif  // PARTITIONUTILS_GROUP_READ_H_
