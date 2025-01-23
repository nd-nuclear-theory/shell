/****************************************************************
  read_wavefunctions.cpp

  adapted read functions from group_read.cpp by Patrick J. Fasano,  University of Notre Dame

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

  return info;
}

std::vector<uint16_t> GenerateSPBinPartitionTable(const MFDnSMWFInfo& info)
{
  std::vector<uint16_t> sp_bin_partition_table{};
  sp_bin_partition_table.reserve(info.num_proton_states + info.num_neutron_states);

  {
    auto partition_it = info.partitioning.proton_partitions.cbegin();
    auto partition_end = info.partitioning.proton_partitions.cend();
    for (int i = 1; i <= info.num_proton_states; ++i)
    {
      if ((std::next(partition_it) != partition_end) && (i >= *std::next(partition_it)))
        partition_it = std::next(partition_it);
      sp_bin_partition_table.push_back(*partition_it);
    }
  }

  {
    auto partition_it = info.partitioning.neutron_partitions.cbegin();
    auto partition_end = info.partitioning.neutron_partitions.cend();
    for (int i = info.num_proton_states + 1;
         i <= info.num_proton_states + info.num_neutron_states;
         ++i)
    {
      if ((std::next(partition_it) != partition_end) && (i >= *std::next(partition_it)))
        partition_it = std::next(partition_it);
      sp_bin_partition_table.push_back(*partition_it);
    }
  }

  return sp_bin_partition_table;
}

// Serial reading of Group IDs

void ReadMBGroups(
    std::string filename_pattern,
    const MFDnSMWFInfo& smwf_info,
    std::vector<std::vector<uint16_t> >& groupid_list,
    bool verbose = false
  )
{
  // convenience mode variable
  static constexpr std::ios_base::openmode mode_argument =
      std::ios_base::in | std::ios_base::binary;

  // get sp_bin
  const auto spbin_partition_table = GenerateSPBinPartitionTable(smwf_info);

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


std::vector<double> ReadCoefficients(int state, int num_states, std::string filename){

  mcutils::FileExistCheck(
      filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
    );
  std::vector<double> coeffs;
  auto stream = std::ifstream(filename, std::ios_base::binary);
  float buffer;
  int count =0;
  int stateCount = 0;
  bool flag = false;
  
  // quite unlikely but here is a sanity check to ensure number of bytes in the file 
  // matches the expected number 
  int file_size = stream.tellg();
  stream.seekg(0, std::ios_base::end);
  file_size = int(stream.tellg()) - file_size;
  if (file_size % ((num_states + 2) * sizeof(float)) != 0){ // 2 is the offset number of bytes after each wf
    std::cout<< "Corrupted file : Unexpected size " << std::endl;
    std::exit(1); // Perhaps a different kind of error must be thrown
  }
  stream.seekg(std::ios_base::beg); // reset position back to beginning of the stream
  while(stream.read(reinterpret_cast<char*>(&buffer), sizeof(float))){
    if(count < num_states){
      if(stateCount==state){
          coeffs.push_back(buffer);
          flag = true;
        }
      count++;
    }else{
      if (flag) break;
      stream.read(reinterpret_cast<char*>(&buffer), sizeof(float));
      count = 0; // reset counter
      stateCount++;
    }
  }

  return coeffs;
}

int main(int argc, char* argv[])
{
  const uint16_t numParticles = 6;
  MBGroupsMetadata metadata{};
  const auto smwf_info = ReadMFDnSMWFInfo("mfdn_smwf.info");
  int num_sp_states = smwf_info.num_proton_states + smwf_info.num_neutron_states;
  int two_mj = smwf_info.twoM;
  std::size_t num_states = smwf_info.dimension;
  int state = 0; // 0 for lowest eigen wavefunction aka ground state in mfdn_smwf001
  std::vector<double> coefficients = ReadCoefficients(state, num_states, "mfdn_smwf001");

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
  
  ReadMBGroups("mfdn_MBgroups{:03d}", smwf_info, groupid_list, false);
  fmt::print("number of groups: {:d}\n", groupid_list.size());
  fflush(stdout);
  
  std::vector<std::vector<uint16_t> > mb_state_list(num_states, std::vector<uint16_t>( numParticles,0)); 
  // Create the mj2_sp vector that contains the 2M values of all the single particle states
  std::vector<int> mj2_sp;
  
  int count = 0;

  const auto& proton_subspace = smwf_info.orbital_space.GetSubspace(0);
  for (int index = 0; index < proton_subspace.size(); ++index)
  {
    for(int j = -1 * TwiceValue(proton_subspace.GetState(index).j()) ; j <= TwiceValue(proton_subspace.GetState(index).j()); j+=2)
      mj2_sp.push_back(j);
  }
  int num_proton_sp_states = mj2_sp.size();
  const auto& neutron_subspace = smwf_info.orbital_space.GetSubspace(0);

  for (int index = 0; index < neutron_subspace.size(); ++index)
  {
    for(int j = -1 * TwiceValue(neutron_subspace.GetState(index).j()) ; j <= TwiceValue(neutron_subspace.GetState(index).j()); j+=2)
      mj2_sp.push_back(j);
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
  
/*
  // write to a text file
  std::string out_filename_txt = "out.txt";
  auto stream1 = std::ofstream(out_filename_txt, std::ios_base::out);
  //mcutils::WriteBinary<int>(stream(),bytes);
  
  stream1 << fmt::format("  {:>4d}   {:>4d}  \n", numParticles, num_states) << std::endl;
  for (uint16_t i = 0; i< num_states; i++){
    stream1 << fmt::format("  {:>4d}   {:+16.7e}  \n", fmt::join(mb_state_list[i],"  "), coefficients[i]) <<std::endl;
  }
  std::cout<< "Writing to output file .. " << num_states << " states" << std::endl;
*/
  std::string out_filename_bin = "out5.MBO";
  auto stream2 = std::ofstream(out_filename_bin, std::ios_base::binary);
  std::vector<uint16_t> buffer(numParticles , 0);
  
  mcutils::WriteBinary(stream2, &numParticles, 1);
  mcutils::WriteBinary(stream2, &num_states, 1);

  for(int i = 0; i< num_states; i++){
    for(int j =0; j<numParticles; j++){
      buffer[j] = (mb_state_list[i])[j];
    }

    mcutils::WriteBinary(stream2, buffer.data(), numParticles);
    mcutils::WriteBinary(stream2, &(coefficients[i]), 1);
    //stream2.write(reinterpret_cast<const char*>(&(coefficients[i])),sizeof(double));

  }

  return 0;
}

// }  // namespace
#endif  // PARTITIONUTILS_GROUP_READ_H_
