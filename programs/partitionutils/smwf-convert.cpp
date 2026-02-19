/****************************************************************
  smwf-convert.cpp

  adapted read functions from group_read.cpp by Patrick J. Fasano,  University of Notre Dame

  Shwetha Vittal
  University of Notre Dame

  + 01/07/25 (slv): First version of code to read many body states from mfdn_MBgroupsxxx
  + 01/22/25 (slv): Create ReadCoefficients function to read the amplitudes, of the 
    many body states, from mfdn_smwf001 file
  + 06/30/35 (mac):
    + Support mismatched sp bases sizes between smwf and trwfn.
    + Update argument processing.
****************************************************************/

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

struct TrwfnInfo{
  int Z;
  int N;
  int num_shells;
  int num_sp_states;
  int Nmax;
  int parity;
  int two_Jz;
};

////////////////////////////////////////////////////////////////
// process arguments
/////////////////////////////////////////////////////////////////

struct RunParameters
// Stores simple parameters for run
{
  // filenames
  std::string mode;
  int state_index;
  std::string source_wf_dir;
  std::string output_filename;
  std::string template_filename;

  // default constructor
  RunParameters()
    : state_index(0), mode(""), source_wf_dir(""), output_filename(""), template_filename("")
  {}

};

void PrintUsage(char **argv) {
  std::cout << "Usage: " << argv[0]
            << " mode seq source_wf_dir [template_filename] output_filename"
            << std::endl
            << std::endl
            << "Modes:" << std::endl
            << "  mbstates: print MB basis states" << std::endl
            << "  mbo: generate MBO" << std::endl
            << "  trwfn: generate trwfn for selected state index (1-based sequence number)" << std::endl
            << std::endl
            << "  A template file defining indexing is required for trwfn conversion." << std::endl;
}

void ProcessArguments(int argc, char *argv[], RunParameters& run_parameters)
{
  
  int arg = 1;
  
  // process options
  while (arg < argc && argv[arg][0] == '-')
    {
      std::istringstream parameter_stream(argv[arg++]);

      if (parameter_stream.str() == "--help" || parameter_stream.str() == "-h")
        {
          PrintUsage(argv);
          std::exit(EXIT_SUCCESS);
        }
      else
        {
          PrintUsage(argv);
          std::cerr << "Unrecognized option '" << parameter_stream.str() << "'" << std::endl;
          std::exit(EXIT_FAILURE);
        }
    }
  
  // process fixed arguments
  if (argc-arg < 4)
    {
      PrintUsage(argv);
      std::cerr << "Insufficient arguments" << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // mode
  run_parameters.mode = argv[arg++]; 
  if(run_parameters.mode != "mbstates" && run_parameters.mode != "mbo" && run_parameters.mode != "trwfn"){
    std::cerr << "ERROR: Unrecognized value for mode" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  // state
  std::istringstream state_stream(argv[arg++]);
  std::size_t seq;
  state_stream >> seq;
  if (!state_stream || seq==0)
    {
      std::cerr << "ERROR: Expecting positive numerical value for state index." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  run_parameters.state_index = seq-1;

  // source wf directory
  run_parameters.source_wf_dir = argv[arg++];
  
  // template filename (optional)
  if (run_parameters.mode == "trwfn"){
    run_parameters.template_filename = argv[arg++];
  }

  // output filename
  run_parameters.output_filename = argv[arg++];

}

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

int SetLastMj(const uint16_t num_particles, int num_sp_states, 
              std::vector<int> &mj2_sp, int two_mj, 
              std::vector<int> &next_bin, 
              std::vector<uint16_t> &temp_state, int flag){
                //int flag = 1;
                int delta_mj = two_mj;
                for (uint16_t i = 0; i < num_particles; i++)
                  { //std::cout << mj2_sp[temp_state[i] ] << std::endl; 
                    delta_mj -= mj2_sp[temp_state[i] ]; // because mj2 indices begin from 0 in c++ whereas they begin from 1 in Fortran
                    }
                
                if (delta_mj< 0) {
                  //std::cout<< "Flag from SetLastMj : " << flag << std::endl;
                  flag = 1;
                  return flag;}
                
                int iLast = temp_state[num_particles-1] + delta_mj/2;
                  //std::cout << "iLast = " << iLast << std::endl;
                
                if (iLast< next_bin[temp_state[num_particles -1]]){
                  //std::cout << "next_bin[temp_state[num_particles -1] ] " << next_bin[temp_state[num_particles -1] ] << std::endl;
                  temp_state[num_particles-1]= iLast;
                  //fmt::print("New temp_state  {:>4d}\n", fmt::join(temp_state," "));
                  flag = 0;
                }
                else flag = -1; 
                //std::cout<< "Flag from SetLastMj : " << flag << std::endl;
                return flag;
              }

int IncrementMj(const uint16_t num_particles, int num_sp_states, 
                  std::vector<int> &mj2_sp, int two_mj, 
                  std::vector<int> &next_bin, 
                  std::vector<uint16_t> &mbgroup,
                  std::vector<uint16_t> &mbstate, int flag){
                  
                  for(int i = num_particles-2; i >= 0; i--){ // index i goes from 4 -> 0 when there are 6 particles
                    if (mbstate[i] < (next_bin[mbgroup[i]] -1)){
                      mbstate[i] += 1;
                      for(int j = i+1; j < num_particles; j++){
                        //std::cout<< "mbgroup[j] " << mbgroup[j] << std::endl;
                        //std::cout<< "mbstate[j-1] " << mbstate[j-1] << " i "<< i << " j "<< j <<std::endl;

                        if ((mbstate[j-1] + 1) >= (next_bin[mbgroup[j]])){
                        
                          //std::cout<< "Flag is 2 "<<std::endl;
                          //fmt::print(" mbgroup {:>4d}\n", fmt::join(mbgroup," "));
                          //fmt::print(" mbstate {:>4d}\n", fmt::join(mbstate," "));
                          flag = 2;
                          break; 
                        }
                        else{
                          mbstate[j] = mbstate[j-1] +1;
                          mbstate[j] = std::max(mbstate[j], mbgroup[j]);
                        }
                      }
                      if (flag==2){
                        flag = 0;
                        //std::cout<< "flag is 2 here----------------------------" << " i " << i << std::endl;
                        continue;
                      }
                      flag = SetLastMj(num_particles, num_sp_states, mj2_sp, two_mj, next_bin, mbstate, flag);

                      if (flag ==1) continue;
                        else {
                          //std::cout<< "Flag from IncrementMj : " << flag << std::endl;
                          return flag;}
                    }
                  }
                  flag = 1;
                  //std::cout<< "Flag from IncrementMj : " << flag << std::endl;
                  return flag;
                }

int MjStatesGen(const uint16_t num_particles, int num_sp_states, 
                  std::vector<int> &mj2_sp, int two_mj, 
                  std::vector<int> &next_bin, 
                  std::vector<uint16_t> &tempVar, int num_states, 
                  std::vector<std::vector<uint16_t> > &mb_state_list,
                  int currentState){
  /****************************************************************
    Functions just like the subroutine mfdn_transitions/src/module_MjStates/MjStatesGen

    num_particles : Total number of particles
    num_sp_states : number of single particle states (it is assumed that there are same number of
                    sp states in both species)
    mj2_sp : contains a list of 2*mj values for each single particle state
    two_mj : the selected 2*mJ
    next_bin : contains the index of the next partition a list of size 1 less than sizeof(mj2_sp)
    tempVar : single GroupID from a list of GroupIDs, the state that marks the beginning of the 
              group(the set of many body states in the partition). It is called ID but it is a 
              vector (of size the num_particles) with sp state indices.
    num_states : Total number of states 
    mb_state_list : List to be updated with the states matching the criteria of having same 2*mj
    currentState : count of number of states matching the criteria of having same 2*mj

    *****************************************************************/

          std::vector<uint16_t> mbstate = tempVar;
          std::vector<uint16_t> mbgroup = mbstate;
          //std::cout<< "MjStatesGen is running .. " << std::endl;
          int flag = 0; 
          flag = SetLastMj(num_particles, num_sp_states, mj2_sp, two_mj, next_bin, mbstate, flag);

          if (flag==0){
            //std::cout<< "Adding a state to mb_state_list---------------------------" << currentState <<std::endl;

            // Sanity check
            int total2mj =0;
            for(int i =0; i< num_particles; i++){
              total2mj +=mj2_sp[mbstate[i]];
            }
            if (total2mj != two_mj) fmt::print(" Wrong MJ {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            
            //fmt::print(" mbstate {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            mb_state_list[currentState] = mbstate;
            currentState += 1;
            }
          flag = 0;
          while(flag==0){
            flag = -1;
            while(flag == -1){
              //std::cout << "flag is -1" << std::endl;
              flag = IncrementMj(num_particles, num_sp_states, mj2_sp, two_mj, next_bin, mbgroup, mbstate, flag); // flag = 1 is the end of loop condition
              //std::cout << "flag is " << flag << std::endl;
            }
            if(flag ==0){
              //std::cout<< "Adding a state to mb_state_list---------------------------" << currentState <<std::endl;
              //std::cout << "flag is 0" << std::endl;
              // Sanity check
            int total2mj =0;
            for(int i =0; i< num_particles; i++){
              total2mj +=mj2_sp[mbstate[i]];
            }
            if (total2mj != two_mj) fmt::print(" Wrong MJ {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            
            //fmt::print(" mbstate {:>4d}------------------------- {:d}\n", fmt::join(mbstate," "), total2mj);
            mb_state_list[currentState] = mbstate;
            currentState += 1;
            }
          }
          return currentState;
        }

std::vector<double> ReadCoefficients(std::string filename_pattern,
                                    const MFDnSMWFInfo& smwf_info,
                                    int state_index, 
                                    std::vector<int>& numStatesPerFile,
                                    bool verbose = true){
  /****************************************************************
  Reads the coefficients of a wavefunction
  
  filename_pattern : "mfdn_smwf{:03d}"
  smwf_info        : consists of wavefunction information from mfdn_smwf.info
  state_index            : 0 for ground state, 1 for next eigen state and so on.
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
      std::exit(EXIT_FAILURE);
    }

    // TO DO (slv): Need to document this 
    int offset = sizeof(float) * (2* state_index +1);
    stream.seekg(std::ios_base::beg + numStatesPerFile[i-1] * sizeof(float) * state_index + offset); // set position back to beginning of the state in the stream
    
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


void ReadTrwfn(
    std::string filename,
    TrwfnInfo &trwfn_info,
    std::vector<std::vector<int16_t> > &sp_state_list_template,
    std::vector<std::vector<uint16_t> > &mb_state_list_template
  ){
  /* Reads list of SP states  from a model trwfn file
    
    filename            : Trwfn filename including the path
    sp_state_list_template  : List of single particle states in trwfn file with columns {n, l, twice_j, twice_mj, species_code}
    mb_state_list_template  : List of many body states in trwfn file with proton and neutron SP indices in ascending order

  */

  // open trwfn file
  fmt::print(" Reading trwfn file .. \n");
  int count =0;

  std::string line;
  int line_count = 0;
  mcutils::FileExistCheck(
      filename, /*exit_on_nonexist=*/true, /*warn_on_overwrite=*/false
    );
  auto stream = std::ifstream(filename, std::ios_base::in);

  for(int i =0; i<4; i++)
  {
    mcutils::GetLine(stream, line, line_count);
    if(i ==0){
      std::istringstream line_stream(line);
      line_stream >> trwfn_info.Z;
    }
    if(i ==1){
      std::istringstream line_stream(line);
      line_stream >> trwfn_info.N;
    }
  }

  {
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> trwfn_info.num_shells;
  }
  
  std::size_t num_sp_states;
  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> trwfn_info.num_sp_states;
  }

  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> trwfn_info.Nmax;
  }
   mcutils::GetLine(stream, line, line_count); // ignored reading number of many-body configuration
  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> trwfn_info.parity;
  }

  {  
    mcutils::GetLine(stream, line, line_count);
    std::istringstream line_stream(line);
    line_stream >> trwfn_info.two_Jz;
  }

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
    for (std::size_t sps_line_count = 0; sps_line_count < trwfn_info.num_sp_states;
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
        sp_state_list_template.push_back(state);
        //spstate_list[count][1]= index;
        count++;
      }
  }
  {
    while(mcutils::GetLine(stream, line, line_count)){
      // set up for parsing
      std::istringstream line_stream(line);
      std::vector<uint16_t> spIndices(trwfn_info.Z + trwfn_info.N, 0);
      for(int i = 0; i < trwfn_info.Z + trwfn_info.N; i++ ){
        line_stream>>spIndices[i];
      }
      mb_state_list_template.push_back(spIndices);
      mcutils::GetLine(stream, line, line_count); // Ignore next line
    }
  }

  if(count != trwfn_info.num_sp_states){
    std::cerr << "ERROR: Missing sp states in trwfn file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
}

void FindAndReplaceSPIndices(std::vector<std::vector<int16_t> > &sp_state_list, 
                            std::vector<std::vector<int16_t> > &sp_state_list_template,
                            std::vector<std::vector<uint16_t> > &mb_state_list,
                            int num_sp_states){
  /*
  Looks up indices of states in sp_states_list from sp_state_list_template (temp for template)
  Replace sp states in mb_state_list with indices found by lookup.
  
  sp_state_list       : List of single particle states in the order that MFDn operates.
                        (n, l, j ,mj, t) 5 columns in the same order as in trwfn
  sp_state_list_template  : List of single particle states in the order of trwfn (Look at docstring for ReadTrwfn)
  mb_state_list       : List of many body states from MFDn
  num_sp_states       : Number of single particle states

  */
  fmt::print(" Replacing SP indices .. \n ");
  std::vector<int> index_list(num_sp_states, 0);
  for(int i =0; i < num_sp_states; i++){
    auto it =  std::find(sp_state_list_template.begin(), sp_state_list_template.end(), sp_state_list[i]);
    if(it != sp_state_list_template.end()){
      index_list[i] = it - sp_state_list_template.begin() + 1; // finding index of the element
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

void SortMBBasisStates(std::vector<std::vector<uint16_t> > &mb_state_list,
                      std::map<std::vector<uint16_t> , std::vector<double> > &mb_states_bigstick,
                      std::vector<std::vector<double> > &coefficients_list,
                      std::vector<int16_t> &phaseFactor,// Not needed anymore
                      uint16_t Z, uint16_t N,
                      int num_eigenvectors){
  /*
  Sorts the single particle state(SPS) indices in the mbstates 
  
  mb_state_list     : List of many body states with SPS indices of trwfn
  mb_states_bigstick       : Map with sorted SPS indices as *key* and coefficients of 
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
    mb_states_bigstick[mbstate] = Coeffs;
  }
}

////////////////////////////////////////////////////////////////
// main program
////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // header
  std::cout << std::endl;
  std::cout << "smwf-convert -- convert MFDn wavefunctions " << std::endl;
  std::cout << std::endl;
  
  // read parameters
  RunParameters run_parameters;
  ProcessArguments(argc, argv, run_parameters);

  int state_index = run_parameters.state_index;
  std::string source_wf_dir = run_parameters.source_wf_dir;
  

  MBGroupsMetadata metadata{};
  const auto smwf_info = ReadMFDnSMWFInfo(source_wf_dir + "/mfdn_smwf.info");
  const uint16_t N = smwf_info.N;
  const uint16_t Z = smwf_info.Z;
  const uint16_t num_particles = N + Z;

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
  std::vector<int> numStatesPerFile = ReadMBGroups(source_wf_dir + "/mfdn_MBgroups{:03d}", smwf_info, groupid_list, false);

  fmt::print("number of groups: {:d}\n", groupid_list.size());
  fflush(stdout);
  
  std::vector<std::vector<uint16_t> > mb_state_list(num_states, std::vector<uint16_t>( num_particles,0)); 
  std::vector<std::vector<int16_t> > sp_state_list; 
  // Create the mj2_sp vector that contains the 2M values of all the single particle states
  std::vector<int> mj2_sp;
  // Create the n_sp, l_sp and twoJ_sp vectors which needs to be printed for mode=2 i.e., smwf --> trwfn
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
      // WARNING 06/30/25 (mac): This variable "state" overloads identically named variable at function scope above.
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
      // WARNING 06/30/25 (mac): This variable "state" overloads identically named variable at function scope above.
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

  int current_num_states = 0;
  for(std::vector<std::vector<uint16_t> >::iterator it = groupid_list.begin(); it != groupid_list.end(); it++ )
  {
    std::vector<uint16_t> groupID = *it;
    // Imitating Fortran subroutine MjStatesGen
    //tempVar = {0, 1, 4, 40, 44, 54}; // test case
    current_num_states = MjStatesGen(num_particles, num_sp_states, mj2_sp, two_mj, next_bin, 
                       groupID, num_states, mb_state_list, current_num_states); 

  }  
  

  // write list of MB states to a text file
  if (run_parameters.mode == "mbstates")  {
    auto output_stream = std::ofstream(run_parameters.output_filename, std::ios_base::out);
    
    output_stream << fmt::format("  {:>4d}   {:>4d}  \n", num_particles, num_states);
    for (int i = 0; i< num_states; i++){
      output_stream << fmt::format("  {:>4d}   \n ", fmt::join(mb_state_list[i],"  "));
    }
    std::cout<< "Writing list of MB states to output file .. " << num_states << " states" << std::endl;
  }


  if (run_parameters.mode == "mbo"){
    std::vector<double> coefficients = ReadCoefficients(source_wf_dir + "/mfdn_smwf{:03d}", smwf_info, state_index, numStatesPerFile); 
    // for(std::vector<double>::iterator it = coefficients.begin(); it !=coefficients.end(); it++)
    //   std::cout<< *it<<std::endl;

    fmt::print("number of states: {:d}\n", coefficients.size());

    auto output_stream = std::ofstream(run_parameters.output_filename, std::ios_base::binary);
    std::vector<uint16_t> buffer(num_particles , 0);
    
    mcutils::WriteBinary(output_stream, &num_particles, 1);
    mcutils::WriteBinary(output_stream, &num_states, 1);

    // Here I need to implement the change of indices
    for(int i = 0; i< num_states; i++){
      for(int j =0; j<num_particles; j++){
        
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
      mcutils::WriteBinary(output_stream, buffer.data(), num_particles);
      mcutils::WriteBinary(output_stream, &(coefficients[i]), 1);
      //output_stream.write(reinterpret_cast<const char*>(&(coefficients[i])),sizeof(double));

    }
  }

  if(run_parameters.mode == "trwfn"){

    // read model trwfn file
    TrwfnInfo trwfn_info;
    std::vector<std::vector<int16_t> > sp_state_list_template;
    std::vector<std::vector<uint16_t> > mb_state_list_template;
    ReadTrwfn(run_parameters.template_filename, trwfn_info, sp_state_list_template, mb_state_list_template);
    
    // check dimensions of single particle bases
    if (trwfn_info.num_sp_states != num_sp_states){
      std::cout << fmt::format("WARNING: Mismatched number of sp states: trwfn template {}, smwf {}", trwfn_info.num_sp_states, num_sp_states) << std::endl;
    }

    // choose state
    std::vector<std::vector<double> > coefficients_list(1, std::vector<double>(1, 0)); // 1-> dimension
    if (state_index <= smwf_info.num_eigenvectors){
      coefficients_list[0] = ReadCoefficients(source_wf_dir + "/mfdn_smwf{:03d}", smwf_info, state_index, numStatesPerFile);
    }
    
    // remap sp basis
    std::map<std::vector<uint16_t> , std::vector<double> > mb_states_bigstick;
    FindAndReplaceSPIndices(sp_state_list, sp_state_list_template, mb_state_list, num_sp_states);
    std::vector<int16_t> phaseFactor;
    SortMBBasisStates(mb_state_list, mb_states_bigstick, coefficients_list, phaseFactor, Z, N, 1);

    // write trwfn header
    auto output_stream = std::ofstream(run_parameters.output_filename, std::ios_base::out);
    std::vector<std::vector<double> > Coeffs;
    output_stream << fmt::format("  {:>4d} ! Z", Z) << std::endl;
    output_stream << fmt::format("  {:>4d} ! N", N) << std::endl;
    output_stream << fmt::format("  ! interaction file") << std::endl;
    output_stream << fmt::format("  {:>4d} ! hw", 0) << std::endl;
    output_stream << fmt::format("  {:>4d} ! # of shell", trwfn_info.num_shells) << std::endl;
    output_stream << fmt::format("  {:>4d} ! total number of p,n s.p. states", trwfn_info.num_sp_states) << std::endl;
    output_stream << fmt::format("  {:>4d} ! Nmax", trwfn_info.Nmax) << std::endl;
    output_stream << fmt::format("  {:>4d} ! # of many-body configurations", num_states) << std::endl;
    output_stream << fmt::format("  {:>4d} ! parity", trwfn_info.parity) << std::endl; // This should be same as smwf_info.parity
    output_stream << fmt::format("  {:>4d} ! 2 x Jz", trwfn_info.two_Jz) << std::endl; // This should be same as smwf_info.two_M
    output_stream << fmt::format("  {:>4d} ! # of eigenstates", 1) << std::endl; 
    output_stream << fmt::format("  {:>4f}   {:>4f}  {:>4f}", smwf_info.energy[state_index], smwf_info.J[state_index], smwf_info.T[state_index]) << std::endl;

    // print sp state listing
    for(int nsp = 0 ; nsp < trwfn_info.num_sp_states; nsp++){
      output_stream << fmt::format("  {:>4d}  ", nsp + 1);
      output_stream << fmt::format(" {:>4d}  ", fmt::join(sp_state_list_template[nsp] , "  ")) << std::endl;
    }

    // print mb states and amplitudes
    int mbstateCount = 0;
    for(int index = 0; index < num_states; index++){
      auto it = mb_states_bigstick.find(mb_state_list_template[index]);
      if(it !=mb_states_bigstick.end()){
        output_stream << fmt::format("  {:>4d}", fmt::join(it->first ,"  ")) << std::endl;
        output_stream << fmt::format("{:>4e}", fmt::join(it->second, "   ")) << std::endl;
        mbstateCount++;
      }
    }
      
    if (mbstateCount != num_states) {
      std::cout  << "Mismatched number of states" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    
  }

  return 0;
}
