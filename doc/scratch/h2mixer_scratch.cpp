/*

h2mixer

out1
Output stream -- out1
  File: out1.dat
  Format: 15099 (text)
  Orbitals: p 6 n 6 (oscillator-like true)
  One-body truncation: p 2.0000 n 2.0000
  Two-body truncation: pp 2.0000 nn 2.0000 pn 2.0000
  Sectors: pp 7 nn 7 pn 7 => total 21
  Matrix elements: pp 32 nn 32 pn 110 => total 174
0x800511b0 0
0x800511b0 0

Closing right away...
out2
Output stream -- out2
  File: out2.dat
  Format: 15099 (text)
  Orbitals: p 6 n 6 (oscillator-like true)
  One-body truncation: p 2.0000 n 2.0000
  Two-body truncation: pp 2.0000 nn 2.0000 pn 2.0000
  Sectors: pp 7 nn 7 pn 7 => total 21
  Matrix elements: pp 32 nn 32 pn 110 => total 174
0x80050d80 0
0x80050d80 0

Closing right away...
Closing output stream out1...
0x800511b0 0
  File: out1.dat
  Format: 15099 (text)
  Orbitals: p 6 n 6 (oscillator-like true)
  One-body truncation: p 2.0000 n 2.0000
  Two-body truncation: pp 2.0000 nn 2.0000 pn 2.0000
  Sectors: pp 7 nn 7 pn 7 => total 21
  Matrix elements: pp 32 nn 32 pn 110 => total 174

Segmentation fault (core dumped)

*/


void OpenTargetStreams(
    const RunParameters& run_parameters,
    const TwoBodyOperatorIndexing& target_indexing,
    std::map<std::string,shell::OutH2Stream>& target_streams
    //std::vector<shell::OutH2Stream>& target_streams
  )
// Open target streams.
{
  for (const auto& id_filename : run_parameters.target_filenames)
    {
      // extract key and value
      const std::string& target_id = id_filename.first;
      const std::string& target_filename = id_filename.second;

      std::cout << target_id << std::endl;

      // create stream
      
      // DEBUGGING: segfault
      //
      // target_streams.insert(
      //     {
      //       target_id,
      //         shell::OutH2Stream(
      //             target_filename,
      //             target_indexing.orbital_space,target_indexing.space,target_indexing.sectors,
      //             run_parameters.output_h2_format
      //           )
      //         }
      //   );

      target_streams[target_id]
        = shell::OutH2Stream(
            target_filename,
            target_indexing.orbital_space,target_indexing.space,target_indexing.sectors,
            run_parameters.output_h2_format
          );

      // write diagnostic
      std::cout << fmt::format("Output stream -- {}",target_id) << std::endl;
      std::cout << target_streams[target_id].DiagnosticStr();
      std::cout << target_streams[target_id].stream_ptr() << " " << target_streams[target_id].stream().is_open()  << std::endl;
      std::cout << target_streams[target_id].stream_ptr() << " " << target_streams[target_id].stream().is_open()  << std::endl;
      std::cout << std::endl;

      // write diagnostic
      std::cout << fmt::format("Closing right away...",target_id) << std::endl;

      target_streams[target_id].Close();
      target_streams[target_id].Close();
    }

  // iterate over output streams -- TEST
  for (auto&& id_stream : target_streams)
    {
      // extract key and value
      const std::string& target_id = id_stream.first;
      shell::OutH2Stream& target_stream = id_stream.second;

      // write diagnostic
      std::cout << fmt::format("Closing output stream {}...",target_id) << std::endl;
      std::cout << target_streams[target_id].stream_ptr() << " " << target_streams[target_id].stream().is_open()  << std::endl;

      std::cout << target_stream.DiagnosticStr();
      std::cout << std::endl;

      // close
      target_stream.Close();

      std::cout << "Done!" << std::endl;

    }
  
}

void CloseTargetStreams(
    std::map<std::string,shell::OutH2Stream>& target_streams
  )
// Close target streams.
{

  // iterate over output streams
  for (auto& id_stream : target_streams)
    {
      // extract key and value
      const std::string& target_id = id_stream.first;
      shell::OutH2Stream& target_stream = id_stream.second;

      // write diagnostic
      std::cout << fmt::format("Closing output stream {}...",target_id) << std::endl;

      // close
      target_stream.Close();

      std::cout << "Done!" << std::endl;

    }
  
}


void GenerateOutputs(
    const typename basis::TwoBodySectorsJJJPN::SectorType& target_sector,
    std::map<std::string,shell::OutH2Stream>& target_streams
)
{

  // iterate over output streams
  for (auto& id_stream : target_streams)
    {
      // extract key and value
      const std::string& target_id = id_stream.first;
      shell::OutH2Stream& target_stream = id_stream.second;

      // write stream matrix
      Eigen::MatrixXd target_matrix=Eigen::MatrixXd::Zero(target_sector.bra_subspace().size(),target_sector.ket_subspace().size());
      target_stream.WriteSector(target_matrix);

      // progress indicator
      std::cout << "." << std::flush;
    }

}





  run_parameters.target_filenames.insert({"out1","test/out1.dat"});
  run_parameters.target_filenames.insert({"out2","test/out2.dat"});


  // direct input stream names 
  std::map<std::string,std::string> direct_input_source_definitions;  // map from source id to input filename
  // target stream names 
  std::map<std::string,std::string> target_filenames;  // map from target id to target filename



      // identity
      if (operator_channel.operator_type==OperatorType::kIdentity)
        operator_matrix = shell::IdentityOperatorMatrixJJJPN(target_sector,run_parameters.A);
      // kinematic
      else if (operator_channel.operator_type==OperatorType::kKinematicVRSqr)
        {
          RadialOperatorData radial_operator_data = radial_operators[{shell::RadialOperatorType::kR,2}];
          bool momentum_space = false;
          operator_matrix = shell::KinematicVTSqrMatrixJJJPN(
              radial_operator_data.ket_orbital_space,
              radial_operator_data.sectors,
              radial_operator_data.matrices,
              momentum_space,
              target_sector,run_parameters.A
          );
        }
      else if (operator_channel.operator_type==OperatorType::kKinematicVR1R2)
        {
          RadialOperatorData radial_operator_data = radial_operators[{shell::RadialOperatorType::kR,1}];
          bool momentum_space = false;
          operator_matrix = shell::KinematicVT1T2MatrixJJJPN(
              radial_operator_data.ket_orbital_space,
              radial_operator_data.sectors,
              radial_operator_data.matrices,
              momentum_space,
              target_sector,run_parameters.A
          );
        }
      else if (operator_channel.operator_type==OperatorType::kKinematicVKSqr)
        {
          RadialOperatorData radial_operator_data = radial_operators[{shell::RadialOperatorType::kK,2}];
          bool momentum_space = true;
          operator_matrix = shell::KinematicVTSqrMatrixJJJPN(
              radial_operator_data.ket_orbital_space,
              radial_operator_data.sectors,
              radial_operator_data.matrices,
              momentum_space,
              target_sector,run_parameters.A
          );
        }
      else if (operator_channel.operator_type==OperatorType::kKinematicVK1K2)
        {
          RadialOperatorData radial_operator_data = radial_operators[{shell::RadialOperatorType::kK,1}];
          bool momentum_space = true;
          operator_matrix = shell::KinematicVT1T2MatrixJJJPN(
              radial_operator_data.ket_orbital_space,
              radial_operator_data.sectors,
              radial_operator_data.matrices,
              momentum_space,
              target_sector,run_parameters.A
          );
        }
      // angular momentum square
      //
      // TODO neaten by mapping id to (family,species)
      else if (operator_channel.operator_type==OperatorType::kAMSqrLp)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrLn)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrL)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kOrbital,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrSp)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrSn)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrS)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kSpin,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrJp)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kP,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrJn)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kN,
            target_sector,run_parameters.A
          );
      else if (operator_channel.operator_type==OperatorType::kAMSqrJ)
        operator_matrix = shell::AngularMomentumMatrixJJJPN(
            shell::AngularMomentumOperatorFamily::kTotal,shell::AngularMomentumOperatorSpecies::kTotal,
            target_sector,run_parameters.A
          );



// enum class OperatorType
// // Possible types of generated operators.
// {
//   // identity
//   kIdentity,
//   // kinematic
//     
//   // angular momentum square
//     kAMSqrLp,kAMSqrLn,kAMSqrL,
//     kAMSqrSp,kAMSqrSn,kAMSqrS,
//     kAMSqrJp,kAMSqrJn,kAMSqrJ
//     };
// const std::array<const char*,3> kOperatorTypeName({"identity","kinematic","am-sqr"});
// std::map<std::string,OperatorType> kOperatorTypeLookup(
//     {
//       {"identity",OperatorType::kIdentity},
//         {"Lp",OperatorType::kAMSqrLp},{"Ln",OperatorType::kAMSqrLn},{"L",OperatorType::kAMSqrL},
//                                                                       {"Sp",OperatorType::kAMSqrSp},{"Sn",OperatorType::kAMSqrSn},{"S",OperatorType::kAMSqrS},
//         {"Jp",OperatorType::kAMSqrJp},{"Jn",OperatorType::kAMSqrJn},{"J",OperatorType::kAMSqrJ}
//                                                                                                                                     }
//   );

          if (operator_channel.id=="Ursqr")
            {
              kinematic_operator_type = shell::KinematicOperatorType::kUTSqr;
              radial_operator_type = shell::RadialOperatorType::kR;
              radial_operator_power = 2;
            }
          else if (operator_channel.id=="Vr1r2")
            {
              kinematic_operator_type = shell::KinematicOperatorType::kVT1T2;
              radial_operator_type = shell::RadialOperatorType::kR;
              radial_operator_power = 1;
            }
          else if (operator_channel.id=="Uksqr")
            {
              kinematic_operator_type = shell::KinematicOperatorType::kUTSqr;
              radial_operator_type = shell::RadialOperatorType::kK;
              radial_operator_power = 2;
            }
          else if (operator_channel.id=="Vk1k2")
            {
              kinematic_operator_type = shell::KinematicOperatorType::kVT1T2;
              radial_operator_type = shell::RadialOperatorType::kK;
              radial_operator_power = 1;
            }


  std::cout
    << fmt::format(
        "{} {} {}",
        char(std::get<0>(radial_operators.begin()->second.labels)),
        std::get<1>(radial_operators.begin()->second.labels),
        radial_operators.begin()->second.sectors.size()
      )
    << std::endl;


std::set<std::string> kKinematicOperatorIdSet(
    {"Ursqr","Vr1r2","Uksqr","Vk1k2"}
  );  // redundant to kKinematicOperatorDefinitions
std::set<std::string> kAMSqrOperatorIdSet(
    {"Lp","Ln","L","Sp","Sn","S","Jp","Jn","J"}
  );
