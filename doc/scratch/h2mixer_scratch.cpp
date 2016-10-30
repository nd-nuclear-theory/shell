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
