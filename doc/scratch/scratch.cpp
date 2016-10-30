  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: version number
  ////////////////////////////////////////////////////////////////

  // Yes, I know using template specialization as a way to distinguish
  // the different implementations is overkill, but isn't it fun?

  template<H2Mode tMode>
  void ReadVersion(std::ifstream stream, InH2Stream& h2_stream);
  template<H2Mode tMode>
  void WriteVersion(std::ofstream stream, OutH2Stream& h2_stream);

  // read version: text mode
  template <>
  void ReadVersion<H2Mode::kText>(std::ofstream stream,OutH2Stream& h2_stream)
  {
    std::string line;
    std::getline(*stream_, line);
    std::stringstream(line) >> header_.format;  // OOPS... requires access to private member
  };

  template <>
  void WriteVersion<H2Mode::kText>(std::ofstream stream,OutH2Stream& h2_stream)
  {
    stream
      << fmt::format("{:10d}",int(h2_format())) << std::endl;
  };

  template <>
  void WriteVersion<H2Mode::kBinary>(std::ofstream stream,OutH2Stream& h2_stream)
  {
    int bytes = kIntegerSize;
    int version_number = int(h2_stream.h2_format());

    // write
    WriteI4(stream,bytes);
    WriteI4(stream,version_number);
    WriteI4(stream,bytes);
  };


        ++line_count;
    std::getline(std::cin,line);
    std::istringstream line_stream(line);
    std::string truncation_rank_code, coupling_code;
    line_stream >> truncation_rank_code
                >> parameters.truncation_cutoff
                >> coupling_code;
    ParsingCheck(line_stream,line_count,line);

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: header
  ////////////////////////////////////////////////////////////////


  template<H2Format tFormat, H2Mode tMode>
  void WriteHeader(std::ofstream stream, OutH2Stream& h2_stream);

  template <>
  void WriteHeader<H2Format::kVersion0,H2Mode::kText>(std::ofstream stream,OutH2Stream& h2_stream)
  {
    // extract parameters
    //
    // assumes oscillator-like weight cutoffs
    // TODO: add tests for oscillator-like orbitals and weight cutoffs
    const int num_types = 2;
    int N1max = int(h2_stream.space().weight_max().one_body[0]);
    int N2max = int(h2_stream.space().weight_max().two_body[0]);
    int size_pp_nn = h2_stream.size_by_type()[0];
    int size_pn = h2_stream.size_by_type()[2];

    stream
      // header line 1: number of particle species
      << fmt::format("{:10d}",num_types) << std::endl
      // header line 2: 1-body basis limit
      << fmt::format("{:10d}",N1max) << std::endl
      // header line 3: 2-body basis limit
      << fmt::format("{:10d}",N2max) << std::endl
      // header line 4: matrix size
      << fmt::format("{:10d} {:10d}",size_pp_nn,size_pn) << std::endl;
    
  };

  template <>
  void WriteHeader<H2Format::kVersion0,H2Mode::kBinary>(std::ofstream stream,OutH2Stream& h2_stream)
  {
    // extract parameters
    //
    // assumes oscillator-like weight cutoffs
    // TODO: add tests for oscillator-like orbitals and weight cutoffs
    const int num_types = 2;
    int N1max = int(h2_stream.space().weight_max().one_body[0]);
    int N2max = int(h2_stream.space().weight_max().two_body[0]);
    int size_pp_nn = h2_stream.size_by_type()[0];
    int size_pn = h2_stream.size_by_type()[2];

    // record size
    int num_header_fields = 5;
    int header_bytes = num_header_fields * kIntegerSize;

    // write
    WriteI4(stream,header_bytes);
    WriteI4(stream,num_types);
    WriteI4(stream,N1max);
    WriteI4(stream,N2max);
    WriteI4(stream,size_pp_nn);
    WriteI4(stream,size_pn);
    WriteI4(stream,header_bytes);
    
  };





  // set up target
  // run_parameters.output_h2_format = 15099;
  run_parameters.truncation_rank = basis::Rank::kTwoBody;
  run_parameters.truncation_cutoff = 2;
  run_parameters.weight_max = basis::WeightMax(run_parameters.truncation_rank,run_parameters.truncation_cutoff);
  // target_indexing.orbital_space = basis::OrbitalSpacePN(0);
  // target_indexing.space = basis::TwoBodySpaceJJJPN(
  //     target_indexing.orbital_space,basis::WeightMax(basis::Rank::kTwoBody,0)
  //   );
  // target_indexing.sectors = basis::TwoBodySectorsJJJPN(target_indexing.space,0,0,0);

  // set up inputs
  input_channels.emplace_back("in1","test/test-tb-2.dat");

  // set up outputs
  target_channels.emplace_back("out1","test/out1.dat");
  target_channels[0].coefficients.insert({"in1",999});
  target_channels.emplace_back("out2","test/out2.dat");
