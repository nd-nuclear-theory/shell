////////////////////////////////////////////////////////////////
// format/mode-specific I/O: v0
////////////////////////////////////////////////////////////////

namespace shell {

  void InH2Stream::ReadHeader_Version0() 
  {
    // header parameters
    int num_types, N1max, N2max, size_pp_nn, size_pn;

    if (h2_mode()==H2Mode::kText)
      {
        // set up line input
        std::string line;

        // header line 1: number of particle types
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> num_types; 
          ParsingCheck(line_stream,line_count_,line);
        }

        // header line 2: 1-body basis limit
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> N1max; 
          ParsingCheck(line_stream,line_count_,line);
        }

        // header line 3: 2-body basis limit
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> N2max; 
          ParsingCheck(line_stream,line_count_,line);
        }
		
        // header line 4: matrix size
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> size_pp_nn >> size_pn; 
          ParsingCheck(line_stream,line_count_,line);
        }
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 5;
        int bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),num_types);
        mcutils::ReadBinary<int>(stream(),N1max);
        mcutils::ReadBinary<int>(stream(),N2max);
        mcutils::ReadBinary<int>(stream(),size_pp_nn);
        mcutils::ReadBinary<int>(stream(),size_pn);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
      }

    // set up indexing
    int J0 = 0;
    int g0 = 0;
    int Tz0 = 0;
    orbital_space_ = basis::OrbitalSpacePN(N1max);
    space_ = basis::TwoBodySpaceJJJPN(orbital_space_,basis::WeightMax(N1max,N2max));
    sectors_ = basis::TwoBodySectorsJJJPN(space_,J0,g0,Tz0);

    // deduce matrix limits by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // validate unused (independently-recalculated) input parameters
    assert(num_types == 2);
    assert(size_pp_nn==size_by_type()[0]);
    assert(size_pn==size_by_type()[2]);
  }

  void OutH2Stream::WriteHeader_Version0()
  {
    // extract parameters
    //
    // assumes oscillator-like weight cutoffs
    // TODO: add tests for oscillator-like orbitals and weight cutoffs
    // assert(orbital_space().is_oscillator_like());
    const int num_types = 2;
    int N1max = int(space().weight_max().one_body[0]);
    int N2max = int(space().weight_max().two_body[0]);
    int size_pp_nn = size_by_type()[0];
    int size_pn = size_by_type()[2];

    if (h2_mode()==H2Mode::kText)
      {
        stream()
          // header line 1: number of particle species
          << fmt::format("{:10d}",num_types) << std::endl
          // header line 2: 1-body basis limit
          << fmt::format("{:10d}",N1max) << std::endl
          // header line 3: 2-body basis limit
          << fmt::format("{:10d}",N2max) << std::endl
          // header line 4: matrix size
          << fmt::format("{:10d} {:10d}",size_pp_nn,size_pn) << std::endl;
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 5;
        int bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),num_types);
        mcutils::WriteBinary<int>(stream(),N1max);
        mcutils::WriteBinary<int>(stream(),N2max);
        mcutils::WriteBinary<int>(stream(),size_pp_nn);
        mcutils::WriteBinary<int>(stream(),size_pn);
        mcutils::WriteBinary<int>(stream(),bytes);
      }
  };

  void InH2Stream::ReadSector_Version0(Eigen::MatrixXd& matrix, bool store)
  {
    // Note: Although H2 format Version0 only suports *diagonal* sectors
    // (operator conserves Tz, J, g), we treat the bra and ket
    // subspaces as distinct below to maintain direct analogy of code
    // with that for later H2 formats.  We only modify this by
    // asserting that the sector is diagonal (and using the ket
    // subspace labels as "the" labels).

    // extract sector
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

    // std::cout << fmt::format(
    //     "sector {} type {} first {} last {}",
    //     sector_index_,basis::kTwoBodySpeciesPNCodeDecimal[int(ket_subspace.two_body_species())],
    //     SectorIsFirstOfType(),SectorIsLastOfType()
    //   )
    //           << std::endl;

    // verify that sector is canonical
    //
    // This is a check that the caller's sector construction
    // followed the specification that only "upper triangle"
    // sectors are stored.
    assert(sector.IsUpperTriangle());

    // allocate matrix
    if (store)
      matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

    // Version0 restrictions -- sector is actually diagonal (in two_body_species, J, g)
    assert(sector.IsDiagonal());
    // const int two_body_species_code = basis::kTwoBodySpeciesPNCodeDecimal[int(ket_subspace.two_body_species())];
    // const int J = sector.ket_subspace().J();
    // const int g = sector.ket_subspace().g();

    // read FORTRAN record beginning delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsFirstOfType())
      {
        int entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }

    // iterate over matrix elements
    for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
      for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
        {

          // diagonal sector: restrict to upper triangle
          if (sector.IsDiagonal())
            if (!(bra_index<=ket_index))
              continue;

          // retrieve states
          const basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
          const basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

          // read input matrix element
	  float input_matrix_element;
          if (h2_mode()==H2Mode::kText)
            {
              // input fields for text mode only
              int input_i1, input_i2, input_i3, input_i4, input_twice_J;
              int input_two_body_species_code;

              // read line
              std::string line;
              ++line_count_;
              std::getline(stream(),line);
              std::istringstream line_stream(line);
              line_stream
                >> input_i1 >> input_i2 >> input_i3 >> input_i4
                >> input_twice_J
                >> input_two_body_species_code
                >> input_matrix_element;
              ParsingCheck(line_stream,line_count_,line);

              // validate input fields against expected values
              bool inputs_as_expected = ( 
                  (input_i1==bra.index1()+1) && (input_i2==bra.index2()+1)
                  && (input_i3==ket.index1()+1) && (input_i4==ket.index2()+1)
                  && (input_twice_J == 2*ket.J())
                  && (input_two_body_species_code==basis::kTwoBodySpeciesPNCodeDecimal[int(ket.two_body_species())])
                );
              if (!inputs_as_expected)
                ParsingError(line_count_,line,"Unexpected matrix element labels in input data");
            }
          else if (h2_mode()==H2Mode::kBinary)
            {
              mcutils::ReadBinary<float>(stream(),input_matrix_element);
            }

          if (store)
            matrix(bra_index,ket_index) = input_matrix_element;
        }

    // read FORTRAN record ending delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsLastOfType())
      {
        int entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }
  }

  void OutH2Stream::WriteSector_Version0(
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
  {
    // Note: Although H2 format Version0 only suports *diagonal* sectors
    // (operator conserves Tz, J, g), we treat the bra and ket
    // subspaces as distinct below to maintain direct analogy of code
    // with that for later H2 formats.  We only modify this by
    // asserting that the sector is diagonal (and using the ket
    // subspace labels as "the" labels).

    // extract sector
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

    // verify that sector is canonical
    //
    // This is a check that the caller's sector construction
    // followed the specification that only "upper triangle"
    // sectors are stored.
    assert(sector.IsUpperTriangle());

    // Version0 restrictions -- sector is actually diagonal (in two_body_species, J, g)
    assert(sector.IsDiagonal());

    // write FORTRAN record beginning delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsFirstOfType())
      {
        int entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      }

    // iterate over matrix elements
    for (int bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
      for (int ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
        {

          // diagonal sector: restrict to upper triangle
          if (sector.IsDiagonal())
            if (!(bra_index<=ket_index))
              continue;

          // retrieve states
          const basis::TwoBodyStateJJJPN bra(bra_subspace,bra_index);
          const basis::TwoBodyStateJJJPN ket(ket_subspace,ket_index);

          // determine matrix element normalization factor
          double conversion_factor = 1.;
          if (conversion_mode == basis::NormalizationConversion::kASToNAS)
            {
              if (bra.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                if (bra.index1()==bra.index2())
                  conversion_factor *= (1/sqrt(2.));
              if (ket.two_body_species()!=basis::TwoBodySpeciesPN::kPN)
                if (ket.index1()==ket.index2())
                  conversion_factor *= (1/sqrt(2.));
            }

          // retrieve matrix element for output
          float output_matrix_element = conversion_factor * matrix(bra_index,ket_index);

          // write line: output matrix element
          int output_i1 = bra.index1()+1;
          int output_i2 = bra.index2()+1;
          int output_i3 = ket.index1()+1;
          int output_i4 = ket.index2()+1;
          int output_twice_J = 2*ket.J();
          int output_two_body_species_code=basis::kTwoBodySpeciesPNCodeDecimal[int(ket.two_body_species())];
          if (h2_mode()==H2Mode::kText)
            {
              stream()
                << fmt::format(
                    "{:3d} {:3d} {:3d} {:3d} {:3d} {:2d} {:+16.7e}",
                    output_i1,output_i2,output_i3,output_i4,
                    output_twice_J,
                    output_two_body_species_code,
                    output_matrix_element
                  )
                << std::endl;
            }
          else if (h2_mode()==H2Mode::kBinary)
            {
              mcutils::WriteBinary<float>(stream(),output_matrix_element);
            }
        }
			
    // write FORTRAN record ending delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsLastOfType())
      {
        int entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      }
  }

}
