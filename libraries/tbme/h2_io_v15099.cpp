////////////////////////////////////////////////////////////////
// format/mode-specific I/O: v15099
////////////////////////////////////////////////////////////////

#include "tbme/h2_io.h"  // include for IDE tools

#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  void InH2Stream::ReadHeader_Version15099()
  {

    // orbital info storage
    std::vector<basis::OrbitalPNInfo> orbitals;

    // header parameters
    int J0, g0, Tz0;
    float wp, wn, wpp, wnn, wpn;
    int twice_Jmax_pp, twice_Jmax_nn, twice_Jmax_pn;
    int size_pp, size_nn, size_pn;

    if (h2_mode()==H2Mode::kText)
      {

        // set up line input
        std::string line;

        // orbital listing header
        std::size_t num_orbitals_p, num_orbitals_n;
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> num_orbitals_p >> num_orbitals_n;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // orbital listing body
        std::size_t num_orbitals = num_orbitals_p + num_orbitals_n;
        std::string orbital_info_str;
        for (std::size_t orbital_line_count=0; orbital_line_count < num_orbitals; ++orbital_line_count)
          {
            mcutils::GetLine(stream(), line, line_count_);
            orbital_info_str.append(line);
            orbital_info_str.append("\n");  // need to restore newline to input line
          }
        std::istringstream orbital_info_stream(orbital_info_str);
        orbitals = basis::ParseOrbitalPNStream(orbital_info_stream,/*standalone=*/false);

        // header line 1: operator properties
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> J0 >> g0 >> Tz0;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 2: 1-body basis limit
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> wp >> wn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 3: 2-body basis limit
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> wpp >> wnn >> wpn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 4: 2-body basis a.m. limit
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> twice_Jmax_pp >> twice_Jmax_nn >> twice_Jmax_pn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 5: matrix size
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> size_pp >> size_nn >> size_pn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }
      }
    else if (h2_mode()==H2Mode::kBinary)
      {

        int num_fields, bytes;

        // orbital listing header
        num_fields=2;
        bytes = num_fields * kIntegerSize;
        int num_orbitals_p, num_orbitals_n;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),num_orbitals_p);
        mcutils::ReadBinary<int>(stream(),num_orbitals_n);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

        // orbital listing body
        /* int num_orbitals = num_orbitals_p + num_orbitals_n; */
        for (basis::OrbitalSpeciesPN orbital_species : {basis::OrbitalSpeciesPN::kP,basis::OrbitalSpeciesPN::kN})
          {
            // select number of orbitals
            std::size_t num_orbitals;
            if (orbital_species==basis::OrbitalSpeciesPN::kP)
              num_orbitals = num_orbitals_p;
            else
              num_orbitals = num_orbitals_n;

            // set up vector storage for quantum numbers
            std::vector<int> orbitals_n;
            std::vector<int> orbitals_l;
            std::vector<int> orbitals_twice_j;
            std::vector<float> orbitals_weight;

            // do vector read for each orbital property
            // read n
            bytes = num_orbitals * kIntegerSize;
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              {
                int n;
                mcutils::ReadBinary<int>(stream(),n);
                orbitals_n.push_back(n);
              }
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            // read l
            bytes = num_orbitals * kIntegerSize;
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              {
                int l;
                mcutils::ReadBinary<int>(stream(),l);
                orbitals_l.push_back(l);
              }
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            // read twice_j
            bytes = num_orbitals * kIntegerSize;
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              {
                int twice_j;
                mcutils::ReadBinary<int>(stream(),twice_j);
                orbitals_twice_j.push_back(twice_j);
              }
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            // read weight
            bytes = num_orbitals * kFloatSize;
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              {
                float weight;
                mcutils::ReadBinary<float>(stream(),weight);
                orbitals_weight.push_back(weight);
              }
            mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

            // store orbital info
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              {
                int n = orbitals_n[orbital_index];
                int l = orbitals_l[orbital_index];
                int twice_j = orbitals_twice_j[orbital_index];
                double weight = orbitals_weight[orbital_index];
                HalfInt j = HalfInt(twice_j,2);

                orbitals.push_back(basis::OrbitalPNInfo(orbital_species,n,l,j,weight));
             }
          }

        // two-body indexing

        // header line 1: operator properties
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),J0);
        mcutils::ReadBinary<int>(stream(),g0);
        mcutils::ReadBinary<int>(stream(),Tz0);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

        // header line 2: 1-body basis limit
        num_fields=2;
        bytes = num_fields * kFloatSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<float>(stream(),wp);
        mcutils::ReadBinary<float>(stream(),wn);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

        // header line 3: 2-body basis limit
        num_fields=3;
        bytes = num_fields * kFloatSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<float>(stream(),wpp);
        mcutils::ReadBinary<float>(stream(),wnn);
        mcutils::ReadBinary<float>(stream(),wpn);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

        // header line 4: 2-body basis a.m. limit
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),twice_Jmax_pp);
        mcutils::ReadBinary<int>(stream(),twice_Jmax_nn);
        mcutils::ReadBinary<int>(stream(),twice_Jmax_pn);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");

        // header line 5: matrix size
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),size_pp);
        mcutils::ReadBinary<int>(stream(),size_nn);
        mcutils::ReadBinary<int>(stream(),size_pn);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
      }

    // set up indexing
    orbital_space_ = basis::OrbitalSpacePN(orbitals);
    // std::cout << orbital_space_.DebugStr() << std::endl;
    space_ = basis::TwoBodySpaceJJJPN(
        orbital_space_,
        basis::WeightMax(wp,wn,wpp,wnn,wpn),
        basis::TwoBodySpaceJJJPNOrdering::kPN
      );
    basis::SectorDirection sector_direction;
    if (Tz0==0)
      {
        // impose upper triangularity on sector selection rule for Tz0=0
        sector_direction = basis::SectorDirection::kCanonical;
      }
    else
      {
        // but need all possible sectors when Tz0!=0
        std::cout << "ERROR: nonzero Tz0 not allowed for version 15099" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    sectors_ = basis::TwoBodySectorsJJJPN(space_,J0,g0,Tz0,sector_direction);

    // reserve vector storage for stream positions
    stream_sector_positions_.reserve(sectors().size());
    // set location of first sector to current stream position
    stream_sector_positions_.push_back(stream().tellg());

    // deduce matrix limits by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // validate unused (independently-recalculated) input parameters
    assert(size_pp==size_by_type()[0]);
    assert(size_nn==size_by_type()[1]);
    assert(size_pn==size_by_type()[2]);
    int Jmax_pp = twice_Jmax_pp / 2;
    int Jmax_nn = twice_Jmax_nn / 2;
    int Jmax_pn = twice_Jmax_pn / 2;
    assert(Jmax_pp==Jmax_by_type()[0]);
    assert(Jmax_nn==Jmax_by_type()[1]);
    assert(Jmax_pn==Jmax_by_type()[2]);
  }

  void OutH2Stream::WriteHeader_Version15099()
  {
    // require that subspaces are in pp/nn/pn order
    assert(space().space_ordering()==basis::TwoBodySpaceJJJPNOrdering::kPN);

    // extract parameters
    if (h2_mode()==H2Mode::kText)
      {
        // dump orbitals
        stream()
          // header line: dimensions
          << fmt::format(
              "{} {}",
              orbital_space().GetSubspace(0).size(),
              orbital_space().GetSubspace(1).size()
            )
          << std::endl
          // orbital listing body
          << basis::OrbitalDefinitionStr(orbital_space().OrbitalInfo(),false);

        // write two-body header info
        stream()
          // header line 1: operator properties
          << fmt::format("{} {} {}",sectors().J0(),sectors().g0(),sectors().Tz0()) << std::endl
          // header line 2: 1-body basis limit
          << fmt::format(
              "{:e} {:e}",
              space().weight_max().one_body[0],space().weight_max().one_body[1]
            )
          << std::endl
          // header line 3: 2-body basis limit
          << fmt::format(
              "{:e} {:e} {:e}",
              space().weight_max().two_body[0],space().weight_max().two_body[1],space().weight_max().two_body[2]
            )
          << std::endl
          // header line 4: 2-body basis a.m. limit
          << fmt::format("{:d} {:d} {:d}",2*Jmax_by_type()[0],2*Jmax_by_type()[1],2*Jmax_by_type()[2]) << std::endl
          // header line 5: matrix size
          << fmt::format("{:d} {:d} {:d}",size_by_type()[0],size_by_type()[1],size_by_type()[2]) << std::endl;
      }
    else if (h2_mode()==H2Mode::kBinary)
      {

        int num_fields, bytes;

        // dump orbitals

        // header: dimensions
        mcutils::WriteBinary<int>(stream(),2*kIntegerSize);
        for (int subspace_index=0; subspace_index < orbital_space().size(); ++subspace_index)
          mcutils::WriteBinary<int>(stream(),orbital_space().GetSubspace(subspace_index).size());
        mcutils::WriteBinary<int>(stream(),2*kIntegerSize);

        // orbital listing body
        for (std::size_t subspace_index=0; subspace_index < orbital_space().size(); ++subspace_index)
          {
            const std::vector<basis::OrbitalPNInfo> orbitals = orbital_space().GetSubspace(subspace_index).OrbitalInfo();
            const std::size_t num_orbitals = orbitals.size();

            // do vector write for each orbital property
            // write n
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),orbitals[orbital_index].n);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            // write l
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),orbitals[orbital_index].l);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            // write twice_j
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),TwiceValue(orbitals[orbital_index].j));
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            // write weight
            mcutils::WriteBinary<int>(stream(),num_orbitals*kFloatSize);
            for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<float>(stream(),orbitals[orbital_index].weight);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kFloatSize);
          }

        // two-body indexing

        // header line 1: operator properties
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),sectors().J0());
        mcutils::WriteBinary<int>(stream(),sectors().g0());
        mcutils::WriteBinary<int>(stream(),sectors().Tz0());
        mcutils::WriteBinary<int>(stream(),bytes);

        // header line 2: 1-body basis limit
        num_fields=2;
        bytes = num_fields * kFloatSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<float>(stream(),space().weight_max().one_body[0]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().one_body[1]);
        mcutils::WriteBinary<int>(stream(),bytes);

        // header line 3: 2-body basis limit
        num_fields=3;
        bytes = num_fields * kFloatSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[0]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[1]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[2]);
        mcutils::WriteBinary<int>(stream(),bytes);

        // header line 4: 2-body basis a.m. limit
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[0]);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[1]);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[2]);
        mcutils::WriteBinary<int>(stream(),bytes);

        // header line 5: matrix size
        num_fields=3;
        bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),size_by_type()[0]);
        mcutils::WriteBinary<int>(stream(),size_by_type()[1]);
        mcutils::WriteBinary<int>(stream(),size_by_type()[2]);
        mcutils::WriteBinary<int>(stream(),bytes);
      }
  };


  void InH2Stream::ReadSector_Version15099(Eigen::MatrixXd& matrix)
  {
    // extract sector
    const auto& sector = sectors().GetSector(sector_index_);
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();

    // verify that sector is canonical
    //
    // This is a check that the caller's sector construction
    // followed the specification that only "upper triangle"
    // sectors are stored.
    assert(sector.IsUpperTriangle());

    // allocate matrix
    matrix = Eigen::MatrixXd::Zero(bra_subspace.size(),ket_subspace.size());

    // read FORTRAN record beginning delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsFirstOfType())
      {
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }

    // iterate over matrix elements
    for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
      for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
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
              int input_i1, input_i2, input_i3, input_i4, input_twice_J_bra, input_twice_J_ket;
              int input_two_body_species_code;

              // read line
              std::string line;
              mcutils::GetLine(stream(),line,line_count_);
              std::istringstream line_stream(line);
              line_stream
                >> input_i1 >> input_i2 >> input_i3 >> input_i4
                >> input_twice_J_bra
                >> input_twice_J_ket
                >> input_two_body_species_code
                >> input_matrix_element;
              mcutils::ParsingCheck(line_stream,line_count_,line);

              // validate input fields against expected values
              bool inputs_as_expected = (
                  (input_i1==bra.index1()+1) && (input_i2==bra.index2()+1)
                  && (input_i3==ket.index1()+1) && (input_i4==ket.index2()+1)
                  && (input_twice_J_bra == 2*bra.J())
                  && (input_twice_J_ket == 2*ket.J())
                  && (input_two_body_species_code==basis::kTwoBodySpeciesPNCodeDecimal[int(ket.two_body_species())])
                );
              if (!inputs_as_expected)
                mcutils::ParsingError(line_count_,line,"Unexpected matrix element labels in input data");
            }
          else if (h2_mode()==H2Mode::kBinary)
            {
              mcutils::ReadBinary<float>(stream(),input_matrix_element);
            }

          matrix(bra_index,ket_index) = input_matrix_element;
        }

    // read FORTRAN record ending delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsLastOfType())
      {
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }
  }

  void InH2Stream::SkipSector_Version15099()
  {
    // extract sector
    const auto& sector = sectors().GetSector(sector_index_);
    const auto& bra_subspace = sector.bra_subspace();
    const auto& ket_subspace = sector.ket_subspace();

    // verify that sector is canonical
    //
    // This is a check that the caller's sector construction
    // followed the specification that only "upper triangle"
    // sectors are stored.
    assert(sector.IsUpperTriangle());

    // read FORTRAN record beginning delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsFirstOfType())
      {
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }

    // calculate number of matrix elements in sector
    std::size_t sector_entries = 0;
    if (sector.IsDiagonal())
      // diagonal sector
      {
        std::size_t dimension = ket_subspace.size();
        sector_entries = dimension*(dimension+1)/2;
      }
    else  // if (sector.IsUpperTriangle())
      // upper triangle sector (but not diagonal)
      {
        std::size_t bra_dimension = bra_subspace.size();
        std::size_t ket_dimension = ket_subspace.size();
        sector_entries = bra_dimension*ket_dimension;
      }

    // skip matrix elements
    if (h2_mode()==H2Mode::kText)
      {
        for (std::size_t index=0; index<sector_entries; ++index)
          // skip sector_entries lines
          {
            std::string line;
            mcutils::GetLine(stream(), line, line_count_);
          }
      }
    else if (h2_mode()==H2Mode::kBinary)
      // simply move file pointer ahead
      {
        stream().seekg(sector_entries*kIntegerSize, std::ios_base::cur);
      }

    // read FORTRAN record ending delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsLastOfType())
      {
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::VerifyBinary<int>(
            stream(),entries*kIntegerSize,
            "Encountered unexpected value in H2 file","record delimiter"
          );
      }
  }

  void OutH2Stream::WriteSector_Version15099(
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
  {
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

    // write FORTRAN record beginning delimiter
    if ((h2_mode()==H2Mode::kBinary) && SectorIsFirstOfType())
      {
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        assert(static_cast<long int>(entries*kIntegerSize) < kMaxRecordLength);
        mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      }

    // iterate over matrix elements
    for (std::size_t bra_index=0; bra_index<bra_subspace.size(); ++bra_index)
      for (std::size_t ket_index=0; ket_index<ket_subspace.size(); ++ket_index)
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
          int output_twice_J_bra = 2*bra.J();
          int output_twice_J_ket = 2*ket.J();
          int output_two_body_species_code=basis::kTwoBodySpeciesPNCodeDecimal[int(ket.two_body_species())];
          if (h2_mode()==H2Mode::kText)
            {
              stream()
                << fmt::format(
                    "{:3d} {:3d} {:3d} {:3d} {:3d} {:3d} {:2d} {:+16.7e}",
                    output_i1,output_i2,output_i3,output_i4,
                    output_twice_J_bra,
                    output_twice_J_ket,
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
        std::size_t entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      }
  }

}

