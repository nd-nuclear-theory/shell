////////////////////////////////////////////////////////////////
// format/mode-specific I/O: v15200
////////////////////////////////////////////////////////////////

#include <cstdint>

#include "tbme/h2_io.h"  // include for IDE tools

#include "mcutils/parsing.h"
#include "mcutils/io.h"

namespace shell {

  void InH2Stream::ReadHeader_Version15200()
  {

    // orbital info storage
    basis::OrbitalPNList orbitals;

    // header parameters
    int J0, g0, Tz0;
    float wp, wn, wpp, wnn, wpn;
    long int size_pp, size_nn, size_pn;

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
        orbitals = basis::ParseOrbitalPNStream(
            orbital_info_stream,
            /*standalone=*/false,
            basis::MFDnOrbitalFormat::kVersion15200
          );

        // header line 1: operator properties
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> J0 >> g0 >> Tz0;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 2: 2-body basis limit
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> wpp >> wpn >> wnn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }

        // header line 3: matrix size
        {
          mcutils::GetLine(stream(), line, line_count_);
          std::istringstream line_stream(line);
          line_stream >> size_pp >> size_pn >> size_nn;
          mcutils::ParsingCheck(line_stream,line_count_,line);
        }
      }
    else if (h2_mode()==H2Mode::kBinary)
      {

        int num_fields, bytes;

        // orbital listing header
        int num_orbitals_p, num_orbitals_n;
        mcutils::ReadBinary<int32_t>(stream(),num_orbitals_p);
        mcutils::ReadBinary<int32_t>(stream(),num_orbitals_n);

        // orbital listing body
        std::size_t num_orbitals = num_orbitals_p + num_orbitals_n;

        // set up vector storage for quantum numbers
        std::vector<int32_t> orbitals_n;
        std::vector<int32_t> orbitals_l;
        std::vector<int32_t> orbitals_twice_j;
        std::vector<int32_t> orbitals_twice_tz;
        std::vector<float> orbitals_weight;

        // do vector read for each orbital property
        // read n
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            int n;
            mcutils::ReadBinary<int32_t>(stream(),n);
            orbitals_n.push_back(n);
          }
        // read l
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            int l;
            mcutils::ReadBinary<int32_t>(stream(),l);
            orbitals_l.push_back(l);
          }
        // read twice_j
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            int twice_j;
            mcutils::ReadBinary<int32_t>(stream(),twice_j);
            orbitals_twice_j.push_back(twice_j);
          }
        // read twice_tz
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            int twice_tz;
            mcutils::ReadBinary<int32_t>(stream(),twice_tz);
            orbitals_twice_tz.push_back(twice_tz);
          }
        // read weight
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            float weight;
            mcutils::ReadBinary<float>(stream(),weight);
            orbitals_weight.push_back(weight);
          }

        // store orbital info
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          {
            int n = orbitals_n[orbital_index];
            int l = orbitals_l[orbital_index];
            HalfInt j = HalfInt(orbitals_twice_j[orbital_index],2);
            HalfInt tz = HalfInt(orbitals_twice_tz[orbital_index],2);
            double weight = orbitals_weight[orbital_index];

            orbitals.push_back(basis::OrbitalPNInfo(n,l,j,tz,weight));
          }

        // two-body indexing

        // header line 1: operator properties
        mcutils::ReadBinary<int32_t>(stream(),J0);
        mcutils::ReadBinary<int32_t>(stream(),g0);
        mcutils::ReadBinary<int32_t>(stream(),Tz0);

        // header line 2: 2-body basis limit
        mcutils::ReadBinary<float>(stream(),wpp);
        mcutils::ReadBinary<float>(stream(),wpn);
        mcutils::ReadBinary<float>(stream(),wnn);

        // header line 3: matrix size
        mcutils::ReadBinary<int64_t>(stream(),size_pp);
        mcutils::ReadBinary<int64_t>(stream(),size_pn);
        mcutils::ReadBinary<int64_t>(stream(),size_nn);
      }

    // set up indexing
    orbital_space_ = basis::OrbitalSpacePN(orbitals);
    // std::cout << orbital_space_.DebugStr() << std::endl;
    space_ = basis::TwoBodySpaceJJJPN(
        orbital_space_,
        basis::WeightMax(wp,wn,wpp,wnn,wpn),
        basis::TwoBodySpaceJJJPNOrdering::kTz
      );
    sectors_ = basis::TwoBodySectorsJJJPN(space_,J0,g0,Tz0);

    // reserve vector storage for stream positions
    stream_sector_positions_.reserve(sectors().size());
    // set location of first sector to current stream position
    stream_sector_positions_.push_back(stream().tellg());

    // deduce matrix limits by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // validate unused (independently-recalculated) input parameters
    assert(size_pp==size_by_type()[0]);  // TODO change if TwoBodySpeciesPN enum changes
    assert(size_pn==size_by_type()[2]);  // TODO change if TwoBodySpeciesPN enum changes
    assert(size_nn==size_by_type()[1]);  // TODO change if TwoBodySpeciesPN enum changes
  }

  void OutH2Stream::WriteHeader_Version15200()
  {
    // require that subspaces are in pp/pn/nn order
    assert(space().space_ordering()==basis::TwoBodySpaceJJJPNOrdering::kTz);

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
          << basis::OrbitalDefinitionStr(
              orbital_space().OrbitalInfo(),
              /*standalone=*/false,
              basis::MFDnOrbitalFormat::kVersion15200
            );

        // write two-body header info
        stream()
          // header line 1: operator properties
          << fmt::format("{} {} {}",sectors().J0(),sectors().g0(),sectors().Tz0()) << std::endl
          // header line 2: 2-body basis limit
          << fmt::format(
              "{:e} {:e} {:e}",
              space().weight_max().two_body[0],space().weight_max().two_body[2],space().weight_max().two_body[1]  // TODO(pjf) change if TwoBodySpeciesPN enum changes
            )
          << std::endl
          // header line 3: matrix size
          << fmt::format("{:d} {:d} {:d}",size_by_type()[0],size_by_type()[2],size_by_type()[1]) << std::endl;  // TODO(pjf) change if TwoBodySpeciesPN enum changes
      }
    else if (h2_mode()==H2Mode::kBinary)
      {

        int num_fields, bytes;

        // dump orbitals

        // header: dimensions
        for (std::size_t subspace_index=0; subspace_index < orbital_space().size(); ++subspace_index)
          mcutils::WriteBinary<int32_t>(stream(),orbital_space().GetSubspace(subspace_index).size());

        // orbital listing body
        const std::vector<basis::OrbitalPNInfo> orbitals = orbital_space().OrbitalInfo();
        const std::size_t num_orbitals = orbitals.size();

        // do vector write for each orbital property
        // write n
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          mcutils::WriteBinary<int32_t>(stream(),orbitals[orbital_index].n);
        // write l
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          mcutils::WriteBinary<int32_t>(stream(),orbitals[orbital_index].l);
        // write twice_j
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          mcutils::WriteBinary<int32_t>(stream(),TwiceValue(orbitals[orbital_index].j));
        // write twice_tz
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          mcutils::WriteBinary<int32_t>(
              stream(),
              TwiceValue(basis::kOrbitalSpeciesPNCodeTz[int(orbitals[orbital_index].orbital_species)])
            );
        // write weight
        for (std::size_t orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
          mcutils::WriteBinary<float>(stream(),orbitals[orbital_index].weight);

        // two-body indexing

        // header line 1: operator properties
        mcutils::WriteBinary<int32_t>(stream(),sectors().J0());
        mcutils::WriteBinary<int32_t>(stream(),sectors().g0());
        mcutils::WriteBinary<int32_t>(stream(),sectors().Tz0());

        // header line 2: 2-body basis limit
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[0]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[2]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[1]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes

        // header line 3: matrix size
        mcutils::WriteBinary<int64_t>(stream(),size_by_type()[0]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes
        mcutils::WriteBinary<int64_t>(stream(),size_by_type()[2]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes
        mcutils::WriteBinary<int64_t>(stream(),size_by_type()[1]);  // TODO(pjf) change if TwoBodySpeciesPN enum changes
      }
  };


  void InH2Stream::ReadSector_Version15200(Eigen::MatrixXd& matrix)
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

    // get offsets to convert index within orbital subspace to index in all orbitals
    std::size_t num_orbitals_p = orbital_space().GetSubspace(0).size();
    int bra_index1_offset = ((bra_subspace.orbital_subspace1().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int bra_index2_offset = ((bra_subspace.orbital_subspace2().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int ket_index1_offset = ((ket_subspace.orbital_subspace1().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int ket_index2_offset = ((ket_subspace.orbital_subspace2().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);

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

              // read line
              std::string line;
              mcutils::GetLine(stream(),line,line_count_);
              std::istringstream line_stream(line);
              line_stream
                >> input_i1 >> input_i2 >> input_i3 >> input_i4
                >> input_twice_J_bra
                >> input_twice_J_ket
                >> input_matrix_element;
              mcutils::ParsingCheck(line_stream,line_count_,line);

              // validate input fields against expected values
              bool inputs_as_expected = (
                  (input_i1==bra.index1()+bra_index1_offset)
                  && (input_i2==bra.index2()+bra_index2_offset)
                  && (input_i3==ket.index1()+ket_index1_offset)
                  && (input_i4==ket.index2()+ket_index2_offset)
                  && (input_twice_J_bra == 2*bra.J())
                  && (input_twice_J_ket == 2*ket.J())
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
  }

  void InH2Stream::SkipSector_Version15200()
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

    // calculate number of matrix elements in sector
    std::size_t sector_entries = GetSectorCount(sector);
    // skip matrix elements
    if (h2_mode()==H2Mode::kText)
      {
        for (std::size_t indx=0; indx<sector_entries; ++indx)
          // skip sector_entries lines
          {
            std::string line;
            mcutils::GetLine(stream(), line, line_count_);
          }
      }
    else if (h2_mode()==H2Mode::kBinary)
      // simply move file pointer ahead
      {
        stream().seekg(sector_entries*kFloatSize, std::ios_base::cur);
      }
  }

  void OutH2Stream::WriteSector_Version15200(
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
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

    // get offsets to convert index within orbital subspace to index in all orbitals
    std::size_t num_orbitals_p = orbital_space().GetSubspace(0).size();
    int bra_index1_offset = ((bra_subspace.orbital_subspace1().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int bra_index2_offset = ((bra_subspace.orbital_subspace2().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int ket_index1_offset = ((ket_subspace.orbital_subspace1().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);
    int ket_index2_offset = ((ket_subspace.orbital_subspace2().orbital_species()==basis::OrbitalSpeciesPN::kP) ? 1 : num_orbitals_p+1);

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
          int output_i1 = bra.index1()+bra_index1_offset;
          int output_i2 = bra.index2()+bra_index2_offset;
          int output_i3 = ket.index1()+ket_index1_offset;
          int output_i4 = ket.index2()+ket_index2_offset;
          int output_twice_J_bra = 2*bra.J();
          int output_twice_J_ket = 2*ket.J();
          if (h2_mode()==H2Mode::kText)
            {
              stream()
                << fmt::format(
                    "{:3d} {:3d} {:3d} {:3d} {:3d} {:3d} {:+16.7e}",
                    output_i1,output_i2,output_i3,output_i4,
                    output_twice_J_bra,
                    output_twice_J_ket,
                    output_matrix_element
                  )
                << std::endl;
            }
          else if (h2_mode()==H2Mode::kBinary)
            {
              mcutils::WriteBinary<float>(stream(),output_matrix_element);
            }
        }
  }
}  // namespace shell
