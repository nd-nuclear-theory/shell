/****************************************************************
  h2_io.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "tbme/h2_io.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "cppformat/format.h"
#include "mcutils/io.h"
#include "mcutils/parsing.h"

namespace shell {

  ////////////////////////////////////////////////////////////////
  // pp/nn/pn matrix size counting
  ////////////////////////////////////////////////////////////////

  void EvaluateJmaxByType(
      const basis::TwoBodySpaceJJJPN& space,
      std::array<int,3>& Jmax_by_type
    )
  // Extract Jmax value by two-body species type for JJJPN space.
  //
  // Arguments:
  //   space (..., input) : space
  //   Jmax_by_type (std::array<int,3>, output) : maximum J values
  {

    // initialize
    Jmax_by_type = std::array<int,3>({0,0,0});

    // scan subspaces for Jmax
    for (int subspace_index=0; subspace_index<space.size(); ++subspace_index)
      {
        // extract subspace
        const basis::TwoBodySubspaceJJJPN& subspace = space.GetSubspace(subspace_index);
        basis::TwoBodySpeciesPN two_body_species = subspace.two_body_species();

        // incorporate sector Jmax
        Jmax_by_type[int(two_body_species)] = std::max(
            Jmax_by_type[int(two_body_species)],
            subspace.J()
          );
      }
  }

  void EvaluateCountsByType(
      const basis::TwoBodySectorsJJJPN& sectors,
      std::array<int,3>& num_sectors_by_type,
      std::array<int,3>& size_by_type
    )
  // Count sectors and matrix elements in the upper triangular portion
  // of a set of pp/nn/pn sectors.
  //
  // Based on basis::UpperTriangularEntries.
  //
  // Lower triangular sectors are ignored.  In diagonal sectors,
  // only upper-triangular entries are counted.
  //
  // Arguments:
  //   sectors (..., input) : container for sectors
  //   num_sectors_by_type (std::array<int,3>, output) : number of sectors
  //   size_by_type (std::array<int,3>, output) : number of matrix elements
  {

    // initialize
    num_sectors_by_type = std::array<int,3>({0,0,0});
    size_by_type = std::array<int,3>({0,0,0});

    // count over sectors
    for (int sector_index=0; sector_index<sectors.size(); ++sector_index)
      {

        // extract sector
        const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors.GetSector(sector_index);
        const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
        const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();

        // characterize pp/nn/pn type for sector
        assert(bra_subspace.two_body_species()==ket_subspace.two_body_species());  // operator Tz=0
        basis::TwoBodySpeciesPN two_body_species = bra_subspace.two_body_species();

        // count sector
        ++num_sectors_by_type[int(two_body_species)];

        // count sector entries
        int sector_entries = 0;
        if (sector.IsDiagonal())
          // diagonal sector
          {
            int dimension = ket_subspace.size();
            sector_entries = dimension*(dimension+1)/2;
          }
        else if (sector.IsUpperTriangle())
          // upper triangle sector (but not diagonal)
          {
            int bra_dimension = bra_subspace.size();
            int ket_dimension = ket_subspace.size();
            sector_entries = bra_dimension*ket_dimension;
          }
        size_by_type[int(two_body_species)] += sector_entries;
      }
  }


  ////////////////////////////////////////////////////////////////
  // file binary I/O
  ////////////////////////////////////////////////////////////////
  // code assumes int and float are both 4-byte

  const int kIntegerSize = 4;
  const int kFloatSize = 4;

  const std::array<const char*,3> kH2ModeDescription({"text","binary","matrix"});
  const std::array<const char*,3> kH2ModeExtension({".dat",".bin",".mat"});

  H2Mode DeducedIOMode(const std::string& filename)
  {
    if (filename.length() < 4 )
      {
	// prevent compare on underlength string
	std::cerr << "H2 file I/O: No extension found (too short) in filename " << filename << std::endl;
	exit(EXIT_FAILURE);
      }
    else if ( ! filename.compare(filename.length()-4,4,".dat") )
      return H2Mode::kText;
    else if ( ! filename.compare(filename.length()-4,4,".mat") )
      return H2Mode::kMatrix;
    else if ( ! filename.compare(filename.length()-4,4,".bin") )
      return H2Mode::kBinary;
    else
      {
	std::cerr << "H2 file I/O: Extension unrecognized in filename " << filename << std::endl;
	exit(EXIT_FAILURE);
      }
  }

  ////////////////////////////////////////////////////////////////
  // H2StreamBase
  ////////////////////////////////////////////////////////////////

  bool H2StreamBase::SectorIsFirstOfType() const
  {
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    int current_type_index = int(sector.ket_subspace().two_body_species());

    // count sectors of prior types (exclusive)
    int prior_type_sectors = 0;
    for (int type_index=0; type_index<current_type_index;++type_index)
      prior_type_sectors+=num_sectors_by_type()[type_index];

    bool is_first_of_type = (sector_index_==prior_type_sectors);
    return is_first_of_type;
  }

  bool H2StreamBase::SectorIsLastOfType() const
  {
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    int current_type_index = int(sector.ket_subspace().two_body_species());

    // count sectors of prior types (inclusive)
    int prior_type_sectors = 0;
    for (int type_index=0; type_index<=current_type_index;++type_index)
      prior_type_sectors+=num_sectors_by_type()[type_index];

    bool is_last_of_type = (sector_index_+1==prior_type_sectors);
    return is_last_of_type;
  }

  std::string H2StreamBase::DiagnosticStr() const
  {

    // TODO: make use of oscillator_like_ check for two-body space when available

    std::string str = fmt::format(
        "  File: {}\n"
        "  Format: {} ({})\n"
        "  Operator properties: J0 {} g0 {} Tz0 {}\n"
        "  Orbitals: p {} n {} (oscillator-like {})\n"
        // "  One-body truncation: p {:.4f} n {:.4f}\n"
        // "  Two-body truncation: pp {:.4f} nn {:.4f} pn {:.4f}\n"
        "  Truncation: p {:.4f} n {:.4f} pp {:.4f} nn {:.4f} pn {:.4f}\n"
        "  Sectors: pp {} nn {} pn {} => total {}\n"
        "  Matrix elements: pp {} nn {} pn {} => total {}\n",
        filename_,
        int(h2_format()),kH2ModeDescription[int(h2_mode())],
        sectors().J0(),sectors().g0(),sectors().Tz0(),
        orbital_space().GetSubspace(0).size(),orbital_space().GetSubspace(1).size(),orbital_space().is_oscillator_like(),
        space().weight_max().one_body[0],space().weight_max().one_body[1],
        space().weight_max().two_body[0],space().weight_max().two_body[1],space().weight_max().two_body[2],
        num_sectors_by_type()[0],num_sectors_by_type()[1],num_sectors_by_type()[2],num_sectors(),
        size_by_type()[0],size_by_type()[1],size_by_type()[2],size()
      );
    return str;
  }

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: version number
  ////////////////////////////////////////////////////////////////

  // read version: text mode
  void InH2Stream::ReadVersion()
  {

    if (h2_mode()==H2Mode::kText)
      {
        // set up line input
        std::string line;

        // line: format code
        {
          ++line_count_;
          std::getline(stream(),line);
          std::istringstream line_stream(line);
          line_stream >> h2_format_; 
          ParsingCheck(line_stream,line_count_,line);
        }
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 1;
        int bytes = num_fields * kIntegerSize;
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
        mcutils::ReadBinary<int>(stream(),h2_format_);
        mcutils::VerifyBinary<int>(stream(),bytes,"Encountered unexpected value in H2 file","record delimiter");
      }

    StreamCheck(bool(stream()),filename_,"Failure while reading H2 file version code");

  };

  void OutH2Stream::WriteVersion()
  {
    if (h2_mode()==H2Mode::kText)
      {
        stream()
          << fmt::format("{:10d}",int(h2_format())) << std::endl;
      }
    else if (h2_mode()==H2Mode::kBinary)
      {
        int num_fields = 1;
        int bytes = num_fields * kIntegerSize;
        mcutils::WriteBinary<int>(stream(),bytes);
        mcutils::WriteBinary<int>(stream(),h2_format());
        mcutils::WriteBinary<int>(stream(),bytes);
      }
    
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file version code");

  }

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: header
  ////////////////////////////////////////////////////////////////


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
    orbital_space_ = basis::OrbitalSpacePN(N1max);
    space_ = basis::TwoBodySpaceJJJPN(orbital_space_,basis::WeightMax(N1max,N2max));
    int J0 = 0;
    int g0 = 0;
    int Tz0 = 0;
    sectors_ = basis::TwoBodySectorsJJJPN(space_,J0,g0,Tz0);

    // store information by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // validate unused input parameters
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

  void OutH2Stream::WriteHeader_Version15099()
  {

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
        // dump orbitals

        // header: dimensions
        mcutils::WriteBinary<int>(stream(),2*kIntegerSize);
        for (int subspace_index=0; subspace_index < orbital_space().size(); ++subspace_index)
          mcutils::WriteBinary<int>(stream(),orbital_space().GetSubspace(subspace_index).size());
        mcutils::WriteBinary<int>(stream(),2*kIntegerSize);

        // orbital listing body
        for (int subspace_index=0; subspace_index < orbital_space().size(); ++subspace_index)
          {
            const std::vector<basis::OrbitalPNInfo> orbitals = orbital_space().GetSubspace(subspace_index).OrbitalInfo();
            const int num_orbitals = orbitals.size();

            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),orbitals[orbital_index].n);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),orbitals[orbital_index].l);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<int>(stream(),TwiceValue(orbitals[orbital_index].j));
            mcutils::WriteBinary<int>(stream(),num_orbitals*kIntegerSize);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kFloatSize);
            for (int orbital_index=0; orbital_index<num_orbitals; ++orbital_index)
              mcutils::WriteBinary<float>(stream(),orbitals[orbital_index].weight);
            mcutils::WriteBinary<int>(stream(),num_orbitals*kFloatSize);
          }         
        
        // two-body indexing

        // header line 1: operator properties
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);
        mcutils::WriteBinary<int>(stream(),sectors().J0());
        mcutils::WriteBinary<int>(stream(),sectors().g0());
        mcutils::WriteBinary<int>(stream(),sectors().Tz0());
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);
        
        // header line 2: 1-body basis limit
        mcutils::WriteBinary<int>(stream(),2*kFloatSize);
        mcutils::WriteBinary<float>(stream(),space().weight_max().one_body[0]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().one_body[1]);
        mcutils::WriteBinary<int>(stream(),2*kFloatSize);

        // header line 3: 2-body basis limit
        mcutils::WriteBinary<int>(stream(),3*kFloatSize);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[0]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[1]);
        mcutils::WriteBinary<float>(stream(),space().weight_max().two_body[2]);
        mcutils::WriteBinary<int>(stream(),3*kFloatSize);

        // header line 4: 2-body basis a.m. limit
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[0]);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[1]);
        mcutils::WriteBinary<int>(stream(),2*Jmax_by_type()[2]);
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);

        // header line 5: matrix size
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);
        mcutils::WriteBinary<int>(stream(),size_by_type()[0]);
        mcutils::WriteBinary<int>(stream(),size_by_type()[1]);
        mcutils::WriteBinary<int>(stream(),size_by_type()[2]);
        mcutils::WriteBinary<int>(stream(),3*kIntegerSize);
      }
  };

  ////////////////////////////////////////////////////////////////
  // format/mode-specific I/O: sector
  ////////////////////////////////////////////////////////////////

  // Note: Underlying iteration logic for sector I/O follows
  // basis::WriteTwoBodyOperatorJJJPN from
  // basis/jjjpnorb_operator.cpp.

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
        int entries = size_by_type()[int(ket_subspace.two_body_species())];
        mcutils::WriteBinary<int>(stream(),entries*kIntegerSize);
      }
  }


  ////////////////////////////////////////////////////////////////
  // InH2Stream
  ////////////////////////////////////////////////////////////////

  InH2Stream::InH2Stream(
      const std::string& filename
    )
    : H2StreamBase(filename), line_count_(0)
  {

    // determine mode
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::in;
    if (h2_mode_==H2Mode::kBinary)
      mode_argument |= std::ios_base::binary;
    stream_ptr_ = new std::ifstream(filename_,mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for input");

    // read format version and header and set up indexing
    ReadVersion();
    if (h2_format()==kVersion0)
      ReadHeader_Version0();
    //else if (h2_format()==kVersion15099))
    //  ReadHeader_Version15099();
    else
      {
        std::cerr << "Unsupported version encountered when opening H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }
  }

  void InH2Stream::Close()
  {
    stream().close();
  };


  void InH2Stream::ReadOrSkipSector(Eigen::MatrixXd& matrix, bool store)
  {
    // read sector
    if (h2_format()==kVersion0)
      ReadSector_Version0(matrix,store);
    // else if (h2_format()==kVersion15099)
    //   ReadSector_Version15099(matrix,store);
    else
      assert(false);  // format version was already checked when writing header

    // validate status
    StreamCheck(bool(stream()),filename_,"Failure while reading H2 file sector");

    // increment to next sector
    sector_index_++;

  }

  void InH2Stream::ReadSector(Eigen::MatrixXd& matrix)
  {
    bool store = true;
    ReadOrSkipSector(matrix,store);
  }

  void InH2Stream::SkipSector()
  {
    // Note that matrix object is still needed as dummy, to use
    // ReadSector machinery, even though matrix is not used
    Eigen::MatrixXd matrix;

    bool store = false;
    ReadOrSkipSector(matrix,store);
  }

  void InH2Stream::SeekToSector(int seek_index)
  {
    assert(seek_index!=basis::kNone);
    assert(seek_index<sectors().size());
    while (sector_index_<seek_index)
      SkipSector();
  }


  ////////////////////////////////////////////////////////////////
  // OutH2Stream
  ////////////////////////////////////////////////////////////////

  OutH2Stream::OutH2Stream(
      const std::string& filename,
      const basis::OrbitalSpacePN& orbital_space,
      const basis::TwoBodySpaceJJJPN& space,
      const basis::TwoBodySectorsJJJPN& sectors,
      H2Format h2_format
    )
    : H2StreamBase(filename)
  {

    // copy in indexing
    orbital_space_ = orbital_space;
    space_ = space;
    sectors_ = sectors;
      
    // store counts by type
    EvaluateJmaxByType(space_,Jmax_by_type_);
    EvaluateCountsByType(sectors_,num_sectors_by_type_,size_by_type_);

    // set format and mode
    h2_format_ = h2_format;
    h2_mode_ = DeducedIOMode(filename_);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::out;
    if (h2_mode_==H2Mode::kBinary)
      mode_argument |= std::ios_base::binary;
    stream_ptr_ = new std::ofstream(filename_,mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening H2 file for output");

    // write version
    WriteVersion();

    // write header
    if (h2_format_==kVersion0)
      WriteHeader_Version0();
    else if (h2_format_==kVersion15099)
      WriteHeader_Version15099();
    else
      {
        std::cerr << "Unsupported version encountered when opening H2 file " << filename_ << std::endl;
        std::exit(EXIT_FAILURE);
      }
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file header");
  }
  
  void OutH2Stream::WriteSector(
      const Eigen::MatrixXd& matrix,
      basis::NormalizationConversion conversion_mode
    )
  {
    // validate matrix dimensions
    const typename basis::TwoBodySectorsJJJPN::SectorType& sector = sectors().GetSector(sector_index_);
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& bra_subspace = sector.bra_subspace();
    const typename basis::TwoBodySectorsJJJPN::SubspaceType& ket_subspace = sector.ket_subspace();
    assert((matrix.rows()==sector.bra_subspace().size())&&(matrix.cols()==sector.ket_subspace().size()));

    // validate output as NAS
    assert(
        (conversion_mode == basis::NormalizationConversion::kNone)
        || (conversion_mode == basis::NormalizationConversion::kASToNAS)
      );

    // write sector
    if (h2_format()==kVersion0)
      WriteSector_Version0(matrix,conversion_mode);
    else if (h2_format()==kVersion15099)
      WriteSector_Version15099(matrix,conversion_mode);
    else
      assert(false);  // format version was already checked when writing header

    // validate status
    StreamCheck(bool(stream()),filename_,"Failure while writing H2 file sector");

    // increment to next sector
    sector_index_++;

  }

  void OutH2Stream::Close()
  {
    stream().close();
  };

  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
} // namespace
