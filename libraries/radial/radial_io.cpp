/****************************************************************
  radial_io.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "eigen3/Eigen/Core"

#include "radial/radial_io.h"
#include "mcutils/parsing.h"
#include "mcutils/eigen.h"

namespace shell {

  InRadialStream::InRadialStream(const std::string& filename)
    : RadialStreamBase(filename), line_count_(0)
  {
    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::in;
    stream_ptr_ = new std::ifstream(filename_, mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening radial operator file for input");

    ReadHeader();
  }

  void InRadialStream::SetToIndexing(
      basis::OrbitalSpaceLJPN& bra_orbital_space__,
      basis::OrbitalSpaceLJPN& ket_orbital_space__,
      basis::OrbitalSectorsLJPN& sectors__
    )
  {
    // subspaces -- simple copy
    bra_orbital_space__ = bra_orbital_space();
    ket_orbital_space__ = ket_orbital_space();

    // sectors -- must reconstruct sectors pointing to these new copies of the subspaces
    sectors__ = basis::OrbitalSectorsLJPN(
        bra_orbital_space__,ket_orbital_space__,
        sectors().l0max(),sectors().Tz0()
      );
  }


  void InRadialStream::Read(basis::MatrixVector& matrices) {
    for (int sector_index=0; sector_index < sectors_.size(); ++sector_index) {
      matrices.push_back(ReadNextSector());
    }
  }

  void InRadialStream::Close() {
    stream().close();
  }

  void InRadialStream::ReadHeader() {
    std::string line;
    char operator_type;
    int l0max, Tz0;
    int num_orbitals_bra, num_orbitals_ket;

    // line 1: version -- but first gobble any comment lines
    while (std::getline(stream(), line), line[0] == '#') {++line_count_;}
    {
      ++line_count_;
      std::istringstream line_stream(line);
      int version;
      line_stream >> version;
      ParsingCheck(line_stream, line_count_, line);
      assert(version == 0);
    }

    // line 2: operator name, l0max, Tz0, number of bra orbitals, number of ket
    // orbitals
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> operator_type
                  >> l0max
                  >> Tz0
                  >> num_orbitals_bra
                  >> num_orbitals_ket;
      ParsingCheck(line_stream, line_count_, line);
    }

    // line 4: blank
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      assert(line.size() == 0);
    }

    // lines 3...num_orbitals_bra+3: bra orbital definitions
    std::vector<basis::OrbitalPNInfo> bra_states;
    int bra_orbital_count = 0;
    while (bra_orbital_count < num_orbitals_bra) {
      getline(stream(), line);
      // count line
      ++line_count_;
      ++bra_orbital_count;

      // set up for parsing
      std::istringstream line_stream(line);
      if (line.size() == 0)
        continue;

      basis::OrbitalPNInfo state;
      line_stream >> state;
      ParsingCheck(line_stream, line_count_, line);

      bra_states.push_back(state);
    }

    // line num_orbitals_bra+4: blank
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      assert(line.size() == 0);
    }

    // lines num_orbitals_bra+4...num_orbitals_bra+4+num_orbitals_ket+3:
    //   ket orbital definitions
    std::vector<basis::OrbitalPNInfo> ket_states;
    int ket_orbital_count = 0;
    while (ket_orbital_count < num_orbitals_ket) {
      getline(stream(), line);
      // count line
      ++line_count_;
      ++ket_orbital_count;

      // set up for parsing
      std::istringstream line_stream(line);
      if (line.size() == 0)
        continue;

      basis::OrbitalPNInfo state;
      line_stream >> state;
      ParsingCheck(line_stream, line_count_, line);

      ket_states.push_back(state);
    }

    // final header line: blank
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      assert(line.size() == 0);
    }

    // set up indexing
    bra_orbital_space_ = basis::OrbitalSpaceLJPN(bra_states);
    ket_orbital_space_ = basis::OrbitalSpaceLJPN(ket_states);
    sectors_ = basis::OrbitalSectorsLJPN(bra_orbital_space_, ket_orbital_space_,
        l0max, Tz0);
    radial_operator_type_ = static_cast<RadialOperatorType>(operator_type);
  }

  Eigen::MatrixXd InRadialStream::ReadNextSector() {
    // get next sector
    const basis::OrbitalSectorsLJPN::SectorType sector = sectors_.GetSector(sector_index_);
    ++sector_index_;

    // construct temporary matrix
    Eigen::MatrixXd sector_matrix(sector.bra_subspace().size(),
                                  sector.ket_subspace().size());

    // Read in one row per line
    std::string line;
    for (int j=0; j < sector.bra_subspace().size(); ++j) {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);

      // read columns from line
      for (int k=0; k < sector.ket_subspace().size(); ++k) {
        line_stream >> sector_matrix(j, k);
      }
      ParsingCheck(line_stream, line_count_, line);
    }

    // enforce blank line at end of sector
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      assert(line.size() == 0);
    }

    // return matrix
    return sector_matrix;
  }

  OutRadialStream::OutRadialStream(const std::string& filename,
      const basis::OrbitalSpaceLJPN& bra_space,
      const basis::OrbitalSpaceLJPN& ket_space,
      const basis::OrbitalSectorsLJPN& sectors,
      const RadialOperatorType radial_operator_type,
      const std::string& format_str,
      bool verbose_mode)
    : RadialStreamBase(filename,bra_space,ket_space,sectors,radial_operator_type), format_str_(format_str), verbose_mode_(verbose_mode)
  {

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::trunc;
    stream_ptr_ = new std::ofstream(filename_, mode_argument);
    StreamCheck(bool(stream()),filename_,"Failure opening radial operator file for output");

    // write header to file
    WriteHeader();
  }

  void OutRadialStream::Close() {
    stream().close();
  }

  void OutRadialStream::Write(const basis::MatrixVector& matrices) {
    // check if this seems like a reasonable MatrixVector
    assert(matrices.size() == sectors_.size());

    // loop over sectors and write
    for (auto&& matrix : matrices) {
      WriteNextSector(matrix);
    }
  }

  void OutRadialStream::WriteHeader() {
    // get orbital info
    std::vector<basis::OrbitalPNInfo> bra_orbitals = bra_orbital_space_.OrbitalInfo();
    std::vector<basis::OrbitalPNInfo> ket_orbitals = ket_orbital_space_.OrbitalInfo();

    // include some header comments
    stream() << "# shell radial matrix elements file" << std::endl;
    stream() << "# version number 0" << std::endl;
    stream() << "# header line:" << std::endl;
    stream() << "# operator l0max Tz0 bra_basis_size ket_basis_size" << std::endl;

    // line 1: version number
    stream() << 0 << std::endl; ++line_count_;

    // line 2: header line
    stream() << " " << static_cast<char>(radial_operator_type_)
             << " " << sectors_.l0max()
             << " " << sectors_.Tz0()
             << " " << bra_orbitals.size()
             << " " << ket_orbitals.size()
             << std::endl; ++line_count_;

    // line 3: blank
    stream() << std::endl; ++line_count_;

    // bra orbital definitions
    if (verbose_mode_) {
      stream() << "# bra space orbitals" << std::endl; ++line_count_;
    }
    for (auto&& state : bra_orbitals) {
      stream() << state << std::endl; ++line_count_;
    }

    // blank line separating bras from kets
    stream() << std::endl; ++line_count_;

    // ket orbital definitions
    if (verbose_mode_) {
      stream() << "# ket space orbitals" << std::endl; ++line_count_;
    }
    for (auto&& state : ket_orbitals) {
      stream() << state << std::endl; ++line_count_;
    }

    // blank line between kets and body of file
    stream() << std::endl; ++line_count_;
  }

  void OutRadialStream::WriteNextSector(const Eigen::MatrixXd& matrix) {
    const int width = 7;
    const int precision = 8;

    // get next sector
    const basis::OrbitalSectorsLJPN::SectorType sector = sectors_.GetSector(sector_index_);
    ++sector_index_;

    // check that this is a valid matrix for this sector
    assert(sector.bra_subspace().size() == matrix.rows());
    assert(sector.ket_subspace().size() == matrix.cols());

    // chop away very small values
    auto chopped_matrix = matrix;
    mcutils::ChopMatrix(chopped_matrix, 1e-14);

    // write labels if in verbose mode
    if (verbose_mode_) {
      stream() << "# bra subspace labels: "
               << "l = " << sector.bra_subspace().l() << " "
               << "2j = " << sector.bra_subspace().j().TwiceValue()
               << std::endl;
      ++line_count_;

      stream() << "# ket subspace labels: "
               << "l = " << sector.ket_subspace().l() << " "
               << "2j = " << sector.ket_subspace().j().TwiceValue()
               << std::endl;
      ++line_count_;
    }

    stream() << mcutils::FormatMatrix(chopped_matrix, format_str_) << std::endl;

    // Write in one row per line
    // std::string line;
    // stream().precision(precision);
    // for (int j=0; j < sector.bra_subspace().size(); ++j) {
    //   // write columns to line
    //   for (int k=0; k < sector.ket_subspace().size(); ++k) {
    //     stream() << " " << std::scientific << std::setw(width+1+precision) << matrix(j, k);
    //   }
    //   stream() << std::endl; ++line_count_;
    // }

    // blank line separating sectors
    stream() << std::endl; ++line_count_;
  }

}  // namespace shell
