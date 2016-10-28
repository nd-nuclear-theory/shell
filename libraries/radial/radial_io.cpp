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

namespace shell {

  std::string RadialOperator::DebugStr() const {
    std::ostringstream os;
    os << "name: " << name << ", l0max: " << l0max << ", Tz0: " << Tz0 << std::endl;
    return os.str();
  }

  InRadialStream::InRadialStream(const std::string& filename)
    : RadialStreamBase(filename), line_count_(0)
  {
    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::in;
    stream_ptr_ = new std::ifstream(filename_, mode_argument);

    ReadHeader();

    for (int sector_index=0; sector_index < sectors_.size(); ++sector_index) {
      ReadNextSector();
    }
  }

  void InRadialStream::Close() {
    (*stream_ptr_).close();
    delete stream_ptr_;
    stream_ptr_ = 0;
  }

  void InRadialStream::ReadHeader() {
    std::string line;
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
      line_stream >> radial_operator_.name
                  >> radial_operator_.l0max
                  >> radial_operator_.Tz0
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
        radial_operator_.l0max, radial_operator_.Tz0);
  }

  void InRadialStream::ReadNextSector() {
    // get next sector
    const basis::BaseSector<basis::OrbitalSubspaceLJPN> sector =
        sectors_.GetSector(sector_index_);
    ++sector_index_;

    // construct temporary matrix
    Eigen::MatrixXd temp_matrix(sector.bra_subspace().size(),
                                sector.ket_subspace().size());

    // Read in one row per line
    std::string line;
    for (int j=0; j < sector.bra_subspace().size(); ++j) {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);

      // read columns from line
      for (int k=0; k < sector.ket_subspace().size(); ++k) {
        line_stream >> temp_matrix(j, k);
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

    // add matrix to vector
    matrices_.push_back(temp_matrix);
  }

  OutRadialStream::OutRadialStream(const std::string& filename,
      const basis::OrbitalSpaceLJPN& bra_space,
      const basis::OrbitalSpaceLJPN& ket_space,
      const RadialOperator& radial_operator)
    : RadialStreamBase(filename)
  {
    // initialize
    bra_orbital_space_ = bra_space;
    ket_orbital_space_ = ket_space;
    radial_operator_ = radial_operator;
    // construct sectors
    sectors_ = basis::OrbitalSectorsLJPN(bra_orbital_space_, ket_orbital_space_,
        radial_operator_.l0max, radial_operator_.Tz0);

    // open stream
    std::ios_base::openmode mode_argument = std::ios_base::trunc;
    stream_ptr_ = new std::ofstream(filename_, mode_argument);
  }

  void OutRadialStream::Close() {
    (*stream_ptr_).close();
    delete stream_ptr_;
    stream_ptr_ = 0;
  }

  void OutRadialStream::Write(const basis::MatrixVector& matrices) {
    // check if this seems like a reasonable MatrixVector
    assert(matrices.size() == sectors_.size());

    // write header to file
    WriteHeader();

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
    stream() << " " << radial_operator_.name
             << " " << radial_operator_.l0max
             << " " << radial_operator_.Tz0
             << " " << bra_orbitals.size()
             << " " << ket_orbitals.size()
             << std::endl; ++line_count_;

    // line 3: blank
    stream() << std::endl; ++line_count_;

    // bra orbital definitions
    for (auto&& state : bra_orbitals) {
      stream() << state << std::endl; ++line_count_;
    }

    // blank line separating bras from kets
    stream() << std::endl; ++line_count_;

    // ket orbital definitions
    for (auto&& state : ket_orbitals) {
      stream() << state << std::endl; ++line_count_;
    }

    // blank line between kets and body of file
    stream() << std::endl; ++line_count_;
  }

  void OutRadialStream::WriteNextSector(const Eigen::MatrixXd& matrix) {
    const int width = 3;
    const int precision = 8;

    // get next sector
    const basis::BaseSector<basis::OrbitalSubspaceLJPN> sector =
        sectors_.GetSector(sector_index_);
    ++sector_index_;

    // check that this is a valid matrix for this sector
    assert(sector.bra_subspace().size() == matrix.rows());
    assert(sector.ket_subspace().size() == matrix.cols());

    // Write in one row per line
    std::string line;
    for (int j=0; j < sector.bra_subspace().size(); ++j) {
      // write columns to line
      for (int k=0; k < sector.ket_subspace().size(); ++k) {
        stream() << " " << std::fixed << std::setw(width+1+precision) << matrix(j, k);
      }
      stream() << std::endl; ++line_count_;
    }

    // blank line separating sectors
    stream() << std::endl; ++line_count_;
  }

}  // namespace shell
