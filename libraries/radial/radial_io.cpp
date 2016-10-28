/****************************************************************
  radial_io.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iostream>
#include <string>
#include <vector>

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

}  // namespace shell
