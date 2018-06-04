/****************************************************************
  obme_io.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "eigen3/Eigen/Core"

#include "mcutils/eigen.h"
#include "mcutils/parsing.h"
#include "obme/obme_operator.h"
#include "obme/obme_io.h"

namespace shell
{
InOBMEStream::InOBMEStream(const std::string& filename)
    : OBMEStreamBase(filename), line_count_(0)
{
  // open stream
  std::ios_base::openmode mode_argument = std::ios_base::in;
  stream_ptr_ = new std::ifstream(filename_, mode_argument);
  StreamCheck(bool(stream()), filename_, "Failure opening radial operator file for input");

  ReadHeader();
}

void InOBMEStream::GetOneBodyOperator(OneBodyOperator& one_body_operator)
{
  one_body_operator.operator_type = operator_type();
  one_body_operator.radial_operator_type = radial_operator_type();
  one_body_operator.radial_operator_power = radial_operator_power();
  SetToIndexing(one_body_operator.bra_orbital_space, one_body_operator.ket_orbital_space, one_body_operator.sectors);
  Read(one_body_operator.matrices);
}

void InOBMEStream::SetToIndexing(basis::OrbitalSpaceLJPN& bra_orbital_space__,
                                 basis::OrbitalSpaceLJPN& ket_orbital_space__,
                                 basis::OrbitalSectorsLJPN& sectors__)
{
  // subspaces -- simple copy
  bra_orbital_space__ = bra_orbital_space();
  ket_orbital_space__ = ket_orbital_space();

  // sectors -- must reconstruct sectors pointing to these new copies of the subspaces
  sectors__ = basis::OrbitalSectorsLJPN(bra_orbital_space__, ket_orbital_space__,
                                        sectors().j0(), sectors().g0(), sectors().Tz0());
}

void InOBMEStream::Read(basis::OperatorBlocks<double>& matrices)
{
  for (int sector_index = 0; sector_index < sectors_.size(); ++sector_index)
  {
    matrices.push_back(ReadNextSector());
  }
}

void InOBMEStream::Close() { stream().close(); }

void InOBMEStream::ReadHeader()
{
  std::string line;
  int version;
  char operator_type, radial_operator_type;
  int l0max, j0, g0, Tz0;
  int num_orbitals_bra, num_orbitals_ket;

  // version -- but first gobble any comment lines
  while (std::getline(stream(), line), line[0] == '#')
  {
    ++line_count_;
  }
  {
    ++line_count_;
    std::istringstream line_stream(line);
    line_stream >> version;
    ParsingCheck(line_stream, line_count_, line);
  }

  // operator info and sector constraints
  if (version == 0)
  {
    // operator name, l0max, Tz0, number of bra orbitals, number of ket orbitals
    // WARNING: only l0max=0 is supported now, with implied g0=0
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> operator_type >> l0max >> Tz0 >> num_orbitals_bra
        >> num_orbitals_ket;
    ParsingCheck(line_stream, line_count_, line);
    radial_operator_type_ = static_cast<RadialOperatorType>(operator_type);
    radial_operator_power_ = l0max;
    assert(l0max==0);
    j0 = 0;
    g0 = 0;
  }
  else if (version == 1)
  {
    // operator type, operator power
    // WARNING: this format only supports radial matrix elements
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> operator_type >> radial_operator_power_;
      ParsingCheck(line_stream, line_count_, line);
      operator_type_ = basis::OneBodyOperatorType::kRadial;
      radial_operator_type_ = static_cast<RadialOperatorType>(operator_type);
    }

    // constraint mode, l0max or j0 and g0, Tz0
    // WARNING: only l0max=0 is supported now, with implied g0=0
    {
      char constraint_mode;
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> constraint_mode;
      ParsingCheck(line_stream, line_count_, line);

      if (constraint_mode == 'R') {
        line_stream >> l0max >> Tz0;
        assert(l0max==0);
        j0 = 0;
        g0 = 0;
      } else if (constraint_mode == 'S') {
        line_stream >> j0 >> g0 >> Tz0;
      }
      ParsingCheck(line_stream, line_count_, line);
    }

    // number of bra orbitals, number of ket orbitals
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> num_orbitals_bra >> num_orbitals_ket;
      ParsingCheck(line_stream, line_count_, line);
    }

  }
  else if (version == 2)
  {
    // operator type, operator power
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> operator_type >> radial_operator_type >> radial_operator_power_;
      ParsingCheck(line_stream, line_count_, line);
      operator_type_ = static_cast<basis::OneBodyOperatorType>(operator_type);
      radial_operator_type_ = static_cast<RadialOperatorType>(radial_operator_type);
    }

    // j0, g0, and Tz0
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> j0 >> g0 >> Tz0;
      ParsingCheck(line_stream, line_count_, line);
    }

    // number of bra orbitals, number of ket orbitals
    {
      ++line_count_;
      std::getline(stream(), line);
      std::istringstream line_stream(line);
      line_stream >> num_orbitals_bra >> num_orbitals_ket;
      ParsingCheck(line_stream, line_count_, line);
    }
  }

  // blank line
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    assert(line.size() == 0);
  }

  // bra orbital definitions
  std::vector<basis::OrbitalPNInfo> bra_states;

  for (int bra_orbital_count = 0; bra_orbital_count < num_orbitals_bra; ++bra_orbital_count)
  {
    // set up for parsing
    getline(stream(), line);
    ++line_count_;
    std::istringstream line_stream(line);
    if (line.size() == 0) continue;

    basis::OrbitalPNInfo state;
    line_stream >> state;
    ParsingCheck(line_stream, line_count_, line);

    bra_states.push_back(state);
  }

  // blank line
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    assert(line.size() == 0);
  }

  // ket orbital definitions
  std::vector<basis::OrbitalPNInfo> ket_states;
  for (int ket_orbital_count = 0; ket_orbital_count < num_orbitals_ket; ++ket_orbital_count)
  {
    getline(stream(), line);
    // count line
    ++line_count_;

    // set up for parsing
    std::istringstream line_stream(line);
    if (line.size() == 0) continue;

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
  sectors_ = basis::OrbitalSectorsLJPN(bra_orbital_space_, ket_orbital_space_, j0, g0, Tz0);
}

Eigen::MatrixXd InOBMEStream::ReadNextSector()
{
  // get next sector
  const basis::OrbitalSectorsLJPN::SectorType sector = sectors_.GetSector(sector_index_);
  ++sector_index_;

  // construct temporary matrix
  Eigen::MatrixXd sector_matrix(sector.bra_subspace().size(), sector.ket_subspace().size());

  // Read in one row per line
  std::string line;
  for (int j = 0; j < sector.bra_subspace().size(); ++j)
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    if (line_stream.peek() == '#') continue;

    // read columns from line
    for (int k = 0; k < sector.ket_subspace().size(); ++k)
    {
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

OutOBMEStream::OutOBMEStream(const std::string& filename,
                             const basis::OrbitalSpaceLJPN& bra_space,
                             const basis::OrbitalSpaceLJPN& ket_space,
                             const basis::OrbitalSectorsLJPN& sectors,
                             const basis::OneBodyOperatorType operator_type,
                             const RadialOperatorType radial_operator_type,
                             int radial_operator_power,
                             const std::string& format_str,
                             bool verbose_mode)
    : OBMEStreamBase(filename, bra_space, ket_space, sectors, operator_type,
                     radial_operator_type, radial_operator_power),
      format_str_(format_str),
      verbose_mode_(verbose_mode)
{
  // open stream
  std::ios_base::openmode mode_argument = std::ios_base::trunc;
  stream_ptr_ = new std::ofstream(filename_, mode_argument);
  StreamCheck(bool(stream()), filename_, "Failure opening radial operator file for output");

  // write header to file
  WriteHeader();
}

void OutOBMEStream::Close() { stream().close(); }

void OutOBMEStream::Write(const basis::OperatorBlocks<double>& matrices)
{
  // check if this seems like a reasonable OperatorBlocks<double>
  assert(matrices.size() == sectors_.size());

  // loop over sectors and write
  for (auto&& matrix : matrices)
  {
    WriteNextSector(matrix);
  }
}

void OutOBMEStream::WriteHeader()
{
  // get orbital info
  std::vector<basis::OrbitalPNInfo> bra_orbitals = bra_orbital_space_.OrbitalInfo();
  std::vector<basis::OrbitalPNInfo> ket_orbitals = ket_orbital_space_.OrbitalInfo();

  // include some header comments
  stream() << "# shell radial matrix elements file" << std::endl;
  stream() << "# version number 2" << std::endl;
  stream() << "# header lines:" << std::endl;
  stream() << "#   type operator power" << std::endl;
  stream() << "#   j0 g0 Tz0" << std::endl;
  stream() << "#   bra_basis_size ket_basis_size" << std::endl;

  // version number
  stream() << 2 << std::endl;
  ++line_count_;

  // operator info
  stream() << " " << static_cast<char>(operator_type())  //
           << " " << static_cast<char>(radial_operator_type())  //
           << " " << radial_operator_power()                    //
           << std::endl;                                        //
  ++line_count_;

  // sectors mode
  stream() << " " << sectors_.j0()   //
           << " " << sectors_.g0()   //
           << " " << sectors_.Tz0()  //
           << std::endl;
  ++line_count_;

  // number of orbitals
  stream() << bra_orbitals.size() << " " << ket_orbitals.size() << std::endl;
  ++line_count_;

  // blank line
  stream() << std::endl;
  ++line_count_;

  // bra orbital definitions
  if (verbose_mode_)
  {
    stream() << "# bra space orbitals" << std::endl;
    ++line_count_;
  }
  for (auto&& state : bra_orbitals)
  {
    stream() << state << std::endl;
    ++line_count_;
  }

  // blank line separating bras from kets
  stream() << std::endl;
  ++line_count_;

  // ket orbital definitions
  if (verbose_mode_)
  {
    stream() << "# ket space orbitals" << std::endl;
    ++line_count_;
  }
  for (auto&& state : ket_orbitals)
  {
    stream() << state << std::endl;
    ++line_count_;
  }

  // blank line between kets and body of file
  stream() << std::endl;
  ++line_count_;
}

void OutOBMEStream::WriteNextSector(const Eigen::MatrixXd& matrix)
{
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
  if (verbose_mode_)
  {
    stream() << "# bra subspace labels:"
             << " l = " << sector.bra_subspace().l()
             << " 2j = " << sector.bra_subspace().j().TwiceValue()
             << " 2Tz = " << sector.bra_subspace().Tz().TwiceValue() << std::endl;
    ++line_count_;

    stream() << "# ket subspace labels:"
             << " l = " << sector.ket_subspace().l()
             << " 2j = " << sector.ket_subspace().j().TwiceValue()
             << " 2Tz = " << sector.ket_subspace().Tz().TwiceValue() << std::endl;
    ++line_count_;
  }

  stream() << mcutils::FormatMatrix(chopped_matrix, format_str_) << std::endl;

  // blank line separating sectors
  stream() << std::endl;
  ++line_count_;
}

}  // namespace shell
