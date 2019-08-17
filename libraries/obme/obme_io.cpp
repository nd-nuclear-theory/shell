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
  mcutils::StreamCheck(bool(stream()), filename_, "Failure opening radial operator file for input");

  ReadHeader();
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
  for (std::size_t sector_index = 0; sector_index < sectors_.size(); ++sector_index)
  {
    matrices.push_back(ReadNextSector());
  }
}

void InOBMEStream::Close() { stream().close(); }

void InOBMEStream::ReadHeader()
{
  std::string line;
  int version;
  char operator_type;
  int j0, g0, Tz0;
  std::size_t num_orbitals_bra, num_orbitals_ket;
  basis::MFDnOrbitalFormat orbital_format;

  // version
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> version;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // operator info and sector constraints
  if (version == 0)
  {
    // operator name, l0max, Tz0, number of bra orbitals, number of ket orbitals
    // WARNING: only l0max=0 is supported now, with implied g0=0

    // legacy fields
    char radial_operator_type;
    int l0max;

    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> radial_operator_type >> l0max >> Tz0 >> num_orbitals_bra
        >> num_orbitals_ket;
    mcutils::ParsingCheck(line_stream, line_count_, line);
    orbital_format = basis::MFDnOrbitalFormat::kVersion15099;
    assert(l0max==0);
    j0 = 0;
    g0 = 0;
  }
  else if (version == 1)
  {
    // operator type, operator power
    // WARNING: this format only supports radial matrix elements
    {
      // legacy fields
      char radial_operator_type;
      int radial_operator_power;  // no longer used
      mcutils::GetLine(stream(), line, line_count_);
      std::istringstream line_stream(line);
      line_stream >> radial_operator_type >> radial_operator_power;
      mcutils::ParsingCheck(line_stream, line_count_, line);
      operator_type_ = basis::OneBodyOperatorType::kRadial;
    }

    // constraint mode, l0max or j0 and g0, Tz0
    // WARNING: only l0max=0 is supported now, with implied g0=0
    {
      char constraint_mode;
      mcutils::GetLine(stream(), line, line_count_);
      std::istringstream line_stream(line);
      line_stream >> constraint_mode;
      mcutils::ParsingCheck(line_stream, line_count_, line);

      if (constraint_mode == 'R') {
        int l0max;
        line_stream >> l0max >> Tz0;
        assert(l0max==0);
        j0 = 0;
        g0 = 0;
      } else if (constraint_mode == 'S') {
        line_stream >> j0 >> g0 >> Tz0;
      }
      mcutils::ParsingCheck(line_stream, line_count_, line);
    }

    // number of bra orbitals, number of ket orbitals
    {
      mcutils::GetLine(stream(), line, line_count_);
      std::istringstream line_stream(line);
      line_stream >> num_orbitals_bra >> num_orbitals_ket;
      mcutils::ParsingCheck(line_stream, line_count_, line);
      orbital_format = basis::MFDnOrbitalFormat::kVersion15099;
    }

  }
  else if (version == 2)
  {
    // operator type, j0, g0, and Tz0
    {
      mcutils::GetLine(stream(), line, line_count_);
      std::istringstream line_stream(line);
      line_stream >> operator_type >> j0 >> g0 >> Tz0;
      operator_type_ = static_cast<basis::OneBodyOperatorType>(operator_type);
      mcutils::ParsingCheck(line_stream, line_count_, line);
    }

    // number of bra orbitals, number of ket orbitals
    {
      mcutils::GetLine(stream(), line, line_count_);
      std::istringstream line_stream(line);
      line_stream >> num_orbitals_bra >> num_orbitals_ket;
      mcutils::ParsingCheck(line_stream, line_count_, line);
      orbital_format = basis::MFDnOrbitalFormat::kVersion15200;
    }
  }

  // bra orbital definitions
  std::string orbital_info_str;
  for (int orbital_line_count=0; orbital_line_count < num_orbitals_bra; ++orbital_line_count)
    {
      mcutils::GetLine(stream(), line, line_count_);
      // older versions (0 and 1) did not store an orbital index
      if ((version==0)||(version==1))
        orbital_info_str.append(std::to_string(orbital_line_count));
      orbital_info_str.append(line);
      orbital_info_str.append("\n");  // need to restore newline to input line
    }
  std::istringstream orbital_info_stream(orbital_info_str);
  basis::OrbitalPNList bra_orbitals = basis::ParseOrbitalPNStream(
      orbital_info_stream,
      /*standalone=*/false,
      orbital_format
    );

  // bra orbital definitions
  orbital_info_str = std::string();
  for (int orbital_line_count=0; orbital_line_count < num_orbitals_ket; ++orbital_line_count)
    {
      mcutils::GetLine(stream(), line, line_count_);
      // older versions (0 and 1) did not store an orbital index
      if ((version==0)||(version==1))
        orbital_info_str.append(std::to_string(orbital_line_count));
      orbital_info_str.append(line);
      orbital_info_str.append("\n");  // need to restore newline to input line
    }
  orbital_info_stream = std::istringstream(orbital_info_str);
  basis::OrbitalPNList ket_orbitals = basis::ParseOrbitalPNStream(
      orbital_info_stream,
      /*standalone=*/false,
      orbital_format
    );


  // set up indexing
  bra_orbital_space_ = basis::OrbitalSpaceLJPN(bra_orbitals);
  ket_orbital_space_ = basis::OrbitalSpaceLJPN(ket_orbitals);
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
  for (std::size_t j = 0; j < sector.bra_subspace().size(); ++j)
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);

    // read columns from line
    for (std::size_t k = 0; k < sector.ket_subspace().size(); ++k)
    {
      line_stream >> sector_matrix(j, k);
    }
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // return matrix
  return sector_matrix;
}

OutOBMEStream::OutOBMEStream(const std::string& filename,
                             const basis::OrbitalSpaceLJPN& bra_space,
                             const basis::OrbitalSpaceLJPN& ket_space,
                             const basis::OrbitalSectorsLJPN& sectors,
                             const basis::OneBodyOperatorType operator_type,
                             const std::string& format_str,
                             bool verbose_mode)
    : OBMEStreamBase(filename, bra_space, ket_space, sectors, operator_type),
      format_str_(format_str),
      verbose_mode_(verbose_mode)
{
  // open stream
  std::ios_base::openmode mode_argument = std::ios_base::trunc;
  stream_ptr_ = new std::ofstream(filename_, mode_argument);
  mcutils::StreamCheck(bool(stream()), filename_, "Failure opening radial operator file for output");

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
  stream() << "#   type j0 g0 Tz0" << std::endl;
  stream() << "#   bra_basis_size ket_basis_size" << std::endl;

  // version number
  stream() << 2 << std::endl;
  ++line_count_;

  // operator info
  stream() << " " << static_cast<char>(operator_type())  //
           << " " << sectors_.j0()   //
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
  stream() << basis::OrbitalDefinitionStr(
      bra_orbitals,
      /*standalone=*/false,
      basis::MFDnOrbitalFormat::kVersion15200
    );
  line_count_ += bra_orbitals.size();

  // blank line separating bras from kets
  stream() << std::endl;
  ++line_count_;

  // ket orbital definitions
  if (verbose_mode_)
  {
    stream() << "# ket space orbitals" << std::endl;
    ++line_count_;
  }
  stream() << basis::OrbitalDefinitionStr(
      ket_orbitals,
      /*standalone=*/false,
      basis::MFDnOrbitalFormat::kVersion15200
    );
  line_count_ += ket_orbitals.size();


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
