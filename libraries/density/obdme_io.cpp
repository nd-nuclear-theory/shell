/****************************************************************
  obdme_io.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <cmath>

#include "eigen3/Eigen/Core"

#include "density/obdme_io.h"
#include "basis/nlj_operator.h"
#include "mcutils/parsing.h"
#include "mcutils/eigen.h"

namespace shell {

constexpr double k_sqrt_4pi = 3.54490770181103205459633496668229;

void InOBDMEStream::GetMultipole(
    int J0,
    basis::OrbitalSectorsLJPN& out_sectors,
    basis::OperatorBlocks<double>& out_matrices
  ) const
{
  out_sectors = sectors(J0);
  out_matrices = matrices(J0);
}

void InOBDMEStream::InitStorage()
{
  // initialize storage
  int num_multipoles = J0_max() - J0_min() + 1;
  sectors_.resize(num_multipoles);
  matrices_.resize(num_multipoles);
  for (int J0 = J0_min(); J0 <= J0_max(); J0++) {
    sectors(J0) = basis::OrbitalSectorsLJPN(orbital_space_, J0, g0(), Tz0());
    basis::SetOperatorToZero(sectors(J0), matrices(J0));
  }
}

/****************************************************************
  Old-style multi-file OBDME input
****************************************************************/
InOBDMEStreamMulti::InOBDMEStreamMulti(
    const std::string& info_filename,
    const std::string& data_filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    HalfInt J_bra, int g_bra, int n_bra,
    HalfInt J_ket, int g_ket, int n_ket
  )
  : info_filename_(info_filename), data_filename_(data_filename),
    InOBDMEStream(orbital_space, J_bra, g_bra, n_bra, J_ket, g_ket, n_ket)
{
  std::ios_base::openmode mode_argument = std::ios_base::in;

  // open info stream
  info_stream_ptr_ = new std::ifstream(info_filename_, mode_argument);
  mcutils::StreamCheck(bool(info_stream()),info_filename_,"Failure opening OBDME info file for input");

  // ingest the file
  ReadInfoHeader();
  ReadInfo();

  // read data
  ReadData();

  // close info stream
  info_stream().close();

}

void InOBDMEStreamMulti::ReadInfoHeader() {
  std::string line;
  line_count_ = 0;

  // line 1: version
  {
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> version_number_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  if (version_number_ == 1405) {
    ReadInfoHeader1405();
  } else if (version_number_ == 1500) {
    ReadInfoHeader1500();
  } else {
    std::cerr << "ERROR: Unknown OBDME info file version" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamMulti::ReadInfoHeader1405() {
  std::string line;

  // make sure we didn't get here by mistake
  assert(version_number_ == 1405);

  // line 2: number of sp orbitals
  {
    std::size_t num_sp_orbitals;
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_sp_orbitals;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_sp_orbitals);
    assert(pn_space.GetSubspace(1).size() == num_sp_orbitals);
  }

  // line 3: max J0 in multipole expansion
  {
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> J0_max_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
    J0_min_ = 0;
  }

  // line 4: number of OBDMEs per type
  {
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_;
    num_neutron_obdme_ = num_proton_obdme_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }
}

void InOBDMEStreamMulti::ReadInfoHeader1500() {
  std::string line;

  // make sure we didn't get here by mistake
  assert(version_number_ == 1500);

  // line 2: max J0 in multipole expansion
  {
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> J0_max_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
    J0_min_ = 0;
  }

  // line 3: number of p and n orbitals
  {
    std::size_t num_proton_orbitals, num_neutron_orbitals;
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_proton_orbitals >> num_neutron_orbitals;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_proton_orbitals);
    assert(pn_space.GetSubspace(1).size() == num_neutron_orbitals);
  }

  // line 4: number of p and n OBDMEs
  {
    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_ >> num_neutron_obdme_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // line 5-6: comment line
  {
    ++line_count_;
    std::getline(info_stream(), line);
    ++line_count_;
    std::getline(info_stream(), line);
  }
}

void InOBDMEStreamMulti::ReadInfo() {
  // switch to the correct parsing function
  if (version_number_ == 1405) {
    ReadInfo1405();
  } else if (version_number_ == 1500) {
    ReadInfo1500();
  } else {
    std::cerr << "ERROR: Unknown OBDME info file version" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamMulti::ReadInfo1405() {
  // make sure we haven't gotten here by mistake
  assert(version_number_ == 1405);

  // loop through lines and put them into obdme_info_
  for (std::size_t i=0; i < num_proton_obdme_; ++i) {
    std::string line;
    int index, na, la, twoja, nb, lb, twojb, J0;

    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> index >> na >> la >> twoja >> nb >> lb >> twojb >> J0;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN::kP, na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN::kP, nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, J0);
  }

  // version 1405 only stores one set of info for both species
  // run through and copy the InfoLines changing only the class
  for (std::size_t i=0; i < num_neutron_obdme_; ++i) {
    InfoLine fake_line = obdme_info_[i];
    std::get<0>(fake_line.bra_labels) = basis::OrbitalSpeciesPN::kN;
    std::get<0>(fake_line.ket_labels) = basis::OrbitalSpeciesPN::kN;
    obdme_info_.push_back(fake_line);
  }
}

void InOBDMEStreamMulti::ReadInfo1500() {
  // make sure we haven't gotten here by mistake
  assert(version_number_ == 1500);

  // loop through lines and put them into obdme_info_
  for (std::size_t i=0; i < num_proton_obdme_ + num_neutron_obdme_; ++i) {
    std::string line;
    int index, ia, na, la, twoja, ib, nb, lb, twojb, J0, cls;

    mcutils::GetLine(info_stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> index >> ia >> na >> la >> twoja >> ib >> nb >> lb >> twojb >> J0 >> cls;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN(cls-1), na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN(cls-1), nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, J0);
  }
}

void InOBDMEStreamMulti::ReadData() {
  int data_line_count = 0;
  double temp_matrix_element = 0;

  // open data stream
  std::ios_base::openmode mode_argument = std::ios_base::in;
  std::ifstream data_stream(data_filename_, mode_argument);
  mcutils::StreamCheck(bool(data_stream), data_filename_, "Failure opening OBDME data file for input");

  // read header
  ReadDataHeader(data_stream, data_line_count);

  // initialize storage (after header read)
  InitStorage();

  // extract data we want for this multipole
  for (const auto& info_line : obdme_info()) {
    data_stream >> temp_matrix_element;

    // get location of matrix element
    std::size_t sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors(info_line.multipole), info_line.bra_labels, info_line.ket_labels
    );
    assert(sector_index != basis::kNone);
    assert(bra_index != basis::kNone);
    assert(ket_index != basis::kNone);

    // multiply by sqrt(4*pi) to fix MFDn bug
    temp_matrix_element *= k_sqrt_4pi;

    // convert to Rose convention (divide by Hat(j_a))
    temp_matrix_element /= Hat(std::get<3>(info_line.bra_labels));

    // store matrix element
    matrices(info_line.multipole)[sector_index](bra_index,ket_index) = temp_matrix_element;
  }

  data_stream.close();
}

void InOBDMEStreamMulti::ReadDataHeader(std::ifstream& data_stream, int& data_line_count)
{
  std::string line;
  int version;

  // line 1: version
  {
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> version;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    assert(version == version_number_);
  }

  // line 2: number of particles per class
  // note: hard coded two classes -- I don't care about neutron matter right now
  {
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> Z_bra_ >> N_bra_;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    Z_ket_ = Z_bra();
    N_ket_ = N_bra();
  }

  // line 3: MFDn quantum numbers for bra
  // note: currently ignored
  {
    int TwiceM_bra__, TwiceJ_bra__, n_bra__, TwiceT_bra__;
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> TwiceM_bra__ >> TwiceJ_bra__ >> n_bra__ >> TwiceT_bra__;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    M_bra_ = HalfInt(TwiceM_bra__,2);
    assert(TwiceJ_bra__ == J_bra().TwiceValue());
    assert(n_bra__ == n_bra());
    T_bra_ = TwiceT_bra__/2.;
  }

  // line 4: MFDn quantum numbers for ket
  // note: currently ignored
  {
    int TwiceM_ket__, TwiceJ_ket__, n_ket__, TwiceT_ket__;
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> TwiceM_ket__ >> TwiceJ_ket__ >> n_ket__ >> TwiceT_ket__;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    M_ket_ = HalfInt(TwiceM_ket__,2);
    assert(TwiceJ_ket__ == J_ket().TwiceValue());
    assert(n_ket__ == n_ket());
    T_ket_ = TwiceT_ket__/2.;
  }

  // line 5: number of OBDMEs
  if (version == 1405) {
    std::size_t num_proton_obdmes, num_neutron_obdmes;
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdmes >> num_neutron_obdmes;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    assert(num_proton_obdmes == num_proton_obdme_);
    assert(num_neutron_obdmes == num_neutron_obdme_);
  } else if (version == 1500) {
    std::size_t num_obdmes;
    mcutils::GetLine(data_stream, line, data_line_count);
    std::istringstream line_stream(line);
    line_stream >> num_obdmes;
    mcutils::ParsingCheck(line_stream, data_line_count, line);
    assert(num_obdmes == (num_proton_obdme_ + num_neutron_obdme_));
  }

}


/****************************************************************
  New-style single-file OBDME input
****************************************************************/
InOBDMEStreamSingle::InOBDMEStreamSingle(
    const std::string& filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    HalfInt J_bra, int g_bra, int n_bra,
    HalfInt J_ket, int g_ket, int n_ket
  )
  : filename_(filename),
    InOBDMEStream(orbital_space, J_bra, g_bra, n_bra, J_ket, g_ket, n_ket)
{
  std::ios_base::openmode mode_argument = std::ios_base::in;

  // open info stream
  stream_ptr_ = new std::ifstream(filename_, mode_argument);
  mcutils::StreamCheck(bool(stream()),filename_,"Failure opening OBDME file for input");

  // ingest the file
  ReadHeader();
  InitStorage();
  ReadData();

  // close info stream
  stream().close();

}

void InOBDMEStreamSingle::ReadHeader() {
  std::string line;
  line_count_ = 0;

  // line 1: version
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> version_number_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  if (version_number_ == 1520) {
    ReadHeader1520();
  } else {
    std::cerr << "ERROR: Unknown OBDME file version " << version_number_ << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamSingle::ReadHeader1520() {
  std::string line;

  // make sure we didn't get here by mistake
  assert(version_number_ == 1520);

  // line 2: bra quantum numbers
  {
    int parity_bra__, TwiceM_bra__, TwiceJ_bra__, n_bra__;
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> Z_bra_ >> N_bra_ >> seq_bra_ >> TwiceJ_bra__ >> TwiceM_bra__
                >> parity_bra__ >> n_bra__ >> T_bra_ >> E_bra_;
    M_bra_ = HalfInt(TwiceM_bra__,2);
    assert(TwiceJ_bra__ == J_bra().TwiceValue());
    assert((parity_bra__==+1)||(parity_bra__==-1));
    assert(g_bra() == (parity_bra__==-1));
    assert(n_bra__ == n_bra());
  }

  // line 3: ket quantum numbers
  {
    int parity_ket__, TwiceM_ket__, TwiceJ_ket__, n_ket__;
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> Z_ket_ >> N_ket_ >> seq_ket_ >> TwiceJ_ket__ >> TwiceM_ket__
                >> parity_ket__ >> n_ket__ >> T_ket_ >> E_ket_;
    M_ket_ = HalfInt(TwiceM_ket__,2);
    assert(TwiceJ_ket__ == J_ket().TwiceValue());
    assert((parity_ket__==+1)||(parity_ket__==-1));
    assert(g_ket() == (parity_ket__==-1));
    assert(n_ket__ == n_ket());
  }

  assert((Tz_bra()-Tz_ket()).IsInteger());

  // line 4: min J0 in multipole expansion
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> J0_min_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // line 5: max J0 in multipole expansion
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> J0_max_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // line 6: number of p and n orbitals
  std::size_t num_orbitals_p, num_orbitals_n;
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_orbitals_p >> num_orbitals_n;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_orbitals_p);
    assert(pn_space.GetSubspace(1).size() == num_orbitals_n);
  }

  // line 7: number of p and n OBDMEs
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_ >> num_neutron_obdme_;
    mcutils::ParsingCheck(line_stream, line_count_, line);
  }

  // orbital listing body
  {
    std::size_t num_orbitals = num_orbitals_p + num_orbitals_n;
    std::string orbital_info_str;
    for (int orbital_line_count=0; orbital_line_count < num_orbitals; ++orbital_line_count)
      {
        mcutils::GetLine(stream(), line, line_count_);
        orbital_info_str.append(line);
        orbital_info_str.append("\n");  // need to restore newline to input line
      }
    std::istringstream orbital_info_stream(orbital_info_str);
    orbital_list_ = basis::ParseOrbitalPNStream(
        orbital_info_stream,
        /*standalone=*/false,
        basis::MFDnOrbitalFormat::kVersion15200
      );
  }
}

void InOBDMEStreamSingle::ReadData() {
  if (version_number_ == 1520) {
    ReadData1520();
  } else {
    std::cerr << "ERROR: Unknown OBDME file version " << version_number_ << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamSingle::ReadData1520() {
  std::string line;
  while (!mcutils::GetLine(stream(), line, line_count_).eof()) {
    int ia, ib, J0;
    double matrix_element;

    std::istringstream line_stream(line);
    line_stream >> ia >> ib >> J0 >> matrix_element;
    mcutils::ParsingCheck(line_stream, line_count_, line);

    // get location of matrix element
    basis::FullOrbitalLabels orbital_a = orbital_list_.at(ia-1);
    basis::FullOrbitalLabels orbital_b = orbital_list_.at(ib-1);
    std::size_t sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors(J0), orbital_a, orbital_b
    );
    assert(sector_index != basis::kNone);
    assert(bra_index != basis::kNone);
    assert(ket_index != basis::kNone);

    // convert to Rose convention (divide by Hat(j_a))
    matrix_element /= Hat(std::get<3>(orbital_a));

    // store matrix element
    matrices(J0)[sector_index](bra_index,ket_index) = matrix_element;
  }

}

}  // namespace shell
