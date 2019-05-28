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
    int j0,
    basis::OrbitalSectorsLJPN& out_sectors,
    basis::OperatorBlocks<double>& out_matrices
  ) const
{
  out_sectors = sectors(j0);
  out_matrices = matrices(j0);
}

void InOBDMEStream::InitStorage()
{
  // initialize storage
  int num_multipoles = j0_max() - j0_min() + 1;
  sectors_.resize(num_multipoles);
  matrices_.resize(num_multipoles);
  for (int j0 = j0_min(); j0 <= j0_max(); j0++) {
    sectors(j0) = basis::OrbitalSectorsLJPN(orbital_space_, j0, g0_, Tz0_);
    basis::SetOperatorToZero(sectors(j0), matrices(j0));
  }
}

/****************************************************************
  Old-style multi-file OBDME input
****************************************************************/
InOBDMEStreamMulti::InOBDMEStreamMulti(
    const std::string& info_filename,
    const std::string& data_filename,
    const basis::OrbitalSpaceLJPN& orbital_space,
    int g0, int Tz0
  )
  : info_filename_(info_filename), data_filename_(data_filename), InOBDMEStream(orbital_space, g0, Tz0)
{
  std::ios_base::openmode mode_argument = std::ios_base::in;

  // open info stream
  info_stream_ptr_ = new std::ifstream(info_filename_, mode_argument);
  StreamCheck(bool(info_stream()),info_filename_,"Failure opening OBDME info file for input");

  // ingest the file
  ReadInfoHeader();
  ReadInfo();

  // read data
  InitStorage();
  ReadData();

  // close info stream
  info_stream().close();

}

void InOBDMEStreamMulti::ReadInfoHeader() {
  std::string line;
  line_count_ = 0;

  // line 1: version
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> version_number_;
    ParsingCheck(line_stream, line_count_, line);
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
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_sp_orbitals;
    ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_sp_orbitals);
    assert(pn_space.GetSubspace(1).size() == num_sp_orbitals);
  }

  // line 3: max j0 in multipole expansion
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> j0_max_;
    ParsingCheck(line_stream, line_count_, line);
    j0_min_ = 0;
  }

  // line 4: number of OBDMEs per type
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_;
    num_neutron_obdme_ = num_proton_obdme_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 5-6: comment line
  {
    ++line_count_;
    std::getline(info_stream(), line);
    ++line_count_;
    std::getline(info_stream(), line);
  }

}

void InOBDMEStreamMulti::ReadInfoHeader1500() {
  std::string line;

  // make sure we didn't get here by mistake
  assert(version_number_ == 1500);

  // line 2: max j0 in multipole expansion
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> j0_max_;
    ParsingCheck(line_stream, line_count_, line);
    j0_min_ = 0;
  }

  // line 3: number of p and n orbitals
  {
    std::size_t num_proton_orbitals, num_neutron_orbitals;
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_orbitals >> num_neutron_orbitals;
    ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_proton_orbitals);
    assert(pn_space.GetSubspace(1).size() == num_neutron_orbitals);
  }

  // line 4: number of p and n OBDMEs
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_ >> num_neutron_obdme_;
    ParsingCheck(line_stream, line_count_, line);
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
    int index, na, la, twoja, nb, lb, twojb, j0;

    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> index >> na >> la >> twoja >> nb >> lb >> twojb >> j0;
    ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN::kP, na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN::kP, nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, j0);
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
    int index, ia, na, la, twoja, ib, nb, lb, twojb, j0, cls;

    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> index >> ia >> na >> la >> twoja >> ib >> nb >> lb >> twojb >> j0 >> cls;
    ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN(cls-1), na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN(cls-1), nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, j0);
  }
}

void InOBDMEStreamMulti::ReadData() {
  int data_line_count = 0;
  double temp_matrix_element = 0;

  // open data stream
  std::ios_base::openmode mode_argument = std::ios_base::in;
  std::ifstream data_stream(data_filename_, mode_argument);
  StreamCheck(bool(data_stream), data_filename_, "Failure opening OBDME data file for input");

  // read header
  ReadDataHeader(data_stream, data_line_count);

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

    // convert to Rose convention (divide by Hat(j_b))
    temp_matrix_element /= Hat(std::get<3>(info_line.ket_labels));

    // store matrix element
    matrices(info_line.multipole)[sector_index](bra_index,ket_index) = temp_matrix_element;
  }

  data_stream.close();
}

void InOBDMEStreamMulti::ReadDataHeader(std::ifstream& data_stream, int& data_line_count) const {
  std::string line;
  std::size_t num_protons, num_neutrons;
  int bra_seq, ket_seq, bra_n, ket_n;
  HalfInt bra_J, ket_J, bra_T, ket_T;
  int version;

  // line 1: version
  {
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> version;
    ParsingCheck(line_stream, data_line_count, line);
    assert(version == version_number_);
  }

  // line 2: number of particles per class
  // note: hard coded two classes -- I don't care about neutron matter right now
  {
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> num_protons >> num_neutrons;
    ParsingCheck(line_stream, data_line_count, line);
  }

  // line 3: MFDn quantum numbers for bra
  // note: currently ignored
  {
    int bra_2J, bra_2T;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> bra_seq >> bra_2J >> bra_n >> bra_2T;
    ParsingCheck(line_stream, data_line_count, line);
    bra_J = HalfInt(bra_2J, 2);
    bra_T = HalfInt(bra_2T, 2);
  }

  // line 4: MFDn quantum numbers for ket
  // note: currently ignored
  {
    int ket_2J, ket_2T;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> ket_seq >> ket_2J >> ket_n >> ket_2T;
    ParsingCheck(line_stream, data_line_count, line);
    ket_J = HalfInt(ket_2J, 2);
    ket_T = HalfInt(ket_2T, 2);
  }

  // line 5: number of OBDMEs
  if (version == 1405) {
    std::size_t num_proton_obdmes, num_neutron_obdmes;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdmes >> num_neutron_obdmes;
    ParsingCheck(line_stream, data_line_count, line);
    assert(num_proton_obdmes == num_proton_obdme_);
    assert(num_neutron_obdmes == num_neutron_obdme_);
  } else if (version == 1500) {
    std::size_t num_obdmes;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> num_obdmes;
    ParsingCheck(line_stream, data_line_count, line);
    assert(num_obdmes == (num_proton_obdme_ + num_neutron_obdme_));
  }

}


/****************************************************************
  New-style single-file OBDME input
****************************************************************/
InOBDMEStreamSingle::InOBDMEStreamSingle(
    const std::string& filename,
    const basis::OrbitalSpaceLJPN& orbital_space
  )
  : filename_(filename), InOBDMEStream(orbital_space)
{
  std::ios_base::openmode mode_argument = std::ios_base::in;

  // open info stream
  stream_ptr_ = new std::ifstream(filename_, mode_argument);
  StreamCheck(bool(stream()),filename_,"Failure opening OBDME file for input");

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
    ParsingCheck(line_stream, line_count_, line);
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
    int parity_bra;
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> Z_bra_ >> N_bra_ >> seq_bra_ >> TwiceJ_bra_ >> TwiceM_bra_
                >> parity_bra >> n_bra_ >> T_bra_ >> E_bra_;
    assert((parity_bra==+1)||(parity_bra==-1));
    g_bra_ = (parity_bra==-1);
  }

  // line 3: ket quantum numbers
  {
    int parity_ket;
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> Z_ket_ >> N_ket_ >> seq_ket_ >> TwiceJ_ket_ >> TwiceM_ket_
                >> parity_ket >> n_ket_ >> T_ket_ >> E_ket_;
    assert((parity_ket==+1)||(parity_ket==-1));
    g_ket_ = (parity_ket==-1);
  }

  g0_ = (g_bra()-g_ket())%2;
  assert((Tz_bra()-Tz_ket()).IsInteger());
  Tz0_ = static_cast<int>(Tz_bra()-Tz_ket());

  // line 4: min j0 in multipole expansion
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> j0_min_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 5: max j0 in multipole expansion
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> j0_max_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 6: number of p and n orbitals
  std::size_t num_orbitals_p, num_orbitals_n;
  {
    mcutils::GetLine(stream(), line, line_count_);
    std::istringstream line_stream(line);
    line_stream >> num_orbitals_p >> num_orbitals_n;
    ParsingCheck(line_stream, line_count_, line);

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
    ParsingCheck(line_stream, line_count_, line);
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

  std::cout << filename_ << std::endl;

  while (!mcutils::GetLine(stream(), line, line_count_).eof()) {
    int ia, ib, j0;
    double matrix_element;

    std::istringstream line_stream(line);
    line_stream >> ia >> ib >> j0 >> matrix_element;
    ParsingCheck(line_stream, line_count_, line);
    std::cout << ia << " " << ib << " " << K << " " << matrix_element << std::endl;

    // get location of matrix element
    basis::FullOrbitalLabels orbital_a = orbital_list_.at(ia-1);
    basis::FullOrbitalLabels orbital_b = orbital_list_.at(ib-1);
    std::size_t sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors(j0), orbital_a, orbital_b
    );
    assert(sector_index != basis::kNone);
    assert(bra_index != basis::kNone);
    assert(ket_index != basis::kNone);

    // convert to Rose convention (divide by Hat(j_b))
    matrix_element /= Hat(std::get<3>(orbital_b));

    // store matrix element
    matrices(j0)[sector_index](bra_index,ket_index) = matrix_element;
  }

}

}  // namespace shell
