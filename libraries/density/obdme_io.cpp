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

void InOBDMEStream::GetMultipole(int K, basis::OrbitalSectorsLJPN& sectors, basis::OperatorBlocks<double>& matrices) const {
  assert((K >= min_K_) && (K <= max_K_));
  sectors = sectors_[K-min_K_];
  matrices = matrices_[K-min_K_];
}

void InOBDMEStream::InitStorage()
{
  // initialize storage
  int num_multipoles = max_K_ + 1;
  sectors_.resize(num_multipoles);
  matrices_.resize(num_multipoles);
  for (int K = min_K_; K <= max_K_; K++) {
    sectors_[K] = basis::OrbitalSectorsLJPN(orbital_space_, K, g0_, Tz0_);
    basis::SetOperatorToZero(sectors_[K-min_K_], matrices_[K-min_K_]);
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
    int num_sp_orbitals;
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

  // line 3: max K in multipole expansion
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> max_K_;
    ParsingCheck(line_stream, line_count_, line);
    min_K_ = 0;
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

  // line 2: max K in multipole expansion
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> max_K_;
    ParsingCheck(line_stream, line_count_, line);
    min_K_ = 0;
  }

  // line 3: number of p and n orbitals
  {
    int num_proton_orbitals, num_neutron_orbitals;
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
  for (int i=0; i < num_proton_obdme_; ++i) {
    std::string line;
    int index, na, la, twoja, nb, lb, twojb, K;

    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> index >> na >> la >> twoja >> nb >> lb >> twojb >> K;
    ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN::kP, na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN::kP, nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, K);
  }

  // version 1405 only stores one set of info for both species
  // run through and copy the InfoLines changing only the class
  for (int i=0; i < num_neutron_obdme_; ++i) {
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
  for (int i=0; i < num_proton_obdme_ + num_neutron_obdme_; ++i) {
    std::string line;
    int index, ia, na, la, twoja, ib, nb, lb, twojb, K, cls;

    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> index >> ia >> na >> la >> twoja >> ib >> nb >> lb >> twojb >> K >> cls;
    ParsingCheck(line_stream, line_count_, line);

    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN(cls-1), na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN(cls-1), nb, lb, HalfInt(twojb,2));

    obdme_info_.emplace_back(orbital_a, orbital_b, K);
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
    int sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors_[info_line.multipole], info_line.bra_labels, info_line.ket_labels
    );
    assert(sector_index != basis::kNone);
    assert(bra_index != basis::kNone);
    assert(ket_index != basis::kNone);

    // store matrix element -- multiply by sqrt(4*pi) to fix MFDn bug
    matrices_[info_line.multipole][sector_index](bra_index,ket_index) = temp_matrix_element * 3.54490770181103205459633496668229;
  }

  data_stream.close();
}

void InOBDMEStreamMulti::ReadDataHeader(std::ifstream& data_stream, int& data_line_count) const {
  std::string line;
  int num_protons, num_neutrons;
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
    int proton_obdmes, neutron_obdmes;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> proton_obdmes >> neutron_obdmes;
    ParsingCheck(line_stream, data_line_count, line);
    assert(proton_obdmes == num_proton_obdme_);
    assert(neutron_obdmes == num_neutron_obdme_);
  } else if (version == 1500) {
    int total_obdmes;
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> total_obdmes;
    ParsingCheck(line_stream, data_line_count, line);
    assert(total_obdmes == (num_proton_obdme_ + num_neutron_obdme_));
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
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> version_number_;
    ParsingCheck(line_stream, line_count_, line);
  }

  if (version_number_ == 1600) {
    ReadHeader1600();
  } else {
    std::cerr << "ERROR: Unknown OBDME file version " << version_number_ << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamSingle::ReadHeader1600() {
  std::string line;

  // make sure we didn't get here by mistake
  assert(version_number_ == 1600);

  // line 2: min K in multipole expansion
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> min_K_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 3: max K in multipole expansion
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> max_K_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 4: parity of multipole expansion
  {
    int parity;
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> parity;
    ParsingCheck(line_stream, line_count_, line);

    g0_ = (parity == -1);
  }

  // hard-coded Tz0
  Tz0_ = 0;

  // line 5: number of p and n orbitals
  {
    int num_proton_orbitals, num_neutron_orbitals;
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_orbitals >> num_neutron_orbitals;
    ParsingCheck(line_stream, line_count_, line);

    // this is a silly way to check that the correct number of orbitals exist
    basis::OrbitalSpacePN pn_space(orbital_space().OrbitalInfo());
    assert(pn_space.GetSubspace(0).size() == num_proton_orbitals);
    assert(pn_space.GetSubspace(1).size() == num_neutron_orbitals);
  }

  // line 6: number of p and n OBDMEs
  {
    ++line_count_;
    std::getline(stream(), line);
    std::istringstream line_stream(line);
    line_stream >> num_proton_obdme_ >> num_neutron_obdme_;
    ParsingCheck(line_stream, line_count_, line);
  }

  // line 7-8: comment line
  {
    ++line_count_;
    std::getline(stream(), line);
    ++line_count_;
    std::getline(stream(), line);
  }
}

void InOBDMEStreamSingle::ReadData() {
  if (version_number_ == 1600) {
    ReadData1600();
  } else {
    std::cerr << "ERROR: Unknown OBDME file version " << version_number_ << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void InOBDMEStreamSingle::ReadData1600() {
  std::string line;
  while (!std::getline(stream(), line).eof()) {
    ++line_count_;
    int index, ia, na, la, twoja, ib, nb, lb, twojb, K, cls;
    double matrix_element;

    std::istringstream line_stream(line);
    line_stream >> index >> ia >> na >> la >> twoja >> ib >> nb >> lb >> twojb >> K >> cls >> matrix_element;
    ParsingCheck(line_stream, line_count_, line);

    // get location of matrix element
    basis::FullOrbitalLabels orbital_a(basis::OrbitalSpeciesPN(cls-1), na, la, HalfInt(twoja,2));
    basis::FullOrbitalLabels orbital_b(basis::OrbitalSpeciesPN(cls-1), nb, lb, HalfInt(twojb,2));
    int sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors_[K-min_K_], orbital_a, orbital_b
    );
    assert(sector_index != basis::kNone);
    assert(bra_index != basis::kNone);
    assert(ket_index != basis::kNone);

    // store matrix element
    matrices_[K-min_K_][sector_index](bra_index,ket_index) = matrix_element;
  }

}

}  // namespace shell
