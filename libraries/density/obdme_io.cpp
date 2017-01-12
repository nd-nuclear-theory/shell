/****************************************************************
  obdme_io.cpp

  Patrick J. Fasano
  University of Notre Dame

****************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "eigen3/Eigen/Core"

#include "density/obdme_io.h"
#include "basis/nlj_operator.h"
#include "mcutils/parsing.h"
#include "mcutils/eigen.h"

namespace shell {

InOBDMEReader::InOBDMEReader(
    const std::string& info_filename,
    const basis::OrbitalSpaceLJPN& orbital_space
  )
  : info_filename_(info_filename), orbital_space_(orbital_space)
{
  std::ios_base::openmode mode_argument = std::ios_base::in;

  // open info stream
  info_stream_ptr_ = new std::ifstream(info_filename_, mode_argument);
  StreamCheck(bool(info_stream()),info_filename_,"Failure opening OBDME info file for input");

  // ingest the file
  ReadInfoHeader();
  ReadInfo();

  // close info stream
  info_stream().close();

}

void InOBDMEReader::SetToIndexing(int order, basis::OrbitalSectorsLJPN& sectors) {
  int g0 = order%2;
  sectors = basis::OrbitalSectorsLJPN(orbital_space_, order, g0, 0);
  // @note hard-coded Tz0=0
}

void InOBDMEReader::ReadInfoHeader() {
  std::string line;
  int line_count_ = 0;

  // line 1: version
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    int version;
    line_stream >> version;
    ParsingCheck(line_stream, line_count_, line);
    assert(version == 1500);
  }

  // line 2: max K in multipole expansion
  {
    ++line_count_;
    std::getline(info_stream(), line);
    std::istringstream line_stream(line);
    line_stream >> max_K_;
    ParsingCheck(line_stream, line_count_, line);
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
    std::getline(info_stream(), line);
  }
}

void InOBDMEReader::ReadInfo() {
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

void InOBDMEReader::ReadMultipole(const std::string& data_filename, int order, basis::MatrixVector& matrices) {
  int data_line_count = 0;
  double temp_matrix_element = 0;

  // open data stream
  std::ios_base::openmode mode_argument = std::ios_base::in;
  std::ifstream data_stream(data_filename, mode_argument);
  StreamCheck(bool(data_stream),data_filename,"Failure opening OBDME data file for input");

  // read header
  ReadDataHeader(data_stream, data_line_count);

  // set up indexing
  basis::OrbitalSectorsLJPN sectors;
  SetToIndexing(order, sectors);

  // initialize matrices
  matrices.resize(sectors.size());
  for (int sector_index=0; sector_index < sectors.size(); ++sector_index) {
    // get next sector
    const basis::OrbitalSectorsLJPN::SectorType sector = sectors.GetSector(sector_index);
    matrices[sector_index].resize(sector.bra_subspace().size(),sector.ket_subspace().size());
  }

  // extract data we want for this multipole
  for (auto& info_line : obdme_info()) {
    data_stream >> temp_matrix_element;

    // skip this element if we don't want it
    if (info_line.multipole != order) { continue; }

    // get location of matrix element
    int sector_index, bra_index, ket_index;
    std::tie(sector_index, bra_index, ket_index) = basis::MatrixElementIndicesLJPN(
      orbital_space(), orbital_space(), sectors, info_line.bra_labels, info_line.ket_labels
    );

    // store matrix element
    matrices[sector_index](bra_index,ket_index) = temp_matrix_element;
  }

  data_stream.close();
}

void InOBDMEReader::ReadDataHeader(std::ifstream& data_stream, int& data_line_count) {
  std::string line;
  int num_protons, num_neutrons;
  int bra_seq, ket_seq, bra_n, ket_n;
  HalfInt bra_J, ket_J, bra_T, ket_T;
  int total_obdmes;

  // line 1: version
  {
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    int version;
    line_stream >> version;
    ParsingCheck(line_stream, data_line_count, line);
    assert(version == 1500);
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

  // line 5: total number of OBDMEs
  {
    ++data_line_count;
    std::getline(data_stream, line);
    std::istringstream line_stream(line);
    line_stream >> total_obdmes;
    ParsingCheck(line_stream, data_line_count, line);
    assert(total_obdmes == (num_proton_obdme_ + num_neutron_obdme_));
  }

}

}  // namespace shell
