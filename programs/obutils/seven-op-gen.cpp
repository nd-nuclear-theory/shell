/******************************************************************************
  seven-op-gen.cpp

  Generate one-body operator matrix files for the seven electroweak
  nuclear multipole operators (ref. Donnelly & Haxton, 1979).

  These matrix files can then be used with obscalc-ob to evaluate
  the many-body matrix element

     <J_f || T_J || J_I> = sum_{alpha,beta} psi_J(alpha,beta) <alpha||T_J||beta>

  where psi_J(alpha,beta) are the one-body density matrix elements (OBDMEs)
  from a many-body calculation (e.g. MFDn), and <alpha||T_J||beta> are the
  single-particle reduced matrix elements computed here.

  Syntax:
    seven-op-gen < input

  Input keywords:

    set-output-file  filename
      Path for the output OBME file.

    set-indexing  orbital_filename
      Read the single-particle space from an orbital file.
      (Use orbital-gen to produce one for a given Nmax.)

    set-operator  operator_name
      Which SevenOperator to generate. One of:
        MJ  DeltaJ  DeltaJP  SigmaJ  SigmaJP  SigmaJPP  OmegaJ  OmegaJP

    set-rank  J
      Operator angular-momentum rank J (integer >= 0).

    set-momentum-transfer  q
      Momentum transfer q (same units as 1/b).

    set-oscillator-length  b
      Harmonic oscillator length parameter b (default 1.0).

  Example input:

    set-output-file  MJ_J1_q05.obme
    set-indexing     orbitals.dat
    set-operator     MJ
    set-rank         1
    set-momentum-transfer  1.0
    set-oscillator-length  1.0


  Victor Duménil
  University of Notre Dame & LPC Caen

******************************************************************************/

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "mcutils/parsing.h"
#include "obme/obme_io.h"
#include "obme/obme_operator.h"
#include "obme/seven_operators.h"

////////////////////////////////////////////////////////////////
// Parameter parsing
////////////////////////////////////////////////////////////////

struct RunParameters {
  // I/O
  std::string output_filename;
  std::string orbital_filename;

  // Operator selection
  shell::SevenOperatorType operator_type = shell::SevenOperatorType::kMJ;
  std::string              operator_name = "MJ";

  // Physics parameters
  int    J = 1;
  double q = 1.0;
  double b = 1.0;
  int    Tz0 = 0;   // 0 = isoscalar (same species); ±1 = isovector (beta decay)
};

const std::unordered_map<std::string, shell::SevenOperatorType> kOperatorNameMap = {
  {"MJ",       shell::SevenOperatorType::kMJ       },
  {"DeltaJ",   shell::SevenOperatorType::kDeltaJ   },
  {"DeltaJP",  shell::SevenOperatorType::kDeltaJP  },
  {"SigmaJ",   shell::SevenOperatorType::kSigmaJ   },
  {"SigmaJP",  shell::SevenOperatorType::kSigmaJP  },
  {"SigmaJPP", shell::SevenOperatorType::kSigmaJPP },
  {"OmegaJ",   shell::SevenOperatorType::kOmegaJ   },
  {"OmegaJP",  shell::SevenOperatorType::kOmegaJP  },
};

void PrintUsage(char** argv)
{
  std::cerr << "Usage: " << argv[0] << " < input_file\n"
            << "See program header for input format.\n";
}

RunParameters ReadParameters()
{
  RunParameters p;
  std::string line;
  int line_count = 0;

  while (mcutils::GetLine(std::cin, line, line_count)) {
    std::istringstream ss(line);
    std::string keyword;
    ss >> keyword;

    if (keyword == "set-output-file") {
      ss >> p.output_filename;
      mcutils::ParsingCheck(ss, line_count, line);
    }
    else if (keyword == "set-indexing") {
      ss >> p.orbital_filename;
      mcutils::ParsingCheck(ss, line_count, line);
      mcutils::FileExistCheck(p.orbital_filename, true, false);
    }
    else if (keyword == "set-operator") {
      ss >> p.operator_name;
      mcutils::ParsingCheck(ss, line_count, line);
      auto it = kOperatorNameMap.find(p.operator_name);
      if (it == kOperatorNameMap.end())
        mcutils::ParsingError(line_count, line,
                              "Unknown operator name. Valid names: "
                              "MJ DeltaJ DeltaJP SigmaJ SigmaJP SigmaJPP "
                              "OmegaJ OmegaJP");
      p.operator_type = it->second;
    }
    else if (keyword == "set-rank") {
      ss >> p.J;
      mcutils::ParsingCheck(ss, line_count, line);
      assert(p.J >= 0);
    }
    else if (keyword == "set-tz0") {
      ss >> p.Tz0;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.Tz0 < -1 || p.Tz0 > 1)
        mcutils::ParsingError(line_count, line, "set-tz0: value must be -1, 0, or +1");
    }
    else if (keyword == "set-momentum-transfer") {
      ss >> p.q;
      mcutils::ParsingCheck(ss, line_count, line);
    }
    else if (keyword == "set-oscillator-length") {
      ss >> p.b;
      mcutils::ParsingCheck(ss, line_count, line);
    }
    else {
      mcutils::ParsingError(line_count, line, "Unrecognized keyword");
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  // header
  std::cout << "\nseven-op-gen -- SevenOperator one-body matrix element generator\n\n";

  RunParameters p = ReadParameters();

  // Validate required parameters
  assert(!p.output_filename.empty() && "set-output-file is required");
  assert(!p.orbital_filename.empty() && "set-indexing is required");

  // Build orbital space from file
  std::ifstream orbital_stream(p.orbital_filename);
  mcutils::StreamCheck(bool(orbital_stream), p.orbital_filename,
                       "Failure opening orbital file");
  std::vector<basis::OrbitalPNInfo> orbitals =
      basis::ParseOrbitalPNStream(orbital_stream, true);
  basis::OrbitalSpaceLJPN space(orbitals);

  // Compute parity change g0 for this operator and rank
  const int J0  = p.J;
  const int g0  = shell::SevenOperatorParityChange(p.operator_type, p.J);
  const int Tz0 = p.Tz0;   

  std::cout << "Operator : " << p.operator_name << "\n"
            << "Rank J   : " << J0  << "\n"
            << "g0       : " << g0  << "\n"
            << "Tz0      : " << Tz0 << "\n"
            << "q        : " << p.q << "\n"
            << "b        : " << p.b << "\n"
            << "Output   : " << p.output_filename << "\n\n";

  // Build sectors consistent with the operator's selection rules
  basis::OrbitalSectorsLJPN sectors(space, space, J0, g0, Tz0);
  std::cout << sectors.DebugStr();

  // Generate single-particle RMEs
  basis::OperatorBlocks<double> matrices;
  shell::SevenOperatorsOneBodyOperator(
      p.operator_type, p.J, p.q, p.b,
      space, sectors, matrices);

  // Write to OBME file (readable by obscalc-ob)
  shell::OutOBMEStream os(
      p.output_filename, space, space, sectors,
      basis::OneBodyOperatorType::kSpherical);
  os.Write(matrices);

  std::cout << "\nDone. Matrix elements written to " << p.output_filename << "\n";
  return 0;
}
