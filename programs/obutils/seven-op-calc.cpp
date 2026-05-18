/******************************************************************************
  seven-op-calc.cpp

  Tabulate the many-body reduced matrix element of a seven electroweak
  multipole operator as a function of momentum transfer q:

      <J_f || T_J(q) || J_I>
          = sum_{alpha,beta} rho_{alpha beta}(J_f, J_I) <alpha || T_J(q) || beta>

  for a mesh of q values.  Single-particle RMEs are computed on-the-fly
  using seven_operators routines (no intermediate OBME files needed).

  This program replaces the two-step workflow
      seven-op-gen | obscalc-ob
  with a single call that loops over q internally.

  Syntax:
    seven-op-calc < input

  Input keywords:

    set-output-file  filename
      Output file for the q-dependent table.

    set-indexing  orbital_filename
      Single-particle space (produced by orbital-gen).

    set-operator  operator_name
      One of: MJ DeltaJ DeltaJP SigmaJ SigmaJP SigmaJPP OmegaJ OmegaJP

    set-rank  J
      Operator angular-momentum rank J (integer >= 0).

    set-tz0  Tz0
      Isospin change: 0 (elastic / EM), +1 (beta-minus), -1 (beta-plus).
      Default: 0.

    set-oscillator-length  b
      HO length parameter b in fm (default 1.0).

    set-nucleon-number  A
      Mass number of the nucleus (for translational-invariance correction).

    set-q-range  q_min  q_max  dq
      Momentum transfer mesh in fm^-1 (both ends inclusive).

    define-densities  Jf gf nf  Ji gi ni  robdme_filename [robdme_info_filename]
      One-body density matrix elements for a single bra-ket pair.
      The same syntax as in obscalc-ob; both single-file and multi-file
      (info + data) ROBDME formats are supported.

  Output format:

    [One-body operator momentum mesh]
    # operator J0 Tz0
    MJ 1 +1
    #   Jf  gf  nf    Ji  gi  ni
       1.0   0   1   0.0   0   1
    #              q               rme
      0.00000000e-00   +1.00000000e-01
      1.00000000e-01   +1.00000000e-01
      ...

  Notes:
    - One section per density pair is written.
    - Sections for which the selection rules forbid the matrix element
      (parity, isospin, triangle) are silently skipped.
    - RMEs are in Edmonds convention (consistent with obscalc-ob / MFDn).

  References:
    [1] Donnelly & Haxton, ADNDT 23 (1979) 103.
    [2] Haxton & Lunardini, CPC 179 (2008) 345.

  Victor Duménil
  University of Notre Dame & LPC Caen

  + 05/15/2026 (vd): Created.

******************************************************************************/

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "am/halfint.h"
#include "am/wigner_gsl.h"
#include "basis/nlj_orbital.h"
#include "basis/proton_neutron.h"
#include "density/obdme_io.h"
#include "fmt/format.h"
#include "mcutils/parsing.h"
#include "obme/ob_observable.h"
#include "obme/obme_operator.h"
#include "obme/seven_operators.h"

#include <chrono>

////////////////////////////////////////////////////////////////
// q-mesh helper
////////////////////////////////////////////////////////////////

// Expand (q_min, q_max, dq) into a vector of q values (both ends inclusive).
std::vector<double> QMesh(double q_min, double q_max, double dq)
{
  std::vector<double> mesh;
  for (double q = q_min; q <= q_max + 0.1 * dq; q += dq)
    mesh.push_back(q);
  return mesh;
}

////////////////////////////////////////////////////////////////
// Parameter storage
////////////////////////////////////////////////////////////////

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

struct RunParameters {
  // I/O
  std::string output_filename;
  std::string orbital_filename;

  // Operator
  shell::SevenOperatorType operator_type = shell::SevenOperatorType::kMJ;
  std::string              operator_name = "MJ";
  int    J   = 0;
  int    Tz0 = 0;

  // Physics
  double b   = 1.0;
  int    A   = 2;

  // q mesh
  double q_min = 0.0;
  double q_max = 1.0;
  double dq    = 0.1;
  bool   q_range_set = false;

  // Orbital space (filled from orbital_filename)
  basis::OrbitalSpaceLJPN space;

  // Density streams (one per bra-ket pair)
  std::vector<std::unique_ptr<shell::InOBDMEStream>> density_streams;
};

////////////////////////////////////////////////////////////////
// Parameter parsing
////////////////////////////////////////////////////////////////

void ReadParameters(RunParameters& p)
{
  std::cout << "START ReadParameters" << std::endl;
  std::string line;
  int line_count = 0;

  //while (mcutils::GetLine(std::cin, line, line_count)) {
  while (std::getline(std::cin, line)) {
    std::cout << "EOF reached" << std::endl;
    std::istringstream ss(line);
    std::string keyword;
    ss >> keyword;

    if (keyword == "set-output-file") {
      ss >> p.output_filename;
      mcutils::ParsingCheck(ss, line_count, line);

    } else if (keyword == "set-indexing") {
      ss >> p.orbital_filename;
      mcutils::ParsingCheck(ss, line_count, line);
      mcutils::FileExistCheck(p.orbital_filename, true, false);

      std::ifstream orbital_stream(p.orbital_filename);
      std::vector<basis::OrbitalPNInfo> orbitals =
          basis::ParseOrbitalPNStream(orbital_stream, true);
      p.space = basis::OrbitalSpaceLJPN(orbitals);

    } else if (keyword == "set-operator") {
      ss >> p.operator_name;
      mcutils::ParsingCheck(ss, line_count, line);
      auto it = kOperatorNameMap.find(p.operator_name);
      if (it == kOperatorNameMap.end())
        mcutils::ParsingError(line_count, line,
            "Unknown operator. Valid: MJ DeltaJ DeltaJP SigmaJ SigmaJP SigmaJPP OmegaJ OmegaJP");
      p.operator_type = it->second;

    } else if (keyword == "set-rank") {
      ss >> p.J;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.J < 0)
        mcutils::ParsingError(line_count, line, "set-rank: J must be >= 0");

    } else if (keyword == "set-tz0") {
      ss >> p.Tz0;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.Tz0 < -1 || p.Tz0 > 1)
        mcutils::ParsingError(line_count, line, "set-tz0: value must be -1, 0, or +1");

    } else if (keyword == "set-oscillator-length") {
      ss >> p.b;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.b <= 0.)
        mcutils::ParsingError(line_count, line, "set-oscillator-length: b must be > 0");

    } else if (keyword == "set-nucleon-number") {
      ss >> p.A;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.A <= 0)
        mcutils::ParsingError(line_count, line, "set-nucleon-number: A must be > 0");

    } else if (keyword == "set-q-range") {
      ss >> p.q_min >> p.q_max >> p.dq;
      mcutils::ParsingCheck(ss, line_count, line);
      if (p.dq <= 0.)
        mcutils::ParsingError(line_count, line, "set-q-range: dq must be > 0");
      if (p.q_min > p.q_max)
        mcutils::ParsingError(line_count, line, "set-q-range: q_min > q_max");
      p.q_range_set = true;

    } else if (keyword == "define-densities") {
      float Jf, Ji;
      int gf, nf, gi, ni;
      std::string robdme_filename, robdme_info_filename = "";
      ss >> Jf >> gf >> nf >> Ji >> gi >> ni >> robdme_filename;
      mcutils::ParsingCheck(ss, line_count, line);
      mcutils::FileExistCheck(robdme_filename, true, false);

      std::cout << "define-densities parsed" << std::endl;
      if (!ss.eof()) {
        ss >> robdme_info_filename;
        mcutils::ParsingCheck(ss, line_count, line);
        mcutils::FileExistCheck(robdme_info_filename, true, false);
        p.density_streams.emplace_back(new shell::InOBDMEStreamMulti(
            robdme_info_filename, robdme_filename, p.space,
            HalfInt(2*Jf, 2), gf, nf, HalfInt(2*Ji, 2), gi, ni));
      } else {
        p.density_streams.emplace_back(new shell::InOBDMEStreamSingle(
            robdme_filename, p.space,
            HalfInt(2*Jf, 2), gf, nf, HalfInt(2*Ji, 2), gi, ni));
      }
      std::cout << "stream count = " << p.density_streams.size() << std::endl;

    } else {
      mcutils::ParsingError(line_count, line, "Unrecognized keyword");
    }
  }
  std::cout << "END ReadParameters" << std::endl;
}

////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  std::cout << "\nseven-op-calc -- seven-operator q-mesh evaluation\n\n";

  RunParameters p;
  ReadParameters(p);

  // Validate required inputs
  assert(!p.output_filename.empty()  && "set-output-file is required");
  assert(!p.orbital_filename.empty() && "set-indexing is required");
  assert(p.q_range_set               && "set-q-range is required");
  assert(!p.density_streams.empty()  && "at least one define-densities is required");

  // Derive operator quantum numbers
  const int g0 = shell::SevenOperatorParityChange(p.operator_type, p.J);

  std::cout << "Operator : " << p.operator_name << "\n"
            << "J        : " << p.J  << "\n"
            << "g0       : " << g0   << "\n"
            << "Tz0      : " << p.Tz0 << "\n"
            << "b [fm]   : " << p.b  << "\n"
            << "A        : " << p.A  << "\n"
            << fmt::format("q range  : [{:.4f}, {:.4f}] step {:.4f}\n",
                           p.q_min, p.q_max, p.dq)
            << "Output   : " << p.output_filename << "\n\n";

  // Build sectors (same for all q)
  basis::OrbitalSectorsLJPN sectors(p.space, p.space, p.J, g0, p.Tz0);

  // Build q mesh
  const std::vector<double> q_mesh = QMesh(p.q_min, p.q_max, p.dq);

  // Open output
  std::ofstream out(p.output_filename);
  mcutils::StreamCheck(bool(out), p.output_filename, "Cannot open output file");

  // ----------------------------------------------------------------
  // Outer loop: density pairs (bra-ket state pairs)
  // ----------------------------------------------------------------
  for (const auto& density_stream : p.density_streams)
  {
    const HalfInt J_bra = density_stream->J_bra();
    const HalfInt J_ket = density_stream->J_ket();
    const int     g_bra = density_stream->g_bra();
    const int     g_ket = density_stream->g_ket();
    const int     n_bra = density_stream->n_bra();
    const int     n_ket = density_stream->n_ket();

    // Build output section in a stringstream;
    // only write it if at least one q point succeeds.
    std::ostringstream section;
    int n_written = 0;

    // ----------------------------------------------------------------
    // Inner loop: q values
    // ----------------------------------------------------------------
    for (double q : q_mesh)
    {
      /*
      // Build SPME matrices for this q
      basis::OperatorBlocks<double> matrices;
      shell::SevenOperatorsOneBodyOperator(
          p.operator_type, p.J, q, p.b, p.A,
          p.space, sectors, matrices);

      // Contract with OBDMEs
      double rme = shell::CalculateOneBodyObservableMatrixElement(
          p.space, sectors, matrices, density_stream);*/
      
      std::cout << "Starting q = " << q << std::endl;

      auto t1 = std::chrono::high_resolution_clock::now();

      basis::OperatorBlocks<double> matrices;
      shell::SevenOperatorsOneBodyOperator(
          p.operator_type, p.J, q, p.b, p.A,
          p.space, sectors, matrices);

      auto t2 = std::chrono::high_resolution_clock::now();

      double rme = shell::CalculateOneBodyObservableMatrixElement(
          p.space, sectors, matrices, density_stream);

      auto t3 = std::chrono::high_resolution_clock::now();

      std::chrono::duration<double> dt_op = t2 - t1;
      std::chrono::duration<double> dt_me = t3 - t2;
      std::cout
          << "  operator build : " << dt_op.count() << " s\n"
          << "  contraction    : " << dt_me.count() << " s\n"
          << "  rme            : " << rme << "\n";

      if (std::isnan(rme)) continue;  // selection rules forbid this pair

      // Write header on first successful q point
      if (n_written == 0)
      {
        section << "[One-body operator momentum mesh]\n";
        section << fmt::format("# operator J0 Tz0\n");
        section << fmt::format("{:s} {:d} {:+d}\n",
                               p.operator_name, p.J, p.Tz0);
        section << fmt::format("# {:>4} {:>3} {:>3}  {:>4} {:>3} {:>3}\n",
                               "Jf", "gf", "nf", "Ji", "gi", "ni");
        section << fmt::format("  {:>4.1f} {:>3d} {:>3d}  {:>4.1f} {:>3d} {:>3d}\n",
                               float(J_bra), g_bra, n_bra,
                               float(J_ket), g_ket, n_ket);
        section << fmt::format("# {:>18s}   {:>20s}\n", "q [fm^-1]", "rme");
      }

      section << fmt::format("  {:18.8e}   {:+20.12e}\n", q, rme);
      ++n_written;
    }  // q loop

    if (n_written > 0) {
      section << "\n";
      out << section.str() << std::flush;
      std::cout << fmt::format(
          "  wrote {:d} q points for <({:.1f},{:d},{:d})| O |({:.1f},{:d},{:d})>\n",
          n_written,
          float(J_bra), g_bra, n_bra,
          float(J_ket), g_ket, n_ket);
    } else {
      std::cerr << fmt::format(
          "WARN: selection rules forbid ({:.1f},{:d},{:d})->({:.1f},{:d},{:d}) "
          "for {:s} J={:d} Tz0={:+d}\n",
          float(J_bra), g_bra, n_bra,
          float(J_ket), g_ket, n_ket,
          p.operator_name, p.J, p.Tz0);
    }
  }  // density loop
  
  std::cout << "\nDone. Results written to " << p.output_filename << "\n";
  return 0;
}
