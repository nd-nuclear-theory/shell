/****************************************************************
  ob_observable.cpp

  Mark A. Caprio and Patrick J. Fasano
  University of Notre Dame

****************************************************************/


#include <Eigen/Core>

#include "am/halfint.h"
#include "am/rme.h"
#include "basis/operator.h"
#include "density/obdme_io.h"

namespace shell
{

  double CalculateOneBodyObservableMatrixElement(
      const basis::OrbitalSpaceLJPN& space,
      const basis::OrbitalSectorsLJPN& sectors,
      const basis::OperatorBlocks<double>& blocks,
      const std::unique_ptr<shell::InOBDMEStream>& density_stream
    )
  {
    // return value for disallowed matrix element
    constexpr double double_NaN = std::numeric_limits<double>::quiet_NaN();

    // convenience dereference of density stream pointer
    const shell::InOBDMEStream& obdme_s = *density_stream;

    // convenience quantum numbers
    HalfInt J_bra = obdme_s.J_bra(), J_ket = obdme_s.J_ket();
    HalfInt M_bra = obdme_s.M_bra(), M_ket = obdme_s.M_ket();
    assert(IsInteger(M_bra-M_ket));
    int M0 = int(M_bra-M_ket);
    int J0 = sectors.J0();
    int g0 = sectors.g0();
    int Tz0 = sectors.Tz0();

    // check for parity and isospin-projection; return NaN if disallowed
    if (obdme_s.g0() != g0) return double_NaN;
    if (obdme_s.Tz0() != Tz0) return double_NaN;

    // check for triangularity; return NaN if triangle-disallowed
    if (!am::AllowedTriangle(J_bra, J0, J_ket)) return double_NaN;

    // check for Clebsch zero; return NaN if accidental zero
    double cg_coeff = am::Wigner3J(J_bra, J0, J_ket, -M_bra, M0, M_ket);
    if (std::abs(cg_coeff) < 1e-8) return double_NaN;

    // output NaN if obdmes missing
    if ((J0 < obdme_s.J0_min()) || (J0 > obdme_s.J0_max())) return double_NaN;

    // get necessary density sectors
    basis::OrbitalSectorsLJPN density_sectors;
    basis::OperatorBlocks<double> density_blocks;
    obdme_s.GetMultipole(sectors.J0(), density_sectors, density_blocks);

    // loop and sum over \sum_{a,b} rho_{ab} T_{ab}
    double value = 0.;
    for (std::size_t subspace_index_a=0; subspace_index_a<space.size(); ++subspace_index_a)
      {
        for (std::size_t subspace_index_b=0; subspace_index_b<space.size(); ++subspace_index_b)
          {
            const auto& subspace_a = space.GetSubspace(subspace_index_a);
            const auto& subspace_b = space.GetSubspace(subspace_index_b);
            auto sector_index =
              sectors.LookUpSectorIndex(subspace_index_a, subspace_index_b);
            if (sector_index == basis::kNone) continue;

            for (std::size_t state_index_a = 0; state_index_a < subspace_a.size(); ++state_index_a) {
              for (std::size_t state_index_b = 0; state_index_b < subspace_b.size(); ++state_index_b) {
                value += blocks[sector_index](state_index_a, state_index_b)
                  * density_blocks[sector_index](state_index_a, state_index_b);
              }
            }
          }
      }
    // convert to Edmonds convention
    value *= Hat(J_bra);
    // store value for return
    return value;
  }

}  // namespace shell
