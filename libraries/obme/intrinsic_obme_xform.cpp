/****************************************************************
  intrinsic_obme_xform.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <Eigen/Core>

#include "am/am.h"
#include "am/wigner_gsl.h"
#include "fmt/format.h"
#include "moshinsky/moshinsky_bracket.h"
#include "obme/obme_operator.h"
#include "obme/intrinsic_obme_xform.h"

namespace shell
{

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNSubspace
  ////////////////////////////////////////////////////////////////

  OneBodyOperatorDeltaNSubspace::OneBodyOperatorDeltaNSubspace(int J0, int g0, int Delta_N, int N1max, int N2max)
  : BaseSubspace{{Delta_N}}, J0_{J0}, g0_{g0}, N1max_{N1max}, N2max_{N2max}
  {

    // validate subspace labels
    assert(ValidLabels());

    // iterate over state labels (n1,l1,j1,n2,l2,j2)

    // iterate over N1 (with N2=N1+Delta_N determined); apply 
    for (int N1=0; N1<=N1max; ++N1)
      {

        // deduce N2
        int N2 = N1-Delta_N;

        // apply oscillator truncation constraints
        //   - square truncation N1,N2<=N1max
        //   - triangular truncation N1+N2<=N2max
        
        if (N2>N1max)
          continue;
        if ((N1+N2)>N2max)
          continue;

        // iterate over (j1,j2) within given oscillatorshells (N1,N2)
        for (HalfInt j1 = HalfInt(1,2); j1 <= N1 + HalfInt(1,2); ++j1)
          for (HalfInt j2 = HalfInt(1,2); j2 <= N2 + HalfInt(1,2); ++j2)
            {
              if (!am::AllowedTriangle(j2, J0, j1))
                  continue;

              // recover derived quantum numbers (n,l) from (N,j)
              int l1 = (TwiceValue(j1)-1)/2 + (N1+(TwiceValue(j1)-1)/2)%2;
              int n1 = (N1-l1)/2;
              int l2 = (TwiceValue(j2)-1)/2 + (N2+(TwiceValue(j2)-1)/2)%2;;
              int n2 = (N2-l2)/2;
              
              PushStateLabels(StateLabelsType(n1, l1, j1, n2, l2, j2));
            }
        
      }
  }

  bool OneBodyOperatorDeltaNSubspace::ValidLabels() const
  {

    bool valid = true;

    // // nonnegativity
    valid &= (Delta_N()-g0())%2==0;
    
    return valid;
  }

  std::string OneBodyOperatorDeltaNSubspace::LabelStr() const
  {
    return fmt::format(
        "[{}]",
        Delta_N()
      );
  }

  std::string OneBodyOperatorDeltaNSubspace::DebugStr() const
  {

    std::string str;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        StateType state(*this,state_index);

        str += fmt::format(
            "  index {:3d}"
            "    n1 {:2d} l1 {:2d} j1 {:4s} N1 {:2d}"
            "    n2 {:2d} l2 {:2d} j2 {:4s} N2 {:2d}"
            "\n",
            state_index,
            state.n1(), state.l1(), state.j1().Str(), state.N1(),
            state.n2(), state.l2(), state.j2().Str(), state.N2()
          );
      }

    return str;

  }

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNState
  ////////////////////////////////////////////////////////////////

  std::string OneBodyOperatorDeltaNState::LabelStr() const
  {
    return fmt::format(
        "[{} {} {} {} {} {} : index {}]",
        n1(), l1(), j1().Str(), n2(), l2(), j2().Str(), index()
      );
  }

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNSpace
  ////////////////////////////////////////////////////////////////
  
  OneBodyOperatorDeltaNSpace::OneBodyOperatorDeltaNSpace(int J0, int g0, int Delta_N_max, int N1max, int N2max)
    : J0_{J0}, g0_{g0}, Delta_N_max_{Delta_N_max}, N1max_{N1max}, N2max_{N2max}
  {

    // validate operator parity and Delta_N_max
    assert((Delta_N_max-g0)%2==0);

    // enumerate subspaces
    for (int Delta_N=-Delta_N_max; Delta_N<=Delta_N_max; Delta_N+=2)
      {
        SubspaceType subspace(J0, g0, Delta_N, N1max, N2max);
        PushSubspace(subspace);
      }
  }
  
  std::string OneBodyOperatorDeltaNSpace::DebugStr() const
  {

    std::string str;

    str += fmt::format(
        "J0 {} g0 {} Delta_N_max {} N1max {} N2max {}\n",
        J0(), g0(), Delta_N_max(), N1max(), N2max()
      );

    for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
      {
        const SubspaceType& subspace = GetSubspace(subspace_index);

        str += fmt::format(
            "  index {:3d}  dim {:4d}  labels {}\n",
            subspace_index, subspace.dimension(), subspace.LabelStr()
          );
      }

    return str;

  }

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNSectors
  ////////////////////////////////////////////////////////////////

  OneBodyOperatorDeltaNSectors::OneBodyOperatorDeltaNSectors(
      const OneBodyOperatorDeltaNSpace& space
    )
    : BaseSectors(space)
  {
    // The M matrix (13) of Navratil (2021) is block diagonal in Delta_N, as
    // established by oscillator quanta conservation of the constituent
    // Moshinsky brackets (see remarks in intrinsic_obme_xform.h).  Since each
    // subspace of OneBodyOperatorDeltaNSpace already corresponds to a single
    // value of Delta_N, only diagonal sectors (bra subspace = ket subspace)
    // are populated.
    
    for (std::size_t subspace_index=0; subspace_index<space.size(); ++subspace_index)
      PushSector(subspace_index,subspace_index);
  }

  
}  // namespace shell
