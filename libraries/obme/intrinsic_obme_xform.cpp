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

  /*OneBodyOperatorDeltaNSubspace::OneBodyOperatorDeltaNSubspace(int J0, int g0, int Delta_N, int N1max, int N2max)
  : BaseSubspace{{Delta_N}}, J0_{J0}, g0_{g0}, N1max_{N1max}, N2max_{N2max}
  {

    // validate subspace labels
    assert(ValidLabels());

    // iterate over state labels

    for (int N_tot = 0; N_tot <= N2max; ++N_tot)
      for (int N1=0; N1 <= N1max; ++N1)
        {
          int N2 = N1+ Delta_N;
          if ((N2<0) || (N2>N1max))
            continue;
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
  }*/

  OneBodyOperatorDeltaNSubspace::OneBodyOperatorDeltaNSubspace(int J0, int g0, int Delta_N, int N1max, int N2max)
  : BaseSubspace{{Delta_N}}, J0_{J0}, g0_{g0}, N1max_{N1max}, N2max_{N2max}
  {

    // validate subspace labels
    assert(ValidLabels());

    // iterate over state labels (n1,l1,j1,n2,l2,j2)
    //
    // DEBUGGING: Must use arguments J0, g0, Delta_N, N1max, N2max (or data
    // members) here, instead of (hidden) accessors, since labels() is not yet
    // established during base class construction.

    // iterate over N1 (with N2=N1+Delta_N determined); apply square
    // truncation N1,N2<=N1max and triangular truncation N1+N2<=N2max
    for (int N1=0; N1<=N1max; ++N1)
      {
        int N2 = N1+Delta_N;
        if (N2>N1max)
          continue;
        if ((N1+N2)>N2max)
          continue;

        // iterate over l1 (same parity as N1, since n1=(N1-l1)/2 integer)
        for (int l1=N1%2; l1<=N1; l1+=2)
          {
            int n1 = (N1-l1)/2;

            // iterate over l2 (same parity as N2)
            for (int l2=N2%2; l2<=N2; l2+=2)
              {
                int n2 = (N2-l2)/2;

                // enforce orbital parity constraint l1+l2~g0
                if (((l1+l2)%2)!=g0)
                  continue;

                // iterate over j1=l1+/-1/2 (only j1=1/2 for l1=0)
                HalfInt j1_min = (l1==0) ? HalfInt(1,2) : HalfInt(2*l1-1,2);
                for (HalfInt j1=j1_min; j1<=HalfInt(2*l1+1,2); j1=j1+1)
                  {

                    // iterate over j2=l2+/-1/2 (only j2=1/2 for l2=0)
                    HalfInt j2_min = (l2==0) ? HalfInt(1,2) : HalfInt(2*l2-1,2);
                    for (HalfInt j2=j2_min; j2<=HalfInt(2*l2+1,2); j2=j2+1)
                      {

                        // enforce multipole triangularity (j1,j2,J0)
                        if (!am::AllowedTriangle(j1,j2,J0))
                          continue;

                        PushStateLabels(StateLabelsType(n1,l1,j1,n2,l2,j2));
                      }
                  }
              }
          }
      }

  }

  bool OneBodyOperatorDeltaNSubspace::ValidLabels() const
  {

    bool valid = true;

    // // nonnegativity
    // valid &= l()>=0;
    // 
    // // nonemptiness in given Nmax truncation
    // valid &= l()<=Nmax();

    // TODO -- Nothing to check regarding Delta_N?  Or check not larger than
    // N1max?  Oh, check Delta_N~g0.
    
    return valid;
  }

  std::string OneBodyOperatorDeltaNSubspace::LabelStr() const
  {
    return fmt::format(
        "[{} {} {}]",
        J0(), g0(), Delta_N()
      );
  }

  std::string OneBodyOperatorDeltaNSubspace::DebugStr() const
  {

    std::string str;

    for (std::size_t state_index=0; state_index<size(); ++state_index)
      {
        StateType state(*this,state_index);

        str += fmt::format(
            "  index {:3d}  n1 {:2d} l1 {:2d} j1 {:4s} N1 {:2d}"
            "  n2 {:2d} l2 {:2d} j2 {:4s} N2 {:2d}\n",
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

    assert((Delta_N_max-J0)%2==0);
    for (int Delta_N=-Delta_N_max; Delta_N<=Delta_N_max; Delta_N+=2)
      {
        SubspaceType subspace(J0, g0, Delta_N, N1max, N2max);
        PushSubspace(subspace);
      }
  }
  
  std::string OneBodyOperatorDeltaNSpace::DebugStr() const
  {

    // std::ostringstream os;
    // 
    // const int width = 3;
    // 
    // os << " Nmax " << Nmax()
    //    << std::endl;
    // 
    // for (std::size_t subspace_index=0; subspace_index<size(); ++subspace_index)
    //   {
    //     const SubspaceType& subspace = GetSubspace(subspace_index);
    // 
    //     os
    //       << " " << "index"
    //       << " " << std::setw(width) << subspace_index
    //       << " " << "dim"
    //       << " " << std::setw(width) << subspace.dimension()
    //       << " " << "labels"
    //       << " " << subspace.LabelStr()
    //       << " " << std::endl;
    //   }
    // 
    // return os.str();

    // TODO

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
    // for (std::size_t bra_subspace_index=0; bra_subspace_index<space.size(); ++bra_subspace_index)
    //   for (std::size_t ket_subspace_index=0; ket_subspace_index<space.size(); ++ket_subspace_index)
    //     {
    //       // enforce canonical ordering
    //       if (
    //           (sector_direction == basis::SectorDirection::kCanonical)
    //           && !(bra_subspace_index<=ket_subspace_index)
    //         )
    //         continue;
    // 
    //       // retrieve subspaces
    //       const SubspaceType& bra_subspace = space.GetSubspace(bra_subspace_index);
    //       const SubspaceType& ket_subspace = space.GetSubspace(ket_subspace_index);
    // 
    //       // verify angular momentum, isosopin, and parity selection rules
    //       bool allowed = true;
    //       allowed &= am::AllowedTriangle(ket_subspace.l(),L0,bra_subspace.l());
    //       allowed &= ((ket_subspace.g()+g0+bra_subspace.g())%2==0);
    // 
    //       // push sector
    //       if (allowed)
    //         PushSector(bra_subspace_index,ket_subspace_index);
    //     }

    // TODO

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
