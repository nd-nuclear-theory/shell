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

    // iterate over state labels

    // for (int n = 0; n <= (Nmax-l)/2; ++n)
    //   // DEBUGGING: Must use arguments Nmax and l (or data members Nmax_ and l_)
    //   // here, instead of (hidden) accessors Nmax() and l().
    //   PushStateLabels(StateLabelsType(n));

    // TODO

  }

  bool OneBodyOperatorDeltaNSubspace::ValidLabels() const
  {

    bool valid = true;

    // // nonnegativity
    // valid &= l()>=0;
    // 
    // // nonemptiness in given Nmax truncation
    // valid &= l()<=Nmax();

    // TODO
    
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

    // // TODO (mac): reimplement in terms of fmt library
    // std::ostringstream os;
    // 
    // const int width = 3;
    // 
    // for (std::size_t state_index=0; state_index<size(); ++state_index)
    //   {
    //     StateType state(*this,state_index);
    // 
    //     os
    //       << " " << "index"
    //       << " " << std::setw(width) << state_index
    //       << " " << "n"
    //       << " " << std::setw(width) << state.n()
    //       << " " << "=> " << "N"
    //       << " " << std::setw(width) << state.N()
    //       << std::endl;
    //   }
    // 
    // return os.str();

    // TODO

  }

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNState
  ////////////////////////////////////////////////////////////////

  std::string OneBodyOperatorDeltaNState::LabelStr() const
  {
    // // TODO (mac): reimplement in terms of fmt library
    // std::ostringstream os;
    // 
    // const int width = 0;  // for now, no fixed width
    // 
    // os << "["
    //    << " " << "l"
    //    << " " << std::setw(width) << l()
    //    << " " << "index"
    //    << " " << std::setw(width) << index()
    //    << " :"
    //    << " " << "n"
    //    << " " << std::setw(width) << n()
    //    << " " << "]";
    // 
    // return os.str();

    // TODO
  }

  ////////////////////////////////////////////////////////////////
  // OneBodyOperatorDeltaNSpace
  ////////////////////////////////////////////////////////////////
  
  OneBodyOperatorDeltaNSpace::OneBodyOperatorDeltaNSpace(int J0, int g0, int Delta_N_max, int N1max, int N2max)
    : J0_{J0}, g0_{g0}, Delta_N_max_{Delta_N_max}, N1max_{N1max}, N2max_{N2max}
  {

    // // iterate over l
    // for (int l=0; l<=Nmax_; ++l)
    //   {
    //     SubspaceType subspace(l,Nmax);
    //     PushSubspace(subspace);
    //   }

    // TODO
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
  }

  
}  // namespace shell
