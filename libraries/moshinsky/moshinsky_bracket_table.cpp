/******************************************************************************
  moshinsky_bracket_table.cpp -- tabulate Moshinsky brackets

  Tabulation is in form

    Nr lr Nc lc N1 l1 N2 l2 L bracket

    Note these are by oscillator labels N.  In contrast, the
    underlying Moshinsky bracket function accepts radial labels n.

    Ordering is by (N,L) sector:

      for N=0..Nmax
        for L=0..N

    Then ordering within (N,L) sector is:

      for bra (relative-cm)
        for ket (1-2)

    Expanding the basis ordering, this becomes

      // bra
      for (Nr,lr)  // [i.e., Nr=0..N, L=Nr%2..Nr (step 2)]
        for (Nc,lc)
          [subject to Nr+Nc=N and triangle(lr,lc,L)]

          // ket
          for (N1,l1)
            for (N2,l2)
              [subject to N1+N2=N and triangle(l1,l2,L)]


  The implementation uses the (N1,l1,N2,l2;L) basis wrappers defined
  in shell_indexing_nl, rather than directly working with these
  labels.

  M. A. Caprio
  University of Notre Dame

  02/12/16 (mac): Created from code in moshinsky_bracket_test.
  07/04/16 (mac): Comment and namespace updates in restructuring of shell package.
  07/18/23 (pjf): Rewrite for updated indexing.

******************************************************************************/

#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include "moshinsky/moshinsky_bracket.h"
#include "basis/lsjt_scheme.h"


int main(int argc, char **argv)
{

  auto relcm_space = basis::RelativeCMSpaceLSJTN(4);
  auto twobody_space = basis::TwoBodySpaceLSJTN(basis::Rank::kOneBody, 2);

  // iterate over sectors
  for (const auto& relcm_subspace : relcm_space)
  {
    if (relcm_subspace.J() != 2) continue;
    if (relcm_subspace.T() != 0) continue;
    if (relcm_subspace.g() != 0) continue;
    if (!twobody_space.ContainsSubspace(relcm_subspace.labels())) continue;
    const auto& twobody_subspace = twobody_space.LookUpSubspace(relcm_subspace.labels());
    fmt::print(
      "{:}\n",
      relcm_subspace.LabelStr()
    );
    for (const auto& relcm_state : relcm_subspace)
      for(const auto& twobody_state : twobody_subspace)
      {
        double bracket = moshinsky::MoshinskyBracket(
            relcm_state.nr(), relcm_state.lr(), relcm_state.nc(), relcm_state.lc(),
            twobody_state.n1(), twobody_state.l1(), twobody_state.n2(), twobody_state.l2(),
            relcm_subspace.L()
          );

        fmt::print(
            "    {:}  {:} -> {:18.8e}\n",
            relcm_state.labels(), twobody_state.labels(), bracket
          );
      }
  }
  // termination
  return 0;
}
