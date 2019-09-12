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

  2/12/16 (mac): Created from code in moshinsky_bracket_test.
  7/4//16 (mac): Comment and namespace updates in restructuring of shell package.

  Note: Currently broken due to deprecation of shell_indexing_nl.

******************************************************************************/

#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>


#include "shell/moshinsky_bracket.h"
#include "shell/shell_indexing_nl.h"


int main(int argc, char **argv)
{

  int N_max = 6;

  // iterate over sectors
  for (int N = 0; N <= N_max; N++)
    for (int L = 0; L <= N; ++L)
      {

        //    // TwoBodySpaceNL enumeration from shell_indexing_nl.cpp:
        //
        //    for (int N1 = 0; N1 <= N; ++N1)
        //      for (int l1 = (N1 % 2); l1 <= N1; l1 += 2)
        //        {
        //          int N2 = N - N1;
        //          int p = (N2 + l1 + L) % 2;
        //          int l2_min = abs(l1 - L) + p;
        //          int l2_max = min(N2, l1 + L - p);
        //          for (int l2 = l2_min; l2 <= l2_max; l2 += 2)
        //            {
        //              state.a1 = SPOrbitalNl(N1,l1);
        //              state.a2 = SPOrbitalNl(N2,l2);
        //              state.L = L;
        //              state_set.push_back(state);
        //            }
        //        }


        // define (N,L) space -- used for both bra and ket iteration
        TwoBodyStateSetNl states = TwoBodySpaceNL(N,L);
        std::size_t two_body_dim = states.size();

        // tabulate brackets within space
        for (std::size_t i_bra = 0; i_bra < two_body_dim; ++i_bra)
          for (std::size_t i_ket = 0; i_ket < two_body_dim; ++i_ket)
            {
              TwoBodyStateNl bra = states[i_bra];
              TwoBodyStateNl ket = states[i_ket];
              double bracket = shell::MoshinskyBracket(bra,ket);

              std::cout << std::setw(4) << bra.a1.GetN() << std::setw(4) << bra.a1.Getl()
                        << std::setw(4) << bra.a2.GetN() << std::setw(4) << bra.a2.Getl()
                        << std::setw(4) << ket.a1.GetN() << std::setw(4) << ket.a1.Getl()
                        << std::setw(4) << ket.a2.GetN() << std::setw(4) << ket.a2.Getl()
                        << std::setw(4) << L
                        << std::scientific << std::setw(18) << std::setprecision(8)  << bracket
                        << std::endl;

            }

          }

  // termination
  return 0;
}
