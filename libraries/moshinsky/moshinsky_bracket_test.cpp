/******************************************************************************

  moshinsky_bracket_test.cpp

  M. A. Caprio
  University of Notre Dame

  10/9/16 (pjf): Rename mcpp -> mcutils.
  05/18/26 (mac): Rewrite orthonormality tests for new indexing.

******************************************************************************/


#include <cmath>
#include <iostream>
#include <ostream>
#include <iomanip>
#include <algorithm>

#include <fmt/format.h>

#include "am/am.h"
#include "mcutils/arithmetic.h"
#include "mcutils/profiling.h"

#include "moshinsky/moshinsky_bracket.h"

typedef std::tuple<int, int, int, int> TwoBodySpatialLabels;
void EnumerateTwoBodyStates(int N, int L, std::vector<TwoBodySpatialLabels>& states)
// Enumerate spatial labels (n1,l1,n2,l2) for given N and L.
//
// Based on old TwoBodySpaceNL enumeration (from shell_indexing_nl.cpp).
{
  states.clear();
  for (int N1 = 0; N1 <= N; ++N1)
    for (int l1 = (N1 % 2); l1 <= N1; l1 += 2)
      {
        int N2 = N - N1;
        for (int l2 = (N2 % 2); l2 <= N2; l2 += 2)
          {
            if (!am::AllowedTriangle(l1, l2, L))
              continue;
            int n1 = (N1-l1)/2;
            int n2 = (N2-l2)/2;
            states.push_back(TwoBodySpatialLabels(n1, l1, n2, l2));
          }
      }
}

void SpotChecks()
{
  ////////////////////////////////
  // spot checks
  ////////////////////////////////

  std::cout << "Spot checks" << std::endl;
  std::cout << std::setprecision(8) << std::fixed;

  // seed test -- from CM test
  // < 1 3 0 2 ; 4 | 0 2 0 5 ; 4 >
  // expect -0.39086801 (TTB p. 3)
  std::cout << moshinsky::MoshinskyBracket(1,3,0,2,0,2,0,5,4) << std::endl;


  // case generic -- from CM test
  // < 5 2 0 0 ; 2 | 3 0 2 2 ; 2 >
  // expect 0.22332586 (TTB p. 119)
  std::cout << moshinsky::MoshinskyBracket(5,2,0,0,3,0,2,2,2) << std::endl;

  // trace_moshinsky = true;

  // < 2 1 1 1 ; 1 | 2 1 1 1 ; 1>
  // expect 0.27500000 (TTB p. 94)
  std::cout << moshinsky::MoshinskyBracket(2,1,1,1,2,1,1,1,1) << std::endl;
  // < 1 2 0 4 ; 2 | 2 1 1 1 ; 2>
  // expect -0.16431679 (TTB p. 94)
  std::cout << moshinsky::MoshinskyBracket(1,2,0,4,2,1,1,1,2) << std::endl;

  std::cout << "****" << std::endl;

  // <0 1 1 0 ; 1 | 0 0 1 1 ; 1>
  // expect -0.45643548 (TTB p. 47)
  moshinsky::trace_moshinsky = true;
  std::cout << moshinsky::MoshinskyBracket(0,1,1,0,0,0,1,1,1) << std::endl;
  moshinsky::trace_moshinsky = false;

  std::cout << "****" << std::endl;
};


void SpotChecks_Generalized()
{
  ////////////////////////////////
  // spot checks
  ////////////////////////////////

  std::cout << "Spot checks for Generalized Moshinsky bracket" << std::endl;
  double d = 1.0;
  std::cout << "with d : " << d << std::endl;
  std::cout << std::setprecision(8) << std::fixed;

  // seed test -- from CM test
  // < 1 3 0 2 ; 4 | 0 2 0 5 ; 4 >
  // expect -0.39086801 (TTB p. 3)
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(1,3,0,2,0,2,0,5,4,d) << std::endl;


  // case generic -- from CM test
  // < 5 2 0 0 ; 2 | 3 0 2 2 ; 2 >
  // expect 0.22332586 (TTB p. 119)
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(5,2,0,0,3,0,2,2,2,d) << std::endl;

  // trace_moshinsky = true;

  // < 2 1 1 1 ; 1 | 2 1 1 1 ; 1>
  // expect 0.27500000 (TTB p. 94)
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(2,1,1,1,2,1,1,1,1,d) << std::endl;
  // < 1 2 0 4 ; 2 | 2 1 1 1 ; 2>
  // expect -0.16431679 (TTB p. 94)
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(1,2,0,4,2,1,1,1,2,d) << std::endl;

  std::cout << "****" << std::endl;

  // <0 1 1 0 ; 1 | 0 0 1 1 ; 1>
  // expect -0.45643548 (TTB p. 47)
  moshinsky::trace_moshinsky = true;
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(0,1,1,0,0,0,1,1,1,d) << std::endl;
  moshinsky::trace_moshinsky = false;

  std::cout << "****" << std::endl;

  // Real values checked against 'Tables of Transformation Brackets 
  // for Nuclear Shell-Model Calculations', Moshinsky & Brody, 2nd Ed. (1967)
  // the values come from Daniel Lay 's code 'com_free_operators'
  // < n l, N L ; Lambda | n1 l1, n2 l2 ; Lambda >_d
  // -> (n,l,N,L,n1,l1,n2,l2,Lambda,d)
  std::cout << "Spot checks for Generalized Moshinsky bracket" << std::endl;
  std::cout << "from Daniel Lay 's code 'com_free_operators'" << std::endl;
  // < 0 0, 0 0 ; 0 | 0 0, 0,0 ; 0 >_1
  // expect 1.0 
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(0,0,0,0,0,0,0,0,0,1) << std::endl;

  // < 2 1, 1 1 ; 1 | 2 1, 1 1 ; 1 >_1
  // expect 0.27500000
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(2,1,1,1,2,1,1,1,1,1) << std::endl;

  // < 6 4, 0 0 ; 4 | 3 2, 2 4 ; 4 >_1
  // expect 0.15256312866340879 
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(6,4,0,0,3,2,2,4,4,1) << std::endl;

  // < 6 4, 0 0 ; 4 | 3 2, 2 4 ; 4 >_0.5
  // expect 9.5244368e-02
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(6,4,0,0,3,2,2,4,4,0.5) << std::endl;

  // < 6 4, 0 0 ; 4 | 3 2, 2 4 ; 4 >_0.3
  // expect 3.878176854944e-02
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(6,4,0,0,3,2,2,4,4,0.3) << std::endl;

  // < 8 2, 0 0 ; 4 | 3 3, 2 5 ; 4 >_2.3
  // expect 0.0 ; Efros2021 incorrect
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(8,2,0,0,3,3,2,5,4,2.3) << std::endl;

  // < 7 4, 0 0 ; 4 | 3 3, 2 5 ; 4 >_2.3
  // expect -4.99503829245e-02 ; Efros2021 fast version incorrect
  std::cout << moshinsky::TrlifajGeneralizedMoshinskyBracket(7,4,0,0,3,3,2,5,4,2.3) << std::endl;

};


void OrthonormalityChecksForKet()
// Perform exhaustive orthonormality checks (over the kets).
{
  std::cout << "Orthonormality checks (ket)" << std::endl;
  double d = 1.0;
  std::cout << "with value of d : " << d << std::endl;


  bool verbose = false;
  int N_max = 10;
  for (int N = 0; N <= N_max; N++)
    for (int L = (N % 2); L <= N; ++L)
      {
        // total (N,L) for two-body space
        std::cout << "N " << std::setw(2) << N << ", "
                  << "L " << std::setw(3) << L << ": ";

        // enumerate spatial labels (n1,l1,n2,l2)
        std::vector<TwoBodySpatialLabels> states;
        EnumerateTwoBodyStates(N, L, states);
        std::size_t two_body_dim = states.size();
        std::cout << "dim " << std::setw(5) << two_body_dim << std::endl;

        // enumerate
        if (verbose)
          for (std::size_t ket_index = 0; ket_index < two_body_dim; ++ket_index)
            {
              int n1, l1, n2, l2;
              std::tie(n1, l1, n2, l2) = states[ket_index];
              fmt::print("{} {} {} {}; {}\n", n1, l1, n2, l2, L);
            }
        

        int diag_count = 0, off_diag_count = 0;
        double max_diag_error = 0., max_off_diag_error = 0.;

        mcutils::SteadyTimer t;
        t.Start();
        for (std::size_t ket_index = 0; ket_index < two_body_dim; ++ket_index)
          // for each ket
          {
            for (std::size_t ket_index_prime = ket_index; ket_index_prime < two_body_dim; ++ket_index_prime)
              // for each ket primed
              {
                int n1, l1, n2, l2;
                int n1_prime, l1_prime, n2_prime, l2_prime;
                std::tie(n1, l1, n2, l2) = states[ket_index];
                std::tie(n1_prime, l1_prime, n2_prime, l2_prime) = states[ket_index_prime];

                // evaluate norm sum
                double norm_sum = 0;
                for (std::size_t bra_index = 0; bra_index < two_body_dim; ++bra_index)
                  {
                    int n1_dot, l1_dot, n2_dot, l2_dot;
                    std::tie(n1_dot, l1_dot, n2_dot, l2_dot) = states[bra_index];
                
                    double factor1 = moshinsky::MoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L);
                    double factor2 =  moshinsky::MoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1_prime, l1_prime, n2_prime, l2_prime, L);

                    //double factor1 = moshinsky::GeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L, d);
                    //double factor2 =  moshinsky::GeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1_prime, l1_prime, n2_prime, l2_prime, L, d);

                    //double factor1 = moshinsky::TrlifajGeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L, d); 
                    //double factor2 =  moshinsky::TrlifajGeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1_prime, l1_prime, n2_prime, l2_prime, L, d) ; 

                    norm_sum += factor1 * factor2;

                    if (verbose)
                      fmt::print("  < {} {} {} {}; {} | {} {} {} {} ; {}> < {} {} {} {}; {} | {} {} {} {} ; {}>:  {} {} \n",
                                 n1_dot, l1_dot, n2_dot, l2_dot, L, n1, l1, n2, l2, L,
                                 n1_dot, l1_dot, n2_dot, l2_dot, L, n1_prime, l1_prime, n2_prime, l2_prime, L,
                                 factor1, factor2
                        );
                  }

                if (verbose)
                  fmt::print("overlap < {} {} {} {}; {} | {} {} {} {} ; {}> = {} {}\n", n1, l1, n2, l2, L, n1_prime, l1_prime, n2_prime, l2_prime, L, norm_sum, ket_index == ket_index_prime ? "*" : " ");

                // process norm sum
                if (ket_index == ket_index_prime)
                  // diagonal entry
                  {
                    ++diag_count;
                    max_diag_error = std::max(max_diag_error,fabs(norm_sum-1));
                  }
                else
                  // off-diagonal entry
                  {
                    ++off_diag_count;
                    max_off_diag_error = std::max(max_diag_error,fabs(norm_sum));
                  }

              }

          }
        t.Stop();

        std::cout << std::setprecision(4) << std::scientific;
        std::cout << "  " << "diagonal: " << "entries " << diag_count << ", max error " << max_diag_error << std::endl;
        std::cout << "  " << "off-diag: " << "entries " << off_diag_count << ", max error " << max_off_diag_error << std::endl;
        std::cout << std::setprecision(4) << std::fixed;
        std::cout << "     " << "time " << t.ElapsedTime() << std::endl;
        std::cout << std::endl;

      }
}


void OrthonormalityChecksForBra()
// Perform normalization checks over bras.
//
// Can optionally be restricted to CMF bras.
{
  std::cout << "Normalization checks (bra)" << std::endl;
  double d = 0.6;
  std::cout << "with value of d : " << d << std::endl;

  bool verbose = false;
  bool restrict_to_cmf = true;
  fmt::print("restrict_to_cmf {}\n", restrict_to_cmf);
 
  int N_max = 10;
  for (int N = 0; N <= N_max; N++)
    for (int L = (N % 2); L <= N; ++L)
      {
        // total (N,L) for two-body space
        fmt::print("N {:3d} L {:3d}\n", N, L);

        // enumerate spatial labels (n1,l1,n2,l2)
        std::vector<TwoBodySpatialLabels> states;
        EnumerateTwoBodyStates(N, L, states);
        std::size_t two_body_dim = states.size();
        std::cout << "dim " << std::setw(5) << two_body_dim << std::endl;

        int diag_count = 0, off_diag_count = 0;
        double max_diag_error = 0., max_off_diag_error = 0.;

        mcutils::SteadyTimer t;
        t.Start();
        for (std::size_t state_index_dot = 0; state_index_dot < two_body_dim; ++state_index_dot)
          // for each ket
          {
            for (std::size_t state_index_dot_prime = state_index_dot; state_index_dot_prime < two_body_dim; ++state_index_dot_prime)
              // for each ket primed
              {
                int n1_dot, l1_dot, n2_dot, l2_dot;
                int n1_dot_prime, l1_dot_prime, n2_dot_prime, l2_dot_prime;
                std::tie(n1_dot, l1_dot, n2_dot, l2_dot) = states[state_index_dot];
                std::tie(n1_dot_prime, l1_dot_prime, n2_dot_prime, l2_dot_prime) = states[state_index_dot_prime];

                // restrict to CMF
                if (restrict_to_cmf)
                  if (! (n1_dot==0 && l1_dot==0 && n1_dot_prime==0 && l1_dot_prime==0))
                    continue;
                
                // evaluate norm sum
                double norm_sum = 0;
                for (std::size_t state_index = 0; state_index < two_body_dim; ++state_index)
                  {
                    int n1, l1, n2, l2;
                    std::tie(n1, l1, n2, l2) = states[state_index];

                    //double factor1 = moshinsky::MoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L);
                    //double factor2 =  moshinsky::MoshinskyBracket(n1_dot_prime, l1_dot_prime, n2_dot_prime, l2_dot_prime, n1, l1, n2, l2, L);
                    
                    //double factor1 = moshinsky::GeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L, d);
                    //double factor2 =  moshinsky::GeneralizedMoshinskyBracket(n1_dot_prime, l1_dot_prime, n2_dot_prime, l2_dot_prime, n1, l1, n2, l2, L, d);
                    double factor1 = moshinsky::TrlifajGeneralizedMoshinskyBracket(n1_dot, l1_dot, n2_dot, l2_dot, n1, l1, n2, l2, L, d); 
                    double factor2 =  moshinsky::TrlifajGeneralizedMoshinskyBracket(n1_dot_prime, l1_dot_prime, n2_dot_prime, l2_dot_prime, n1, l1, n2, l2, L, d) ; 
                    norm_sum += factor1 * factor2;

                  }

                // process norm sum
                if (state_index_dot == state_index_dot_prime)
                  // diagonal entry
                  {
                    ++diag_count;
                    max_diag_error = std::max(max_diag_error,fabs(norm_sum-1));
                  }
                else
                  // off-diagonal entry
                  {
                    ++off_diag_count;
                    max_off_diag_error = std::max(max_diag_error,fabs(norm_sum));
                  }

              }

          }
        t.Stop();

        std::cout << std::setprecision(4) << std::scientific;
        std::cout << "  " << "diagonal: " << "entries " << diag_count << ", max error " << max_diag_error << std::endl;
        std::cout << "  " << "off-diag: " << "entries " << off_diag_count << ", max error " << max_off_diag_error << std::endl;
        std::cout << std::setprecision(4) << std::fixed;
        std::cout << "     " << "time " << t.ElapsedTime() << std::endl;
        std::cout << std::endl;

      }
}

int main(int argc, char **argv)
{

  SpotChecks();
  SpotChecks_Generalized();
  OrthonormalityChecksForKet();
  OrthonormalityChecksForBra();

  // termination
  return EXIT_SUCCESS;
}
