/******************************************************************************

  statecount.cpp -- count J-coupled states (at two-body and three-body level)

  Was written to assist in dimension studies for transforming
  three-body interactions to general radial bases.

  Syntax: statecount

  M. A. Caprio, University of Notre Dame.

  10/20/13 (mac): Originated.
  4/25/15 (mac): Source file reformatted.

******************************************************************************/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>


#include <mcpp/profiling.h>
#include <am/angular_momentum.h>

#include <shell/shell_indexing_nlj.h>


using namespace shell;

////////////////////////////////////////////////////////////////
// utility functions
////////////////////////////////////////////////////////////////

// OrbitalDimension(N,j) returns the number of orbitals of am j through shell N
int OrbitalDimension (int N, HalfInt j)
{
  return IValue(N - j + HalfInt(3,2));
}

////////////////////////////////////////////////////////////////
// two-body counting functions
////////////////////////////////////////////////////////////////

long TwoBodyDimensionOneBodyCutAMCentric (int N1b, HalfInt J)
{
  // do j-centric count 
  // approximates strict canonical ordering by j ordering

  long state_count = 0;

  // iterate over ja <= jb 
  for (HalfInt ja = HalfInt(1,2); ja <= Shelljmax(N1b); ++ja)
    for (HalfInt jb = ja; jb <= Shelljmax(N1b); ++jb)
      if ( AllowedTriangle(ja,jb,J) )
	state_count += OrbitalDimension(N1b,ja) * OrbitalDimension(N1b,jb);
  return state_count;
}

long TwoBodyDimensionOneBodyCut (int N1b, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b 
  for (SPOrbitalNlj a; a.GetN() <= N1b; ++a)
    for (SPOrbitalNlj b = a; b.GetN() <= N1b; ++b)
      if ( AllowedTriangle(a.Getj(),b.Getj(),J) )
	state_count++;
  return state_count;
}

long TwoBodyDimensionTwoBodyCut (int N2b, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b 
  for (SPOrbitalNlj a; a.GetN() <= N2b; ++a)
    for (SPOrbitalNlj b = a; a.GetN() + b.GetN() <= N2b; ++b)
      if ( AllowedTriangle(a.Getj(),b.Getj(),J) )
	state_count++;
  return state_count;
}

////////////////////////////////////////////////////////////////
// three-body counting functions
////////////////////////////////////////////////////////////////

long ThreeBodyDimensionOneBodyCutAMCentric (int N1b, HalfInt J1, HalfInt J)
{
  // do j-centric count 
  // approximates strict canonical ordering by j ordering

  long state_count = 0;

  // iterate over ja <= jb <= jc
  for (HalfInt ja = HalfInt(1,2); ja <= Shelljmax(N1b); ++ja)
    for (HalfInt jb = ja; jb <= Shelljmax(N1b); ++jb)
      for (HalfInt jc = jb; jc <= Shelljmax(N1b); ++jc)
	if ( AllowedTriangle(ja,jb,J1) && AllowedTriangle(J1,jc,J) )
	  state_count += OrbitalDimension(N1b,ja) * OrbitalDimension(N1b,jb) * OrbitalDimension(N1b,jc);
  return state_count;
}

long ThreeBodyDimensionOneBodyCutAMCentric (int N1b, HalfInt J)
{
  long state_count = 0;

  // iterate over 0 <= J1 <= 2 jmax(N)
  for (HalfInt J1 = 0; J1 <= TwiceValue(Shelljmax(N1b)); ++J1)
    state_count += ThreeBodyDimensionOneBodyCutAMCentric (N1b, J1, J);

  return state_count;
}
			
long ThreeBodyDimensionOneBodyCut (int N1b, HalfInt J1, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b <= c
  for (SPOrbitalNlj a; a.GetN() <= N1b; ++a)
    for (SPOrbitalNlj b = a; b.GetN() <= N1b; ++b)
      for (SPOrbitalNlj c = b; c.GetN() <= N1b; ++c)
	if ( AllowedTriangle(a.Getj(),b.Getj(),J1) && AllowedTriangle(J1,c.Getj(),J) )
	  state_count++;
  return state_count;
}

int ThreeBodyDimensionOneBodyCut (int N1b, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b <= c
  for (SPOrbitalNlj a; a.GetN() <= N1b; ++a)
    for (SPOrbitalNlj b = a; b.GetN() <= N1b; ++b)
      for (SPOrbitalNlj c = b; c.GetN() <= N1b; ++c)
	// add in number of allowed intermediate J1 values
	state_count += std::max(
				-std::max(
					  abs(IValue(a.Getj()-b.Getj())),
					  abs(IValue(J-c.Getj()))
					  )
				+std::min(
					  IValue(a.Getj()+b.Getj()),
					  IValue(J+c.Getj())
					  )
				+ 1,
				0
				);
				
  return state_count;
}


long ThreeBodyDimensionThreeBodyCut (int N3b, HalfInt J1, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b <= c
  for (SPOrbitalNlj a; a.GetN() <= N3b; ++a)
    for (SPOrbitalNlj b = a; a.GetN() + b.GetN() <= N3b; ++b)
      for (SPOrbitalNlj c = b; a.GetN() + b.GetN() + c.GetN() <= N3b; ++c)
	if ( AllowedTriangle(a.Getj(),b.Getj(),J1) && AllowedTriangle(J1,c.Getj(),J) )
	  state_count++;
  return state_count;
}

long ThreeBodyDimensionThreeBodyCut (int N3b, HalfInt J)
{
  long state_count = 0;

  // iterate over a <= b <= c
  for (SPOrbitalNlj a; a.GetN() <= N3b; ++a)
    for (SPOrbitalNlj b = a; a.GetN() + b.GetN() <= N3b; ++b)
      for (SPOrbitalNlj c = b; a.GetN() + b.GetN() + c.GetN() <= N3b; ++c)
	// add in number of allowed intermediate J1 values
	state_count += std::max(
				-std::max(
					  abs(IValue(a.Getj()-b.Getj())),
					  abs(IValue(J-c.Getj()))
					  )
				+std::min(
					  IValue(a.Getj()+b.Getj()),
					  IValue(J+c.Getj())
					  )
				+ 1,
				0
				);

  return state_count;
}

////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  ////////////////////////////////////////////////////////////////
  // initialization
  ////////////////////////////////////////////////////////////////

  // header
  std::cout << std::endl;
  std::cout << "statecount" << std::endl;
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // two-body tabulation
  ////////////////////////////////////////////////////////////////

  // N J ~dim(N1b) dim(N1b) dim(N2b)

  // run parameters
  const int N2b_max = 40;

  // start timing
  Timer total_time;
  total_time.Start();

  // two-body table
  for (int N = 0; N <= N2b_max; ++N)
    for (HalfInt J = 0; J <= Shelljmax(N) + Shelljmax(N); ++J)
      {
	std::cout << std::setw(2) << N 
		  << " " << std::setw(2) << J
		  << " " << std::setw(6) << TwoBodyDimensionOneBodyCutAMCentric(N, J) 
		  << " " << std::setw(6) << TwoBodyDimensionOneBodyCut(N, J) 
		  << " " << std::setw(6) << TwoBodyDimensionTwoBodyCut(N, J) 
		  << std::endl;
      }


  // end timing
  total_time.Stop();
  std::cout << std::endl;
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // three-body tabulation by J1
  ////////////////////////////////////////////////////////////////

  // N J J1 ~dim(N1b) dim(N1b) dim(N3b)

  // run parameters
  const int N3b_max_J1 = 0;

  // start timing
  total_time.Start();

  // three-body table
  for (int N = 0; N <= N3b_max_J1; ++N)
    for (HalfInt J = HalfInt(1,2); J <= Shelljmax(N) + Shelljmax(N) + Shelljmax(N); ++J)
      for (HalfInt J1 = 0; J1 <= Shelljmax(N) + Shelljmax(N); ++J1)
	{
	  std::cout << std::setw(2) << N 
		    << " " << std::setw(4) << DValue(J)
		    << " " << std::setw(4) << IValue(J1)
		    << " " << std::setw(10) << ThreeBodyDimensionOneBodyCutAMCentric(N, J1, J) 
		    << " " << std::setw(10) << ThreeBodyDimensionOneBodyCut(N, J1, J) 
		    << " " << std::setw(10) << ThreeBodyDimensionThreeBodyCut(N, J1, J) 
		    << std::endl;
	}


  // end timing
  total_time.Stop();
  std::cout << std::endl;
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;


  ////////////////////////////////////////////////////////////////
  // three-body tabulation
  ////////////////////////////////////////////////////////////////

  // N J ~dim(N1b) dim(N1b) dim(N3b)

  // run parameters
  const int N3b_max = 0;

  // start timing
  total_time.Start();

  // three-body table
  for (int N = 0; N <= N3b_max; ++N)
    for (HalfInt J = HalfInt(1,2); J <= Shelljmax(N) + Shelljmax(N) + Shelljmax(N); ++J)
      {
	std::cout << std::setw(2) << N 
		  << " " << std::setw(4) << DValue(J)
		  << " " << std::setw(10) << ThreeBodyDimensionOneBodyCutAMCentric(N, J) 
		  << " " << std::setw(10) << ThreeBodyDimensionOneBodyCut(N, J) 
		  << " " << std::setw(10) << ThreeBodyDimensionThreeBodyCut(N, J) 
		  << std::endl;
      }


  // end timing
  total_time.Stop();
  std::cout << std::endl;
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////
	
  // termination
  return 0;
}
