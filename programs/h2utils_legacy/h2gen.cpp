/******************************************************************************

  h2gen.cpp -- generation of r^2, r1.r2, k^2, and k1.k2 TBMEs in MFDn
  H2 format

  Syntax: h2gen

  Standard input: 

    N1b N2b
    radial_basename_p radial_basename_n 
    beta_p beta_n
    output_basename

  Mark A. Caprio
  University of Notre Dame

  3/13/12 (mac): Originated.  
  4/4/12 (mac): Phase patch.  
  1/30/13 (mac): Updated to mfdn_h2 class-based I/O.  
  2/14/13 (mac): Updated to pn bases.  
  4/25/15 (mac): Reformat source file.
  10/11/16 (mac): Integrate into shell project.

******************************************************************************/


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "am/wigner_gsl.h"
#include "legacy/shell_radial_nl.h"
#include "mcutils/profiling.h"
#include "mcutils/parsing.h"
//#include "tbme/h2_io.h"
#include "legacy/mfdn_h2.h"

////////////////////////////////////////////////////////////////
// parameter input from stdin
////////////////////////////////////////////////////////////////

// RunParameters -- stores simple parameters for run

struct RunParameters
{
  // output truncation
  int N1b, N2b;

  // radial bases
  std::string radial_basename_p, radial_basename_n;
  double beta_p, beta_n;

  // output file
  std::string output_basename; 
};


void ReadInput(RunParameters& run_parameters)
{
  // debugging note: input line isstringstream is localized on
  // each read as workaround, since .str() fails to properly
  // reset the input buffer to the new string (gcc 4.3.4)

	
  std::string line;

  {
    getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> run_parameters.N1b >> run_parameters.N2b;
    ParsingCheck (line_stream, 1, line);
  }

  // radial matrix element file basename
  {
    getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> run_parameters.radial_basename_p >> run_parameters.radial_basename_n;
    ParsingCheck (line_stream, 2, line);
  }

  // radial matrix element dilations
  {
    getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> run_parameters.beta_p >> run_parameters.beta_n;
    ParsingCheck (line_stream, 3, line);
  }

  // output matrix file basename
  {
    getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> run_parameters.output_basename;
    ParsingCheck (line_stream, 4, line);
  }
}

////////////////////////////////////////////////////////////////
// radial matrix setup
////////////////////////////////////////////////////////////////

// R2K2RadialMatrices -- container structure for all radial matrices

struct R2K2RadialMatrices 
{
  // RadialMatrices r2, r1, k2, k1;
  legacy::RadialMatrices r2_p, r1_p, k2_p, k1_p;
  legacy::RadialMatrices r2_n, r1_n, k2_n, k1_n;

  void Initialize (const std::string& radial_basename_p, const std::string& radial_basename_n);
};

// initializer -- read all radial matrices, given basenames

void R2K2RadialMatrices::Initialize (const std::string& radial_basename_p, const std::string& radial_basename_n) 
{

  // read proton radial matrices

  r2_p.Initialize(radial_basename_p + "-r2.dat");
  r1_p.Initialize(radial_basename_p + "-r1.dat");
  k2_p.Initialize(radial_basename_p + "-k2.dat");
  k1_p.Initialize(radial_basename_p + "-k1.dat");

  // read neutron radial matrices
  r2_n.Initialize(radial_basename_n + "-r2.dat");
  r1_n.Initialize(radial_basename_n + "-r1.dat");
  k2_n.Initialize(radial_basename_n + "-k2.dat");
  k1_n.Initialize(radial_basename_n + "-k1.dat");
}

////////////////////////////////////////////////////////////////
// basis and matrix setup
////////////////////////////////////////////////////////////////

// OutputMatrices -- container structure for all output two-body matrices, 
//   and associated streams, including a copy of their basis

struct OutputMatrices 
{
  legacy::TwoBodyBasisNljTzJP basis; 
  legacy::TwoBodyMatrixNljTzJP r2, r1r2, k2, k1k2;
  legacy::OutMFDnH2Stream os_r2, os_r1r2, os_k2, os_k1k2;

  void Initialize (int N1b, int N2b, const std::string& output_basename);
};

// initializer -- initialize basis, matrices, and output streams

void OutputMatrices::Initialize (int N1b, int N2b, const std::string& output_basename)
{
  // basis initialization
  basis.InitializeN1bN2b(N1b, N2b);
	
  // output matrix definition
  r2.SetBasis(basis);
  r1r2.SetBasis(basis);
  k2.SetBasis(basis);
  k1k2.SetBasis(basis);

  // output file header 
  //   common to all output files

  legacy::MFDnH2Header output_header;
  output_header.Initialize(basis);
  os_r2.Open(output_basename + "-r2.bin",output_header);
  os_r2.PrintDiagnostic();

  // open r1r2
  os_r1r2.Open(output_basename + "-r1r2.bin",output_header);
  os_r1r2.PrintDiagnostic();

  // open k2
  os_k2.Open(output_basename + "-k2.bin",output_header);
  os_k2.PrintDiagnostic();

  // open k1k2
  os_k1k2.Open(output_basename + "-k1k2.bin",output_header);
  os_k1k2.PrintDiagnostic();
}

////////////////////////////////////////////////////////////////
// sector generation
////////////////////////////////////////////////////////////////

// ShellScalarOBME gives <b|t|a>

double ShellScalarOBME (const legacy::RadialMatrices& radial_matrix, const legacy::SPOrbitalNlj& b, const legacy::SPOrbitalNlj& a)
{

  int na = a.Getn();
  int nb = b.Getn();
  int la = a.Getl();
  int lb = b.Getl();
  HalfInt ja = a.Getj();
  HalfInt jb = b.Getj();

  double matrix_element = 0.;

  if ( (lb == la) && (jb == ja) )
    matrix_element += radial_matrix.GetMatrixElement(lb, la, nb, na);

  return matrix_element;
}

// ShellScalarOBMEUpgrade uses <b|t|a> to compute <cd|V_T|ab>_AS 
//    for pp/nn states, or (cd|V_T|ab)_pn for pn states
//    applies radial basis dilations

double ShellScalarOBMEUpgrade (const legacy::RadialMatrices& radial_matrix_p, const legacy::RadialMatrices& radial_matrix_n,
			       double scale_p, double scale_n, bool momentum_space, 
			       legacy::TwoSpeciesStateType state_type, int J, const legacy::TwoBodyStateNlj& s2, const legacy::TwoBodyStateNlj& s1)
{

  legacy::SPOrbitalNlj a = s1.a1;
  legacy::SPOrbitalNlj b = s1.a2;
  legacy::SPOrbitalNlj c = s2.a1;
  legacy::SPOrbitalNlj d = s2.a2;

  double matrix_element = 0.;

  if (state_type == legacy::kPP)
    {
      if (d == b)
	matrix_element += scale_p * ShellScalarOBME(radial_matrix_p, c, a);
      if (c == a)
	matrix_element += scale_p * ShellScalarOBME(radial_matrix_p, d, b);
      int phase = - ParitySign(J - a.Getj() - b.Getj());
      if (d == a)
	matrix_element += phase * scale_p * ShellScalarOBME(radial_matrix_p, c, b);
      if (c == b)
	matrix_element += phase * scale_p * ShellScalarOBME(radial_matrix_p, d, a);
    }
  else if (state_type == legacy::kNN)
    {
      if (d == b)
	matrix_element += scale_n * ShellScalarOBME(radial_matrix_n, c, a);
      if (c == a)
	matrix_element += scale_n * ShellScalarOBME(radial_matrix_n, d, b);
      int phase = - ParitySign(J - a.Getj() - b.Getj());
      if (d == a)
	matrix_element += phase * scale_n * ShellScalarOBME(radial_matrix_n, c, b);
      if (c == b)
	matrix_element += phase * scale_n * ShellScalarOBME(radial_matrix_n, d, a);
    }
  else if (state_type == legacy::kPN)
    {
      if (d == b)
	matrix_element += scale_p * ShellScalarOBME(radial_matrix_p, c, a);
      if (c == a)
	matrix_element += scale_n * ShellScalarOBME(radial_matrix_n, d, b);
    }

  int momentum_sign = +1;
  return momentum_sign * matrix_element;
	
}

// ShellVectorOBRME gives <b||T||a>

double ShellVectorOBRME (const legacy::RadialMatrices& radial_matrix, const legacy::SPOrbitalNlj& b, const legacy::SPOrbitalNlj& a)
{

  int na = a.Getn();
  int nb = b.Getn();
  int la = a.Getl();
  int lb = b.Getl();
  HalfInt ja = a.Getj();
  HalfInt jb = b.Getj();

  double matrix_element = 0.;

  if ( am::AllowedTriangle(ja,1,jb) && ((la+lb+1)%2==0) )
    {
      // std::cout << "in OBRME " << ParitySign(jb-ja+1) << " "
      //      << Hat(ja) << " " << am::ClebschGordan(ja,HalfInt(1,2),1,0,jb,HalfInt(1,2))
      //      << " " << radial_matrix.GetMatrixElement(lb, la, nb, na) << std::endl;
      matrix_element += ParitySign(jb-ja+1) * Hat(ja) * am::ClebschGordan(ja,HalfInt(1,2),1,0,jb,HalfInt(1,2))
	* radial_matrix.GetMatrixElement(lb, la, nb, na);
    }

  return matrix_element;
}

// ShellVectorDotTBMEProduct uses <a||T||b> to compute (cd|T.T|ab) for direct product state
//    Note: The species to be used for "particle 1" and "particle 2" are 
//    determined by the calling function, which chooses the radial matrices.

double ShellVectorDotTBMEProduct (const legacy::RadialMatrices& radial_matrix_1, const legacy::RadialMatrices& radial_matrix_2, 
				  int J, 
				  const legacy::SPOrbitalNlj& c, const legacy::SPOrbitalNlj& d, const legacy::SPOrbitalNlj& a, const legacy::SPOrbitalNlj& b)
{

  double matrix_element = 0.;

  if ( am::AllowedTriangle(a.Getj(),c.Getj(),1) && am::AllowedTriangle(b.Getj(),d.Getj(),1))
    {
      matrix_element += ParitySign(d.Getj() + a.Getj() + J) * am::Wigner6J(c.Getj(),d.Getj(),J,b.Getj(),a.Getj(),1)
	* ShellVectorOBRME(radial_matrix_1,c,a) * ShellVectorOBRME(radial_matrix_2,d,b);
    }
	
  return matrix_element;
}

// ShellVectorDotTBMEUpgrade uses (cd|T.T|ab) to compute <cd|T.T|ab>_AS 
//   for pp/nn states, or passes through (cd|T.T|ab)_pn for pn states
//   applies radial basis dilations

double ShellVectorDotTBMEUpgrade (const legacy::RadialMatrices& radial_matrix_p, const legacy::RadialMatrices& radial_matrix_n, 
				  double scale_p, double scale_n, bool momentum_space, 
				  legacy::TwoSpeciesStateType state_type, int J, 
				  const legacy::TwoBodyStateNlj& s2, const legacy::TwoBodyStateNlj& s1)
{

  legacy::SPOrbitalNlj a = s1.a1;
  legacy::SPOrbitalNlj b = s1.a2;
  legacy::SPOrbitalNlj c = s2.a1;
  legacy::SPOrbitalNlj d = s2.a2;

  // bool debug = (extra_sign == +1);

  double matrix_element;

  // if (debug)
  // {
  // 	std::cout << c << d << a << b << J << " " << state_type << std::endl;
  // }

  if (state_type == legacy::kPP)
    {
      matrix_element = scale_p * scale_p * ShellVectorDotTBMEProduct (radial_matrix_p, radial_matrix_p, J, c, d, a, b);
      // DEBUG: std::cout << "direct " << matrix_element << std::endl;
      int phase = - ParitySign(J - a.Getj() - b.Getj());
      matrix_element += phase * scale_p * scale_p * ShellVectorDotTBMEProduct (radial_matrix_p, radial_matrix_p, J, c, d, b, a);
      // DEBUG: std::cout << "total " << matrix_element << std::endl;
    }
  else if (state_type == legacy::kNN)
    {
      matrix_element = scale_n * scale_n * ShellVectorDotTBMEProduct (radial_matrix_n, radial_matrix_n, J, c, d, a, b);
      int phase = - ParitySign(J - a.Getj() - b.Getj());
      matrix_element += phase * scale_n *scale_n * ShellVectorDotTBMEProduct (radial_matrix_n, radial_matrix_n, J, c, d, b, a);
    }
  else if (state_type == legacy::kPN)
    {
      matrix_element = scale_p * scale_n * ShellVectorDotTBMEProduct (radial_matrix_p, radial_matrix_n, J, c, d, a, b);
    }

  int momentum_sign = +1;
  if (momentum_space)
    momentum_sign *= ParitySign( ( c.Getl() + d.Getl() - a.Getl() - b.Getl() )/2 );

  return momentum_sign * matrix_element;
}


void GeneratorSector(const R2K2RadialMatrices& r2k2_radial_matrices,
		     double beta_p, double beta_n,
		     OutputMatrices& output_matrices, 
		     const legacy::SectorNljTzJP& sector
		     )
{
  // std::cout << "sector " << sector.GetStateType() << " " << sector.GetJ() << " " << sector.GetGrade() << ": " << output_matrices.r2.GetDimension(sector) << std::endl;
  std::cout << "." << std::flush;

  // allocate storage for sector
  output_matrices.r2.Initialize(sector);
  output_matrices.r1r2.Initialize(sector);
  output_matrices.k2.Initialize(sector);
  output_matrices.k1k2.Initialize(sector);


  // recover sector properties
  const legacy::TwoSpeciesStateType state_type = sector.GetStateType();
  const int J = sector.GetJ();
  const int g = sector.GetGrade();
  const int dimension = output_matrices.basis.GetDimension(state_type,J,g);

  // for canonical pairs of states in destination two-body space (in lexicographical order)
  for (int k1 = 0; k1 < dimension; ++k1)
    for (int k2 = k1; k2 < dimension; ++k2)
      {
	// identify states
	legacy::TwoBodyStateNlj s1 = output_matrices.basis.GetState(state_type,J,g,k1);
	legacy::TwoBodyStateNlj s2 = output_matrices.basis.GetState(state_type,J,g,k2);
			
	// calculate matrix elements
	const bool kCoordinate = false;
	const bool kMomentum = true;
	double matrix_element_r2 = ShellScalarOBMEUpgrade(
							  r2k2_radial_matrices.r2_p, r2k2_radial_matrices.r2_n, (beta_p*beta_p), (beta_n*beta_n), kCoordinate, state_type, J, s2, s1
							  );
	double matrix_element_r1r2 = ShellVectorDotTBMEUpgrade(
							       r2k2_radial_matrices.r1_p, r2k2_radial_matrices.r1_n, beta_p, beta_n, kCoordinate, state_type, J, s2, s1
							       );
	double matrix_element_k2 = ShellScalarOBMEUpgrade(
							  r2k2_radial_matrices.k2_p, r2k2_radial_matrices.k2_n, 1/(beta_p*beta_p), 1/(beta_n*beta_n), kMomentum, state_type, J, s2, s1
							  );
	double matrix_element_k1k2 = ShellVectorDotTBMEUpgrade(
							       r2k2_radial_matrices.k1_p, r2k2_radial_matrices.k1_n, 1/beta_p, 1/beta_n, kMomentum, state_type, J, s2, s1
							       );

	// store matrix elements
	output_matrices.r2.SetMatrixElementUNAS(state_type, s1, s2, matrix_element_r2); 
	output_matrices.r1r2.SetMatrixElementUNAS(state_type, s1, s2, matrix_element_r1r2); 
	output_matrices.k2.SetMatrixElementUNAS(state_type, s1, s2, matrix_element_k2); 
	output_matrices.k1k2.SetMatrixElementUNAS(state_type, s1, s2, matrix_element_k1k2); 
      }

  // output results for sector
  output_matrices.os_r2.WriteSector(output_matrices.r2);
  output_matrices.os_r1r2.WriteSector(output_matrices.r1r2);
  output_matrices.os_k2.WriteSector(output_matrices.k2);
  output_matrices.os_k1k2.WriteSector(output_matrices.k1k2);

  // free storage for sector
  output_matrices.r2.Free(sector);
  output_matrices.r1r2.Free(sector);
  output_matrices.k2.Free(sector);
  output_matrices.k1k2.Free(sector);

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
  std::cout << "h2gen" << std::endl;
  std::cout << std::endl;

  // start timing
  Timer total_time;
  total_time.Start();

  // read input parameters
  RunParameters run_parameters;
  ReadInput(run_parameters); 

  // I/O and calculation initialization

  R2K2RadialMatrices r2k2_radial_matrices;
  r2k2_radial_matrices.Initialize(run_parameters.radial_basename_p, run_parameters.radial_basename_n);

  OutputMatrices output_matrices;
  output_matrices.Initialize(run_parameters.N1b, run_parameters.N2b, run_parameters.output_basename);

  ////////////////////////////////////////////////////////////////
  // evaluation
  ////////////////////////////////////////////////////////////////

  for (legacy::SectorNljTzJP sector(output_matrices.basis); sector.InRange(); ++sector)
    {
      GeneratorSector(r2k2_radial_matrices, run_parameters.beta_p, run_parameters.beta_n, output_matrices, sector);
      std::cout << "." << std::flush;
    }
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // end timing
  total_time.Stop();
  std::cout << std::endl;
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;


  // termination
  return 0;
}