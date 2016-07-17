/******************************************************************************

  h2xform.cpp -- MFDn TBME change of radial basis

  Syntax: 

    h2xform 

  Environment:
 
    OMP_NUM_THREADS: sets OpenMP parallelization if OpenMP enabled on compilation

  Standard input:
    input_basename
    [Then lines of the form:]
       output_basename xform_basename_p xform_basename_n N1b_cut N2b_cut N1b N2b

  Outline for seqential transformation of sectors:

    for each (J,Tz,g) sector defined in maximum of input and output truncations
      if given sector is defined in input truncation
        allocate sector (input)
        read sector (input)
      for each output xform
        if given sector is defined in output truncation
          allocate sector
          transform sector (or leave as zeros if no input sector)
          write sector
          free sector (output)
      if given sector is defined in input truncation
        free sector (input)

  Created by M. A. Caprio, University of Notre Dame.
  4/18/11 (mac): Originated.
  4/25/11 (mac): Required headers updated.
  8/15/11 (mac): Output of multiple transformation files.
  8/27/11 (mac): Additional instrumentation.
  3/15/12 (mac): Update to use sector iterators.
  3/20/12 (mac): Instruction input syntax update.  Supplant shell_xform with shell_radial.
  8/11/12 (mac): Allow input file truncation smaller than output file truncation, provided Ncut 
                 is accommodated in input truncation.
  1/30/13 (mac): Update to mfdn_h2 class-based I/O.
  2/15/13 (mac): Updated to pn bases.
  12/18/14 (mac): Add OpenMP parallelization of main xform loop.
  4/25/15 (mac): Reformat source file.
  9/23/15 (mac): Fix missing if on free sector.

******************************************************************************/

#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <mcpp/profiling.h>

#include <shell/shell_radial_nl.h>
#include <shell/mfdn_h2.h>


using namespace shell;


////////////////////////////////////////////////////////////////
// radial matrix setup
////////////////////////////////////////////////////////////////

// XformRadialMatrices -- container structure for all radial matrices

struct XformRadialMatrices 
{
  RadialMatrices xform_p, xform_n;

  void Initialize (const std::string& radial_basename_p, const std::string& radial_basename_n);
};

// initializer -- read all radial matrices, given basenames

void XformRadialMatrices::Initialize (const std::string& radial_basename_p, const std::string& radial_basename_n) 
{
  xform_p.Initialize(radial_basename_p + ".dat");
  xform_n.Initialize(radial_basename_n + ".dat");
}


////////////////////////////////////////////////////////////////
// sector transformation
////////////////////////////////////////////////////////////////

void TwoBodyMatrixSectorTransform (const TwoBodyMatrixNljTzJP& source_matrix, TwoBodyMatrixNljTzJP& destination_matrix, 
				   const XformRadialMatrices& xform_radial_matrices, 
				   int N1b_cut, int N2b_cut,
				   const SectorNljTzJP& sector
				   )
{
  // recover sector properties
  const TwoSpeciesStateType state_type = sector.GetStateType();
  const int J = sector.GetJ();
  const int g = sector.GetGrade();

  // Notation: unprimed variables are for untransformed
  // (oscillator) states, primed variables are for transformed
  // states

  // extract dimension of subspace
  const int dimension = destination_matrix.GetTwoBodyBasis().GetDimension(state_type,J,g);

  // for canonical pairs of states (in lexicographical order)

  // OMP: 0uter loop parallelization is "embarassingly parallel"
  // parallelization of target matrix elements involving
  // differen bra-ket pairs.  This avoids thread overhead of
  // parallelizing the reduction operation in the inner loops.
  // Optimal approach in tests based on a simple computational
  // load (i.e., no memory access) in inner loop.

#pragma omp parallel for collapse(2) 
  for (int k1p = 0; k1p < dimension; ++k1p)
    for (int k2p = k1p; k2p < dimension; ++k2p)
      {
	// identify target matrix element
	TwoBodyStateNlj s1p = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k1p);
	TwoBodyStateNlj s2p = destination_matrix.GetTwoBodyBasis().GetState(state_type,J,g,k2p);

	// DBG: std::cout << "  Evaluating:" 
	// DBG:      << " type " << state_type << " J " << J
	// DBG:      << " " << s1p.a1 << s1p.a2 << s2p.a1 << s2p.a2 
	// DBG:      << std::endl;

	// define shorthands for n and l quantum numbers (for more readable coding below)
	int n11p = s1p.a1.Getn();
	int n12p = s1p.a2.Getn();
	int n21p = s2p.a1.Getn();
	int n22p = s2p.a2.Getn();

	int l11p = s1p.a1.Getl();
	int l12p = s1p.a2.Getl();
	int l21p = s2p.a1.Getl();
	int l22p = s2p.a2.Getl();

	// determine summation limits under N1b truncation for source matrix elements
	//    n = (N-l)/2 and rely on integer division to give floor of (N1b - l) /2
	//    but if N1b_cut < l this cannot be relied upon, since integer division truncation with 
	//    negative argument is initially implementation defined with late standard towards zero
	int n11_max = (N1b_cut >= l11p) ? (N1b_cut - l11p) / 2 : -1;
	int n12_max = (N1b_cut >= l12p) ? (N1b_cut - l12p) / 2 : -1;
	int n21_max = (N1b_cut >= l21p) ? (N1b_cut - l21p) / 2 : -1;
	int n22_max = (N1b_cut >= l22p) ? (N1b_cut - l22p) / 2 : -1;

	// carry out sum

	// Note: Sum is explicitly over radial quantum numbers, rather than over all two-body 
	// states in unprimed basis, since this allows much greater a priori restriction 
	// of the sum, by imposition of the Kroenecker delta on l and j quantum numbers.
	// (Whether or not these reduced loops are actually warranted for efficiency is not
	// thoroughly explored and may be expected to depend on basis size.)  The for loops
	// directly implement the one-body cutoff on the sum over the unprimed basis (i.e., Ncut 
	// as defined in csbasis).  But the two-body cutoff is then implemented with a short-circuit
	// text inside the loop.

	double matrix_element = 0.;
	for (int n11 = 0; n11 <= n11_max; ++n11)
	  for (int n12 = 0; n12 <= n12_max; ++n12)
	    for (int n21 = 0; n21 <= n21_max; ++n21)
	      for (int n22 = 0; n22 <= n22_max; ++n22)
		{
		  // look up orbitals for source matrix element
		  int N11 = 2*n11 + l11p;
		  int N12	= 2*n12 + l12p;
		  int N21	= 2*n21 + l21p;
		  int N22	= 2*n22 + l22p;
		  SPOrbitalNlj a11(N11,s1p.a1.Getj());
		  SPOrbitalNlj a12(N12,s1p.a2.Getj());
		  SPOrbitalNlj a21(N21,s2p.a1.Getj());
		  SPOrbitalNlj a22(N22,s2p.a2.Getj());
							
		  // build states for bracket
		  TwoBodyStateNlj s1(a11,a12,J);
		  TwoBodyStateNlj s2(a21,a22,J);
		  // DBG: std::cout << "    Adding: " << s1.a1 << s1.a2 << s2.a1 << s2.a2 << std::endl;

		  // impose 2-body cutoff on unprimed states
		  if ( 
		      ( (s1.a1.GetN() + s1.a2.GetN()) > N2b_cut)
		      ||
		      ( (s2.a1.GetN() + s2.a2.GetN()) > N2b_cut)
		       )
		    continue;

		  // skip symmetry forbidden states
		  if ( !SymmetryAllowedState(state_type,s1) ||  !SymmetryAllowedState(state_type,s2) )
		    continue;

		  // revise bracket to canonical order
		  int phase = CanonicalizeTwoBodyBracket(state_type, s1, s2, 0);
		  // DBG: std::cout << "      ====> " << s1.a1 << s1.a2 << s2.a1 << s2.a2 << std::endl;

		  // accumulate source matrix element to target matrix element

		  // Note (2013): Transformation brackets in radial xform file are stored  
		  // with unprimed (HO) basis as bra (row index) and primed (new) basis as 
		  // ket (column).
							

		  // The patern for pn bases is...
		  // 	double coefficient 
		  // 		= transformation1.GetMatrixElement(l11p,l11p,n11,n11p)
		  // 		* transformation2.GetMatrixElement(l12p,l12p,n12,n12p)
		  // 		* transformation1.GetMatrixElement(l21p,l21p,n21,n21p)
		  // 		* transformation2.GetMatrixElement(l22p,l22p,n22,n22p);

		  double coefficient;
		  if (state_type == kPP)
		    {
		      coefficient
			= xform_radial_matrices.xform_p.GetMatrixElement(l11p,l11p,n11,n11p)
			* xform_radial_matrices.xform_p.GetMatrixElement(l12p,l12p,n12,n12p)
			* xform_radial_matrices.xform_p.GetMatrixElement(l21p,l21p,n21,n21p)
			* xform_radial_matrices.xform_p.GetMatrixElement(l22p,l22p,n22,n22p);
		    }
		  else if (state_type == kNN)
		    {
		      coefficient
			= xform_radial_matrices.xform_n.GetMatrixElement(l11p,l11p,n11,n11p)
			* xform_radial_matrices.xform_n.GetMatrixElement(l12p,l12p,n12,n12p)
			* xform_radial_matrices.xform_n.GetMatrixElement(l21p,l21p,n21,n21p)
			* xform_radial_matrices.xform_n.GetMatrixElement(l22p,l22p,n22,n22p);
		    }
		  else if (state_type == kPN)
		    {
		      coefficient
			= xform_radial_matrices.xform_p.GetMatrixElement(l11p,l11p,n11,n11p)
			* xform_radial_matrices.xform_n.GetMatrixElement(l12p,l12p,n12,n12p)
			* xform_radial_matrices.xform_p.GetMatrixElement(l21p,l21p,n21,n21p)
			* xform_radial_matrices.xform_n.GetMatrixElement(l22p,l22p,n22,n22p);
		    }

		  // DBG: std::cout << "      Contribution:" 
		  // DBG:      << " phase " << phase 
		  // DBG:      << " coefficient " << coefficient 
		  // DBG:      << " NAS " << source_matrix.GetMatrixElement(state_type, s1, s2) 
		  // DBG:      << " UNAS " << source_matrix.GetMatrixElementUNAS(state_type, s1, s2)
		  // DBG:      << " total " << phase * coefficient * source_matrix.GetMatrixElementUNAS(state_type, s1, s2)
		  // DBG:      << std::endl;

		  matrix_element += phase * coefficient * source_matrix.GetMatrixElementUNAS(state_type, s1, s2); 
		}
			
	// save resultant matrix element

	// OMP: The "critical" designation for saving the matrix element is out of caution since
	// STL containers are supposedly not thread-safe.
#pragma omp critical
	destination_matrix.SetMatrixElementUNAS(state_type, s1p, s2p, matrix_element); 
			
	// DBG: std::cout << "    Resultant:" 
	// DBG:      << " calculated UNAS " << matrix_element 
	// DBG:      << " readback NAS " << destination_matrix.GetMatrixElement(state_type, s1p, s2p) 
	// DBG:      << std::endl;
			
      }
}

////////////////////////////////////////////////////////////////
// transformation control
////////////////////////////////////////////////////////////////

// Transformation descriptor holds all data necessary for a transformation 
// plus the output stream for the transformation


class TransformationDescriptor
{
public:
  TransformationDescriptor (const std::string& os_basename, 
			    const std::string& radial_basename_p, const std::string& radial_basename_n, 
			    int N1b_cut, int N2b_cut, 
			    int N1b, int N2b);
  void Initialize (int N1b_in, int N2b_in);
  void TransformSector (const TwoBodyMatrixNljTzJP& input_matrix, const SectorNljTzJP& sector);

  // truncation definition retrieval
  int GetN1b () const
  {return N1b_;};

  int GetN2b () const
  {return N2b_;};


private:
  // source stream information
  MFDnH2Header is_header_;

  // parameters
  std::string os_basename_;
  OutMFDnH2Stream os_;
	
  std::string radial_basename_p_, radial_basename_n_;
  int N1b_cut_, N2b_cut_;
  int N1b_, N2b_;

  // data
  XformRadialMatrices transformations_;
  TwoBodyBasisNljTzJP output_basis_;
  TwoBodyMatrixNljTzJP output_matrix_;

};

TransformationDescriptor::TransformationDescriptor(const std::string& os_basename, 
						   const std::string& radial_basename_p, const std::string& radial_basename_n, 
						   int N1b_cut, int N2b_cut, 
						   int N1b, int N2b)
{

  // Note: is_header_ must be defined later

  // set parameters
  os_basename_ = os_basename;
  radial_basename_p_ = radial_basename_p;
  radial_basename_n_ = radial_basename_n;
  N1b_cut_ = N1b_cut;
  N2b_cut_ = N2b_cut;
  N1b_ = N1b;
  N2b_ = N2b;

}

void TransformationDescriptor::Initialize(int N1b_in, int N2b_in)
{
  // validation of cutoff
  if ( ( N1b_in < N1b_cut_ ) || ( N2b_in < N2b_cut_ ) )
    {
      std::cerr << "input file truncation does not support requested oscillator summation cut " << std::endl;
      std::exit(EXIT_FAILURE);
    }

  // read transformation coefficients
  transformations_.Initialize(radial_basename_p_, radial_basename_n_);

  // transformed matrix setup
  output_basis_.InitializeN1bN2b(N1b_, N2b_);
  output_matrix_.SetBasis(output_basis_);

  // initialize output stream
  MFDnH2Header output_header;
  output_header.Initialize(output_basis_);
  os_.Open(os_basename_+".bin", output_header);
  os_.PrintDiagnostic();
}

void TransformationDescriptor::TransformSector(const TwoBodyMatrixNljTzJP& input_matrix, const SectorNljTzJP& sector)
{

  // short circuit for nonexistent output sector
  if (!output_matrix_.HasSector(sector))
    return;

  // allocate memory for sector
  output_matrix_.Initialize(sector);

  // generate matrix elements (if input exisits)
  if (input_matrix.HasSector(sector))
    TwoBodyMatrixSectorTransform(input_matrix, output_matrix_, transformations_, N1b_cut_, N2b_cut_, sector);

  // write sector
  os_.WriteSector(output_matrix_);

  // free memory for sector
  output_matrix_.Free(sector);
}

////////////////////////////////////////////////////////////////
// parameter input from stdin
////////////////////////////////////////////////////////////////

void ReadParameters(std::string& is_basename, std::vector<TransformationDescriptor>& transformations)
{

  std::string line;
  int line_count = 0;

  {
    ++line_count;
    std::getline(std::cin, line);
    std::istringstream line_stream(line);
    line_stream >> is_basename;
    ParsingCheck (line_stream, line_count, line);
  }


  while ( std::getline(std::cin, line) )
    {
      // count line
      ++line_count;

      // set up for parsing
      std::istringstream line_stream(line);

      // read transformation parameters
      std::string os_basename; 
      std::string radial_basename_p, radial_basename_n;
      int N1b_cut, N2b_cut;
      int N1b, N2b;
      std::istringstream ls(line);
      line_stream >> os_basename >> radial_basename_p >> radial_basename_n >> N1b_cut >> N2b_cut >> N1b >> N2b;
      ParsingCheck (line_stream, line_count, line);

      // generate transformation descriptor
      transformations.push_back(TransformationDescriptor(os_basename, radial_basename_p, radial_basename_n, N1b_cut, N2b_cut, N1b, N2b));

    }
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
  std::cout << "h2xform" << std::endl;
  std::cout << std::endl;

  // start timing
  Timer total_time;
  total_time.Start();

  ////////////////////////////////////////////////////////////////
  // transformation parameter setup 
  ////////////////////////////////////////////////////////////////
	
  // read transformation parameters
  std::string is_basename;
  std::vector<TransformationDescriptor> transformations;
  ReadParameters(is_basename, transformations);

  // initialize master sector iteration truncation
  int N1b_master(0), N2b_master(0);
  // update master truncation to encompass output truncations
  for (std::vector<TransformationDescriptor>::iterator it = transformations.begin(); it != transformations.end(); ++it)
    {
      N1b_master = std::max( N1b_master, it->GetN1b() );
      N2b_master = std::max( N2b_master, it->GetN2b() );
    }
		

  std::cout << "Maximum output truncation" << " N1b " << N1b_master << " N2b " << N2b_master << std::endl;

  ////////////////////////////////////////////////////////////////
  // input initialization
  ////////////////////////////////////////////////////////////////

  // input stream initialization
  InMFDnH2Stream is;
  MFDnH2Header is_header;
  is.Open (is_basename + ".bin", is_header);
  is.PrintDiagnostic ();

  // input basis/matrix initialization
  TwoBodyBasisNljTzJP input_basis;
  input_basis.InitializeN1bN2b(is_header.N1b, is_header.N2b);
  TwoBodyMatrixNljTzJP input_matrix(input_basis);

  // process input stream truncation information
  int N1b_in = is_header.N1b;
  int N2b_in = is_header.N2b;
  N1b_master = std::max( N1b_master, N1b_in );
  N2b_master = std::max( N2b_master, N2b_in );

  ////////////////////////////////////////////////////////////////
  // output initialization
  ////////////////////////////////////////////////////////////////

  for (std::vector<TransformationDescriptor>::iterator it = transformations.begin(); it != transformations.end(); ++it)
    it->Initialize(N1b_in, N2b_in);

  ////////////////////////////////////////////////////////////////
  // sector-by-sector transformation
  ////////////////////////////////////////////////////////////////


  for (SectorNljTzJP sector(N1b_master,N2b_master); sector.InRange(); ++sector)
    {
      if (input_matrix.HasSector(sector))
	{
	  input_matrix.Initialize(sector);
	  is.ReadSector(input_matrix);
	}
      for (std::vector<TransformationDescriptor>::iterator it = transformations.begin(); it != transformations.end(); ++it)
	it->TransformSector(input_matrix, sector);
      if (input_matrix.HasSector(sector))
	input_matrix.Free(sector);
      std::cout << "." << std::flush;
    }
	
  std::cout << std::endl;

  ////////////////////////////////////////////////////////////////
  // termination
  ////////////////////////////////////////////////////////////////

  // end timing
  total_time.Stop();
  std::cout << "(Total time: " << total_time.ElapsedTime() << ")" << std::endl;
  std::cout << std::endl;


  // termination
  return 0;
}
