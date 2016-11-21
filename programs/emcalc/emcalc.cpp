/*******************emcalc.cpp****************************

emcalc.cpp was first written for the calculation of the electromagnetic multiple moments
using the one body densities of medium heavy nuclei (A~40) obtained with the generalized seniority 
model. 

The code was later adapted to work with the one body densities obtained with the no-core shell model.
In the first version of the program the uncoupled one body densities were used. The uncloupled densities 
had the form <\Psi| c_{a,ma,mta}^{dagger} c_{b,mb,mtb} | \Psi>, where a=(nalajata), ma is the z projection of ja, and
mta is the z projection of isospin ta.

This format with the uncoupled densities was changed in the most recent versions of mfdn (version 14 beta 06 
as far as I know). The output is now the coupled one body densities <\Psi || c_{a}^{dagger} c_{b}^{\tilde} || \Psi>,
thus the calculation is simplified a little bit, because a clebsch gordan coefficient is not needed in the sum anymore, 
and the number of elements that need to be summed is reduced because we don't need to sum over m anymore.

In the old format the code was reading in the orbital information file mfdn.obdme.info, which contained all the information 
for the creation and annihilation operators c_{a,ma,mta}^{dagger} c_{b,mb,mtb}, and either the file containing the static one body 
densities mfdn.statobdme.seqXXX.2JXX.pX.nXX.2TXX or the transition densities mfdn.obdme.seqXXX.2JXX.pX.nXX.2TXX.seqXXX.2JXX.pX.nXX.2TXX.

  

A program that calculates the reduced matrix elements for E/M transitions                    
for both protons (tz=1) and neutrons(tz=-1). For the evaluation of the                
radial integrals (R_ab) given in Suhonen 6.32 an input file is used. The program needs           
to be given the paths to the orbital information file mfdn.obdme.info, the path to 
the obdme file called either mfdn.obdme.seq002.2J03.p0.n01.2T01.seq004.2J01.p0.n01.2T01
or mfdn.statobdme.seq006.2J05.p1.n01.2T01 and it also needs the path to the radial integral
file radial-me-BASIS-r*.dat, where BASIS is a label indicating the type of basis function, such as
HO (Harmonic Oscillator) or CSFN (Coulomb Sturmian with first node) and * is 0 1 or 2. The program is 
fully updated from an older version written in 05/23/13. Major changes include a change to the 
input and output files as well the inclusion of several functions that process the input file 
internally to produce the result. The user only needs to modify the input file as he/she wishes
to get the desired result. The input and output files are described below.     
                                                                                                      
Valentino Constantinou 07/31/2013         

Parity included in output 08/13/13          

M2 calculation added on 11/23/13   

Addition on 01/26/14: There was a problem with emcalc not working for the 9Be runs. The problem 
was located in the function InitializeObdmeHeaderInformation(string&, ObdmeHeaderInformation&) that
reads the first 2 lines of the file containing the one body densities, in other words from files like
mfdn.obdme.seq012.2J05.p1.n01.2T01.seq010.2J03.p0.n03.2T01. The problem was that a stringstream was used
twice in the same function without being cleared after the first usage. Therefore the values of ji, jf, pi,
pf, etc were read wrong. The problem was resolved with the use of the function clear() that clears the 
stringstream.      

04/03/15: The intension is to upgrade the code to work with the new format of reduced one body densities. The old code
          worked with the uncoupled one densities. Also the function NumberOperator() will be added to conveniently calculate
          the sum of the one body densities..                                
**************************************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
//#include "radial_integrals.h"
#include "shell_radial.h"
#include "halfint.h"
#include "wigner_gsl.h"
#include "angular_momentum.h"
#include "mfdn_file_processing.h"

using namespace std;  
  

double number_operator(PostProcessingParameters& input_file_parameters)
{ 
   double protons = 0, neutrons = 0;
   vector <OrbitalInformation> orbitals;
   vector <double> densities;

   ReadOrbitalInformationFile(orbitals, input_file_parameters.path_to_orbital_information_file);  
   ReadOneBodyDensitiesFile(densities,input_file_parameters.path_to_obdme_file);

   int half_size = (densities.size())/2;

   //proton sector

   for(int i=0; i < half_size; i++)
      {
        if( (orbitals[i].na == orbitals[i].nb) && (orbitals[i].la == orbitals[i].lb) && 
            (orbitals[i].ja == orbitals[i].jb) && (orbitals[i].k == 0)  )
           {
              protons = protons + Hat(orbitals[i].ja) * densities[i];
           }
      }

   //neutron sector

   for(int i=half_size; i < 2*half_size; i++)
      {
        if( (orbitals[i-half_size].na == orbitals[i-half_size].nb) && (orbitals[i-half_size].la == orbitals[i-half_size].lb) && 
            (orbitals[i-half_size].ja == orbitals[i-half_size].jb) && (orbitals[i-half_size].k == 0)  )
           {
              neutrons = neutrons + Hat(orbitals[i-half_size].ja) * densities[i];
           }
      }  
     
      //cout << "The sum is: " << protons + neutrons << endl;

      return (protons+neutrons);
}                                                             
                                                                                           
/*************RadialIntegrals ReturnRadialIntegrals (string& radial_integrals_partial_path, char& type,
				       int& lambda)

It returns a RadialIntegrals class. It takes a string which is the full path to the radial integrals file
a char with the type of transition and lambda. It then internally constructs the full path, opens the file
and initialize the RadialIntegral class which is used by the program. An obvious extension in the future 
is from the program to take the type of basis (HO or CS) and construct the appropriate RadialIntegral class

*********************************************************************************************************/

shell::RadialMatrices ReturnRadialIntegrals (string& radial_integrals_partial_path, char& type,
				       int& lambda)
{
  string transition_type;                     // string which will be added in proton_integrals
  ostringstream convert_lambda_to_string;     // stream used for the conversion of integer lambda to char
  string radial_integrals_full_path;
  shell::RadialMatrices radial_matrix;

  if(type=='E')
    {
      convert_lambda_to_string << lambda;      
      transition_type = convert_lambda_to_string.str(); 

      radial_integrals_full_path = radial_integrals_partial_path + "-r" + transition_type + ".dat";
    }
  else if(type=='M')
    {
      convert_lambda_to_string << (lambda-1);      
      transition_type = convert_lambda_to_string.str(); 

      radial_integrals_full_path = radial_integrals_partial_path + "-r" + transition_type + ".dat";
    }
  else if ((type != 'E') || (type !='M'))
    {
      cerr << "Not an accepted transition. Please check the first line of your input file \n"
	"The first letter should be M or E" << endl;
      exit(EXIT_FAILURE);
    }    
  
  ifstream open_integrals_file(radial_integrals_full_path.c_str(),ios::in); 

  if( !open_integrals_file )  
    {
      cerr << "The radial integrals file cannon be opened. Please check the file emcalc.in \n "
	"and make sure that the path is correct" << endl;
      exit(EXIT_FAILURE);
    }
  else
    {      
      radial_matrix.Initialize(open_integrals_file);     
    }

  return radial_matrix;
}

/*************double OscillatorLength(double& scale_length, double hw)**************

A function that calculates the oscillator length which can be different for protons
and neutrons 

***********************************************************************************/

double OscillatorLength(double& scale_length, double& hw)
{
  double hbar_c = 197.327;                         //MeV*fm
  double mass_nucleon_c = 938.92;                  //MeV, nucleon mass
  double pi = 3.141593;
  double oscillator_length = scale_length * (hbar_c/sqrt(mass_nucleon_c * hw));

  return oscillator_length;
}

/*******CalculateElectromagneticMatrixElements(PostProcessingParameters& input_file_parameters)     
                                                                                                         
The function takes a structure PostProcessingParameters as an input.
The structure contains all the variables described above, i.e paths to input
files, parameters etc.

The function sends its output to the file emcalc.out which has             
the following format:                                                                               
                                                                                                         
If the transition is electric                                                                           
                                                                                                         
Elambda Ji  pi  ni  Jf  pf  nf  ME_p  ME_n  RME_p  RME_n 

where ME_p is the matrix element for protons, RME_p is the reduced matrix element for protons,
ni/nf is the state label for the initial/final state and pi/pf is the parity. For example if there
are two initial states with J=3/2, the first has a label ni=1 and the second has a label ni = 2                                                               
                                                                                                    
If the transition is magnetic                                                                            
                                                                                                         
Mlambda Ji  pi  ni  Jf  pf  nf dl_p dl_n ds_p ds_n ME_total RME_total                                     
                                                                                                         
The definition of dl_p, dl_n etc is given in page 128 Suhonen. ME_total is the total
matrix element for both protons and neutrons and RME_total is the total reduced matrix element.

If the Clebsch Gordan coefficient needed to calculate the reduced matrix elements is zero,
the reduced matrix element is not printed (an empty space "" is printed) because otherwise the division with 
zero will print an inf as the output!
                                                             
*******************************************************************************************************/ 

void CalculateElectromagneticMatrixElements (PostProcessingParameters& input_file_parameters)
{  
  FILE * output_file;
  output_file = fopen("emcalc.out","w");
  int lambda_value = input_file_parameters.lambda;       

  /*-----------Debug on 01/26/14------------

  cout << "Lambda from input file: " << input_file_parameters.lambda << endl;
  cout << "Path to obdme file from input file: " << input_file_parameters.path_to_obdme_file << endl;
  cout << "Observable type from input file: " << input_file_parameters.observable_type << endl;*/		                    

  // (04/03/15) The struct header_info is passed to InitializeObdmeHeaderInformation along with the path
  // to the density file and so it gets initialized and returned. Then by acting on header_info we can
  // extract the information we want.. 


  ObdmeHeaderInformation header_info;

  InitializeObdmeHeaderInformation (input_file_parameters.path_to_obdme_file, header_info); 

  /*-----------Debug on 14/01/26------------

  cout << "Header Info contents " << endl << endl;

  cout << "magic number: " << header_info.magic_number << " number of elements: " << header_info.number_of_elements << " parity_i: " << header_info.parity_i << " parity_f " << header_info.parity_f <<endl;

  cout << "ji: " << header_info.ji << " jf: " << header_info.jf << " ti: " << header_info.ti << " tf: " <<
       header_info.tf << " mm: " << header_info.mm_basis << endl;*/ 

  
  if( input_file_parameters.observable_type == 'E' &&
      ( (header_info.parity_i == header_info.parity_f && lambda_value%2 ==1) || 
	(header_info.parity_i != header_info.parity_f && lambda_value%2 ==0)  ))
    { 
      fprintf(stdout,"This is a transition that is not allowed because the parity condition is not satisfied \n");	
    
      fprintf(output_file," ");	       
      
      fclose(output_file);      
    }
  else if( input_file_parameters.observable_type == 'M' &&
	   ( (header_info.parity_i == header_info.parity_f && lambda_value%2 ==0) ||
	     (header_info.parity_i != header_info.parity_f && lambda_value%2 ==1) ) )
    {
      fprintf(stdout,"This is a transition that is not allowed because the parity condition is not satisfied \n");	
      
      fprintf(output_file," ");	       
      
      fclose(output_file);      
    }    
  else if( ( input_file_parameters.observable_type == 'M' || input_file_parameters.observable_type == 'E' ) 
	   && AllowedTriangle(header_info.ji,header_info.jf,HalfInt(lambda_value)) ==0 )
    {
      fprintf(stdout,"This is a transition that is not allowed because the triangular condition between Ji"
	      ", Jf and lambda is not satisfied\n");	
    
      fprintf(output_file," ");	       
      
      fclose(output_file);      
    }
  else if ( AllowedTriangle(header_info.ji,header_info.jf,HalfInt(lambda_value)) == 0 || 
	    TwiceValue(header_info.ji-header_info.mm_basis) %2 !=0 ||
	    TwiceValue(header_info.jf-header_info.mm_basis) %2 !=0 )
    {
      fprintf(stdout, "There is a corrupted J value in the file which contains the obdme values. "
	      "The J and T values are calculated as doubles by mfdn\n");    

      fprintf(output_file," ");	       
      
      fclose(output_file);      
    }  
  else if( (input_file_parameters.observable_type=='E') && 
	   (AllowedTriangle(header_info.ji,header_info.jf,HalfInt(lambda_value)) == 1) )
    {
      vector <OrbitalInformation> orbital_file;

      ReadOrbitalInformationFile (orbital_file, input_file_parameters.path_to_orbital_information_file);



      shell::RadialMatrices radial_matrix_protons = ReturnRadialIntegrals (input_file_parameters.path_to_radial_integrals_file_protons,
								     input_file_parameters.observable_type,
								     input_file_parameters.lambda);

      shell::RadialMatrices radial_matrix_neutrons = ReturnRadialIntegrals (input_file_parameters.path_to_radial_integrals_file_neutrons,
								     input_file_parameters.observable_type,
								     input_file_parameters.lambda);





      vector <double> one_body_densities;
      ReadOneBodyDensitiesFile (one_body_densities, input_file_parameters.path_to_obdme_file);

      double boscil_p = OscillatorLength(input_file_parameters.beta_p,input_file_parameters.hw);  
      double boscil_n = OscillatorLength(input_file_parameters.beta_n,input_file_parameters.hw); 
      double pi = 3.141593;

      double reduced_many_body_matrix_element_protons = 0, reduced_many_body_matrix_element_neutrons = 0;

      double single_particle_me_coefficient, single_particle_me_p, single_particle_me_n, radial_integral_p, radial_integral_n;  
    
      int file_size = header_info.number_of_elements;

      // protons sector

      for (int i=0; i < (file_size/2); i++)
	{
	  if ( ((orbital_file[i].la+orbital_file[i].lb+lambda_value)%2==0) && 
	       (AllowedTriangle(orbital_file[i].ja,orbital_file[i].jb,HalfInt(lambda_value)) == 1) && 
               (orbital_file[i].k == lambda_value) 
             ) 
	    {    
	      radial_integral_p = radial_matrix_protons.GetMatrixElement(orbital_file[i].la,orbital_file[i].lb,orbital_file[i].na,orbital_file[i].nb);	     

	      single_particle_me_coefficient = (ParitySign(orbital_file[i].ja+HalfInt(lambda_value)-HalfInt(1,2))/sqrt(4*pi)) * Hat(orbital_file[i].jb) 
		* Hat(orbital_file[i].ja) * ClebschGordan(orbital_file[i].ja,HalfInt(1,2),orbital_file[i].jb,-HalfInt(1,2),HalfInt(lambda_value),HalfInt(0));

	      single_particle_me_p = single_particle_me_coefficient*radial_integral_p;
	      
	      reduced_many_body_matrix_element_protons = reduced_many_body_matrix_element_protons + pow(boscil_p,lambda_value) * single_particle_me_p * (one_body_densities[i]);
            }
	}	
      for (int i = (file_size/2); i < file_size; i++)
        {
	  if ( ((orbital_file[i-(file_size/2)].la+orbital_file[i-(file_size/2)].lb+lambda_value)%2==0) && 
	       (AllowedTriangle(orbital_file[i-(file_size/2)].ja,orbital_file[i-(file_size/2)].jb,HalfInt(lambda_value)) == 1) && (orbital_file[i-(file_size/2)].k == lambda_value) ) 
	    {             
	      radial_integral_n = radial_matrix_neutrons.GetMatrixElement(orbital_file[i-(file_size/2)].la,orbital_file[i-(file_size/2)].lb,orbital_file[i-(file_size/2)].na,orbital_file[i-(file_size/2)].nb);     
                
	      single_particle_me_coefficient= (ParitySign(orbital_file[i-(file_size/2)].ja+HalfInt(lambda_value)-HalfInt(1,2))/sqrt(4*pi)) * Hat(orbital_file[i-(file_size/2)].jb) 
		* Hat(orbital_file[i-(file_size/2)].ja) * ClebschGordan(orbital_file[i-(file_size/2)].ja,HalfInt(1,2),orbital_file[i-(file_size/2)].jb,-HalfInt(1,2),HalfInt(lambda_value),HalfInt(0));

	      single_particle_me_n = single_particle_me_coefficient*radial_integral_n;
	      
	      reduced_many_body_matrix_element_neutrons = reduced_many_body_matrix_element_neutrons + pow(boscil_n,lambda_value) * single_particle_me_n * (one_body_densities[i]);
            }
	}

      fprintf(stdout,"%c%d  %.1f %2d %2d  %.1f %2d %2d %10f %10f \n",input_file_parameters.observable_type,lambda_value,
	      header_info.ji.DValue(), header_info.parity_i, header_info.state_label_i, header_info.jf.DValue(), header_info.parity_f, 
	      header_info.state_label_f, reduced_many_body_matrix_element_protons/Hat(lambda_value), reduced_many_body_matrix_element_neutrons/Hat(lambda_value));
    
      fprintf(output_file,"%c%d  %.1f %2d %2d  %.1f %2d %2d %10f %10f \n",input_file_parameters.observable_type,lambda_value,
	      header_info.ji.DValue(), header_info.parity_i, header_info.state_label_i, header_info.jf.DValue(), header_info.parity_f, 
	      header_info.state_label_f, reduced_many_body_matrix_element_protons/Hat(lambda_value), reduced_many_body_matrix_element_neutrons/Hat(lambda_value));
      
      fclose(output_file);
    }	  
  else if( (input_file_parameters.observable_type == 'M') && 
	   (AllowedTriangle(header_info.ji,header_info.jf,HalfInt(lambda_value)) ==1) )
    { 
      vector <OrbitalInformation> orbital_file;
      ReadOrbitalInformationFile (orbital_file, input_file_parameters.path_to_orbital_information_file);

      shell::RadialMatrices radial_matrix_protons = ReturnRadialIntegrals (input_file_parameters.path_to_radial_integrals_file_protons,
								     input_file_parameters.observable_type,
								     input_file_parameters.lambda);
      shell::RadialMatrices radial_matrix_neutrons = ReturnRadialIntegrals (input_file_parameters.path_to_radial_integrals_file_neutrons,
								     input_file_parameters.observable_type,
								     input_file_parameters.lambda);

      vector <double> one_body_densities;
      ReadOneBodyDensitiesFile (one_body_densities, input_file_parameters.path_to_obdme_file);
      double pi = 3.141593;  
   
      double boscil_p = OscillatorLength(input_file_parameters.beta_p,input_file_parameters.hw);  
      double boscil_n = OscillatorLength(input_file_parameters.beta_n,input_file_parameters.hw); 

      double dl_p = 0, dl_n = 0, ds_p = 0, ds_n = 0, total_reduced_matrix_element; 
     
      double single_particle_matrix_element_coefficient, single_particle_matrix_element_p, single_particle_matrix_element_n;

      double kappa, radial_integral_p, radial_integral_n;
 
      int file_size = header_info.number_of_elements; 
   
      // proton sector      

      for (int i=0; i < file_size/2; i++)
	{
	  if( (orbital_file[i].la+orbital_file[i].lb+lambda_value)%2!=0 && 
              (AllowedTriangle(orbital_file[i].ja,orbital_file[i].jb,HalfInt(lambda_value)) == 1) &&
              (orbital_file[i].k == lambda_value)
            )  
	    {          
	      //cout << "Indexes (lb la) are " << orbital_file[i].lb << " " << orbital_file[i].la << endl;
	      //Angular momentum indexes print to check what kind of input file we need for the radial integrals

 	      radial_integral_p = radial_matrix_protons.GetMatrixElement(orbital_file[i].la,orbital_file[i].lb,orbital_file[i].na,orbital_file[i].nb); 	      	     

	      single_particle_matrix_element_coefficient = (ParitySign(orbital_file[i].ja+HalfInt(lambda_value)-HalfInt(1,2))/sqrt(4*pi)) 
		* Hat(orbital_file[i].jb) * Hat(orbital_file[i].ja) 
		* ClebschGordan(orbital_file[i].ja,HalfInt(1,2),orbital_file[i].jb,-HalfInt(1,2),HalfInt(lambda_value),HalfInt(0)); 

	      kappa = ParitySign(orbital_file[i].ja+HalfInt(orbital_file[i].la)+HalfInt(1,2))*DValue(orbital_file[i].ja+HalfInt(1,2))
		+ParitySign(orbital_file[i].jb+HalfInt(orbital_file[i].lb)+HalfInt(1,2))*DValue(orbital_file[i].jb+HalfInt(1,2));
	      
	      single_particle_matrix_element_p = single_particle_matrix_element_coefficient * radial_integral_p * (lambda_value-kappa); 	      
	      
	      dl_p = dl_p +  pow(boscil_p,lambda_value-1) * (1+kappa/(lambda_value+1)) * single_particle_matrix_element_p * (one_body_densities[i]);

	      ds_p = ds_p + (-0.5) *  pow(boscil_p,lambda_value-1)* single_particle_matrix_element_p * (one_body_densities[i]);	           
		
	    }
	}
      for (int i=file_size/2; i < file_size; i++)
	{
	  if( (orbital_file[i-file_size/2].la+orbital_file[i-file_size/2].lb+lambda_value)%2!=0 && 
              (AllowedTriangle(orbital_file[i-file_size/2].ja,orbital_file[i-file_size/2].jb,HalfInt(lambda_value)) == 1) &&
              (orbital_file[i-file_size/2].k == lambda_value)
            )  
	    {          
	      //cout << "Indexes (lb la) are " << orbital_file[i].lb << " " << orbital_file[i].la << endl;
	      //Angular momentum indexes print to check what kind of input file we need for the radial integrals

 	      radial_integral_n = radial_matrix_neutrons.GetMatrixElement(orbital_file[i-file_size/2].la,orbital_file[i-file_size/2].lb,orbital_file[i-file_size/2].na,orbital_file[i-file_size/2].nb); 	     

	      single_particle_matrix_element_coefficient = (ParitySign(orbital_file[i-file_size/2].ja+HalfInt(lambda_value)-HalfInt(1,2))/sqrt(4*pi)) 
		* Hat(orbital_file[i-file_size/2].jb) * Hat(orbital_file[i-file_size/2].ja) 
		* ClebschGordan(orbital_file[i-file_size/2].ja,HalfInt(1,2),orbital_file[i-file_size/2].jb,-HalfInt(1,2),HalfInt(lambda_value),HalfInt(0)); 

	      kappa = ParitySign(orbital_file[i-file_size/2].ja+HalfInt(orbital_file[i-file_size/2].la)+HalfInt(1,2))*DValue(orbital_file[i-file_size/2].ja+HalfInt(1,2))
		+ParitySign(orbital_file[i-file_size/2].jb+HalfInt(orbital_file[i-file_size/2].lb)+HalfInt(1,2))*DValue(orbital_file[i-file_size/2].jb+HalfInt(1,2));

	      single_particle_matrix_element_n = single_particle_matrix_element_coefficient * radial_integral_n * (lambda_value-kappa); 	      
	      
	      dl_n = dl_n +  pow(boscil_n,lambda_value-1) * (1+kappa/(lambda_value+1)) * single_particle_matrix_element_n * (one_body_densities[i]);

	      ds_n = ds_n + (-0.5) *  pow(boscil_n,lambda_value-1)* single_particle_matrix_element_n * (one_body_densities[i]);	       
		
	    }
	}  

      total_reduced_matrix_element = dl_p * 1 + dl_n * 0 + ds_p * 5.586 + ds_n *(-3.826);         

      fprintf(stdout,"%c%d  %.1f %2d %2d  %.1f %2d %2d %10f %10f %10f %10f %10f \n",input_file_parameters.observable_type,lambda_value,
	      header_info.ji.DValue(), header_info.parity_i, header_info.state_label_i, header_info.jf.DValue(), header_info.parity_f,
	      header_info.state_label_f, dl_p/Hat(lambda_value), dl_n/Hat(lambda_value), ds_p/Hat(lambda_value), ds_n/Hat(lambda_value), total_reduced_matrix_element/Hat(lambda_value));	
    
       fprintf(output_file,"%c%d  %.1f %2d %2d  %.1f %2d %2d %10f %10f %10f %10f %10f \n",input_file_parameters.observable_type,lambda_value,
	       header_info.ji.DValue(), header_info.parity_i, header_info.state_label_i, header_info.jf.DValue(), header_info.parity_f,
	       header_info.state_label_f, dl_p/Hat(lambda_value), dl_n/Hat(lambda_value), ds_p/Hat(lambda_value), ds_n/Hat(lambda_value), total_reduced_matrix_element/Hat(lambda_value));	       
      
      fclose(output_file);
    } 
}

/********************************************************************************************************************************

Main function: The only thing needed in this new version is to declare a struct PostProcessingParameters struct, initialize 
it using the function ReadPostProcessingParameters and then pass the struct into the CalculateElectromagneticMatrixElements function

**********************************************************************************************************************************/  

int main()
{  
  PostProcessingParameters input;
  ReadPostProcessingParameters(input);  

  CalculateElectromagneticMatrixElements(input); 


  // Debugging the upgraded version on 04/05/15..

  // cout << number_operator(input) << endl;
    
  /*
  cout << input.observable_type << " " << input.lambda << " " << input.hw << " " << input.beta_p << " " << input.beta_n << endl;
  cout << input.path_to_orbital_information_file << endl;
  cout << input.path_to_obdme_file << endl;
  cout << input.path_to_radial_integrals_file_protons << endl;
  cout << input.path_to_radial_integrals_file_neutrons << endl;
  */
  
  /*
  vector <OrbitalInformation> orbitals;
  ReadOrbitalInformationFile(orbitals, input.path_to_orbital_information_file);

  vector <double> densities;
  ReadOneBodyDensitiesFile(densities,input.path_to_obdme_file);
  */

  /*for(int i=0; i<orbitals.size(); i++)
     {
       cout << orbitals[i].na << " " << orbitals[i].la << " " << orbitals[i].ja << " " << orbitals[i].nb << " " << orbitals[i].lb << " " <<
               orbitals[i].jb << " " << orbitals[i].k  << " " << orbitals[i].dp << endl;   
     }*/
  
  /*for(int i = 0; i< densities.size(); i++)
     {
       cout << densities[i] << endl;
     }*/
  
  /*
  ObdmeHeaderInformation header;
  InitializeObdmeHeaderInformation (input.path_to_obdme_file, header);

  cout << header.number_of_elements << endl;
  cout << header.parity_i << " " << header.parity_f << endl;
  cout << header.state_label_i << " " << header.state_label_f << endl;
  cout << header.ji << " " << header.jf << endl;
  cout << header.ti << " " << header.tf << endl;
  cout << header.mm_basis << endl;
  */
  
  return 0;
}


