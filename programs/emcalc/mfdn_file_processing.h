/**************************************
  mfdn_file_processing.h
 
  Utility functions which extract information from the mfdn
  one body info and one body density output files.
                                  
  Created by Valentino Constantinou on 04/03/15.
  Code extracted from emcalc originally used for the calculation 
  of electromagnetic observables. The code is also useful for the 
  calculation of the natural orbitals:

  -- Renamed struct FilePathsAndCalculationParameters to PostProcessingParameters,
     ReadInputFile to ReadPostProcessingParameters, ReadObdmeFile to ReadOneBodyDensitiesFile  


  

**************************************/

#ifndef mfdn_file_processing_h
#define mfdn_file_processing_h


#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include "halfint.h"

using namespace std;

/****** struct PostProcessingParameters ************                                   
                                                                                                  
 observable_type: char variable that can take the value E or M 
 for electric or magnetic moments         
                                                                                                         
 lambda: Type of transition, dipole for lambda = 1 for instance                              
                                                                                                                   
 hw: double variable which gives the oscillator length                                       
                                                                                                           
 beta_p, beta_n: Allows the selection of different length scales for protons or neutrons. The      
 length has to do with the length parameter of the basis                                     
                                                                                                         
 path_to_orbital_information_file: string variable which gives the full path (including the file name) 
 to the file mfdn.rppobdme.info. Example /home/valentino/mfdn.rppobdme.info                                    
                                                                                                         
 path_to_obdme_file: string with the full path to the obdme file. 
 Example /home/valentino/mfdn.statrobdme.seq001.2J02.p0.n01.2T00                                   
                                                                                                         
 path_to_radial_integrals_file_protons/neutrons: string that gives the path to the files that contain 
 the radial integrals R_ab^{/lambda} as given by Suhonen 6.32. One can choose different radial 
 integrals for protons or neutrons. The integrals can be either calculated
 in the CS basis or the HO basis. The file shall be given as radial-me-HO
 (i.e not radial-me-HO-r2.dat). The file is then constructed internally as 
 radial-me-HO-r2.dat. An example of a path is /home/valentino/radial-me-HO     

***************************************/                                                                                                    

struct PostProcessingParameters
{
  char observable_type;                                 
  int lambda;                               
  double hw;
  double beta_p, beta_n;
  string path_to_orbital_information_file;
  string path_to_obdme_file;
  string path_to_radial_integrals_file_protons;
  string path_to_radial_integrals_file_neutrons;
};

/***************struct OrbitalInformation************************

04/03/15: The struct was upgraded to read data from the new file mfdn.rppobdme.info. 
          In the new format the isospin tz is no longer included and neither does the magnetic 
          projection m. Instead two new numbers are included. k which denotes the coupling of 
          [ca^{dagger} cb^{tilde}]_{k} and dp which denotes the parity.
          zero denotes positive parity.. For example we need k=2 to calculate E2 values.. 

          After the information is stored it is later used to sum over orbitals and calculate the 
          E1, E2 or M1 values we want.. Also we can calculate the number operator sum..

******************************************************************/      

struct OrbitalInformation
{
  int na, la, nb, lb, k, dp;
  HalfInt ja, jb;
};

/************** void ReadPostProcessingParameters ( PostProcessingParameters& input_information)***********

Takes a structure PostProcessingParameters as an input and it initializes all of its parameters
by reading the following input file line by line. Example Input:   
                                                                  
E     1                                                                                
25                                                                            
1     1                                                                         
/home/valentino/mfdn.rppobdme.info                                               
/home/valentino/mfdn.statrobdme.seq001.2J02.p0.n01.2T00                  
/home/valentino/radial-me-HO          /home/valentino/radial-me-HO                                    
         
*****************************************************************************************************/                                              

/************struct ObdmeHeaderInformation********************

It holds all the information contained in the first two lines 
of the obdme file which is the output of mfdn
 TODO: Instead of looking into the density file perhaps its better to read the filename
       and extract the information from there in the future..
**************************************************************/

struct ObdmeHeaderInformation
{
  int number_of_elements, parity_i, 
    parity_f, state_label_i, state_label_f;

// (04/03/15) Please note that ji is actually the angular momentum of the many body state Ji! Same for Jf, 
// Ti, Tf and M..   

  HalfInt ji, jf, ti, tf, mm_basis;  
};


void ReadPostProcessingParameters(PostProcessingParameters&);

/************void ReadOrbitalInformationFile( vector <OrbitalInformation>& 
orbital_data_vector, string& input_file_path )***********   

Takes a vector of OrbitalInformation structs as an input since we need one OrbitalInformation struct per line of input 
of the mfdn.rppobdme.info. It also takes a string which is the path to the mfdn.rppobdme.info file. Eventually the data are stored 
in the vector orbital_data_vector which is a vector of OrbitalInformation
                                                                                                         
*********************************/

void ReadOrbitalInformationFile(vector <OrbitalInformation>&, string&);

/*********void ReadObdmeFile(vector <double>& obdme_vector, string& obdme_file_path)**********

A function that initializes the vector of doubles that will hold the values of the rppobdme. The string
is the path to the obdme file 

***********************************************************************************************/

void ReadOneBodyDensitiesFile(vector <double>&,string&);

/**************void InitializeObdmeHeaderInformation (string& file_path, 
				     ObdmeHeaderInformation& header)

The function opens the obdme file and initializes the struct 
ObdmeHeaderInformation

************************************************************************/

void InitializeObdmeHeaderInformation (string&, ObdmeHeaderInformation&);

#endif
