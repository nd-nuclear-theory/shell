
#include "no/mfdn_file_processing.h"

void ReadPostProcessingParameters ( PostProcessingParameters& input_information )
{
  string input_line;

  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.observable_type 
		>> input_information.lambda;
  }

  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.hw;
  }

  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.beta_p 
		>> input_information.beta_n;
  }
  
  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.path_to_orbital_information_file;
  }

  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.path_to_obdme_file;
  }
  
  {
    getline(cin,input_line);
    istringstream line_stream(input_line);
    line_stream >> input_information.path_to_radial_integrals_file_protons 
		>> input_information.path_to_radial_integrals_file_neutrons;
  }  
}


void ReadOrbitalInformationFile( vector <OrbitalInformation>& 
				      orbital_data_vector, string& input_file_path )
{  
  string dataline; //dataline temporarily holds each line of the input 
  int counter = 0;

  ifstream orbital_file(input_file_path.c_str(),ios::in);

  if(!orbital_file)  
    {
      cerr << "The orbital information file mfdn.obdme.info cannon be opened. "
	      "Please check the file em_calc.in and make sure that the path"
	      " to the orbital information file (mfdn.obdme.info) is "
	      "correct. Also make sure that the file exists." 
	   << endl;
      exit(EXIT_FAILURE);
    }  
  else
    {
      for(int i = 0; i<6;i++) //Getting rid of the header in the mfdn.obdme.info.
	{
	  getline(orbital_file,dataline);
	}
  
      while (getline (orbital_file,dataline))
	{
	  orbital_data_vector.push_back(OrbitalInformation());
	  istringstream buffer(dataline);
	  int index; // Note that index is not stored. And as far as the info file is concerned its just an index.. 
	  int twice_ja, twice_jb; 

	  buffer >> index 
		 >> orbital_data_vector[counter].na 
		 >> orbital_data_vector[counter].la 
		 >> twice_ja >> orbital_data_vector[counter].nb 
		 >> orbital_data_vector[counter].lb >> twice_jb 
		 >> orbital_data_vector[counter].k >> orbital_data_vector[counter].dp;
      
	  orbital_data_vector[counter].ja = HalfInt(twice_ja,2);
	  orbital_data_vector[counter].jb = HalfInt(twice_jb,2);
      
	  counter++;            
	}
    }
  
  orbital_file.close();      
}


void ReadOneBodyDensitiesFile( vector <double>& obdme_vector, string& obdme_file_path )
{
  ifstream read_densities(obdme_file_path.c_str(),ios::in);  
  string header_lines; 
  
  //The information in the first 2 lines are read in the
  //ObdmeFileInformation struct
  
  if(!read_densities)  
    {
      cerr << "The obdme file cannon be opened. Please check the file emcalc.in"
	" and make sure that the path is correct or that the file is not empty" 
	   << endl;
      exit(EXIT_FAILURE);
    }
  else
    {
      double density_value; 
      double sqrt_four_pi = 3.544907702; // Note: When the issue with the sqrt of 4Pi is fixed this line should be removed..

      for(int i = 0; i<5;i++) //Getting rid of the header..
	{
	  getline(read_densities,header_lines);
	}
      
      while(read_densities >> density_value)
	{
	  obdme_vector.push_back(density_value*sqrt_four_pi);
	}  
    }

  read_densities.close();
}


void InitializeObdmeHeaderInformation (string& file_path, 
				     ObdmeHeaderInformation& header)
{
  ifstream open_file(file_path.c_str(),ios::in);
  stringstream buffer;
  string dataline;
  int twice_ji, twice_jf, twice_ti, twice_tf, twice_mm, proton_elements, neutron_elements;

  if(!open_file)
    {
      cerr << "The obdme cannot be opened. Please check"
	" the input file" << endl;
      exit(EXIT_FAILURE);
    }

  {
    // (04/03/15): Getting rid of first two lines..
    for(int i = 0; i<2;i++) //Getting rid of the first two lines..
	{
	  getline(open_file,dataline);
	}

    getline(open_file, dataline);
    buffer << dataline;
    
    buffer >> twice_mm >> twice_ji >> header.parity_i >>
           header.state_label_i >> twice_ti;

    header.ji = HalfInt(twice_ji,2);
    header.ti = HalfInt(twice_ti,2);
    header.mm_basis = HalfInt(twice_mm,2);
    buffer.str("");
    buffer.clear();

    /*--------Debug on 01/26/14------------

    cout << "File content" << endl;	
    cout << header.magic_number << " " << header.number_of_elements << " " <<  header.mm_basis << " " << header.ji << " " << header.parity_i << " " << header.state_label_i << endl;*/
  }

  {
    getline(open_file, dataline);
    buffer << dataline;
    
    buffer >>  twice_mm >> twice_jf >> header.parity_f >> 
           header.state_label_f >> twice_tf;

    header.jf = HalfInt(twice_jf,2);
    header.tf = HalfInt(twice_tf,2);
    buffer.str("");
    buffer.clear();

    /*-------------Debug on 01/26/14------------------    
    cout << "File content" << endl;	
    cout << header.ti << " " << header.jf << " " <<  header.parity_f << " " << header.state_label_f << " " << header.tf << endl;*/
  }

  {
    getline(open_file, dataline);
    buffer << dataline;
    
    buffer >> proton_elements >> neutron_elements;

    header.number_of_elements = proton_elements + neutron_elements;
    buffer.str("");
    buffer.clear();
  }
}



