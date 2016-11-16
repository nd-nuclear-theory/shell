#include "radial_integrals.h"

using namespace std;

void RadialIntegrals::Initialize(istream& opened_file)
{
  string dataline;
  double integral;

  getline(opened_file,dataline);
  istringstream buffer(dataline);
  buffer >> lmax_ >> dl_ >> nmax_ >> mmax_; 
  
  std::vector<std::vector<std::vector <double> > > three_dimentional_array;
  std::vector<std::vector<double> > two_dimentional_array;
  std::vector<double> one_dimentional_array;

  for (int la=0; la<=lmax_;la++)
    {
      for (int lb=(la+dl_%2);lb<=min(la+dl_,lmax_);lb+=2)
	{	  
	  for (int i=0; i < nmax_;i++)
	    {
	      for (int j=0; j< nmax_;j++)
		{
		  opened_file >> integral;
		  one_dimentional_array.push_back(integral); 	  
		}	    
	      two_dimentional_array.push_back(one_dimentional_array); 
	      one_dimentional_array.clear();
	    }
	  three_dimentional_array.push_back(two_dimentional_array);
	  two_dimentional_array.clear();
	}
      matrix_elements_.push_back(three_dimentional_array);
      three_dimentional_array.clear();
    }
}

double RadialIntegrals::GetElement(int& na,int& la,int& nb,int& lb)
{
  double the_value;
  int index=(abs(la-lb)-dl_%2)/2;

  if( (la >= 0) && (lb >= 0) && (la <= lmax_) && (lb <= lmax_) && (na >= 0) && (nb >= 0) && (na <= nmax_) && (nb <= nmax_)
     && ((la+lb+dl_)%2 == 0) && (abs(la-lb) <= dl_) && nmax_ == mmax_ )
    {
      if(la>lb)
	the_value = matrix_elements_[lb][index][nb][na];
      else 
	the_value = matrix_elements_[la][index][na][nb];
    }
  else
    {
      cerr << "The radial integral is invalid. Please check indices" << endl;
      exit(EXIT_FAILURE);     
    }
  
  return the_value;
} 

void RadialIntegrals::Free()
{
  matrix_elements_.clear();
}
