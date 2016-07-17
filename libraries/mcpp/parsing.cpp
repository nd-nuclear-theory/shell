/****************************************************************
  parsing.cpp

  Mark A. Caprio
  University of Notre Dame

****************************************************************/

#include "parsing.h"

#include <cstdlib>

void OpenCheck(bool success, const std::string& filename)
{
  if (!success)
    {
      std::cerr << std::endl;
      std::cerr << "Failed to open file: " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }
}

void ParsingError(const std::string& message, int line_count, const std::string& line)
{
  std::cerr << std::endl;
  std::cerr << message << std::endl;
  std::cerr << "Input line " << line_count << ": " << line  << std::endl;
  std::exit(EXIT_FAILURE);
}
	
void ParsingCheck(std::istringstream& line_stream, int line_count, const std::string& line)
{
  if (!line_stream)
    ParsingError("Failed parsing line (missing or incorrect arguments)", line_count, line);
}
