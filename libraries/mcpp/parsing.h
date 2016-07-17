/****************************************************************
  parsing.cpp

  Generic input parsing and error message utilities.
`
  Mark A. Caprio
  University of Notre Dame

  3/10/16 (mac): Extracted from mfdn_h2.

****************************************************************/

#ifndef PARSING_H_
#define PARSING_H_

#include <iostream>
#include <sstream>
#include <string>

void OpenCheck(bool success, const std::string& filename);
// Terminates with error message on stream open failure.
//
// Arguments:
//   success (bool) : meant to be used as automatic conversion from string state
//   filename (string) : file name to use in error messasge
//
// Example:
//   OpenCheck(in_stream,in_stream_name);

void ParsingCheck(std::istringstream& line_stream, int line_count, const std::string& line);
// Provide error message upon parsing failure for line of input.
//
// Limitations: Would ideally also support arguments to give filename
// and optional supplementary information on expected content.
//
// Arguments:
//   line_stream (istringstream) : string stream from which parsing was attempted
//   line_count (int) : line count for error message
//   line (string) : line text for error message
//
// Example:
//
//    // scan input file
//    std::string line;
//    int line_count = 0;
//    while ( std::getline(in_stream, line) )
//    {
//      // count line
//      ++line_count;
//
//      // set up for parsing
//      std::istringstream line_stream(line);
//
//      // parse line
//      int a, b, ...;
//      line_stream >> a >> b >> ...;
//      ParsingCheck(line_stream, line_count, line);
//  
//      // do stuff with input
//      ...
//    }

#endif
