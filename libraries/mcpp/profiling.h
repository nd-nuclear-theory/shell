/****************************************************************

  profiling.h

  Simple timing with clock().

  Created by Mark A. Caprio, 2/23/11.
                                  
****************************************************************/

#ifndef PROFILING_H_
#define PROFILING_H_

#include <ctime>

////////////////////////////////////////////////////////////////
// simple clock timing
////////////////////////////////////////////////////////////////

class Timer{

public:
	////////////////////////////////
	// constructors
	////////////////////////////////

	// default constructor: initializes timer values to zero
	Timer() : start_time_(0), end_time_(0) {}; 

	////////////////////////////////
	// stopwatch
	////////////////////////////////

	void Start() {start_time_ = clock();};
	void Stop() {end_time_ = clock();};

	double ElapsedTime() const {
		return (static_cast<double>(end_time_ - start_time_)) / CLOCKS_PER_SEC;
	};

	////////////////////////////////
	// 
	////////////////////////////////

private:

	clock_t start_time_;
	clock_t end_time_;

};

#endif
