/******************************************************************************
  
  Created by M. A. Caprio, 2/23/11.
  Patch include path 7/4/16.

******************************************************************************/

#include "memoizer.h"
#include "profiling.h"

#include <iostream>
#include "am/halfint.h"

using namespace std;

HalfInt f(int i)
{

	static Memoizer<int,HalfInt> m;

	HalfInt v;
	if (m.Seek(i))
	{
		return m.GetValue();
	}
	else
	{
		cout << "(" << i << ")";
		return m.SetValue(i,HalfInt(i,2));
	}

}

int g(int i)
{

	static Memoizer<int,int> m;

	return MEMOIZE(m, i, (cout << "(" << i << ")", i*i));
	
}

int factorial_disabled(int i)
{

	static Memoizer<int,int> m;
	m.EnableCaching(false);

	if (i==0) 
		return 1;
	else
		return MEMOIZE(m, i, (cout << "(" << i << ")",i*factorial_disabled(i-1)) );
	
}

int factorial_cached(int i)
{

	static Memoizer<int,int> m;

	if (i==0) 
		return 1;
	else
		return MEMOIZE(m, i, (cout << "(" << i << ")", i*factorial_cached(i-1)) );
	
}

int factorial_dump(int i)
{

	static Memoizer<int,int> m;
	
	int v;

	if (i==0) 
		v = 1;
	else
		v = MEMOIZE(m, i, (i*factorial_dump(i-1)) );

	if (i==10)
		cout << m;

	return v;
}


int factorial_plain_timing(int i)
{

	if (i==0) 
		return 1;
	else
		return (i*factorial_plain_timing(i-1));
	
}

int factorial_disabled_timing(int i)
{

	static Memoizer<int,int> m(false);

	if (i==0) 
		return 1;
	else
		return MEMOIZE(m, i, (i*factorial_disabled_timing(i-1)) );
	
}

int factorial_cached_timing(int i)
{

	static Memoizer<int,int> m;

	if (i==0) 
		return 1;
	else
		return MEMOIZE(m, i, (i*factorial_cached_timing(i-1)) );
	
}

int main(int argc, char **argv)
{

	static Memoizer<vector<int>,int> m;

	// function defined using if structure

	cout << f(1) << endl;
	cout << f(2) << endl;
	cout << f(3) << endl;
	cout << f(1) << endl;
	cout << f(2) << endl;

	cout << "****" << endl;

	cout << g(1) << g(2) << g(3) << g(1) << endl;

	cout << "****" << endl;

	cout << "No caching..." << endl;
	for (int i = 1; i<=5; ++i)
		cout << factorial_disabled(i) << endl;
	cout << factorial_disabled(10) << endl;
	cout << factorial_disabled(5) << endl;
	cout << "****" << endl;

	cout << "With caching..." << endl;
	for (int i = 1; i<=5; ++i)
		cout << factorial_cached(i) << endl;
	cout << factorial_cached(10) << endl;
	cout << factorial_cached(5) << endl;
	cout << "****" << endl;

	cout << "Cache dump..." << endl;
	factorial_dump(10);
	cout << "****" << endl;

	Timer t;
	int x;
	int n_max = 10000;

	t.Start();
	x=0;
	for (int i = 1; i<=n_max; ++i)
		x+= factorial_plain_timing(i);
	t.Stop();
	cout << "Time without caching: " << t.ElapsedTime() << ",  result " << x << endl;

	t.Start();
	x=0;
	for (int i = 1; i<=n_max; ++i)
		x+= factorial_disabled_timing(i);
	t.Stop();
	cout << "Time with caching disabled: " << t.ElapsedTime() << ",  result " << x << endl;

	t.Start();
	x=0;
	for (int i = 1; i<=n_max; ++i)
		x+= factorial_cached_timing(i);
	t.Stop();
	cout << "Time with caching: " << t.ElapsedTime() << ",  result " << x  << endl;
	
	// On mac 2/23/11:
	// for n_max = 10000;
	// Time without caching: 0.469,  result -125961703
	// Time with caching disabled: 35.109,  result -125961703
	// Time with caching: 0.047,  result -125961703

	// termination
	return 0;
}

		
