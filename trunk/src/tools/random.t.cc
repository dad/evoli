#include <iostream>
#include <cmath>
#include "random.hh"

using namespace std;

void test( bool testPassed )
{
	if ( testPassed )
		cout << "test passed" << endl;
	else
		cout << "test failed" << endl;
}

int main()
{
	// the first 21 random integers for seed 4357
	uint randints[21] = { 3510405877U, 4290933890U, 2191955339U,
			564929546U,  152112058U, 4262624192U, 2687398418U,
			268830360U, 1763988213U,  578848526U, 4212814465U,
			3596577449U, 4146913070U, 950422373U, 1908844540U,
			1452005258U, 3029421110U, 142578355U, 1583761762U,
			1816660702U, 2530498888U };
	Random::seed( 4357U );

	bool testPassed = true;
	for( int i=0; i<21; i++ )
		testPassed = testPassed && ( Random::rint() == randints[i] );
	test( testPassed );

	testPassed = true;
	for( int i=0; i<100000; i++ )
	{
		double x = Random::runif();
		testPassed = testPassed && ( x >= 0. && x < 1. );
	}
	test( testPassed );

	// The first 10 random normal deviates for random seed 23
	double randnorms[10] = { 0.345669, 1.56275, -0.118777, -0.858446, 1.66391,
				-0.326751, -0.428157, 0.00572826, -0.598339, -0.907825 };

	Random::seed( 23U );
	testPassed = true;
	for( int i=0; i<10; i++ )
		testPassed = testPassed && ( fabs( Random::rnorm() - randnorms[i] ) < 1e-5 ) ;
	test( testPassed );

	Random::seed( 23U );
	testPassed = true;
	for( int i=0; i<10; i++ )
		testPassed = testPassed && ( fabs( Random::rnorm( 3.1, 2.4 ) - 2.4*randnorms[i] - 3.1 ) < 1e-4 ) ;
	test( testPassed );

}

