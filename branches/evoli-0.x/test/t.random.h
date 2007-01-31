/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


#ifndef _T_RANDOM_H__
#define _T_RANDOM_H__
#include "cutee.h"
#include "random.hh"

#include <cmath>

struct TEST_CLASS( random )
{
	void TEST_FUNCTION( rint )
	{
		// the first 21 random integers for seed 4357
		uint randints[21] = { 3510405877U, 4290933890U, 2191955339U,
			564929546U,  152112058U, 4262624192U, 2687398418U,
			268830360U, 1763988213U,  578848526U, 4212814465U,
			3596577449U, 4146913070U, 950422373U, 1908844540U,
			1452005258U, 3029421110U, 142578355U, 1583761762U,
			1816660702U, 2530498888U };
		Random::seed( 4357U );

		bool rint_test_passed = true;
		for( int i=0; i<21; i++ )
			rint_test_passed = rint_test_passed && ( Random::rint() == randints[i] );
		TEST_ASSERT( rint_test_passed );

		bool rint_lim_test_passed = true;
		for( int i=0; i<1000; i++ )
		{
			int k = Random::rint(100);
			rint_lim_test_passed = rint_lim_test_passed && ( k >=0 && k < 100 );
		}
		TEST_ASSERT( rint_lim_test_passed );
		return;
	}

	void TEST_FUNCTION( runif )
	{
		bool runif_test_passed = true;
		for( int i=0; i<100000; i++ )
		{
			double x = Random::runif();
			runif_test_passed = runif_test_passed && ( x >= 0. && x < 1. );
		}
		TEST_ASSERT( runif_test_passed );
		return;
	}

	void TEST_FUNCTION( rnorm )
	{
		// The first 10 standard normal random deviates for random seed 23
		double randnorms[10] = { 0.345669, 1.56275, -0.118777, -0.858446, 1.66391,
				-0.326751, -0.428157, 0.00572826, -0.598339, -0.907825 };

		Random::seed( 23U );
		bool rstdnorm_test_passed = true;
		for( int i=0; i<10; i++ )
			rstdnorm_test_passed = rstdnorm_test_passed && ( fabs( Random::rnorm() - randnorms[i] ) < 1e-5 ) ;
		TEST_ASSERT( rstdnorm_test_passed );

		Random::seed( 23U );
		bool rnorm_test_passed = true;
		for( int i=0; i<10; i++ )
			rnorm_test_passed = rnorm_test_passed && ( fabs( Random::rnorm( 3.1, 2.4 ) - 2.4*randnorms[i] - 3.1 ) < 1e-4 ) ;
		TEST_ASSERT( rnorm_test_passed );
		return;
	}

	void TEST_FUNCTION( rpois )
	{
		// The first 10 poisson random deviates for random seed 45 and mean 2.7.
		uint randpois[10] = { 5, 4, 0, 1, 6, 2, 2, 4, 4, 1 };
		Random::seed( 45U );
		bool rpois_test_passed = true;
		for( int i=0; i<10; i++ )
			rpois_test_passed = rpois_test_passed && ( Random::rpois( 2.7 ) == randpois[i] );
		TEST_ASSERT( rpois_test_passed );
	}

	void TEST_FUNCTION( randint )
	{
		double distr[10] = { .1, .2, .2, .25, .26, .3, .6, .7, .8, 1. };
		int counts[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int n=100000;
		for( int i=0; i<n; i++ )
			counts[Random::randintFromDistr( distr, 10 )] += 1;
		bool randint_test_passed
				= fabs( counts[0]/static_cast<double>(n) - distr[0] ) < .003;
		for( int i=1; i<10; i++ )
			randint_test_passed = randint_test_passed && ( fabs( counts[i]/static_cast<double>(n) - distr[i] + distr[i-1] ) < .003 );
		TEST_ASSERT( randint_test_passed );
	}

	void TEST_FUNCTION( dpois  )
	{
		// first 10 Poission probabilities for mean 2:
		double dpois[10] = { 0.1353352832, 0.2706705665, 0.2706705665,
			0.1804470443, 0.0902235222, 0.0360894089, 0.0120298030,
			0.0034370866, 0.0008592716, 0.0001909493 };
		bool dpois_test_passed = true;
		for( int i=0; i<10; i++ )
			dpois_test_passed = dpois_test_passed && ( fabs( Random::dpois( i, 2 ) - dpois[i] ) < 1e-8 );
		TEST_ASSERT( dpois_test_passed );
	}

	void TEST_FUNCTION( dbinom )
	{
		// Binomial probabilities for n = 10 and p = .13:
		double dbinom[11] = { 2.484234e-01, 3.712074e-01, 2.496050e-01,
			9.945945e-02, 2.600808e-02, 4.663517e-03, 5.807061e-04,
			4.958410e-05, 2.778420e-06, 9.225914e-08, 1.378585e-09 };
		bool dbinom_test_passed = true;
		for( int i=0; i<11; i++ )
		{
			dbinom_test_passed = dbinom_test_passed && ( fabs( Random::dbinom( i, 10, .13 ) - dbinom[i] )/dbinom[i] < 1e-6 );
		}
		TEST_ASSERT( dbinom_test_passed );
	}
};


#endif //_T_RANDOM_H__
