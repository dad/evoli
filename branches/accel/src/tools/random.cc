/*
This file is part of the evoli project.
Copyright (C) 2006 Claus Wilke <cwilke@mail.utexas.edu>

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

This file contains wrapper code around the implemtation by 
Shawn Cokus of the Mersenne Twister RNG by Makoto Matsumoto 
and Takuji Nishimura. See http://www.math.keio.ac.jp/~matumoto/emt.html
for more on the algorithm.

Part of this code was taken from the Octave modules for the Mersenne
Twister (MT199337) Random Number Generator, Copyright (C) 1998,
1999 Dirk Eddelbuettel <edd@debian.org>
*/

#include "random.hh"
#include <climits>
#include <cmath>
#include <cassert>


// define DRAND48 to use drand48() instead of the Mersenne twister
// useful for debugging purposes
// #define DRAND48
#ifdef DRAND48
	#include <cstdlib>
#endif



typedef unsigned long uint32;

extern "C" {                    // include as C
	#include "cokus.h"
}

/**
Convenience function that turns a random integer into a random double on [0,1).
*/
inline double randomMTunit(void) {
#ifdef DRAND48
	return drand48();
#else
	// multiply by 1/0xFFFFFFFF
	// code taken from mt199379-1.c
	return static_cast<double>( randomMT() ) * 2.328306437080794e-10;
#endif
}


/**
Routine that generates random standard normal variates, using the Box-Muller method.
(G. E. P. Box and M. E. Muller, Ann. Math. Statist. 29:610-611, 1958).
*/
static double snorm( void )
{
	static double gauss2; // second normal random variate
	static double have_gauss2 = false; // do we have gauss2 already?
	double pi = 3.14159265358979323846264338327950288419716939937510;

	if ( have_gauss2 )
	{
		have_gauss2 = false;
		return gauss2;
	}
	double uv = randomMTunit(); // first uniform random variate
	double x = sqrt(-2*log(1-uv)); // need 1-uv1 to avoid zero in log
	uv = randomMTunit(); // second uniform random variate
	double y = 2*pi*uv;
	gauss2 = x*sin(y); // store second variate for future use
	have_gauss2 = true;
	return x*cos(y);
}


// The following is taken from http://www.herine.net/stat/papers/dbinom.pdf
// the Loader algorithm for fast, accurate computation of binomial probabilities.

// NTYPE is the type used for the n and x arguments.
//For 32-bit integers, the maximum n is 2^31-1=2147483647.
//If larger n is required, NTYPE must be double.

typedef int NTYPE;
#define PI2 6.283185307179586476925286
#define S0 0.083333333333333333333 // 1/12
#define S1 0.00277777777777777777778 // 1/360
#define S2 0.00079365079365079365079365 // 1/1260
#define S3 0.000595238095238095238095238 // 1/1680
#define S4 0.0008417508417508417508417508 // 1/1188
static double sfe[16] = {
	0, 0.081061466795327258219670264,
	0.041340695955409294093822081, 0.0276779256849983391487892927,
	0.020790672103765093111522771, 0.0166446911898211921631948653,
	0.013876128823070747998745727, 0.0118967099458917700950557241,
	0.010411265261972096497478567, 0.0092554621827127329177286366,
	0.008330563433362871256469318, 0.0075736754879518407949720242,
	0.006942840107209529865664152, 0.0064089941880042070684396310,
	0.005951370112758847735624416, 0.0055547335519628013710386899};

// stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
inline double stirlerr(NTYPE n) {
	double nn;
	if (n<16) return(sfe[(int)n]);
	nn = (double)n;
	nn = nn*nn;
	if (n>500) return((S0-S1/nn)/n);
	if (n>80) return((S0-(S1-S2/nn)/nn)/n);
	if (n>35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
	return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}

// Evaluate the deviance term
//bd0(x,np) = x log(x/np) + np - x

inline double bd0(NTYPE x, double np) {
	double ej, s, s1, v;
	int j;
	if (fabs(x-np)<0.1*(x+np)) {
		s = (x-np)*(x-np)/(x+np);
		v = (x-np)/(x+np);
		ej = 2*x*v;
		for (j=1; ;j++) {
			ej *= v*v;
			s1 = s+ej/(2*j+1);
			if (s1==s) return(s1);
			s = s1;
		}
	}
	return(x*log(x/np)+np-x);
}

double Random::dbinom(NTYPE x, NTYPE n, double p) {
	double lc;
	if (p==0.0) return( (x==0) ? 1.0 : 0.0);
	if (p==1.0) return( (x==n) ? 1.0 : 0.0);
	if (x==0) return(exp(n*log(1-p)));
	if (x==n) return(exp(n*log(p)));
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x)- bd0(x,n*p) - bd0(n-x,n*(1.0-p));
	return(exp(lc)*sqrt(n/(PI2*x*(n-x))));
}

double Random::dpois(NTYPE x, double lb) {
	if (lb==0) return( (x==0) ? 1.0 : 0.0);
	if (x==0) return(exp(-lb));
	return(exp(-stirlerr(x)-bd0(x,lb))/sqrt(PI2*x));
}

// end Loader algorithm


uint Random::rint()
{
	return randomMT();
}

uint Random::rint( uint max )
{
	return static_cast<uint>( max*randomMTunit() );
}

double Random::runif()
{
	return randomMTunit();
}


double Random::rnorm( double mean, double sd )
{
	return sd*snorm() + mean;
}

uint Random::rpois( double mean )
{
	// Draw from a Poisson Dist with mean
	// if cannot calculate, returns UINT_MAX
	// Uses Rejection Method
	assert( mean >= 0 );
	uint k = 0;
	double a = exp( -mean );
	double u = randomMTunit();
	assert( a > 0 ); // a should really be always larger than zero
	if( a<=0 ) return UINT_MAX; // however, if it is not, return UINT_MAX
	while( u>=a )
	{
		u *= randomMTunit();
		++k;
	}
	return k;
}

uint Random::randintFromDistr( double *distr, int n )
{
        double randVar = runif();
        int result, resultMin, resultMax;

        resultMin = 0;
        resultMax = n-1;
        do
        {
                result = ( resultMax - resultMin ) / 2 + resultMin;
                if ( 0 == result )
                {
                        if ( randVar < distr[0] ) break;
                        else resultMin = 1;
                }
                else if ( result == n - 1 )
                {
                        if ( randVar < distr[result-1] )
                                resultMax = result-1;
                        else break;
                }
                else if ( randVar < distr[result] )
                {
                        if ( randVar >= distr[result-1] )
                        {
                                break;
                        }
                        else resultMax = result - 1;
                }
                else if ( randVar < distr[result+1] )
                {
                        result++;
                        break;
                }
                else resultMin = result + 1;
        }
        while ( 1 );

        return result;
}




void Random::seed( uint i )
{
#ifdef DRAND48
	srand48( i );
#else
	seedMT( i );
#endif
}
