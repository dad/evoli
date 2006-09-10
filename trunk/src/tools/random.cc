#include "random.hh"
#include <climits>
#include <cmath>
#include <cassert>

// Wrapper code around the implemtation by Shawn Cokus of the Mersenne
// Twister RNG by Makoto Matsumoto and Takuji Nishimura. See
//    http://www.math.keio.ac.jp/~matumoto/emt.html
// for more on the algorithm.
//
// Copyright (C) 2006 Claus O. Wilke, <cwilke@mail.utexas.edu>
//
// Part of this code was taken from the Octave modules for the Mersenne
// Twister (MT199337) Random Number Generator, Copyright (C) 1998,
// 1999 Dirk Eddelbuettel <edd@debian.org>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.



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
Routine that transforms U(0,UINT_MAX) into N(0,1), taken from 
the GNU GPL'ed randlib library by Brown, Lovato, Russell and Venier 
which is available from   ftp://odin.mdacc.tmc.edu/pub/source
*/
static double snorm(void)
/*
**********************************************************************

     (STANDARD-)  NORMAL  DISTRIBUTION

**********************************************************************

     FOR DETAILS SEE:
               AHRENS, J.H. AND DIETER, U.
               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
               SAMPLING FROM THE NORMAL DISTRIBUTION.
               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.

     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

     Modified by Dirk Eddelbuettel <edd@debian.org> on 1 Aug 1999
     to use randomMTunit() instead of RANF

**********************************************************************
*/
//     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
//     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
{
  static double a[32] = {
    0.0,3.917609E-2,7.841241E-2,0.11777,0.1573107,0.1970991,0.2372021,0.2776904,
    0.3186394,0.36013,0.4022501,0.4450965,0.4887764,0.5334097,0.5791322,
    0.626099,0.6744898,0.7245144,0.7764218,0.8305109,0.8871466,0.9467818,
    1.00999,1.077516,1.150349,1.229859,1.318011,1.417797,1.534121,1.67594,
    1.862732,2.153875
  };
  static double d[31] = {
    0.0,0.0,0.0,0.0,0.0,0.2636843,0.2425085,0.2255674,0.2116342,0.1999243,
    0.1899108,0.1812252,0.1736014,0.1668419,0.1607967,0.1553497,0.1504094,
    0.1459026,0.14177,0.1379632,0.1344418,0.1311722,0.128126,0.1252791,
    0.1226109,0.1201036,0.1177417,0.1155119,0.1134023,0.1114027,0.1095039
  };
  static double t[31] = {
    7.673828E-4,2.30687E-3,3.860618E-3,5.438454E-3,7.0507E-3,8.708396E-3,
    1.042357E-2,1.220953E-2,1.408125E-2,1.605579E-2,1.81529E-2,2.039573E-2,
    2.281177E-2,2.543407E-2,2.830296E-2,3.146822E-2,3.499233E-2,3.895483E-2,
    4.345878E-2,4.864035E-2,5.468334E-2,6.184222E-2,7.047983E-2,8.113195E-2,
    9.462444E-2,0.1123001,0.136498,0.1716886,0.2276241,0.330498,0.5847031
  };
  static double h[31] = {
    3.920617E-2,3.932705E-2,3.951E-2,3.975703E-2,4.007093E-2,4.045533E-2,
    4.091481E-2,4.145507E-2,4.208311E-2,4.280748E-2,4.363863E-2,4.458932E-2,
    4.567523E-2,4.691571E-2,4.833487E-2,4.996298E-2,5.183859E-2,5.401138E-2,
    5.654656E-2,5.95313E-2,6.308489E-2,6.737503E-2,7.264544E-2,7.926471E-2,
    8.781922E-2,9.930398E-2,0.11556,0.1404344,0.1836142,0.2790016,0.7010474
  };
  static long i;
  static double snorm,u,s,ustar,aa,w,y,tt;
  u = randomMTunit();
  s = 0.0;
  if(u > 0.5) s = 1.0;
  u += (u-s);
  u = 32.0*u;
  i = (long) (u);
  if(i == 32) i = 31;
  if(i == 0) goto S100;
  /*
    START CENTER
  */
  ustar = u-(double)i;
  aa = *(a+i-1);
 S40:
  if(ustar <= *(t+i-1)) goto S60;
  w = (ustar-*(t+i-1))**(h+i-1);
 S50:
  /*
    EXIT   (BOTH CASES)
  */
  y = aa+w;
  snorm = y;
  if(s == 1.0) snorm = -y;
  return snorm;
 S60:
  /*
    CENTER CONTINUED
  */
  u = randomMTunit();
  w = u*(*(a+i)-aa);
  tt = (0.5*w+aa)*w;
  goto S80;
 S70:
  tt = u;
  ustar = randomMTunit();
 S80:
  if(ustar > tt) goto S50;
  u = randomMTunit();
  if(ustar >= u) goto S70;
  ustar = randomMTunit();
  goto S40;
 S100:
  /*
    START TAIL
  */
  i = 6;
  aa = *(a+31);
  goto S120;
 S110:
  aa += *(d+i-1);
  i += 1;
 S120:
  u += u;
  if(u < 1.0) goto S110;
  u -= 1.0;
 S140:
  w = u**(d+i-1);
  tt = (0.5*w+aa)*w;
  goto S160;
 S150:
  tt = u;
 S160:
  ustar = randomMTunit();
  if(ustar > tt) goto S50;
  u = randomMTunit();
  if(ustar >= u) goto S150;
  u = randomMTunit();
  goto S140;
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
