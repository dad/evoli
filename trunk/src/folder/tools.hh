#ifndef TOOLS_HH
#define TOOLS_HH

#include <cmath>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <vector>

using namespace std;

static const char* tab = "\t";

typedef unsigned int uint16;


/**
 * This function encapsulates the random number generator.
 **/
double myRand();

/**
 * Poissonian distribution with mean 'mean'.
 **/
int poissonian ( double mean );

/**
 * This function chooses randomly from a probability distribution
 * defined in the array 'bins'.
 **/
int chooseAtRandom( double *bins, int noOfBins );

/**
 * This function calculates mean and variance from the sums s1=sum_i^n x_i and s2=sum_i^n x_i^2.
 **/

pair<double, double> meanvar( double s1, double s2, int n );

/**
 * Computes the mean of v.
 **/
template <typename T> double mean(const vector<T>& v);

/**
 * Computes the variance of v.
 **/
template <typename T> double variance(const vector<T>& v);

/**
 * Computes the mean and variance of v.
 **/
template <typename T> pair<double,double> meanvar(const vector<T>& v);


//************************************************************
//
//			 Implementation
//
//************************************************************

template <class T> inline double mean(const vector<T>& v) {
	double sum = 0.0;
	uint16 count = 0;
	for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
		sum += (double)(*it);
		count++;
	}
	return sum/count;
}

template <class T> inline double variance(const vector<T>& v) {
	double sum = 0.0;
	double sumsq = 0.0;
	uint16 count = 0;
	for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
		double val = *it;
		sum += val;
		sumsq += val*val;
		count++;
	}
	// Sample variance is (1/(N-1)) sum_i (x_i - E(x))^2
	return (sumsq - (sum*sum)/count)/(count-1);
}

/**
 * Compute mean and variance.
 **/
template <class T> inline pair<double,double> meanvar(const vector<T>& v) {
	double sum = 0.0;
	double sumsq = 0.0;
	uint16 count = 0;
	for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
		double val = *it;
		sum += val;
		sumsq += val*val;
		count++;
	}
	// Sample variance is (1/(N-1)) sum_i (x_i - E(x))^2
	return pair<double,double>(sum/count, (sumsq - (sum*sum)/count)/(count-1));
}

inline double myRand()
{
        return drand48();
}


inline int poissonian ( double mean )
{
        // Draw from a Poisson Dist with mean
        // if cannot calculate, returns UINT_MAX
        // Uses Rejection Method
        int k=0;
        double a=exp(-mean);
        double u=myRand();
        if( a <=0 ) return UINT_MAX; // cannot calculate, so return UINT_MAX
        while( u>=a )
        {
                u*=myRand();
                ++k;
        }
        return k;
}


inline int chooseAtRandom( double *bins, int noOfBins )
{
        double randVar = myRand();
        int result, resultMin, resultMax;

        resultMin = 0;
        resultMax = noOfBins-1;
        do
        {
                result = ( resultMax - resultMin ) / 2 + resultMin;
                if ( 0 == result )
                {
                        if ( randVar < bins[0] ) break;
                        else resultMin = 1;
                }
                else if ( result == noOfBins - 1 )
                {
                        if ( randVar < bins[result-1] )
                                resultMax = result-1;
                        else break;
                }
                else if ( randVar < bins[result] )
                {
                        if ( randVar >= bins[result-1] )
                        {
                                break;
                        }
                        else resultMax = result - 1;
                }
                else if ( randVar < bins[result+1] )
                {
                        result++;
                        break;
                }
                else resultMin = result + 1;
        }
        while ( 1 );

        return result;
}


inline pair<double, double> meanvar( double s1, double s2, int n )
{
        pair<double, double> p;
        double x = (double) n;
        p.first = s1 / x;
        p.second = ( x * s2 - s1 * s1 )/( x * ( x - 1 ) );
        return p;
}

inline string itoa(int value, int base) {
	enum { kMaxDigits = 35 };
	string buf;
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.
	// Check that the base is valid
	if (base < 2 || base > 16) return buf;
	int quotient = value;
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ std::abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += '-';
	reverse( buf.begin(), buf.end() );
	return buf;
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

inline double dbinom(NTYPE x, NTYPE n, double p) {
	double lc;
	if (p==0.0) return( (x==0) ? 1.0 : 0.0);
	if (p==1.0) return( (x==n) ? 1.0 : 0.0);
	if (x==0) return(exp(n*log(1-p)));
	if (x==n) return(exp(n*log(p)));
	lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x)- bd0(x,n*p) - bd0(n-x,n*(1.0-p));
	return(exp(lc)*sqrt(n/(PI2*x*(n-x))));
}
inline double dpois(NTYPE x, double lb) {
	if (lb==0) return( (x==0) ? 1.0 : 0.0);
	if (x==0) return(exp(-lb));
	return(exp(-stirlerr(x)-bd0(x,lb))/sqrt(PI2*x));
}

// end Loader algorithm


#endif
