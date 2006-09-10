#ifndef TOOLS_HH
#define TOOLS_HH

#include <cmath>
#include <climits>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

static const char* tab = "\t";

typedef unsigned int uint16;



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



inline pair<double, double> meanvar( double s1, double s2, int n )
{
        pair<double, double> p;
        double x = (double) n;
        p.first = s1 / x;
        p.second = ( x * s2 - s1 * s1 )/( x * ( x - 1 ) );
        return p;
}

/**
 * Converts an integer into a character string.
 * Numbers must be no longer than 35 characters, and base must be between 2 (binary) and 16 (hexadecimal).
 * 
 * @param value Integer to convert into a string.
 * @param base Numerical base for integer conversion (e.g. 10=decimal, 2=binary).
 * @return A string representation of this integer.
 **/
inline string itoa(int value, int base) {
	enum { kMaxDigits = 35 };
	string buf;
	buf.reserve( kMaxDigits ); // Pre-allocate enough space.
	// Check that the base is valid
	if (base < 2 || base > 16) return buf;
	int quotient = value;
	// Translating number to string with base:
	do {
		buf += "0123456789abcdef"[ abs( quotient % base ) ];
		quotient /= base;
	} while ( quotient );
	// Append the negative sign for base 10
	if ( value < 0 && base == 10) buf += '-';
	reverse( buf.begin(), buf.end() );
	return buf;
}


#endif
