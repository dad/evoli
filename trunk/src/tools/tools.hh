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
typedef unsigned int uint;

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

/** \brief A class to keep track of a running sum or average.
*/
class Accumulator {
private:
	double m_running_sum;
	double m_running_sum_squared;
	int m_count;
public:
	Accumulator() { m_running_sum=0.0; m_running_sum_squared=0.0; m_count=0; }
	Accumulator(double sum, double sum_squared, int count) { m_running_sum = sum; m_running_sum_squared = sum_squared; m_count = count; }

	/**
	@return The running average over all values added
	*/
	double value() const { return m_running_sum/m_count; }

	/**
	Adds an individual value to the running sum.
	@param val The value to be added.
	*/
	void add(double val) { m_running_sum += val; m_running_sum_squared += val*val; m_count++; }

	/**
	Sets the running sum and the number of values added to zero.
	*/
	void reset()
	{
		m_running_sum = 0.0;
		m_running_sum_squared = 0.0;
		m_count = 0;
	}

	/**
	@return The number of values added so far.
	*/
	int count() const { return m_count; }

	/**
	@return The running sum.
	*/
	double sum() const { return m_running_sum; }

	/**
	@return The mean of the summands.
	*/
	double mean() const { return value(); }

	/**
	@return The standard deviation of the summands.
	*/
	double stdev() const
	{
		pair<double,double> mv = meanvar(m_running_sum, m_running_sum_squared, m_count);
		return sqrt(mv.second);
	}

	/**
	@return The standard error of the summands.
	*/
	double stderror() const
	{
		pair<double,double> mv = meanvar(m_running_sum, m_running_sum_squared, m_count);
		return sqrt(mv.second/static_cast<double>(m_count));
	}

	/**
	Operator version of the function \ref value().
	*/
	operator double() { return value(); }

	void operator+=(double x) { add(x); }
};




//************************************************************
//
//			 Implementation
//
//************************************************************

template <class T> inline double mean(const vector<T>& v) {
	double sum = 0.0;
	uint count = 0;
	for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
		sum += (double)(*it);
		count++;
	}
	return sum/count;
}

template <class T> inline double variance(const vector<T>& v) {
	double sum = 0.0;
	double sumsq = 0.0;
	uint count = 0;
	for (typename vector<T>::const_iterator it = v.begin(); it != v.end(); it++) {
		double val = *it;
		sum += val;
		sumsq += val*val;
		count++;
	}
	// Sample variance is (1/(N-1)) sum_i (x_i - E(x))^2
	return (sumsq - (sum*sum)/count)/(count-1);
}

template <class T> inline pair<double,double> meanvar(const vector<T>& v) {
	double sum = 0.0;
	double sumsq = 0.0;
	uint count = 0;
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
