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


#ifndef CODON_HH
#define CODON_HH

#include "random.hh"

#include <iostream>
#include <cassert>

using namespace std;

class CodonUtil
{
public:
	/**
	* Outputs the given codon in three letter format to the given ostream.
	**/
	static void printCodon( ostream &s, int codon )
	{
		assert( codon < 64 );

		for ( int i=0; i<3; i++ )
		{
			int x = (codon >> (2*(2-i))) & 3;
			switch( x )
			{
			case 0:
				s << 'A';
				break;
			case 1:
				s << 'C';
				break;
			case 2:
				s << 'G';
				break;
			case 3:
				s << 'U';
				break;
			}
		}
	}


	/**
	* Translates a base (A, C, G, U/T) into its numeric code (0, 1, 2, 3).
	**/
	static int baseToInt( char c )
	{
		switch( c )
		{
		case 'A':
		case 'a':
			return 0;
		case 'C':
		case 'c':
			return 1;
		case 'G':
		case 'g':
			return 2;
		default:
			return 3;
		}
	}

	/**
	* Translates a numeric code (0, 1, 2, 3) into the corresponding base (A, C, G, U).
	**/
	static int intToBase( int c )
	{
		switch( c )
		{
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
				return 'U';
		default:
			return '?';
		}
	}

	/**
	* Tests whether a mutations from base 1 to base 2 is a transition (A<->G, C<->U) or not.
	* The integers are mapped to bases in this way:
	*  A: 0, C: 1, G: 2, U: 3
	**/
	static bool isTransition( int l1, int l2 )
	{
		if ( l1 == 0 && l2 == 2 )
			return true;
		if ( l1 == 1 && l2 == 3 )
			return true;
		if ( l1 == 2 && l2 == 0 )
			return true;
		if ( l1 == 3 && l2 == 1 )
			return true;
		return false;
	}

	/**
	* Tests whether codon co can be turned into codon ct in one nucleotide change.
	**/
	static bool isPointMutant( int co, int ct ) {
			if ( co == ct )
				return false;

			int lo1, lo2, lo3;
			int lt1, lt2, lt3;

			CodonUtil::codonToLetters( lo1, lo2, lo3, co );
			CodonUtil::codonToLetters( lt1, lt2, lt3, ct );

			bool res = (( lo1 != lt1 && lo2 == lt2 && lo3 == lt3 ) ||
					( lo2 != lt2 && lo1 == lt1 && lo3 == lt3 ) ||
					( lo3 != lt3 && lo1 == lt1 && lo2 == lt2 ));
			//cout << lo1 << lo2 << lo3 << tab << lt1 << lt2 << lt3 << tab << res << endl;
			return res;
		}


	/**
	* Combines the three given letters into a single codon.
	* The integers are mapped to bases in this way:
	*  A: 0, C: 1, G: 2, U: 3
	**/
	static int lettersToCodon( int l1, int l2, int l3 )
	{
		return 16*l1 + 4*l2 + l3;
	}

	/**
	* Combines the three given letters into a single codon.
	* Takes characters A, C, G, U/T instead of integer values
	**/
	static int lettersToCodon( char c1, char c2, char c3 )
	{
		int l1 = baseToInt( c1 );
		int l2 = baseToInt( c2 );
		int l3 = baseToInt( c3 );

		return 16*l1 + 4*l2 + l3;
	}
	/**
	* Separates a given codon into its three letters.
	**/
	static void codonToLetters( int &l1, int &l2, int &l3, int codon )
	{
		l1 = (codon >> 4) & 3;
		l2 = (codon >> 2) & 3;
		l3 = codon & 3;
	}

	/**
	* Randomly replaces the bases in the codon with other bases. The probability
	* that a single base is changed is given in the parameter U. All alternative
	* bases are equally likely to be substituted.
	**/
	static int mutateCodon( double U, int codon )
	{
		int result = 0;
		for ( int i=0; i<3; i++ )
		{
			result = result << 2;
			int x = (codon >> (2*(2-i))) & 3; // separate base
			if ( Random::runif() < U )
			{ // do we have to mutate this base?
				x = ( x + Random::rint( 3 ) + 1 ) & 3; // mutate to a random other base
			}
			result += x; // add to the result
		}
		return result;
	}
};



#endif









