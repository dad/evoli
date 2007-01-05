/*
This file is part of the E.voli project.  Copyright (C) 2004, 2005,
2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>, Allan Drummond
<drummond@alumni.princeton.edu>

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


#ifndef _CODON_HH__
#define _CODON_HH__

#include "nucleotide-sequence.hh"

#include <iostream>
#include <cassert>

using namespace std;

class Codon : public NucleotideSequence {

public:
	/**
	Default constructor. Constructs a codon of all A's.
	**/
	Codon() : NucleotideSequence("AAA") {}
	/**
	Extracts a codon from a sequence at a given position.
	@param s The sequence from which the codon is extracted.
	@param start The position in the sequence from where the codon should be extracted.
	*/
	Codon(const NucleotideSequence& s, unsigned int start) : NucleotideSequence(s, start, 3) {
	}

	/**
	Builds a codon from a C++ string of length 3.
	*/
	Codon(const string&s) : NucleotideSequence(s) {
		assert( size() == 3 );
	}

	/**
	Builds a codon from a C string of length 3.
	*/
	Codon(const char*s) : NucleotideSequence(s) {
		assert( size() == 3 );
	}

	/**
	Builds a codon from an index, as explained in \ref codonToIndex().
	*/
	Codon( int index, bool RNA = true ) : NucleotideSequence( "AAA" ) {
		*this = indexToCodon( index, RNA );
	}

	/**
	* Translates a base (A, C, G, U/T) into its numeric code (0, 1, 2, 3).
	**/
	static int baseToInt( char c );
        /**
	* Translates a numeric code (0, 1, 2, 3) into the corresponding base (A, C, G, U/T).
	* @param c The numeric code representing the base
	* @param RNA A bool indicating whether we are interested in DNA or RNA codons (i.e., U or T).
	* @return The correct base. 
	**/
	static char intToBase( int c, bool RNA=true );

	/**
	* Retrieve an integer index for a codon.
	*
	* We represent each base by two binary digits:
	*   A: 00, C: 01, G: 10, U/T: 11
	* Then, each codon is an integer between 0 and 63. For example,
	* GUG = 101110 = 46.
	* @return The integer value corresponding to the codon.
	**/
	static int codonToIndex( const Codon& codon );

	/**
	* Turn the integral representation of a codon back into a codon object.
	* @param index The index, coded as described in \ref codonToIndex().
	* @param RNA A bool indicating whether we are interested in DNA or RNA codons (i.e., U or T).
	* @return The codon corresponding to the given index.
	**/
	static Codon indexToCodon( int index, bool RNA = true );

	/**
	Returns the index of the codon, as explained in \ref codonToIndex().
	*/
	int getIndex() {
		return codonToIndex( *this ); }


	/**
	 * Transcribes DNA pseudo-codon into RNA codon.
	 * @return A Codon, with instances of T in the original Codon replaced by U.
	 **/
	Codon transcribe() const;
};

#endif // _CODON_HH__









