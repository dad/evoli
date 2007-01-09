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


#ifndef NUCLEOTIDE_SEQUENCE_HH
#define NUCLEOTIDE_SEQUENCE_HH

#include <vector>
#include <string>
#include <cassert>
#include "sequence.hh"

using namespace std;

/** \brief A class holding a protein-coding DNA sequence.
*/
class NucleotideSequence : public Sequence {
public:
	NucleotideSequence(const string& s) : Sequence(s) {}
	NucleotideSequence(unsigned int length) : Sequence (length) {}
	NucleotideSequence(unsigned int length, char val) : Sequence (length, val) {}
	NucleotideSequence(const NucleotideSequence& s, unsigned int start, unsigned int length) : Sequence(s,start,length) {}

	/**
	 *
	 **/
	double getGCFraction() const {
		int n_gc = 0;
		const Sequence& s = *this;
		for (unsigned int i=0; i<size(); i++) {
			if (s[i] == 'G' || s[i] == 'C') {
				n_gc++;
			}
		}
		return n_gc/double(length());
	}

};

#endif // NUCLEOTIDE_SEQUENCE_HH
