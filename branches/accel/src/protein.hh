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


#ifndef PROTEIN_HH
#define PROTEIN_HH

#include <vector>
#include <string>
#include <cassert>
#include "sequence.hh"
#include "codon.hh"

using namespace std;

typedef pair<int, int> Contact;

/** \brief A class holding a protein sequence.
*/
class Protein : public Sequence {
protected:
public:
	Protein(unsigned int length);
	Protein(unsigned int length, char val);
	Protein(const string& s) : Sequence(s) {}
	~Protein() {}

	/**
	Converts the protein sequence into a string useful for output etc.
	@return A string object containing the protein sequence, using the single-letter amino-acid alphabet.
	*/
	string toString() const { return *this; }

	/**
	Creates a random protein.
	@param length Length of the desired protein, in residues
	@return The random protein.
	*/
	static Protein createRandom(unsigned int length);

};


#endif //PROTEIN_HH
