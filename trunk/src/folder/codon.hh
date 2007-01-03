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

#include "sequence.hh"

#include <iostream>
#include <cassert>

using namespace std;

class Codon : public Sequence {
public:
	Codon() : Sequence("AAA") {}
	Codon(const Sequence& s, unsigned int start) : Sequence(s, start, 3) {
	}
	Codon(const string&s) : Sequence(s) {}
	Codon(const char*s) : Sequence(s) {}

	/**
	 * Transcribes DNA pseudo-codon into RNA codon.
	 * @return A Codon, with instances of T in the original Codon replaced by U.
	 **/
	Codon transcribe() const;
};

#endif // _CODON_HH__









