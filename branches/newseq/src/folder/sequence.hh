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


#ifndef SEQUENCE_HH
#define SEQUENCE_HH

#include <string>
#include <iostream>

using namespace std;


/** \brief Basic class to store a genetic sequence of any type (DNA, protein).
**/

//typedef string Sequence;

class Sequence : public string {
public:
	Sequence(const string& s) : string(s) {}
	Sequence(const char* s) : string(s) {}
	Sequence(const Sequence& s) : string(s) {}
	Sequence(unsigned int length, char val) : string(length, val) {}
	Sequence(unsigned int length) : string(length, 'A') {}
	Sequence(const Sequence& s, unsigned int start, unsigned int length) : string(s.substr(start,length)) {
		//cout << endl << "wakka: " << start << " " << length << " " << *this << endl;
	}
	
	const string& toString() { return *this; }
};

#endif //SEQUENCE_HH
