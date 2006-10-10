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

#include <vector>

using namespace std;


/** \brief Basic class to store a genetic sequence of any type (DNA, protein).
**/

class Sequence {
private:
	Sequence();
protected:
	vector<int> m_sequence;

public:
	/**
	Copy constructor.
	*/
	Sequence( const Sequence &s ) {
		m_sequence = s.m_sequence;
	}

	/**
	Copy operator.
	*/
	Sequence & operator=( const Sequence &s ) {
		m_sequence = s.m_sequence;
		return m_sequence;
	}

	/**
	Another constructor. Can construct a sequence from a vector of ints.
	*/
	Sequence(const vector<int>& v) {
		m_sequence = v;
	}

	/**
	Another constructor. This one creates an empty sequence of a given length.
	*/
	Sequence(const int length) : m_sequence(length) {}
 
	virtual ~Sequence(void) {}

	/**
	@return The length of the sequence.
	*/
	virtual uint length() const {
		return m_sequence.size();
	}

	int& operator[](const int index) {
		return m_sequence[index];
	}

	int operator[](const int index) const {
		return m_sequence[index];
	}

	typedef vector<int>::iterator iterator;
	typedef vector<int>::const_iterator const_iterator;

	iterator begin() { return m_sequence.begin(); }
	iterator end() { return m_sequence.end(); }

	const_iterator begin() const { return m_sequence.begin(); }
	const_iterator end() const { return m_sequence.end(); }


	/**
	Comparison operator. Two sequences are equal if they have the same
	length and exactly the same contents.
	@return A boolean indicating whether the two sequences are identical.
	*/
	bool operator==(const Sequence& s) const {
		bool identical = (length() == s.length());
		const_iterator qit = begin();
		const_iterator pit = s.begin();
		for (; qit != end() && pit != s.end() && identical; qit++, pit++) {
			identical = (*pit == *qit);
		}
		return identical;
	}
};

#endif //SEQUENCE_HH
