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

#include "codon.hh"

int Codon::baseToInt( char c )
{
	switch( c )
	{
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'U':
	case 'T':
		return 3;
	default:
		assert( false );
	}
	return 0;
}


char Codon::intToBase( int c, bool RNA )
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
		if ( RNA ) return 'U';
		else return 'T';
	default:
		assert( false ); // should never get here.
	}
	return 0;
}

int Codon::codonToIndex( const Codon& codon )
{
	int l1 = baseToInt( codon[0] );
	int l2 = baseToInt( codon[1] );
	int l3 = baseToInt( codon[2] );

	return 16*l1 + 4*l2 + l3;
}

Codon Codon::indexToCodon( int index, bool RNA )
{
	assert (index >= 0 && index < 64);	

	/* Do we want indices >63 to indicate DNA codons?
	if ( index > 63 )
	{
		RNA = false;
	}
	*/

	char codon[4];
	codon[0] = intToBase( (index >> 4) & 3, RNA );
	codon[1] = intToBase( (index >> 2) & 3, RNA );
	codon[2] = intToBase( index & 3, RNA );
	codon[3] = 0;
	return Codon( codon );
}

Codon Codon::transcribe() const {
	Codon res = *this;
	for (int i=0; i<3; i++) {
		if (res[i] == 'T')
			res[i] = 'U';
	}
	return res;
}

