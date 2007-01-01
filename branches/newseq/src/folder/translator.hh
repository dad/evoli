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

#ifndef TRANSLATOR_HH
#define TRANSLATOR_HH

#include <vector>
#include "protein.hh"

class Translator
{
private:
	double m_mutation_prob;

	Translator( const Translator & );
	Translator& operator=( const Translator & );
public:
	Translator( double mutation_prob);
	Translator();

	bool translateErrorFree( const CodingRNA &g, Protein& residue_sequence ) const;
	bool translate( const CodingRNA &g, Protein& residue_sequence ) const;

	int translateWeighted( const CodingRNA &g, Protein& residue_sequence, const vector<vector<pair<double, int> > >& weights,
						   const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);
	int translateRelativeWeighted( const CodingRNA &g, Protein& residue_sequence, const double relative_gene_weight,
								   const vector<vector<pair<double, int> > >& weights, const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);
};



#endif
