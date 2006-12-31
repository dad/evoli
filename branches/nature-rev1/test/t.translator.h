/*
This file is part of the E.voli project.
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


#ifndef _T_TRANSLATOR_H__
#define _T_TRANSLATOR_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "protein.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include <fstream>
#include <cmath>
#include <memory>

class TranslationTester : public ErrorproneTranslation { // Get access to protected fields.
public:
	TranslationTester(Folder* folder) : ErrorproneTranslation(folder, 25, 574, -5, 1, 6, 0.0114, 59, 104.5) {
	}
	//int Translator::translateRelativeWeighted( const Gene &g, Protein& residue_sequence, const double error_weight,
	//									   const vector<vector<pair<double, int> > >& weights, const double* prefCodons, 
	//									   const double nonPrefCodonPenalty, bool& truncated)

	const vector<vector<pair<double, int> > >& getWeights() {
		return m_cum_weight_matrix;
	}
	const double* getPrefCodons() {
		return m_codon_cost;
	}
};

struct TEST_CLASS( translator_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( native_sequence_truncated )
	{
		CompactLatticeFolder folder(side_length);
		Gene g("AUGCGUUAAGGG");
		Translator t(0);
		TranslationTester t_tester(&folder);
		Protein p(g.codonLength());
		bool truncated = false;
		int numErrors = -10;
		numErrors = t.translateRelativeWeighted(g, p, 50, t_tester.getWeights(), t_tester.getPrefCodons(), 6, truncated);
		TEST_ASSERT( numErrors == 0 );
		TEST_ASSERT( truncated );
		string prot_str = p.toString();
		TEST_ASSERT( prot_str[0] == 'M' );
		return;
	}
};

#endif
