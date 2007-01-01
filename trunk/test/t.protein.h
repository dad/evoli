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


#ifndef _T_PROTEIN_H__
#define _T_PROTEIN_H__
#include "cutee.h"
#include "protein.hh"
#include "gene-util.hh"
#include "compact-lattice-folder.hh"
#include <iostream>
using namespace std;

struct TEST_CLASS( protein_gene_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( length_test )
	{
		Gene g = Gene::createRandom(gene_length);
		TEST_ASSERT( g.length()==gene_length );
		return;
	}
	void TEST_FUNCTION( reverse_translate )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		Protein p2 = GeneUtil::reverseTranslate(p).translate();
		TEST_ASSERT( p == p2 );
		return;
	}
	void TEST_FUNCTION( gene_string_gene )
	{
		Gene g = Gene::createRandom(gene_length);
		string str(g.toString());
		Gene g2(str);
		TEST_ASSERT( g == g2 );
		return;
	}
	void TEST_FUNCTION( translate_known )
	{
		Gene g("ATGTGGGGG");
		Protein p = g.translate();
		string str(p.toString());
		TEST_ASSERT( str == "MWG" );
		return;
	}
	void TEST_FUNCTION( full_length_true )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		TEST_ASSERT( g.encodesFullLength() );
		return;
	}
	void TEST_FUNCTION( full_length_false )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		int codon = CodonUtil::lettersToCodon('U','A','G');
		g[4] = codon;
		TEST_ASSERT( !g.encodesFullLength() );
		return;
	}
	void TEST_FUNCTION( full_length_protein )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		TEST_ASSERT( g.codonLength() == p.length() );
		return;
	}
	void TEST_FUNCTION( protein_from_string )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		string s = p.toString();
		Protein p2(s);
		TEST_ASSERT( p2 == p );
		return;
	}

	void TEST_FUNCTION( protein_from_strings_equality ) {
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p1 = g.translate();
		Protein p2(p1.toString());
		TEST_ASSERT( p1 == p2 );
	}
	void TEST_FUNCTION( sequence_for_structure )
	{
		CompactLatticeFolder folder(side_length);
		double max_dg = -1;
		int sid = 574;
		Gene g = GeneUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( folder.fold(p) );
		TEST_ASSERT( fi->getDeltaG() <= max_dg );
		TEST_ASSERT( fi->getStructure() == (StructureID)sid );
		return;
	}
};


#endif //_T_PROTEIN_H__
