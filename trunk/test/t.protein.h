#ifndef _T_PROTEIN_H__
#define _T_PROTEIN_H__
#include "cutee.h"
#include "protein.hh"
#include "gene-util.hh"
#include "compact-lattice-folder.hh"

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
		CompactLatticeFolder* folder = new CompactLatticeFolder(side_length);
		double max_dg = -5;
		double sid = 574;
		Gene g = GeneUtil::getSequenceForStructure(*folder, gene_length, max_dg, sid);
		Protein p = g.translate();
		FoldInfo fi = folder->fold(p);
		TEST_ASSERT( fi.getFreeEnergy() <= max_dg );
		TEST_ASSERT( fi.getStructure() == (StructureID)sid );
		delete folder;
		return;
	}
};


#endif //_T_PROTEIN_H__
