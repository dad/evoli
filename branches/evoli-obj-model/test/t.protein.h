#ifndef _T_PROTEIN_H__
#define _T_PROTEIN_H__
#include "cutee.h"
#include "protein.hh"

struct TEST_CLASS( genotype_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;
	void TEST_FUNCTION( length_test )
	{
		Gene g = Gene::createRandom(gene_length);
		TEST_ASSERT( g.length()==gene_length );
		return;
	}
	void TEST_FUNCTION( translate_test )
	{

		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		Protein p2 = p.reverseTranslate().translate();
		TEST_ASSERT( p == p2 );
		return;
	}
};


#endif //_T_GENOTYPE_H__

