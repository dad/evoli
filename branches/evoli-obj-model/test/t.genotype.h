#ifndef _T_GENOTYPE_H__
#define _T_GENOTYPE_H__
#include "cutee.h"
#include "genotype.hh"

struct TEST_CLASS( genotype_basic )
{
	void TEST_FUNCTION( length_test )
	{
		unsigned int len = 75;
		Genotype g(len);
		TEST_ASSERT( g.size()==len );
		return;
	}
};


#endif //_T_GENOTYPE_H__

