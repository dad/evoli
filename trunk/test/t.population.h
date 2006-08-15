#ifndef _T_POPULATION_H__
#define _T_POPULATION_H__
#include "cutee.h"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include <iostream>

struct TEST_CLASS( population )
{
	const static int side_length = 4;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( init )
	{
		int N = 20;
		double U = .001;
		Population pop( N );
		CompactLatticeFolder folder( 4 );
		folder.enumerateStructures();
		ProteinFreeEnergyFitness fe( &folder );
		Gene g( "AAAAAAAAGAGUCCUACCACCCUUGACCUCAUGUCCUGUGCAGAUAAU" );
		
		pop.init( g, &fe, U );
		// without evolution, mean fitness should equal the fitness
		// of the incoming gene, which should be 1.62403
		TEST_ASSERT( abs(fe.getFitness( g )-1.62403) < 1e-4 );
		TEST_ASSERT( abs(pop.getAveFitness()-1.62403) < 1e-4 );
	}
};


#endif // _T_POPULATION_H__
