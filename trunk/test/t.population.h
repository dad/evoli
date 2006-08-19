#ifndef _T_POPULATION_H__
#define _T_POPULATION_H__
#include "cutee.h"
#include "compact-lattice-folder.hh"
#include "decoy-contact-folder.hh"
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
		// of the incoming gene, which should be 0.691722
		TEST_ASSERT( abs(fe.getFitness( g )-0.691722) < 1e-4 );
		TEST_ASSERT( abs(pop.getAveFitness()-0.691722) < 1e-4 );
	}

	void TEST_FUNCTION( decoy_folder )
	{
		int N = 100;
		double U = .001;
		Population pop( N );

		int protein_length = 300;
		vector<DecoyContactStructure*> structs;
		ifstream fin("test/data/williams_contact_maps/maps.txt");
		TEST_ASSERT( fin.good() );
		if (!fin.good()) // if we can't read the contact maps file, bail out
			return;
		readContactMapsFromFile(fin, "test/data/williams_contact_maps/", structs);

		double log_nconf = 160 * log(10);
		DecoyContactFolder folder( protein_length, log_nconf, structs);
		ProteinFreeEnergyFitness fe( &folder );
		Gene g;
		FoldInfo fi;
		do {
			g = Gene::createRandomNoStops(protein_length*3);
			Protein p = g.translate();
			fi = folder.fold(p);
		} while (fi.getStructure() != (StructureID)0);

		double dG = fi.getFreeEnergy();
		double fitness = exp(-dG/0.6)/(1+exp(-dG/0.6));
		
		pop.init( g, &fe, U );
		// without evolution, mean fitness should equal the fitness
		// of the incoming gene, which should be 0.691722
		TEST_ASSERT( pop.getAveFitness() > fitness );
	}
};


#endif // _T_POPULATION_H__
