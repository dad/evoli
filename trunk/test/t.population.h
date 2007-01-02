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


#ifndef _T_POPULATION_H__
#define _T_POPULATION_H__
#include "cutee.h"

#include "random.hh"
#include "compact-lattice-folder.hh"
#include "decoy-contact-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "mutator.hh"

struct TEST_CLASS( population )
{
	const static int side_length = 4;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( init )
	{
		int N = 20;
		double U = .001;
		CompactLatticeFolder folder( 4 );
		ProteinFreeEnergyFitness fe( &folder );
		Population<Gene, ProteinFreeEnergyFitness, SimpleMutator>  pop( N );
		SimpleMutator mut(U);
		Gene g( "AAAAAAAAGAGUCCUACCACCCUUGACCUCAUGUCCUGUGCAGAUAAU" );
		
		pop.init( g, &fe, &mut );
		// without evolution, mean fitness should equal the fitness
		// of the incoming gene, which should be 0.691722
		TEST_ASSERT( fabs(fe.getFitness( g )-0.691722) < 1e-4 );
		TEST_ASSERT( fabs(pop.getAveFitness()-0.691722) < 1e-4 );
	}

	void TEST_FUNCTION( num_generations )
	{
		int N = 20;
		double U = .001;
		CompactLatticeFolder folder( 4 );
		ProteinFreeEnergyFitness fe( &folder );
		Population<Gene, ProteinFreeEnergyFitness, SimpleMutator>  pop( N );
		SimpleMutator mut(U);
		Gene g( "AAAAAAAAGAGUCCUACCACCCUUGACCUCAUGUCCUGUGCAGAUAAU" );
		
		pop.init( g, &fe, &mut );
		int num_generations = 102;
		for ( int j=0; j<num_generations; j++ ) {
			pop.evolve();
		}
		TEST_ASSERT( pop.getNumGenerations() == num_generations );
	}

	void TEST_FUNCTION( lattice_folder )
	{
		int N = 50;
		int size = 4;
		int length = size*size;
		double U = 0.01;
		Random::seed( 37 ); // initialize random number generator
		CompactLatticeFolder b(size);
		ProteinFreeEnergyFitness fe( &b );
		Population<Gene, ProteinFreeEnergyFitness, SimpleMutator>  pop( N );
		SimpleMutator mut(U);
		Gene g( "GGGAAGUGCGUCCAGCAGAGUUGGGUAUGGGAGGGAUCUAAGUUAAAG" );
		//Gene g = Gene::createRandom( length*3 );
		bool randomOK = true;
		TEST_ASSERT( randomOK = ( g == Gene( "GGGAAGUGCGUCCAGCAGAGUUGGGUAUGGGAGGGAUCUAAGUUAAAG" ) ) );
		//std::cout << g << std::endl;
		if ( !randomOK )
			cout << "Test failures in function lattice_folder likely due to differences in random number generator" << endl;
		pop.init( g, &fe, &mut );
		GenebankAnalyzer<Gene> analyzer(pop.getGenebank());

		int equil_time = 100;
		int window_size = 30;
		for ( int i=0; ; i++ )
		{
			for ( int j=0; j<100; j++ )
			{
				pop.evolve();
			}
			analyzer.prepareCoalescenceCalcs(pop.begin(), pop.end(), pop.getNumGenerations());
			if ( analyzer.calcCoalescenceTime() > equil_time + window_size )
				break;
		}
//		p.printGenebank( cout );
		double ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop;
		vector<bool> is_optimal;
		for( int i=0; i<64; i++ )
			is_optimal.push_back( false );

		analyzer.analyzeDnDs( window_size, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, is_optimal );

		std::cout << ave_dn << " " << ave_ds << " " << ave_N << " " << ave_S << " " << ave_f << " " << ave_fop << std::endl;

		TEST_ASSERT( fabs( ave_dn - 7) < 1e-4 );
		TEST_ASSERT( fabs( ave_ds - 1) < 1e-4 );
		TEST_ASSERT( fabs( ave_N - 37.7667) < 1e-4 );
		TEST_ASSERT( fabs( ave_S - 10.2333) < 1e-4 );
		TEST_ASSERT( fabs( ave_f - 0.998525) < 1e-4 );
		TEST_ASSERT( fabs( ave_fop - 0) < 1e-10 );
	}

	void TEST_FUNCTION( decoy_folder )
	{
		int protein_length = 300;
		int N = 100;
		double U = 1.0/(N*protein_length);
		Population<Gene, ProteinFreeEnergyFitness, SimpleMutator>  pop( N );
		SimpleMutator mut(U);
		string fname = "test/data/rand_contact_maps/maps.txt";
		ifstream fin(fname.c_str());
		string dir = "test/data/rand_contact_maps/";
		double log_nconf = 10 * log(10.0);
		DecoyContactFolder folder( protein_length, log_nconf, fin, dir);
		TEST_ASSERT( folder.good() );
		if (!folder.good() ) {
			cout << "# Couldn't initialize folder with " << fname << endl;
			return;
		}
		ProteinFreeEnergyFitness fe( &folder );
		Gene g;
		auto_ptr<FoldInfo> fi;
		StructureID target_struct = (StructureID)0;
		do {
			g = Gene::createRandomNoStops(protein_length*3);
			Protein p = g.translate();
			fi = auto_ptr<FoldInfo>( folder.fold(p) );
		} while (fi->getStructure() != target_struct);

		double fitness = fe.getFitness(g);
		
		pop.init( g, &fe, &mut );
		// with evolution, mean fitness should be greater than the fitness
		// of the incoming gene.
		
		for (int i=0, j=0; i<100; i++, j++) {
			pop.evolve();
			/*if (j>=100){
				cout << i << "\t" << pop.getAveFitness() << endl;
				j = 0;
				}*/
		}
		
		TEST_ASSERT( pop.getAveFitness() > fitness );
	}
};


#endif // _T_POPULATION_H__
