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


/** \page struct-driver struct-driver
The program \c struct-drivers simulates evolution with under a model
of constrained structure.  If a structure matches the target structure,
the n the fitness of its genotype is 1.  Otherwise, the fitness is zero.

   ./struct-driver 25 1000 599 -5 -100 0.00001 100000 50000 50 111
                   0  1    2    3  4   5       6      7     8  9
\endverbatim
(The numbers under the call count the parameters.)

The parameters are, in order (the zeroth parameter is the program name):
-#  %Protein length in amino acids
-#  %Population size
-#  Structure ID -- index of target structure
-#  Maximum delta free energy of folding
-#  Minimum delta free energy of folding
-#  Mutation rate per nucleotide
-#  Window time -- number of generations for which evolutionary data is being collected
-#  Equilibration time -- evolution time in generations before collection of evolutionary data begins
-#  Number of replicates
-#  %Random number seed -- best to use an odd number

*/

#include "folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>

using namespace std;

class Parameters {
public:
	unsigned int protein_length;
	double free_energy_cutoff;
	double free_energy_minimum;
	double u;
	int repetitions;
	int window_size;
	int equilibration_time;
	int random_seed;
	int N;
	mutable int structure_ID;
	string run_id;
	bool valid;

	Parameters( int ac, char **av ) {
		
		if ( ac < 11 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <prot length> <pop size> <structure id> <free energy cutoff> ";
			cout << "<free energy minimum> <mutation rate> <window time> <equilibration time> <repetitions> <random seed> <run ID>" << endl;
			return;
		}

		int i = 1;
		protein_length = atoi( av[i++] );
		N = atoi( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		u = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );
		if ( ac == 12 ) {
			run_id = av[i++];
		}
		else {
			run_id = itoa(random_seed, 10);
		}

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p );

int getStructureID( Folder &b, const Gene &g );
bool analyzeReplica( ProteinStructureFitness *fe, const Parameters &p, ostream &s );
void evolutionExperiment( const Parameters &p, ProteinStructureFitness& fe);

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   structure_ID: " << p.structure_ID << endl;
	s << "#   free-energy max: " << p.free_energy_cutoff << endl;
	s << "#   free-energy min: " << p.free_energy_minimum << endl;
	s << "#   per-site mutation rate u: " << p.u << endl;
	s << "#   window size tau: " << p.window_size << endl;
	s << "#   equilibration time: " << p.equilibration_time << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   run ID: " << p.run_id << endl;
	s << "#" << endl;
	return s;
}


StructureID getStructureID( Folder &b, const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( b.fold(p) );
		return fi->getStructure();
	}
	else
		return (StructureID)-1;
}

bool runAndAnalyzeReplica( ProteinStructureFitness *fe, const Parameters &p, ostream &s ) {
	
	// initialize the population
	Population pop( p.N );

	Folder& folder = *(fe->getFolder());
	
	// Find a sequence.
	Gene g = GeneUtil::getSequenceForStructure(folder, p.protein_length * 3, p.free_energy_cutoff, p.structure_ID );
	s << "# Starting genotype: " << g << endl;
	
	// Fill the population with the genotype that we found above
	pop.init( g, fe, p.u );

	int loop_length = 100;
	int n_folded = 0;
	for ( int i = 0; ; i++ ) {
		
		for ( int j = 0; j<loop_length; j++ ) {
			pop.evolve();
		}
		
		pop.prepareCoalescenceCalcs();
		if ( pop.calcCoalescenceTime() > p.equilibration_time + p.window_size )
			break;
		
		cout << "t=" << i*loop_length + 1 << "; " << flush;
		
	}
	
	cout << endl;

	pop.printGenebank( s );

}

void evolutionExperiment( const Parameters &p, ProteinStructureFitness& fe ) {
	
	// set random seed
	Random::seed( p.random_seed );

	int count = 0;
	for ( int i = 0; i < p.repetitions; i++ ) {
		
		string filename;
		stringstream repfname;
		repfname << "struct-run-gb-rep" << i << "-id" << p.run_id << ".dat";
		filename = repfname.str();
		ofstream gen_file( filename.c_str(), ios::out );
		gen_file << p;

		runAndAnalyzeReplica( &fe, p, gen_file );
		
		cout << "[" << i+1 << "/" << p.repetitions << "] " << endl;
		
	}
	
	cout << endl;
	
}

int main( int ac, char **av)
{
	Parameters p( ac, av );

	if ( ! p.valid ) {
		exit( 1 );
	}

	// seed the random number generator
	Random::seed( p.random_seed );

	// initialize the protein folder
	int side_length = (int)(sqrt(float(p.protein_length)));
	assert(side_length*side_length == p.protein_length);
	CompactLatticeFolder folder(side_length);

	// Print the command-line parameters
	cout << p;
	
	// Construct FitnessEvaluator
	ProteinStructureFitness* fe = new ProteinStructureFitness( &folder, p.structure_ID, p.free_energy_cutoff );

	if ( ! fe ) {
		cerr << "Problem creating the fitness evaluator.  Exiting..." << endl;
		exit(1);
	}

	cout << setprecision(4);
	evolutionExperiment( p, *fe );
	//evolutionTest( p, *fe );
	delete fe;

	return 0;
}




