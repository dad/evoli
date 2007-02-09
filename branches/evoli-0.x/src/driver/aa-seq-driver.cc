/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>, Matt Cowperthwaite <mattccowp@mac.com>

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


/** \page aaseq-driver aaseq-driver
The program \c aa-seq-driver is used to simulate evoltion with a constrained amino
acid sequence.  The fitness can be set to scale with the hamming distance from the 
target amino-acid sequence.  The fourth argument is a scaling factor, which must be
between 0 & 1.

The fitness is (1-s)^d, where s is the scaling factor and d is the hamming distance
between the two amino acids sequences.  If s == 1, then anything that does not perfectly
match has fitness 0.

The program is called with the following command-line parameters:
\verbatim
   ./aaseq-driver 25 1000 0.00001 0.0   100000 50000 50 111
                  0  1    2       3     4      5     6  7
\endverbatim
(The numbers under the call count the parameters.)

The parameters are, in order (the zeroth parameter is the program name):
-#  %Protein length in amino acids
-#  %Population size
-#  Mutation rate per nucleotide
-#  Scaling factor for distance
-#  Window time -- number of generations for which evolutionary data is being collected
-#  Equilibration time -- evolution time in generations before collection of evolutionary data begins
-#  Number of replicates
-#  %Random number seed -- best to use an odd number
-#  Run ID (optional)

*/

#include "compact-lattice-folder.hh"
#include "folder.hh"
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
	double u, dist_scale_factor;
	int repetitions;
	int window_size;
	int equilibration_time;
	int random_seed;
	int N;
	string run_id;
	bool valid;

	Parameters( int ac, char **av ) {
		
		if ( ac < 9 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <prot length> <pop size> <mutation rate> <scale factor>";
			cout << " <window time> <equilibration time> <repetitions> <random seed> <run ID>" << endl;
			return;
		}

		int i = 1;
		protein_length = atoi( av[i++] );
		N = atoi( av[i++] );
		u = atof( av[i++] );
		dist_scale_factor = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );
		if ( ac == 10 ){
			run_id = av[i++];
		}
		else{
			run_id = itoa( random_seed, 10 );
		}

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p );

int getStructureID( Folder &b, const Gene &g );
bool analyzeReplica( const Parameters &p, ostream &s );
void evolutionExperiment( const Parameters &p );

ostream & operator<<( ostream &s, const Parameters &p ) {
	
	s << "# Parameters:" << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   per-site mutation rate u: " << p.u << endl;
	s << "#   scaling factor s: " << p.dist_scale_factor << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   window size tau: " << p.window_size << endl;
	s << "#   equilibration time: " << p.equilibration_time << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   run ID: " << p.run_id << endl;
	s << "#" << endl;
	return s;
	
}


StructureID getStructureID ( Folder &b, const Gene &g ) {
	
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( b.fold(p) );
		return fi->getStructure();
	}
	else
		return (StructureID)-1;
		
}


void runAndAnalyzeReplica ( const Parameters &p, ostream &s ) {
	
	// initialize the population
	Population pop( p.N );
	
	// Make a random target sequence
	Protein targ = Protein::createRandom( p.protein_length );
	s << "# Target sequence: " << targ << endl;
	
	// Get a fitness evaluator with the random target
	AASequenceFitness* fe = new AASequenceFitness( targ, p.dist_scale_factor );

	// Find a sequence.
	Gene g = GeneUtil::reverseTranslate( targ );
	s << "# Starting genotype: " << g << endl;
	
	// Fill the population with the genotype that we found above
	pop.init( g, fe, p.u );

	int loop_length = 100;
	int n_folded = 0;
	for ( int i = 0; ; i++ ) 	{
		for ( int j = 0; j < loop_length; j++ ) {
			pop.evolve();
		}
		pop.prepareCoalescenceCalcs();
		if ( pop.calcCoalescenceTime() > p.equilibration_time + p.window_size )
			break;
		cout << "t=" << i*loop_length+1 << "; " << flush;
	}
	cout << endl;

	pop.printGenebank( s );
	
	delete fe;
	
}


void evolutionExperiment( const Parameters &p )
{
	// set random seed
	Random::seed( p.random_seed );
	
	int count = 0;
	for ( int i = 0; i < p.repetitions; i++ )
	{
		string filename;
		stringstream repfname ;
		repfname << "aaseq_run-gb-rep" << i << "-id" << p.run_id << ".dat";
		filename = repfname.str();
		ofstream gen_file( filename.c_str(), ios::out );
		gen_file << p;

		runAndAnalyzeReplica( p, gen_file );
		
		cout << "[" << i+1 << "/" << p.repetitions << "] " << endl;
	}
	cout << endl;
	
}

int main( int ac, char **av)
{
  Parameters p( ac, av );
  
  if ( ! p.valid ) {
    exit(1);
  }
  
  // seed the random number generator
  Random::seed( p.random_seed );
  
  cout << p;
  
  cout << setprecision(4);
  evolutionExperiment( p );
  
  return 0;

}



