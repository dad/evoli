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


/** \page neut-driver neut-driver
The program \c neut-driver is used to run simulations of neutrally evolving sequences.

The program is called with the following command-line parameters:
\verbatim
   ./neut-driver 25 1000 0.00001 100000 50000 50 111
             0   1  2    5       6      7     8  9 
\endverbatim
(The numbers under the call count the parameters.)

The parameters are, in order (the zeroth parameter is the program name):
-#  %Protein length in amino acids
-#  %Population size
-#  Mutation rate per nucleotide
-#  Window time -- number of generations for which evolutionary data is being collected
-#  Equilibration time -- evolution time in generations before collection of evolutionary data begins
-#  Number of replicates
-#  %Random number seed -- best to use an odd number
-#  Run ID (optional)

*/

#include "compact-lattice-folder.hh"
#include "neut-driver.hh"

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   per-site mutation rate u: " << p.u << endl;
	s << "#   population size N: " << p.N << endl;
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


void runAndAnalyzeReplica( NeutralFitness *fe, const Parameters &p, ostream &s )
{
	// initialize the population
	Population pop( p.N );

	// Use any old random sequence.
	Gene g = Gene::createRandom( p.protein_length * 3 );
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
		
}

void evolutionExperiment( const Parameters &p, NeutralFitness& fe )
{
	// set random seed
	Random::seed( p.random_seed );
	
	int count = 0;
	for ( int i=0; i < p.repetitions; i++ )
	{
		string filename;
		stringstream repfname ;
		repfname << "neutral-run-gb-rep" << i << "-id" << p.run_id << ".dat";
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

  cout << p;
  
  // fitness evaluator for neutral evolution
  NeutralFitness* fe = new NeutralFitness( );
  
  cout << setprecision(4);
  evolutionExperiment( p, *fe );
  
  //evolutionTest( p, *fe );
  delete fe;
  
  return 0;

}




