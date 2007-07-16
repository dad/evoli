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


/** \page get-weights get-weights
The program \c get-weights is needed to determine the correct weights for the program \c tr-driver. The \ref tr-driver needs parameters error rate and weights. We want to choose these parameters such that, given a particular structure and codon adaptation cost (which affects translational accuracy), a random gene encoding a protein that folds to that structure is accurately translated a particular fraction of the time (say, 85% of the time).  These parameters can be obtained by running the get-weights program as follows:
\verbatim
   ./get-weights 3 -5 111 stable-sequence.txt 10000 10 0.85 5
               0 1  2   3                   4     5  6    7 8
\endverbatim
-#  Codon adaptation cost -- average accuracy ratio between optimal and non-optimal codons
-#  Maximum free energy of folding
-#  %Random number seed -- best to use an odd number
-#  Seed gene sequence file containing nucleotide sequence that folds into the target structure with less than maximum free energy
-#  Equilibration time -- number of generations to equilibrate seed sequence by mutational drift
-#  Window time -- number of generations for which to measure random sequence statistics
-#  Target translational accuracy -- fraction of time a random gene is accurately translated
-#  Number of replicates -- averaged data will also be provided

If the sequence file name does not represent a valid file, the program
will attempt to interpret the name as a sequence identifier
(SequenceID, e.g. an integer between 0 and 1080 in the case of the 5x5
compact lattice model).
*/

#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "folder-util.hh"
#include "random.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <ctime>
#include <iomanip>

class Stats
{
public:
	double sum;
	double sumSquared;
	uint samples;
 
  Stats() {
	reset();
  }

	void addValue(double val)
	{
		sum += val;
		sumSquared += val*val;
		samples++;
	}

  void operator+=(double val) {
	addValue(val);
  }

	void reset(void)
	{
		sum = 0.0;
		sumSquared = 0.0;
		samples = 0;
	}

	double getMean(void) { return sum/samples; }
	double getVariance(void) { return sumSquared/samples - (sum/samples)*(sum/samples); }
	double getSampleVariance(void) { return sumSquared/(samples-1) - (sum/samples)*(sum/samples)*(double)samples/(samples-1); }
	double getStandardDeviation(void) { return sqrt(getVariance()); }
	double getStandardError(void) { return sqrt(getVariance()/samples); }
	double getZscore(double val) { return (val - getMean())/getStandardDeviation(); }
};

class Parameters {
public:
	double ca_cost;
	mutable int structure_ID;
	double free_energy_cutoff;
	int random_seed;
	string genotype_file_name;
	int num_to_fold;
	int num_equil;
	double target_accuracy;
	int reps;
	bool valid;


	Parameters( int ac, char **av ) {
		if ( ac != 9 )	{
			valid = false;
			return;
		}

		int i = 1;
		ca_cost = atof( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		random_seed = atoi( av[i++] );
		genotype_file_name = av[i++];
		num_to_fold = atoi( av[i++] );
		num_equil = atoi( av[i++] );
		target_accuracy = atof( av[i++] );
		reps = atoi( av[i++] );
		valid = true;
	}
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   genotype file | structure ID: " << p.genotype_file_name << endl;
	s << "#   num. gens to equilibrate: " << p.num_equil << endl;
	s << "#   num. gens to measure: " << p.num_to_fold << endl;
	s << "#   target accuracy: " << p.target_accuracy << endl;
	s << "#   reps: " << p.reps << endl;
	s << "#" << endl;
	return s;
}

StructureID getStructureID( Folder *b, const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( b->fold(p) );
		cout << fi->getStructure() << " " << fi->getDeltaG() << flush << endl;
		return fi->getStructure();
	}
	else
		return (StructureID)-1;
}

void getWeightsExperiment(Parameters& p)
{
	// Seed random number generator
	long seconds = (long)time(NULL);
	Random::seed(p.random_seed);

	// Get protein and extract structure ID
	Folder* folder;
	Gene seed_gene;
	StructureID structure_ID;
	ifstream fin(p.genotype_file_name.c_str());
	if (fin.good()) {
		fin >> seed_gene;
		int side_length = (int)(sqrt(float(seed_gene.codonLength())));
		// initialize the protein folder
		folder = new CompactLatticeFolder(side_length);
		structure_ID = getStructureID(folder, seed_gene);
	}
	else {
		// If can't find the file, interpret filename as structure ID
		cout << "# Invalid filename; interpreting " << p.genotype_file_name << " as structure ID and assuming length 25." << endl;
		StructureID struct_id = (StructureID)(atoi(p.genotype_file_name.c_str()));
		folder = new CompactLatticeFolder(5);
		seed_gene = FolderUtil::getSequenceForStructure(*folder, 75, p.free_energy_cutoff, struct_id);
		cout << "# Seed gene " << seed_gene << endl;
		cout << "# Translates to " << seed_gene.translate() << endl;
		auto_ptr<FoldInfo> fi( folder->fold(seed_gene.translate()) );
		structure_ID = fi->getStructure();
	}
	fin.close();

	if ( structure_ID < 0 )	{
		cerr << "# Input sequence does not translate!" << endl;
		exit(1);
	}
	else {
		cout << "# Structure ID " << structure_ID << endl;
	}

	ErrorproneTranslation* ept = new ErrorproneTranslation( folder, seed_gene.codonLength(), structure_ID, p.free_energy_cutoff, 1, p.ca_cost, 0.1, 0.1, 0.1 );
	double m_error_rate;
	double m_accuracy_weight;
	double m_error_weight;
	Stats er, aw, ew;
	cout << "error rate\taccuracy weight\terror weight" << endl;
	for (int ni=0; ni<p.reps; ni++) {
	  ept->getWeightsForTargetAccuracy(seed_gene, p.target_accuracy, m_error_rate, m_accuracy_weight, m_error_weight, p.num_equil, p.num_to_fold);	
	  er += m_error_rate;
	  ew += m_error_weight;
	  aw += m_accuracy_weight;
	  cout << m_error_rate << tab << m_accuracy_weight << tab << m_error_weight << endl;
	}
	cout << "# Averages:" << endl;
	cout << er.getMean() << tab << aw.getMean() << tab << ew.getMean() << endl;

	// Clean up.
	delete folder;
}

int main( int ac, char **av)
{
	Parameters p( ac, av );
	if (p.valid) {
	  cout << p;
		getWeightsExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <ca cost> <free energy cutoff> <random seed> <gene file name | structure ID> <num. to equil> <num. to measure> <target accuracy> <reps>" << endl;
	}
}
