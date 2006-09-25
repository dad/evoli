#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"
#include "random.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <ctime>
#include <iomanip>

typedef unsigned int uint16;

class Stats
{
public:
	double sum;
	double sumSquared;
	uint16 samples;
 
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
	s << "#   genotype file: " << p.genotype_file_name << endl;
	s << "#   num. gens to equilibrate: " << p.num_equil << endl;
	s << "#   num. gens to measure: " << p.num_to_fold << endl;
	s << "#   target accuracy: " << p.target_accuracy << endl;
	s << "#   reps: " << p.reps << endl;
	s << "#" << endl;
	return s;
}

StructureID getStructureID( Folder &b, const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		return b.fold(p).getStructure();
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
	Gene seed_gene;
	ifstream fin(p.genotype_file_name.c_str());
	if (fin.good()) {
	  fin >> seed_gene;
	}
	fin.close();

	int side_length = (int)(sqrt(seed_gene.codonLength()));
	// initialize the protein folder
	CompactLatticeFolder folder(side_length);

	int structure_ID = getStructureID(folder, seed_gene);


	if ( structure_ID < 0 )	{
		cerr << "# Input sequence does not translate!" << endl;
		exit(1);
	}
	ErrorproneTranslation* ept = new ErrorproneTranslation();
	ept->init( &folder, seed_gene.codonLength(), structure_ID, p.free_energy_cutoff, 1, p.ca_cost, 0.1, 0.1, 0.1 );
	double m_error_rate;
	double m_accuracy_weight;
	double m_error_weight;
	Stats er, aw, ew;
	for (int ni=0; ni<p.reps; ni++) {
	  ept->getWeightsForTargetAccuracy(seed_gene, p.target_accuracy, m_error_rate, m_accuracy_weight, m_error_weight, p.num_equil, p.num_to_fold);	
	  er += m_error_rate;
	  ew += m_error_weight;
	  aw += m_accuracy_weight;
	  cout << m_error_rate << tab << m_accuracy_weight << tab << m_error_weight << endl;
	}
	cout << "# Averages:" << endl;
	cout << er.getMean() << tab << aw.getMean() << tab << ew.getMean() << endl;
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
		cout << "\t" << av[0] << " <ca cost> <free energy cutoff> <random seed> <gene file name> <num. to equil> <num. to measure> <target accuracy> <reps>" << endl;
	}
}




