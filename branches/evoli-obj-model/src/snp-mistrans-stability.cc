#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

const int size = 5;

class Parameters {
public:
	string eval_type;
	int num_to_fold;
	double free_energy_cutoff;
	double free_energy_minimum;
	double ca_cost;
	double error_rate;
	double error_weight;
	double accuracy_weight;
	int random_seed;
	mutable int structure_ID;
	string genotype_file;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 12 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <num to fold> <eval type> <ca cost> <error rate> <accuracy weight> <error weight> <structure id> <free energy cutoff> <free energy minimum> <random seed> <genotype file>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		num_to_fold = atoi( av[i++] );
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		random_seed = atoi( av[i++] );
		genotype_file = av[i++];

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   num. to fold: " << p.num_to_fold << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   genotype file: " << p.genotype_file << endl;
	s << "#" << endl;
	return s;
}

struct RunRecord
{
	double cost;
	string cost_id;
	int runNumber;
	double fitness;
	double s;
	int nonsynonymous;
	Genotype gene;
};

void evolStabExperiment(Parameters& p) {
	// Goal: read in each evolved genotype and compute its dG distribution

	// size of the lattice proteins is hardcoded
	const int size = 5;

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	srand48(p.random_seed);


	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation();
		ept->init( &folder, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation();
		afe->init( &folder, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation();
		rob->init( &folder, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate );
		fe = rob;
	}
	else if (p.eval_type == "con") {
		StabilityConstraint* sc = new StabilityConstraint();
		sc->init( &folder, p.structure_ID, p.free_energy_cutoff, p.free_energy_minimum, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = sc;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(0);
	}


	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(p.genotype_file.c_str());
	char buf[1000];
	fin.getline(buf,100);
	fin.getline(buf,100);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.cost_id >> rec.runNumber >> rec.fitness >> rec.s >> rec.nonsynonymous >> rec.gene;
		if (rec.gene.size() > 0) {
			rec.cost = pow(10.0,atof(rec.cost_id.c_str()));
			runResults.push_back(rec);
			//cout << rec.cost_id << tab << rec.runNumber << endl;
		}
	}
	fin.close();

	cout << "tr\trun\tns\ttrial\texpr\ts\tddg\tdg" << endl;
	for (vector<RunRecord>::iterator it = runResults.begin(); it != runResults.end(); it++)	{
		RunRecord& rec = *it;
		if (rec.cost >= 0.01 && rec.cost <= 100.0) {
			pair<double, int> fdata = GenotypeUtil::translateAndFold( folder, rec.gene );
			vector<double> stabilities;
			stabilities.reserve(p.num_to_fold);
			fe->stabilityOutcomes(rec.gene, p.num_to_fold, stabilities);
			for (int i=0; i<p.num_to_fold; i++) {
				cout << rec.cost_id << tab << rec.runNumber << tab << rec.nonsynonymous << tab << i << tab << rec.cost << tab 
					 << rec.s << tab << (stabilities[i]-fdata.first) << tab << stabilities[i] << endl;
			}
		}
	}
}


int main( int ac, char **av) {
	Parameters p( ac, av );
	if (p.valid) {
		evolStabExperiment(p);
	}
}




