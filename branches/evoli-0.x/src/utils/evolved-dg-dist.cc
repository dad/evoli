#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

const int size = 5;

class Parameters {
public:
	string eval_type;
	int protein_length;
	int N;
	double ca_cost;
	double error_rate;
	double accuracy_weight;
	double error_weight;
	mutable int structure_ID;
	double free_energy_cutoff;
	double free_energy_minimum;
	int random_seed;
	string genotype_file_name;
	int num_to_fold;
	bool valid;


	Parameters( int ac, char **av ) {
		if ( ac != 11 )	{
			valid = false;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		protein_length = atoi( av[i++] );
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		//free_energy_cutoff = atof( av[i++] );
		//free_energy_minimum = atof( av[i++] );
		random_seed = atoi( av[i++] );
		genotype_file_name = av[i++];
		num_to_fold = atoi( av[i++] );

		valid = true;
	}
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   protein_length: " << p.protein_length << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	//s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	//s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   genotype file: " << p.genotype_file_name << endl;
	s << "#" << endl;
	return s;
}

struct RunRecord
{
	string orf;
	string run_id;
	double cost;
	int rep;
	Gene gene;
};

ostream& operator<<(ostream& os, const RunRecord& rec) {
	os << rec.orf << tab << rec.cost << tab << rec.run_id << tab << "rep " << rec.rep << tab << rec.gene << endl;
	return os;
}

void evolStabExperiment(Parameters& p) {
	// Goal: read in each evolved genotype and compute its dG distribution

	// size of the lattice proteins is hardcoded
	const int size = 5;

	// initialize the protein folder
	CompactLatticeFolder folder(size);

	Random::seed(p.random_seed);


	double free_energy_cutoff = -5;
	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation( &folder, p.protein_length, p.structure_ID, free_energy_cutoff, 0.0, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation( &folder, p.protein_length, p.structure_ID, free_energy_cutoff, 0.0, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation( &folder, p.protein_length, p.structure_ID, free_energy_cutoff, 0.0, p.ca_cost, p.error_rate );
		fe = rob;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(1);
	}

	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(p.genotype_file_name.c_str());
	char buf[1000];
	// Skip two lines
	fin.getline(buf,100);
	fin.getline(buf,100);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.orf >> rec.run_id >> rec.cost >> rec.rep >> rec.gene;
		if (rec.gene.length() > 0) {
			runResults.push_back(rec);
			//cout << rec;
		}
	}
	fin.close();

	cout << "orf\trun_id\ttr\trep\ttrial\tddg\tdg" << endl;
	for (vector<RunRecord>::iterator it = runResults.begin(); it != runResults.end(); it++)	{
		RunRecord& rec = *it;
		fe->setMisfoldingCost(rec.cost);
		Protein prot = rec.gene.translate();
		auto_ptr<FoldInfo> fi(folder.fold(prot));
		double dG = fi->getFreeEnergy();
		vector<double> stabilities;
		//cout << prot << tab << dG << endl;
		stabilities.reserve(p.num_to_fold);
		fe->stabilityOutcomes(rec.gene, p.num_to_fold, stabilities);
		for (int i=0; i<p.num_to_fold; i++) {
			cout << rec.orf << tab << rec.run_id << tab << rec.cost << tab << rec.rep << tab << i << tab << (stabilities[i]-dG) << tab << stabilities[i] << endl;
		}
	}
}


int main( int ac, char **av) {
	// evolved-ddg-dist tr 6 0.0114 59.0 104.5 599 0 ..\analysis\trs599ca6-genes.txt > ..\analysis\trs599ca6-ddg.txt
	Parameters p( ac, av );
	if (p.valid) {
		evolStabExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <eval type> <prot length> <ca cost> <error rate> <accuracy weight> <error weight> <structure ID> <random seed> <gene file name> <num. to fold>" << endl;
	}
}




