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
	string output_file_name;
	bool valid;


	Parameters( int ac, char **av ) {
		if ( ac != 12 )	{
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
		output_file_name = av[i++];

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
	s << "#   output file: " << p.output_file_name << endl;
	s << "#" << endl;
	return s;
}

struct RunRecord
{
	//orf	tr	rep	fitness	s	ns	gene
	string orf;
	double cost;
	int rep;
	double fitness; // fitness
	double s;       // fitness advantage
	bool ns;        // is SNP nonsynonymous?
	Gene gene;
};

ostream& operator<<(ostream& os, const RunRecord& rec) {
	os << rec.orf << tab << rec.cost << tab << "rep " << rec.rep << tab << rec.ns << tab << rec.s << tab << rec.fitness << tab << rec.gene << endl;
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
	fin.getline(buf,100);
	fin.getline(buf,100);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.orf >> rec.cost >> rec.rep >> rec.fitness >> rec.s >> rec.ns >> rec.gene;
		if (rec.gene.codonLength() == p.protein_length) {
			runResults.push_back(rec);
			//cout << rec.cost_id << tab << rec.runNumber << endl;
		}
	}
	fin.close();

	ofstream fout(p.output_file_name.c_str());

	fout << "orf\ttr\trep\tns\ttrial\ts\tddg\tdg" << endl;
	for (vector<RunRecord>::iterator it = runResults.begin(); it != runResults.end(); it++)	{
		RunRecord& rec = *it;
		if (rec.cost <0.4) {
			continue;
		}
		fe->setMisfoldingCost(rec.cost);
		Protein prot = rec.gene.translate();
		auto_ptr<FoldInfo> fi(folder.fold(prot));
		double dG = fi->getDeltaG();
		vector<double> stabilities;
		stabilities.reserve(p.num_to_fold);
		fe->stabilityOutcomes(rec.gene, p.num_to_fold, stabilities);
		for (int i=0; i<p.num_to_fold; i++) {
			fout << rec.orf << tab << rec.cost << tab << rec.rep << tab << rec.ns << tab << i << tab
				 << rec.s << tab << (stabilities[i]-dG) << tab << stabilities[i] << endl;
		}
	}
	fout.close();
}


int main( int ac, char **av) {
	Parameters p( ac, av );
	if (p.valid) {
		evolStabExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <eval type> <prot length> <ca cost> <error rate> <accuracy weight> <error weight> <structure ID> <random seed> <SNP file name> <num. to fold>" << endl;
	}
}




