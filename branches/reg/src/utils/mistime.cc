#include "folder.hh"
#include "decoy-contact-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "protein.hh"
#include "codon.hh"
#include "upstream.hh"
#include "tools.hh"
#include "gene-util.hh"
#include "expressible-gene.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

typedef Population<ExpressibleGene, ExpressibleGeneFitness, ExpressibleGenePolymerase> ExprPop;

class Parameters {
public:
	unsigned int protein_length;
	unsigned int N;
	double mutation_rate;
	mutable int structure_ID;
	double free_energy_cutoff;
	int random_seed;
	string genotype_file_name;
	string binding_energy_file_name;
	unsigned int reps;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 10 )	{
			valid = false;
			return;
		}

		int i = 1;
		protein_length = atoi( av[i++] );
		N = atoi( av[i++] );
		mutation_rate = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		random_seed = atoi( av[i++] );
		genotype_file_name = av[i++];
		binding_energy_file_name = av[i++];
		reps = atoi( av[i++] );

		valid = true;
	}
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   protein_length: " << p.protein_length << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   mutation rate: " << p.mutation_rate << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   genotype file: " << p.genotype_file_name << endl;
	s << "#   binding energy file: " << p.binding_energy_file_name << endl;
	s << "#   reps: " << p.reps << endl;
	s << "#" << endl;
	return s;
}

ExpressibleGene readGene(ifstream& s, const Parameters& params) {
	string gene_str, promoter_str;
	s >> promoter_str;
	s >> gene_str;
	CodingDNA dna(gene_str);
	Promoter prom(promoter_str);
	//cout << prom << endl;
	//cout << dna << endl;
	return ExpressibleGene(prom, dna);
}

void printPopulation(ExprPop* pop) {
	ExprPop::iterator it = pop->begin();
	int i=0;
	for (; it != pop->end(); it++) {
		ExpressibleGene g = (*it)->getOrganism();
		CodingDNA dna = g.getCodingDNA();
		ExpressibleGeneFitness* fe = pop->getFitnessEvaluator();
		double fitness = fe->getFitness(dna);
		cout << i << " " << fitness << endl;
		i++;
	}
}

/**
 * Return true if any sequence in the population encodes a folded protein.
 **/
bool anySequenceFolds(ExprPop* pop, ProteinFolder* folder, StructureID structure_ID, double deltag_cutoff) {
	ExprPop::iterator it = pop->begin();
	bool res = false;
	int i=0;
	Translator trans(0);
	for (; it != pop->end() && !res; it++) {
		ExpressibleGene g = (*it)->getOrganism();
		CodingDNA dna = g.getCodingDNA();
		double fitness = pop->getFitnessEvaluator()->getFitness(g);
		if (fitness == 1.0) {
			Protein p(dna.codonLength());
			trans.translateErrorFree(dna.transcribe(), p);
			auto_ptr<FoldInfo> fi(folder->fold(p));
			//cout << i << " " << fitness << " " << fi->getDeltaG() << " " << fi->getStructure() << endl;
			res = res || (fi->getDeltaG() <= deltag_cutoff && fi->getStructure() == structure_ID);
			i++;
		}
	}
	return res;
}

void mistimeExperiment(Parameters& params) {
	BindingInteraction bi(params.binding_energy_file_name);
	//cout << setw(2) << setprecision(2);
	//cout << bi << endl;
	TranscriptionFactor tf(Protein::createRandom(100));
	
	ifstream gene_fin(params.genotype_file_name.c_str());
	ExpressibleGene seed_gene = readGene( gene_fin, params );
	cout << "# " << seed_gene << endl;
	CodingDNA gene = seed_gene.getCodingDNA();
	gene_fin.close();

	string fname = "../../test/data/p450-structs/maps.txt";
	ifstream fin(fname.c_str());
	string dir = "../../test/data/p450-structs/";
	double log_nconf = 10 * log(10);
	DecoyContactFolder folder( params.protein_length, log_nconf, fin, dir);

	if (!folder.good()) {
		cout << "# Couldn't initialize folder with " << fname << endl;
		return;
	}
	cout <<"# Initialized folder" << endl;
	//ProteinStructureFitness fe( &folder, params.structure_ID, params.free_energy_cutoff );
	ProteinFreeEnergyFitness fe( &folder ); //, params.structure_ID, params.free_energy_cutoff );
	cout << "# Initialized PSF" << endl;
	ExpressibleGeneFitness egf(&fe, &tf, &bi);
	cout <<"# Initialized FitnessEvaluator" << endl;

	cout << "# Avg. mutations per population per generation N*U*L = " << (params.N*params.mutation_rate*gene.length()) << endl;
	
	Protein p = gene.translate();
	auto_ptr<DecoyFoldInfo> fi(folder.fold(p));
	cout << "# seed prot " << p << endl;
	cout << "# dG = " << fi->getDeltaG() << " " << fi->getStructure() << endl;
 
	if (fi->getStructure() != params.structure_ID) {
		cerr << "Seed gene folded into structure " << fi->getStructure() << " not " << params.structure_ID << endl;
	}

	ostream& out = cout; // write to file
	out << "rep\ttmis" << endl;
	for (unsigned int rep=0; rep < params.reps; rep++ ) {
		// Now seed a population and let it evolve until all genes encoding a folded protein are lost from the population.
		ExprPop pop(params.N);
		ExpressibleGenePolymerase poly(params.mutation_rate);
		pop.init(seed_gene, &egf, &poly);
		bool seq_folds = fi->getDeltaG() <= params.free_energy_cutoff;
		unsigned int num_generations = 0;
		unsigned int num_gens_per_check = 10;
		while (seq_folds) {
			for (int gen=0; gen<num_gens_per_check; gen++) {
				pop.evolve();
				num_generations++;
			}
			seq_folds = anySequenceFolds(&pop, &folder, params.structure_ID, params.free_energy_cutoff);
			//cout << "gen: " << num_generations << endl;
		}
		out << rep << tab << num_generations << endl;
	}
}

int main(int ac, char** av) {
	Parameters p( ac, av );
	if (p.valid) {
		Random::seed(p.random_seed);
		cout << p;
		mistimeExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <prot length> <pop size> <mutation rate> <structure id> <free energy cutoff> "
			 << "<random seed> <gene file name> <binding energy file name> <reps>" << endl;
	}
}
