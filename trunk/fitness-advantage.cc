#include "genotype.hh"
#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "translator.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

const int size = 5;
const char* tab = "\t";

class Parameters {
public:
	int N;
	double ca_cost;
	double transl_error_rate;
	double free_energy_cutoff;
	double free_energy_minimum;
	int random_seed;
	int use_aa_diffs; // Use differing accuracy values for amino acids?
	int structure_ID;
	string file_name;

	bool valid;

	Parameters() {
		valid = false;
	}

	Parameters(int ac, char**av, const string& usage) {
		if ( ac != 9 )	{
			cout << usage << endl;
			valid = false;
		}
		else {
			int i = 1;
			N = atoi( av[i++] );
			ca_cost = atof( av[i++] );
			transl_error_rate = atof( av[i++] );
			free_energy_cutoff = atof( av[i++] );
			free_energy_minimum = atof( av[i++] );
			random_seed = atoi( av[i++] );
			use_aa_diffs = atoi( av[i++] );
			file_name = av[i++];
			valid = true;
		}
	}

};

int getStructureID( ProteinFolder &b, const Genotype &g ) {
	int l = b.getProteinLength();
	Translator t( 0, l );
	int *seq = new int[l];

	if ( t.translateErrorFree( g, seq ) ) {
		b.foldProtein( seq );
		return b.getLastFoldedProteinStructureID();
	}
	else
		return -1;
	delete [] seq;
}

/*double fixation_probability(int N, double s) {
	double res = 0.0;
	if (s == 0.0) {
		// completely neutral evolution
		res = 1.0/N;
	}
	else if (s <= -1.0) {
		// mutant has zero fitness.
		res = 0.0;
	}
	else {
		// Sella and Hirsh probability of fixation
		double delta = log(1.0+s);
		double num = 1.0 - exp(-2.0 * delta);
		double denom = 1.0 - exp(-2.0 * N * delta);
		if (denom != 0.0) {
			res = num/denom;
		}
	}
	return res;
}*/

bool isPointMutant( int co, int ct ) {
	if ( co == ct )
		return false;

	int lo1, lo2, lo3;
	int lt1, lt2, lt3;

	CodonUtil::codonToLetters( lo1, lo2, lo3, co );
	CodonUtil::codonToLetters( lt1, lt2, lt3, ct );

	bool res = (( lo1 != lt1 && lo2 == lt2 && lo3 == lt3 ) ||
			( lo2 != lt2 && lo1 == lt1 && lo3 == lt3 ) ||
			( lo3 != lt3 && lo1 == lt1 && lo2 == lt2 ));
	//cout << lo1 << lo2 << lo3 << tab << lt1 << lt2 << lt3 << tab << res << endl;
	return res;
}


/**
 * Calculates the neutrality of the given genotype.
 **/
void calcMutantFitnesses( FitnessEvaluator& fe, const Genotype &gorig, vector<double>& fitnesses ) {
	Genotype g = gorig;
	int num_codons = g.size();
	// go through all positions in the gene
	for ( int i=0; i<num_codons; i++ ) {
		// go through all possible point mutations
		// (avoid operator %, which can be very slow)
		int codon = g[i];
		for (int to_codon=0; to_codon<64; to_codon++) {
			if (isPointMutant(codon, to_codon)) {
				g[i] = to_codon;
				double fitness = fe.getFitness(g);
				fitnesses.push_back(fitness);
			}
		}
		// Reset the genotype
		g[i] = gorig[i];
	}
}



int main( int ac, char **av) {
	string prog_name = av[0];
	stringstream usage;
	usage << "Start program like this:\n  " << prog_name  << " <pop size> <ca cost> <transl. error rate> <free energy cutoff> <free energy minimum> <random seed> <use AA diffs> <file>";
	Parameters p( ac, av, usage.str() );

	if (!p.valid) {
		return 0;
	}

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	bool fe_initialized = false;
	ErrorproneTranslation* pfe = NULL;

	Genotype g;
	double tr = -1.0;
	double oldtr = tr;
	int rep;

	// Read file of genotypes
	ifstream fin;
	fin.open(p.file_name.c_str(), ifstream::in);
	char buf[1000];
	fin.getline(buf, 1000);
	fin.getline(buf, 1000);
	ostream& os = cout;
	os << "tr\trep\ts\tpfix" << endl;
	while (!fin.eof()) {
		oldtr = tr;
		fin >> tr >> rep >> g;
		//if (rep >= 5) {
		//	continue;
		//}
		os << "# " << tr << tab << rep << tab << g << endl;
		if (oldtr != tr) {
			// initialize fitness evaluator
			p.structure_ID = getStructureID( folder, g );
			delete pfe;
			pfe = new ErrorproneTranslation();
			pfe->init( &folder, p.structure_ID, p.free_energy_cutoff, tr, p.ca_cost, p.transl_error_rate );
			if (p.use_aa_diffs == 1) {
				cout << "# setting codon costs" << endl;
				pfe->setCodonCosts();
				pfe->printCodonCosts(cout);
			}
			fe_initialized = true;
		}
		double base_fitness = pfe->getFitness(g);
		// Loop over all nucleotide mutants and output fitness advantage and probability of fixation.
		vector<double> fitnesses;
		calcMutantFitnesses(*pfe, g, fitnesses);
		for (vector<double>::iterator it=fitnesses.begin(); it != fitnesses.end(); it++) {
			double s = *it/base_fitness - 1.0;
			double pfix = fixation_probability(p.N, s);
			os << tr << tab << rep << tab << s << tab << pfix << endl;
		}
	}
	delete pfe;
	fin.close();
}




