#include "genotype.hh"
#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "translator.hh"
#include "codon.hh"
#include "genotype-util.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>

const int size = 5;
const char* tab = "\t";

class Parameters {
public:
	string eval_type;
	int aa_diffs;
	double free_energy_cutoff;
	double free_energy_minimum;
	double ca_cost;
	double base_fraction_accurate;
	int num_to_fold;
	string file_name;
	Genotype g;
	const static string usage;
	bool valid;

	Parameters(int ac, char**av) {
		valid = false;
		if ( ac != 9 )	{
			cout << "Start program like this:" << endl;
			cout << av[0] << Parameters::usage << endl;
		}
		else {
			int i = 1;
			eval_type = av[i++];
			ca_cost = atof( av[i++] );
			aa_diffs = atoi( av[i++] );
			base_fraction_accurate = atof( av[i++] );
			free_energy_cutoff = atof( av[i++] );
			free_energy_minimum = atof( av[i++] );
			num_to_fold = atoi( av[i++] );
			file_name = av[i];

			// Read genotype file
			ifstream fin;
			fin.open(file_name.c_str(), ifstream::in);
			fin >> g;
			fin.close();
			valid = true;
		}
	}
};
const string Parameters::usage = " <eval type> <ca cost> <aa diffs> <base accuracy> <free energy cutoff> <free energy minimum> <num to fold> <genotype file>";

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   fitness evaluation type: " << p.eval_type << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   use aa cost differences: " << p.aa_diffs << endl;
	s << "#   baseline frac. accurate: " << p.base_fraction_accurate << endl;
	s << "#   free energy cutoff: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   num. proteins to attempt to fold: " << p.num_to_fold << endl;
	s << "#   genotype file: " << p.file_name << endl;
	s << "#" << endl;
	return s;
}

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


int main( int ac, char **av) {
	Parameters p( ac, av );

	if (!p.valid) {
		return 1;
	}

	srand48(0);

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	// find initial structure
	int structure_ID = getStructureID( folder, p.g );
	if ( structure_ID < 0 )	{
		cerr << "Input sequence does not translate!" << endl;
		return 1;
	}

	double transl_error_rate = 1 - pow(p.base_fraction_accurate, (double)p.g.size());

	double cost = 1.0;
	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation();
		afe->init( &folder, p.g, p.free_energy_cutoff, cost, p.ca_cost, transl_error_rate );
		fe = afe;
	}
	else if (p.eval_type == "con") {
		StabilityConstraint* sc = new StabilityConstraint();
		sc->init( &folder, structure_ID, p.free_energy_cutoff, p.free_energy_minimum, cost, p.ca_cost, transl_error_rate );
		fe = sc;
	}
	else if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation();
		ept->init( &folder, structure_ID, p.free_energy_cutoff, cost, p.ca_cost, transl_error_rate );
		fe = ept;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation();
		rob->init( &folder, structure_ID, p.free_energy_cutoff, cost, p.ca_cost, p.base_fraction_accurate );
		fe = rob;
	}
	if (p.aa_diffs > 0) {
		fe->setCodonCosts();
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(0);
	}

	fe->setTargetAccuracyOfRandomGenes(p.g, p.base_fraction_accurate, 5000, 5000);
	double error_rate = fe->getErrorRate();
	cout << p;
	cout << "# Error rate = " << error_rate << endl;
	vector<bool> isOptimal = fe->getOptimalCodons(true);

	cout << "acc\trob\ttrunc\tfold\tfacc\tfrob\tftrunc\tffold\tcfacc\tcfrob\tcftrunc\tcffold\tncalcfold" << endl;

	Genotype& g = p.g;

	unsigned int num_codons = g.size();
	double cfacc, cfrob, cftrunc, cffold;

	int numRobust = 0;
	int numTruncated = 0;
	int numFolded = 0;
	int numAccurate = 0;

	int step_count = 0;
	int step_max = 1000000;

	double dg = 0;

	int to_codon = -1, randpos = -1, from_codon = -1;

	while (step_count < step_max ) {
		if (fe->getFolded(g)) {
			dg = fe->getLastFreeEnergy();
			fe->countOutcomes(g, p.num_to_fold, numAccurate, numRobust, numTruncated, numFolded);
			int num_folded_in_calculation = folder.getNumFolded();
			fe->calcOutcomes( g, cfacc, cfrob, cftrunc, cffold);
			num_folded_in_calculation = folder.getNumFolded() - num_folded_in_calculation;

			cout << numAccurate << tab << numRobust << tab << numTruncated << tab << numFolded << tab << setprecision(4)
				 << cfacc << tab << cfrob << tab << cftrunc << tab << cffold << tab
				 << numAccurate/(double)p.num_to_fold << tab << numRobust/(double)(p.num_to_fold-numAccurate) << tab
				 << numTruncated/(double)p.num_to_fold << tab << numFolded/(double)p.num_to_fold << tab << num_folded_in_calculation << endl;
			step_count++;
		}
		else {
			if (randpos >= 0) {
				g[randpos] = from_codon;
			}
		}
		double rand = myRand();
		randpos = (int)(rand * num_codons+1)-1;
		// go through all possible point mutations
		from_codon = g[randpos];
		rand = myRand();
		to_codon = from_codon;
		do {
			rand = myRand();
			to_codon = (int)(64*rand+1)-1;
		} while (to_codon == from_codon);

		g[randpos] = to_codon;
	}
	delete fe;
}




