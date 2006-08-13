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
#include <list>

const int size = 5;


class Parameters {
public:
	string eval_type;
	double free_energy_cutoff;
	double free_energy_minimum;
	double ca_cost;
	double error_rate;
	double error_weight;
	double accuracy_weight;
	int random_seed;
	double min_frac_folded;
	string file_name;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 11 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval type> <ca cost> <error rate> <accuracy weight> <error weight> <free energy cutoff> <free energy minimum> <random seed> <min. ffold> <gene file>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		random_seed = atoi( av[i++] );
		min_frac_folded = atof( av[i++] );
		file_name = av[i++];

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   min. fraction folded: " << p.min_frac_folded << endl;
	s << "#   gene file name: " << p.file_name << endl;
	s << "#" << endl;
	return s;
}


struct Move {
	int to;
	int from;
	int pos;

	Move(int to, int from, int pos) : to(to), from(from), pos(pos) {
	}
};

struct Outcomes {
	double ffold;
	double facc;
	double frob;
	double ftrunc;

	void setZero() {
		ffold = 0.0;
		facc = 0.0;
		frob = 0.0;
		ftrunc=0.0;
	}
};

Move getMove(const Genotype& g) {
	double rand = myRand();
	int randpos = (int)(rand * g.size()+1)-1;
	int from_codon = g[randpos];
	int to_codon = from_codon;
	do {
		rand = myRand();
		to_codon = (int)(64*rand+1)-1;
	} while (to_codon == from_codon || GeneticCodeUtil::geneticCode[to_codon]<0);
	return Move(to_codon, from_codon, randpos);
}

struct RunRecord
{
	string tr;
	double cost;
	int runNumber;
	Genotype gene;
};

ostream& operator<<(ostream& os, const RunRecord& rec) {
	os << rec.cost << tab << rec.runNumber << tab << rec.gene << endl;
	return os;
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

Genotype optimizeFracFolded(const Genotype& gorig) {
	Genotype g = gorig;
	return g;
}

bool acceptMove(const Move& move, ErrorproneTranslation* fe, Genotype& g, Outcomes& outcomes, double min_frac_folded) {
	bool accept = false;
	g[move.pos] = move.to;
	accept = fe->getFolded(g);
	if (!accept) {
		outcomes.setZero();
	}
	else {
		fe->calcOutcomes( g, outcomes.facc, outcomes.frob, outcomes.ftrunc, outcomes.ffold);
		accept = (outcomes.ffold >= min_frac_folded);
	}
	g[move.pos] = move.from;
	return accept;
}


int main( int ac, char **av) {
	Parameters p( ac, av );

	if (!p.valid) {
		return 1;
	}

	srand48(0);

	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(p.file_name.c_str());
	char buf[500];
	// Skip two lines
	fin.getline(buf,500);
	fin.getline(buf,500);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.tr >> rec.runNumber >> rec.gene;
		rec.cost = pow(10.0, atof(rec.tr.c_str()));
		if (rec.gene.size() > 1 && rec.cost >= 100) {
			runResults.push_back(rec);
			//cout << rec;
		}
	}
	fin.close();

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	vector<RunRecord>::iterator it = runResults.begin();
	int structure_ID = -1;
	structure_ID = getStructureID( folder, runResults[0].gene );
	if ( structure_ID < 0 )	{
		cerr << "Input sequence does not translate!" << endl;
	}

	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation();
		ept->init( &folder, structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation();
		afe->init( &folder, structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation();
		rob->init( &folder, structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate );
		fe = rob;
	}
	else if (p.eval_type == "con") {
		StabilityConstraint* sc = new StabilityConstraint();
		sc->init( &folder, structure_ID, p.free_energy_cutoff, p.free_energy_minimum, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = sc;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(0);
	}

	vector<bool> isOptimal = fe->getOptimalCodons(true);

	cout << "tr\trun\tstep\tfacc\tfrob\tftrunc\tffold\tfop\tdg" << endl;
	for (; it != runResults.end(); it++) {
		RunRecord& rec = *it;
		//if (rec.cost < 100) {
		//	continue;
		//}
		Genotype& g = rec.gene;

		Outcomes outcomes;
		fe->calcOutcomes( g, outcomes.facc, outcomes.frob, outcomes.ftrunc, outcomes.ffold);
		if (outcomes.ffold < p.min_frac_folded) {
			continue;
		}

		double fop = 0.0, dg = 0.0;

		int step_count = 0;
		int step_max = 1000;
		Move move(-1, -1, -1);


		while (step_count < step_max ) {
			move = getMove(g);
			//cout << move.pos << tab << move.to << endl;
			if (acceptMove(move, fe, g, outcomes, p.min_frac_folded)) {
				// Do move
				g[move.pos] = move.to;
				fop = GenotypeUtil::calcFop(g, isOptimal);
				dg = fe->getLastFreeEnergy();
				step_count++;
				cout << rec.tr << tab << rec.runNumber << tab << step_count << tab << outcomes.facc << tab << outcomes.frob << tab << outcomes.ftrunc << tab << outcomes.ffold
					 << tab << fop << tab << dg << endl;
			}
		}
	}
	delete fe;
	return 0;
}




