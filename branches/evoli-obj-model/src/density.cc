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
	mutable int structure_ID;
	double bias;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 11 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval type> <ca cost> <error rate> <accuracy weight> <error weight> <structure id> <free energy cutoff> <free energy minimum> <random seed> <bias>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		random_seed = atoi( av[i++] );
		bias = atof( av[i++] );

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   bias: " << p.bias << endl;
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

double calcAcceptanceProbability(const Move& move, Genotype& g, ErrorproneTranslation* fe, const double bias, Outcomes& outcomes) {
	double prob = 1.0;
	g[move.pos] = move.to;
	bool folded = fe->getFolded(g);
	if (!folded) {
		prob = 0.0;
		outcomes.setZero();
	}
	else {
		fe->calcOutcomes( g, outcomes.facc, outcomes.frob, outcomes.ftrunc, outcomes.ffold);
		prob = exp(-bias*(1.0-outcomes.ffold));
	}
	g[move.pos] = move.from;
	if (prob > 0.0) {
		prob = max(prob,1e-4);
	}
	return prob;
}

int main( int ac, char **av) {
	Parameters p( ac, av );
	if (!p.valid) {
		exit(1);
	}
	cout << p;

	srand48(p.random_seed);

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

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

	Genotype g = GenotypeUtil::getSequenceForStructure(folder, p.free_energy_cutoff, p.structure_ID);

	vector<bool> isOptimal = fe->getOptimalCodons(true);

	Outcomes outcomes;
	double fop = 0.0, dg = 0.0;

	int step_count = 0;
	int step_max = 1000000;
	double last_ffold = 0.0;
	Move move(-1, -1, -1);

	cout << "pos\ttocodon\tfacc\tfrob\tftrunc\tffold\tfop\tdg" << endl;
	while (step_count < step_max ) {
		move = getMove(g);
		double prob = calcAcceptanceProbability(move, g, fe, p.bias, outcomes);
		//if (myRand() < prob) {
		if (outcomes.ffold>last_ffold || (prob>0.0 && myRand()<1/p.bias)) {
			// Do move
			g[move.pos] = move.to;
			fop = GenotypeUtil::calcFop(g, isOptimal);
			dg = fe->getLastFreeEnergy();
			step_count++;
			cout << move.pos << tab << move.to << tab << outcomes.facc << tab << outcomes.frob << tab << outcomes.ftrunc << tab << outcomes.ffold << tab << fop << tab << dg << endl;
			last_ffold = outcomes.ffold;
		}
	}
	delete fe;
}




