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
	int step_max;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 12 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval type> <ca cost> <error rate> <accuracy weight> <error weight> <structure id> <free energy cutoff> <free energy minimum> <random seed> <bias> <max steps>" << endl;
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
		step_max = atoi( av[i++] );

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
	s << "#   max steps: " << p.step_max << endl;
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

class Histogram {
	vector<int> m_bins;
	unsigned int m_nbins;
	double m_binmin;
	double m_binmax;
	int m_count;
public:
	Histogram(unsigned int nbins, double binmin, double binmax) {
		m_bins.resize(nbins,0);
		m_binmin = binmin;
		m_binmax = binmax;
		m_nbins = nbins;
		m_count = 0;
	}

	void reset() {
		m_count = 0;
		m_bins.resize(m_nbins,0);
	}

	void add(double x) {
		unsigned int whichbin = (unsigned int)((x-m_binmin)*m_bins.size()/(m_binmax-m_binmin));
		if (whichbin > 0 && whichbin < m_bins.size()) {
			m_bins[whichbin]++;
			m_count++;
		}
	}
	vector<double> levels() {
		vector<double> levs(m_nbins);
		double binwidth = (m_binmax-m_binmin)/m_nbins;
		for (unsigned int i=0; i<m_nbins; i++) {
			levs[i] = m_binmin + i*binwidth;
		}
		return levs;
	}
	const vector<int>& counts() const {
		return m_bins;
	}

	const int count() const {
		return m_count;
	}

	const int operator[](unsigned int i) {
		int res = -1;
		if (i<m_nbins) {
			res = m_bins[i];
		}
		return res;
	}

	void operator+=(double x) {
		add(x);
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

	vector<bool> isOptimal = fe->getOptimalCodons(false);

	Outcomes outcomes;
	double fop = 0.0;

	int step_count = 0;
	Move move(-1, -1, -1);

	long t = clock();

	int nbins = 10000;
	Histogram faccs(nbins, 0.0, 1.0);
	Histogram frobs(nbins, 0.0, 1.0);
	Histogram fops(nbins, 0.0, 1.0);
	Histogram ffolds(nbins, 0.0, 1.0);

	double cond_cutoff = 0.9633;
	Histogram cond_faccs(nbins, 0.0, 1.0);
	Histogram cond_frobs(nbins, 0.0, 1.0);
	Histogram cond_fops(nbins, 0.0, 1.0);
	Histogram cond_ffolds(nbins, 0.0, 1.0);

	cout << "lev\tfacc\tfrob\tffold\tfop\tcfacc\tcfrob\tcffold\tcfop" << endl;
	while (step_count < p.step_max ) {
		move = getMove(g);
		double prob = calcAcceptanceProbability(move, g, fe, p.bias, outcomes);
		//if (myRand() < prob) {
		if (prob>0.0) {
			// Do move
			g[move.pos] = move.to;
			fop = GenotypeUtil::calcFop(g, isOptimal);
			faccs += outcomes.facc;
			frobs += outcomes.frob;
			fops += fop;
			ffolds += outcomes.ffold;

			if (outcomes.ffold>=cond_cutoff) {
				cond_faccs += outcomes.facc;
				cond_frobs += outcomes.frob;
				cond_fops += fop;
				cond_ffolds += outcomes.ffold;
			}
			step_count++;
		}
	}
	delete fe;
	cout << "# Conditional cutoff = " << cond_cutoff << endl;
	cout << "# Enumerating " << step_count << " genes took " << (clock()-t)/(double)CLK_TCK << " seconds" << endl;
	vector<double> levels = faccs.levels();
	for (unsigned int i=0; i<levels.size(); i++) {
		cout << levels[i] << tab << faccs[i] << tab << frobs[i] << tab << ffolds[i] << tab << fops[i] << tab
			 << cond_faccs[i] << tab << cond_frobs[i] << tab << cond_ffolds[i] << tab << cond_fops[i] << endl;
	}
}




