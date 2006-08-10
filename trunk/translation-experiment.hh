#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

class Parameters {
public:
	string eval_type;
	double free_energy_cutoff;
	double free_energy_minimum;
	double u;
	double tr_cost;
	string tr_cost_str;
	double ca_cost;
	double error_rate;
	double error_weight;
	double accuracy_weight;
	int repetitions;
	int window_size;
	int equilibration_time;
	int random_seed;
	int N;
	mutable int structure_ID;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac != 16 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval type> <pop size> <log10 tr cost> <ca cost> <error rate> <accuracy weight> <error weight> <structure id> <free energy cutoff> <free energy minimum> <mutation rate> <window time> <equilibration time> <repetitions> <random seed>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		N = atoi( av[i++] );
		tr_cost_str = av[i++];
		tr_cost = pow(10.0,atof( tr_cost_str.c_str() ));
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		u = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p );

int getStructureID( ProteinFolder &b, const Genotype &g );
bool analyzeReplica( ErrorproneTranslation *fe, const Parameters &p, ostream &s,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop );
void evolutionTest( const Parameters &p, ErrorproneTranslation& fe);
void evolutionExperiment( const Parameters &p, ErrorproneTranslation& fe);

