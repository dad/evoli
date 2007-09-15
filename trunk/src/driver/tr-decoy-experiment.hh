/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


#include "folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"
#include "folder-util.hh"
#include "mutator.hh"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

class Parameters {
public:
	string eval_type;
	unsigned int protein_length;
	string contact_map_dir;
	double log_nconf;
	double free_energy_cutoff;
	double free_energy_minimum;
	double u;
	double tr_cost;
	string tr_cost_str;
	double ca_cost;
	double target_fraction_accurate;
	int repetitions;
	int window_size;
	int equilibration_time;
	int random_seed;
	int N;
	mutable int structure_ID;
	string run_id;
	bool valid;

	Parameters( int ac, char **av ) {
		if ( ac < 17 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval type> <prot length> <contact map dir> <log10 num confs> <pop size> <log10 tr cost> <ca cost> <target facc> <structure id> <free energy cutoff> <free energy minimum> <mutation rate> <window time> <equilibration time> <repetitions> <random seed> <run ID>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		protein_length = atoi( av[i++] );
		contact_map_dir = av[i++];
		log_nconf = atof( av[i++] );
		N = atoi( av[i++] );
		tr_cost_str = av[i++];
		tr_cost = pow(10.0,atof( tr_cost_str.c_str() ));
		ca_cost = atof( av[i++] );
		target_fraction_accurate = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		u = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );
		if (ac==18){
			run_id = av[i++];
		}
		else{
			run_id = itoa(random_seed, 10);
		}

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p );

int getStructureID( Folder &b, const Gene &g );
bool analyzeReplica( ErrorproneTranslation *fe, const Parameters &p, ostream &s,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop );
void evolutionTest( const Parameters &p, ErrorproneTranslation& fe, Polymerase& poly);
void evolutionExperiment( const Parameters &p, ErrorproneTranslation& fe, Polymerase& poly);

