/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>, Matt Cowperthwaite <mattccowp@mac.com>

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

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sstream>
using namespace std;

class Parameters {
public:

	unsigned int protein_length;
	double u;
	int repetitions;
	int window_size;
	int equilibration_time;
	int random_seed;
	int N;
	string run_id;
	bool valid;

	Parameters( int ac, char **av ) {
		
		if ( ac < 8 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <prot length> <pop size> <mutation rate> ";
			cout << " <window time> <equilibration time> <repetitions> <random seed> <run ID>" << endl;
			return;
		}

		int i = 1;
		protein_length = atoi( av[i++] );
		N = atoi( av[i++] );
		u = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );
		if ( ac == 9 ){
			run_id = av[i++];
		}
		else{
			run_id = itoa( random_seed, 10 );
		}

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p );

int getStructureID( Folder &b, const Gene &g );
bool analyzeReplica( NeutralFitness *fe, const Parameters &p, ostream &s );
void evolutionExperiment( const Parameters &p, NeutralFitness& fe );
