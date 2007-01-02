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


#include "random.hh"
#include "folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"
#include "compact-lattice-folder.hh"
#include "mutator.hh"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

bool runAndAnalyzeReplica( ErrorproneTranslation *fe, vector<bool>& is_optimal,
						   double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
						   double &ave_fop );

void evolutionExperiment( ErrorproneTranslation& fe);

int main( int ac, char **av)
{
	cout << "# Starting performance test" << endl;
	time_t overall_time = time(NULL);
	Random::seed(13);
	cout << "# Initializing folder" << endl;
	time_t start_time = time(NULL);
	// initialize the protein folder
	int side_length = 5;
	CompactLatticeFolder folder(side_length);
	time_t duration = time(NULL) - start_time;
	cout << "# Folder initialization: " << duration << " seconds" << endl;
	cout << "# Initializing fitness evaluator" << endl;
	start_time = time(NULL);
	// Choose the FitnessEvaluator based on input parameters (p.eval_type).
	ErrorproneTranslation ept( &folder, side_length*side_length, 599, -5, 100, 6, 0.0114, 59.0, 104.5 );
	cout << "# FE initialization: " << (time(NULL)-start_time) << " seconds" << endl;
	
	cout << "# Running evolution" << endl;
	start_time = time(NULL);
	evolutionExperiment( ept );
	cout << "# Evolution: " << (time(NULL)-start_time) << " seconds" << endl;
	
	cout << "# Performance test took " << (time(NULL)-overall_time) << " seconds" << endl;
	cout << "# Folded " << folder.getNumFolded() << " proteins" << endl;

	return 0;
}


void evolutionExperiment( ErrorproneTranslation& fe)
{
	vector<bool> is_optimal = fe.getOptimalCodons(false);

	double dn_s1, dn_s2, ds_s1, ds_s2, N_s1, N_s2, S_s1, S_s2, f_s1, f_s2, fop_s1, fop_s2;
	double dn, ds, N, S, f, fop;
	dn_s1 = dn_s2 = ds_s1 = ds_s2 = N_s1 = N_s2 = S_s1 = S_s2 = f_s1 = f_s2 = fop_s1 = fop_s2 = 0;

	int count = 0;
	int reps = 1;
	for ( int i=0; i<reps; i++ )
	{
		stringstream repfname;

		if ( runAndAnalyzeReplica( &fe, is_optimal, dn, ds, N, S, f, fop ) )
		{
			count += 1;
			dn_s1 += dn; dn_s2 += dn*dn;
			ds_s1 += ds; ds_s2 += ds*ds;
			N_s1 += N; N_s2 += N*N;
			S_s1 += S; S_s2 += S*S;
			f_s1 += f; f_s2 += f*f;
			fop_s1 += fop; fop_s2 += fop*fop;
		}
		//cout << "[" << i+1 << "/" << reps << "] " << endl;
	}
}

bool runAndAnalyzeReplica( ErrorproneTranslation *fe, vector<bool>& is_optimal,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop )
{
	// initialize the population
	Population<Gene, ErrorproneTranslation, SimpleMutator> pop( 1000 );
	SimpleMutator mut(0.00001);

	Folder& folder = *(fe->getFolder());
	// Find a sequence.
	Gene g = GeneUtil::getSequenceForStructure(folder, 75, -5, 599);
	cout << "nf: " <<  folder.getNumFolded() << endl;
	// Fill the population with the genotype that we found above
	pop.init( g, fe, &mut );
	GenebankAnalyzer<Gene> analyzer(pop.getGenebank());

	int loop_length = 100;
	int n_folded = 0;
	for ( int i=0; ; i++ ) 	{
		for ( int j=0; j<loop_length; j++ ) {
			pop.evolve();
		}
		cout << "nf ev: " << folder.getNumFolded() << endl;
		analyzer.prepareCoalescenceCalcs(pop.begin(), pop.end(), pop.getNumGenerations());
		if ( analyzer.calcCoalescenceTime() > 2000 )
			break;
	}
	cout << "# Evolved for " << pop.getNumGenerations() << " generations" << endl;

	return analyzer.analyzeDnDs( 1000, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, is_optimal );
}
