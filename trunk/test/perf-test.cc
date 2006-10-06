#include "random.hh"
#include "folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "gene-util.hh"
#include "compact-lattice-folder.hh"

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
	Random::seed(11);
	// initialize the protein folder
	int side_length = 5;
	CompactLatticeFolder* folder = new CompactLatticeFolder(5);

	// Choose the FitnessEvaluator based on input parameters (p.eval_type).
	ErrorproneTranslation* ept = new ErrorproneTranslation();
	ept->init( folder, 25, 599, -5, 100, 6, 0.0114, 59.0, 104.5 );

	time_t start_time = time(NULL);
	evolutionExperiment( *ept );
	time_t duration = time(NULL) - start_time;
	cout << "# Performance test took " << duration << " seconds" << endl;
	cout << "# Folded " << folder->getNumFolded() << " proteins" << endl;
	delete ept;
	delete folder;

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
	Population pop( 1000 );

	Folder& folder = *(fe->getFolder());
	// Find a sequence.
	Gene g = GeneUtil::getSequenceForStructure(folder, 75, -5, 599);
	// Fill the population with the genotype that we found above
	pop.init( g, fe, 0.00001 );

	int loop_length = 100;
	int n_folded = 0;
	for ( int i=0; ; i++ ) 	{
		for ( int j=0; j<loop_length; j++ ) {
			pop.evolve();
		}
		pop.prepareCoalescenceCalcs();
		if ( pop.calcCoalescenceTime() > 2000 )
			break;
	}

	return pop.analyzeDnDs( 1000, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, is_optimal );
}
