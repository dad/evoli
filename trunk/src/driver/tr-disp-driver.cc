/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006, 2008 Claus Wilke <cwilke@mail.utexas.edu>,
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
#include "protein.hh"
#include "gene-util.hh"
#include "folder-util.hh"
#include "mutator.hh"
#include "compact-lattice-folder.hh"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

class Parameters {
public:
	string eval_type;
	unsigned int protein_length;
	double free_energy_cutoff;
	double u;
	double tr_cost;
	string tr_cost_str;
	double diff_cost;
	string diff_cost_str;
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
		if ( ac < 15 )	{
			valid = false;
			cout << "Start program like this:" << endl;
			cout << "\t" << av[0] << " <eval. type> <prot length> <pop size> <log10 tr cost> <log10 diff cost> <ca cost> <target frac. accurate> <structure id> <free energy cutoff> <mutation rate> <window time> <equilibration time> <repetitions> <random seed> <run ID>" << endl;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		protein_length = atoi( av[i++] );
		N = atoi( av[i++] );
		tr_cost_str = av[i++];
		tr_cost = pow(10.0,atof( tr_cost_str.c_str() ));
		diff_cost_str = av[i++];
		diff_cost = atof( diff_cost_str.c_str() ); //pow(10.0,atof( diff_cost_str.c_str() ));
		ca_cost = atof( av[i++] );
		target_fraction_accurate = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		u = atof( av[i++] );
		window_size = atoi( av[i++] );
		equilibration_time = atoi( av[i++] );
		repetitions = atoi( av[i++] );
		random_seed = atoi( av[i++] );
		if (ac==16){
			run_id = av[i++];
		}
		else{
			run_id = itoa(random_seed, 10);
		}

		valid = true;
	}
};

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   per-site mutation rate u: " << p.u << endl;
	s << "#   transl. robust. cost factor: " << p.tr_cost << endl;
	s << "#   diff. cost factor: " << p.diff_cost << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   target fraction accurate: " << p.target_fraction_accurate << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   window size tau: " << p.window_size << endl;
	s << "#   equilibration time: " << p.equilibration_time << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   run ID: " << p.run_id << endl;
	s << "#" << endl;
	return s;
}

int getStructureID( Folder &b, const Gene &g );
bool analyzeReplica( ErrorproneTranslation *fe, const Parameters &p, ostream &s,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop );
void evolutionExperiment( const Parameters &p, ErrorproneTranslation* fe, Polymerase& poly);

/** \page tr-fxnloss-driver tr-fxnloss-driver
The program \c tr-fxnloss-driver is used to run actual translational robustness experiments. It takes a number of command-line parameters. A typical call is:
\verbatim
   ./tr-fxnloss-driver 25 1000 0.00 3 0.85 599 -5 0.00001 100000 50000 50 111 LVLRRPCNRINSSMPDIWFLALDKK
                     0 1  2    3    4 5    6   7  8       9      10    11 12  13
\endverbatim
(The numbers under the call count the parameters.)

The parameters are, in order (the zeroth parameter is the program name):
-#  %Protein length in amino acids
-#  %Population size
-#  Base-10 logarithm of the translational cost factor
-#  Base-10 logarithm of the difference cost factor
-#  Codon adaptation cost -- average accuracy ratio between optimal and non-optimal codons
-#  Target translational accuracy -- fraction of time a random gene is accurately translated
-#  Structure ID -- index of structure
-#  Maximum delta free energy of folding
-#  Mutation rate per nucleotide
-#  Window time -- number of generations for which evolutionary data is being collected
-#  Equilibration time -- evolution time in generations before collection of evolutionary data begins
-#  Number of replicates
-#  Random number seed -- best to use an odd number.  n and n+1 yield identical results if n is even.

*/



int main( int ac, char **av)
{
	Parameters p( ac, av );

	if (!p.valid) {
		exit(1);
	}

	// seed the random number generator
	Random::seed(p.random_seed);

	// initialize the protein folder
	int side_length = (int)(sqrt(float(p.protein_length)));
	assert(side_length*side_length == p.protein_length);
	CompactLatticeFolder folder(side_length);

	cout << p;
	// Create Polymerase based on input parameter p.mutation_rate
	double GCtoAT = 69.;
	double ATtoGC = 41.;
	double GCtoTA = 88.;
	double GCtoCG = 52.;
	double ATtoCG = 5.;
	double ATtoTA = 11.;
	//Polymerase poly(p.u, GCtoAT, ATtoGC, GCtoTA, GCtoCG, ATtoCG, ATtoTA );
	Polymerase poly(p.u);

	cout << setprecision(4);
	if (p.eval_type == "tr") {
		ErrorproneTranslation ept( &folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.target_fraction_accurate );
		evolutionExperiment( p, &ept, poly);
	}
	else if (p.eval_type == "disp") {
		DispensabilityErrorproneTranslation dept(&folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.target_fraction_accurate, p.diff_cost );
		evolutionExperiment( p, &dept, poly);
	}

	return 0;
}


StructureID getStructureID( Folder &b, const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( b.fold(p) );
		return fi->getStructure();
	}
	else
		return (StructureID)-1;
}

bool runAndAnalyzeReplica( ErrorproneTranslation *fe, Polymerase *poly, const Parameters &p, ostream &s, vector<bool>& is_optimal,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop )
{
	// initialize the population
	Population<Gene, ErrorproneTranslation, Polymerase> pop( p.N );

	Folder& folder = *(fe->getFolder());
	// Find a sequence.
	CodingDNA g = FolderUtil::getSequenceForStructure(folder, p.protein_length*3, p.free_energy_cutoff, p.structure_ID);
	s << "# Starting genotype: " << g << endl;
	// Fill the population with the genotype that we found above
	pop.init( g, fe, poly );
	// If using AccuracyOnlyTranslation, initialize the evaluator with this sequence.
	/*AccuracyOnlyTranslation *aot = dynamic_cast<AccuracyOnlyTranslation*>(fe);
	if (aot != NULL) {
		Protein target_prot = g.translate();
		aot->setTargetSequence(target_prot);
		}*/

	int loop_length = 100;
	int n_folded = 0;
	for ( int i=0; ; i++ ) 	{
		for ( int j=0; j<loop_length; j++ ) {
			pop.evolve();
		}
		if ( pop.calcCoalescenceTime() > p.equilibration_time + p.window_size )
			break;
		cout << "t=" << i*loop_length+1 << "; " << flush;
	}
	cout << endl;

	pop.printGenebank( s );

	GenebankAnalyzer<Gene> analyzer(pop.getGenebank(), pop.begin(), pop.end(), pop.getNumGenerations(), pop.calcCoalescenceTime() );
	return analyzer.analyzeDnDs( p.window_size, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, is_optimal );
}

void evolutionExperiment( const Parameters &p, ErrorproneTranslation* fe, Polymerase& poly) {
	// set random seed
	Random::seed( p.random_seed );

	vector<bool> is_optimal = fe->getOptimalCodons(true);

	// set up output file
	stringstream fname;
	fname << "run-" << p.eval_type << "_s" << p.structure_ID << "tr" << p.tr_cost_str << "di" << p.diff_cost_str << "ca" << p.ca_cost << "-id" << p.run_id << ".dat";
	string filename = fname.str();

	ofstream data_file( filename.c_str(), ios::out );

	data_file << p;
	data_file << "# <nonsyn. substs.> <syn. substs.> <nonsyn. sites> <syn. sites> <fitness> <fop>" << endl;

	double dn_s1, dn_s2, ds_s1, ds_s2, N_s1, N_s2, S_s1, S_s2, f_s1, f_s2, fop_s1, fop_s2;
	double dn, ds, N, S, f, fop;
	dn_s1 = dn_s2 = ds_s1 = ds_s2 = N_s1 = N_s2 = S_s1 = S_s2 = f_s1 = f_s2 = fop_s1 = fop_s2 = 0;

	int count = 0;
	for ( int i=0; i<p.repetitions; i++ )
	{
		stringstream repfname;
		repfname << "run-" << p.eval_type << "_s" << p.structure_ID << "tr" << p.tr_cost_str << "di" << p.diff_cost_str << "ca" << p.ca_cost << "-gb-rep" << i << "-id" << p.run_id << ".dat";
		filename = repfname.str();
		ofstream gen_file( filename.c_str(), ios::out );
		gen_file << p;

		if ( runAndAnalyzeReplica( fe, &poly, p, gen_file, is_optimal, dn, ds, N, S, f, fop ) )
		{
			count += 1;
			data_file << "# " << dn << tab << ds << tab << N << tab << S << tab << f
				  << tab << fop << endl;
			dn_s1 += dn; dn_s2 += dn*dn;
			ds_s1 += ds; ds_s2 += ds*ds;
			N_s1 += N; N_s2 += N*N;
			S_s1 += S; S_s2 += S*S;
			f_s1 += f; f_s2 += f*f;
			fop_s1 += fop; fop_s2 += fop*fop;
		}
		cout << "[" << i+1 << "/" << p.repetitions << "] " << endl;
	}
	cout << endl;
	Folder* folder = fe->getFolder();
	if ( folder != NULL) {
		cout << "# Folded " << folder->getNumFolded() << " proteins" << endl;
		data_file << "# Folded " << folder->getNumFolded() << " proteins" << endl;
	}


	data_file << "# Summary:\n# <tr cost> <ca cost> <pop. size> <nonsyn. substs.> <var> <syn. substs.> <var> <nonsyn. sites> <var> <syn. sites> <var> <ave. fitness> <var> <ave fop> <var>" << endl;

	data_file << p.tr_cost << tab << p.ca_cost << tab << p.N << tab;
	pair<double, double> stat = meanvar( dn_s1, dn_s2, count );
	data_file << stat.first << tab << stat.second << tab;
	stat = meanvar( ds_s1, ds_s2, count );
	data_file << stat.first << tab << stat.second << tab;
	stat = meanvar( N_s1, N_s2, count );
	data_file << stat.first << tab << stat.second << tab;
	stat = meanvar( S_s1, S_s2, count );
	data_file << stat.first << tab << stat.second << tab;
	stat = meanvar( f_s1, f_s2, count );
	data_file << stat.first << tab << stat.second << tab;
	stat = meanvar( fop_s1, fop_s2, count );
	data_file << stat.first << tab << stat.second << endl;
}

