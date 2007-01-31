/*
This file is part of the E.voli project.
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


#include "translation-experiment.hh"
#include "compact-lattice-folder.hh"
#include <sstream>


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   per-site mutation rate u: " << p.u << endl;
	s << "#   transl. robust. cost factor: " << p.tr_cost << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   window size tau: " << p.window_size << endl;
	s << "#   equilibration time: " << p.equilibration_time << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   run ID: " << p.run_id << endl;
	s << "#" << endl;
	return s;
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

void evolutionTest( const Parameters &p, ErrorproneTranslation& fe) {
	if ( p.structure_ID < 0 )
	{
		cerr << "Input sequence does not translate!" << endl;
		cerr << "Exiting. No data written." << endl;
		exit( - 1 );
	}

	// set random seed
	Random::seed( p.random_seed );

	Population pop( p.N );
	Folder& folder = *(fe.getFolder());
	Gene g = GeneUtil::getSequenceForStructure(folder, p.protein_length, p.free_energy_cutoff, p.structure_ID);
	pop.init(g, &fe, p.u);

	vector<bool> is_optimal = fe.getOptimalCodons();

	cout << "gen\twav\tfop\tdg\tfacc\tfrob\tftrunc\tffold" << endl;
	int gap = 100;

	Accumulator fops, faccs, frobs, ftruncs, ffolds, fitnesses, dgs;

	for ( int i=0; ; i++ )	{
		fitnesses.reset();
		fops.reset();
		faccs.reset();
		frobs.reset();
		ftruncs.reset();
		ffolds.reset();
		dgs.reset();
		for (int k=0; k<p.N; k++) {
			const Gene& g = pop[k];

			bool folded = fe.getFolded(g);
			double dG = 0;
			if (folded) {
				Protein p = g.translate();
				auto_ptr<FoldInfo> fold_data( folder.fold(p) );
				dG = fold_data->getFreeEnergy();
			}
			double facc, frob, ftrunc, ffold;
			fe.calcOutcomes( g, facc, frob, ftrunc, ffold);
			double fop = GeneUtil::calcFop( g, is_optimal );

			fops += fop;
			faccs += facc;
			frobs += frob;
			ftruncs += ftrunc;
			ffolds += ffold;
			dgs += dG;
		}
		double w_av = pop.getAveFitness();
		cout << i*gap << tab << setprecision(6) << w_av << tab << (double)fops << tab << (double)dgs << tab
			 << (double)faccs << tab << (double)frobs << tab << (double)ftruncs << tab << (double)ffolds << endl;

		for ( int j=0; j<gap; j++ )	{
			pop.evolve();
		}
	}
}

bool runAndAnalyzeReplica( ErrorproneTranslation *fe, const Parameters &p, ostream &s, vector<bool>& is_optimal,
		     double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f,
		     double &ave_fop )
{
	// initialize the population
	Population pop( p.N );

	Folder& folder = *(fe->getFolder());
	// Find a sequence.
	Gene g = GeneUtil::getSequenceForStructure(folder, p.protein_length*3, p.free_energy_cutoff, p.structure_ID);
	s << "# Starting genotype: " << g << endl;
	// Fill the population with the genotype that we found above
	pop.init( g, fe, p.u );
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
		pop.prepareCoalescenceCalcs();
		if ( pop.calcCoalescenceTime() > p.equilibration_time + p.window_size )
			break;
		cout << "t=" << i*loop_length+1 << "; " << flush;
	}
	cout << endl;

	pop.printGenebank( s );

	return pop.analyzeDnDs( p.window_size, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, is_optimal );
}

void evolutionExperiment( const Parameters &p, ErrorproneTranslation& fe)
{
	// set random seed
	Random::seed( p.random_seed );

	vector<bool> is_optimal = fe.getOptimalCodons(true);

	// set up output file
	stringstream fname;
	fname << "run-" << p.eval_type << "_s" << p.structure_ID << "tr" << p.tr_cost_str << "ca" << p.ca_cost << "-id" << p.run_id << ".dat";
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
		repfname << "run-" << p.eval_type << "_s" << p.structure_ID << "tr" << p.tr_cost_str << "ca" << p.ca_cost << "-gb-rep" << i << "-id" << p.run_id << ".dat";
		filename = repfname.str();
		ofstream gen_file( filename.c_str(), ios::out );
		gen_file << p;

		if ( runAndAnalyzeReplica( &fe, p, gen_file, is_optimal, dn, ds, N, S, f, fop ) )
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
	Folder* folder = fe.getFolder();
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
