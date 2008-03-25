/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <drummond@alumni.princeton.edu>

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


#ifndef _T_FITNESSEVALUATOR_H__
#define _T_FITNESSEVALUATOR_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "protein.hh"
#include "folder-util.hh"
#include "tools.hh"
#include <fstream>
#include <cmath>

struct TEST_CLASS( fitness_evaluator_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( create_EPT_without_weights )
	{
		CompactLatticeFolder folder(side_length);
		ErrorproneTranslation ept(&folder, gene_length/3, 599, -2.0, 1.0, 6.0, 0.85);
		CodingDNA test_gene = CodingDNA::createRandomNoStops(gene_length);
		Protein p = test_gene.translate();
		auto_ptr<FoldInfo> fi(folder.fold(p));
		double fitness = ept.getFitness(test_gene);
		TEST_ASSERT(fitness >= 0 && fitness <= 1);
	}

	void TEST_FUNCTION( create_EPT )
	{
		CompactLatticeFolder folder(side_length);
		ErrorproneTranslation ept(&folder, gene_length/3, 599, -1, 1, 6, 0.1, 0.1, 0.1 );
		TEST_ASSERT(true);
	}

	void TEST_FUNCTION( test_EPT_approximation_for_accuracy ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -1;
		int sid = 599;
		CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		// Test with pre-discovered weights (generated using ./get-weights 6 -5 11 599 1000 1000 0.85 10)
		ErrorproneTranslation ept(&folder, g.codonLength(), sid, max_dg, 1.0, 6.0, 0.0114735, 57.9439, 102.567);
		TEST_ASSERT_M(ept.getFolded(g), "Generated test gene not folded.");
		Accumulator accuracies;
		// Get the stats
		accumulateStatistics(ept, accuracies, g);

		double std_error = accuracies.stderror();
		double std_dev = accuracies.stdev();
		double mean = accuracies.mean();
		// Really should be using standard error
		double spread = 3*std_dev;
		stringstream ss;
		ss << endl << "Target accuracy = " << target_accuracy << "; actual mean accuracy = " << mean << " +/- " << spread << endl;
		//cout << ss.str() << endl;
		TEST_ASSERT_M( ((target_accuracy <= mean+spread) && (target_accuracy >= mean-spread)), ss.str());
	}

	void TEST_FUNCTION( test_automatic_init_EPT_approximation_for_accuracy ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -1;
		int sid = 599;

		CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		// Test with automatically determined weights.
		ErrorproneTranslation ept(&folder, g.codonLength(), sid, max_dg, 1.0, 6.0, target_accuracy);
		TEST_ASSERT_M(ept.getFolded(g), "Generated test gene not folded.");
		Accumulator accuracies;
		// Get the stats
		accumulateStatistics(ept, accuracies, g);

		double std_error = accuracies.stderror();
		double std_dev = accuracies.stdev();
		double mean = accuracies.mean();
		// Really should be using standard error
		double spread = 3*std_dev;
		stringstream ss;
		ss << endl << "Target accuracy = " << target_accuracy << "; actual mean accuracy = " << mean << " +/- " << spread << endl;
		//cout << ss.str() << endl;
		TEST_ASSERT_M( ((target_accuracy <= mean+spread) && (target_accuracy >= mean-spread)), ss.str());
	}

	void accumulateStatistics(ErrorproneTranslation& ept, Accumulator& accuracies, CodingDNA g) {
		// Evolve while preserving fold for tot_equil steps to
		// equilibrate, then for tot_rand steps, recording weights.
		int num_rand = 20;
		int num_equil = 100;
		int nrand=0, nequil=0;
		while ( nrand < num_rand ) {
			int randpos = Random::rint(g.codonLength());
			// go through all possible point mutations
			Codon from_codon = g.getCodon(randpos);
			Codon to_codon = from_codon;
			do {
				to_codon = GeneticCodeUtil::indexToCodon(Random::rint(64));
			} while (to_codon == from_codon);

			g.setCodon(randpos, to_codon);
			if (ept.getFolded(g)) {
				nequil++;
				if (nequil > num_equil) {
					// We've equilibrated enough; check accuracy.
					double ffold, frob, facc, ftrunc;
					double fitness = ept.calcOutcomes(g, facc, frob, ftrunc, ffold);
					accuracies += facc;
					nrand++;
				}
			}
			else {
				g.setCodon(randpos, from_codon);
			}
		}
	}

	void TEST_FUNCTION( test_functional_loss_EPT ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -1;
		int sid = 599;
		double ca_cost = 5.0;
		double diff_cost = -10.0;
		double tr_cost = 1.0;

		CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		Protein p = g.translate();
		// Test with automatically determined weights.
		FunctionalLossErrorproneTranslation flept(&folder, g.codonLength(), sid, max_dg, tr_cost, ca_cost, target_accuracy, diff_cost, p);
		double ffold, frob, facc, ftrunc;
		flept.calcOutcomes(g, facc, frob, ftrunc, ffold);
		double fitness = flept.getFitness(g);
		// exp{-s[1-(a(1-k/L) + (1-a)r(1-(k+1)/L))]}
		// exp{-s[1-(a + (1-a)r(1-1/L))]}
		double target_fitness = exp(diff_cost*(1.0-(facc + (1-facc)*frob*(1-1.0/25))));
		TEST_ASSERT(abs(fitness - target_fitness) < 1e-6);
	}

	void TEST_FUNCTION( test_functional_loss_EPT_pseudogene ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -1;
		int sid = 599;
		double ca_cost = 5.0;
		double diff_cost = -10.0;
		double tr_cost = 1.0;

		CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		Protein p = g.translate();
		FunctionalLossErrorproneTranslation flept(&folder, g.codonLength(), sid, -5, tr_cost, ca_cost, target_accuracy, diff_cost, p);
		double fitness = flept.getFitness(g);
		// Should be complete knockout, because required -5 stability is too low to get by chance when seeking -1 stability.
		double target_fitness = exp(diff_cost);
		TEST_ASSERT(abs(fitness - target_fitness) < 1e-6);
	}

	void TEST_FUNCTION( test_functional_loss_EPT_get_fitness_virtual ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -1;
		int sid = 699;
		double ca_cost = 5.0;
		double diff_cost = -10.0;
		double tr_cost = 0.0724436;

		Random::seed(27);
		CodingDNA g("UGGCGUAUUCUUGAAAUGGACCGGAUAGACGUCGUTAGAACCGAAAUGAAGCCUUUUAAGAACAAGGAAGUGAAG");
		Protein p = g.translate();
		FunctionalLossErrorproneTranslation flept(&folder, g.codonLength(), sid, -5, tr_cost, ca_cost, target_accuracy, diff_cost, p);
		ErrorproneTranslation* ept = &flept;
		double fitness_fl = flept.getFitness(g);
		double fitness_ept = ept->getFitness(g);
		//cout << fitness_fl << " " << fitness_ept << endl;
		TEST_ASSERT(fitness_fl == fitness_ept);
	}
};

#endif
