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


#ifndef _T_FITNESSEVALUATOR_H__
#define _T_FITNESSEVALUATOR_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "protein.hh"
#include "gene-util.hh"
#include "tools.hh"
#include <fstream>
#include <cmath>

struct TEST_CLASS( fitness_evaluator_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( create_EPT )
	{
		CompactLatticeFolder folder(side_length);
		ErrorproneTranslation ept(&folder, gene_length, 599, -5, 1, 6, 0.1, 0.1, 0.1 );
		TEST_ASSERT(true);
	}

	void TEST_FUNCTION( create_EPT_without_weights )
	{
		CompactLatticeFolder folder(side_length);
		ErrorproneTranslation ept(&folder, gene_length, 599, -2.0, 1.0, 6.0, 0.85);
		Gene test_gene = Gene::createRandomNoStops(gene_length);

		double fitness = ept.getFitness(test_gene);
		TEST_ASSERT(fitness >= 0 && fitness <= 1);
	}

	void TEST_FUNCTION( test_EPT_approximation_for_accuracy ) {
		CompactLatticeFolder folder(side_length);
		double target_accuracy = 0.85;
		double max_dg = -5;
		int sid = 599;
		//ErrorproneTranslation ept(&folder, gene_length, sid, max_dg, 1.0, 6.0, target_accuracy);
		ErrorproneTranslation ept(&folder, gene_length, sid, max_dg, 1.0, 6.0, 0.0114, 59.0, 104.5);
		Gene g = GeneUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);

		TEST_ASSERT_M(ept.getFolded(g), "Generated test gene not folded.");

		Accumulator accuracies;
		// Evolve while preserving fold for tot_equil steps to
		// equilibrate, then for tot_rand steps, recording weights.
		int num_rand = 200;
		int num_equil = 4000;
		int nrand=0, nequil=0;
		while ( nrand < num_rand ) {
			int randpos = Random::rint(g.codonLength());
			// go through all possible point mutations
			int from_codon = g[randpos];
			int to_codon = from_codon;
			do {
				to_codon = Random::rint(64);
			} while (to_codon == from_codon);

			g[randpos] = to_codon;
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
				g[randpos] = from_codon;
			}
		}
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
};

#endif
