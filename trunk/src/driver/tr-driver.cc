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


#include "translation-experiment.hh"
#include "compact-lattice-folder.hh"

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
	CompactLatticeFolder* folder = new CompactLatticeFolder(side_length);

	cout << p;

	// Choose the FitnessEvaluator based on input parameters (p.eval_type).
	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate );
		fe = rob;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(1);
	}

	cout << setprecision(4);
	evolutionExperiment( p, *fe );
	//evolutionTest( p, *fe );
	delete fe;
	delete folder;

	return 0;
}




