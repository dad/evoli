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
#include "folder-util.hh"
#include "decoy-contact-folder.hh"
#include "mutator.hh"

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>

int main() {
	int protein_length = 300;
	double log_nconf = 100.0*log(10.0);
	ifstream fin("test/data/williams_contact_maps/maps.txt");
	DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
	if (!folder.good() )
		return 1;

	double max_dg = 0.0;
	int sid = 0;
	Random::seed(11);
	int nfolded = folder.getNumFolded();
	CodingDNA g = FolderUtil::getSequenceForStructure( folder, 3*protein_length, max_dg, sid);
	cout << "# Found initial sequence: " << folder.fold(g.translate())->getDeltaG() << endl;;

	CodingDNA g2 = g;
	Polymerase poly(0.001);
	// Evolve under low mutation pressure
	uint max_steps = 10000;
	uint fold_steps = 0;
	time_t start_time = time(NULL);
	for (uint step=0; step<max_steps;) {
		g2 = g;
		bool mutated = poly.mutate(g2);
		if (mutated) {
			Protein p = g2.translate();
			auto_ptr<FoldInfo> fi(folder.fold(p));
			fold_steps++;
			if (fi->getDeltaG()<=max_dg) {
				g = g2;
				step++;
			}
		}
	}
	cout << "# Evolution: " << (time(NULL)-start_time) << " seconds (" << fold_steps << " folded)" << endl;
	return 0;
}
