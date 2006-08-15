#ifndef _T_FOLDER_H__
#define _T_FOLDER_H__
#include "cutee.h"
#include "protein-folder.hh"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include <fstream>

struct TEST_CLASS( folder_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;
	const static double myeps = 1e-2;
	void TEST_FUNCTION( init_lattice )
	{
		ProteinFolder* folder = new CompactLatticeFolder(side_length);
		dynamic_cast<CompactLatticeFolder*>(folder)->enumerateStructures();
		Protein p("CSVMQGGKTVFQMPIIERVMQAYNI"); //Gene::createRandomNoStops(gene_length);
		FoldInfo fi = p.fold(*folder);
		TEST_ASSERT(abs(fi.getFreeEnergy()-0.564) < 1e-2);
		TEST_ASSERT(fi.getStructure() == (StructureID)225);
		delete folder;
		return;
	}
	void TEST_FUNCTION( init_decoy )
	{
		vector<DecoyContactStructure*> structs;
		ifstream fin;
		int n_structs = 10;
		for (int i=0; i<n_structs; i++) {
			structs.push_back(new DecoyContactStructure());
			fin.open("");
			structs[i]->read(fin);
			fin.close();
		}

		int length = 100;
		double log_nconf = 160.0*log(10.0);
		ProteinFolder* folder = new DecoyContactFolder(length, log_nconf, structs);
		Gene g = Gene::createRandomNoStops(length);
		Protein p = g.translate();
		FoldInfo fi = p.fold(*folder);
		delete folder;
		for (int i=0; i<n_structs; i++) {
			delete structs[i];
		}

		return;
	}
};


#endif // _T_FOLDER_H__
