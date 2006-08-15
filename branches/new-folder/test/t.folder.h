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
	void TEST_FUNCTION( init_lattice )
	{
		ProteinFolder* lattice_folder = new CompactLatticeFolder(side_length);
		dynamic_cast<CompactLatticeFolder*>(lattice_folder)->enumerateStructures();
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		FoldInfo fi = p.fold(*lattice_folder);
		return;
	}
	void TEST_FUNCTION( init_decoy )
	{
		vector<DecoyContactStructure> structs;
		ifstream fin("");
		for (int i=0; i<10; i++) {
			structs.push_back(DecoyContactStructure());
			fin.open("");
			structs[i].read(fin);
			fin.close();
		}

		int length = 100;
		ProteinFolder* folder = new DecoyContactFolder(length, structs);
		Gene g = Gene::createRandomNoStops(length);
		Protein p = g.translate();
		FoldInfo fi = p.fold(*folder);
		delete folder;
		return;
	}
};


#endif // _T_FOLDER_H__
