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

/** \page  neutral-evolve neutral-evolve
This program operates similarly to \c sequence-generator. However, it does not always start anew from random, but rather continues to evolve the previously found sequence. The command-line parameter "neutral steps" indicates how many steps of neutral evolution should be made between two sequences that are reported.
*/


#include "compact-lattice-folder.hh"
#include "translator.hh"
#include "folder-util.hh"
#include "tools.hh"
#include "random.hh"

#include <fstream>


struct Parameters
{
	string folder_type;
	int protein_length;
	double free_energy_cutoff;
	int repetitions;
	int neutral_steps;
	int random_seed;
	int struct_id;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#   free energy cutoff: " << p.free_energy_cutoff << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   neutral steps: " << p.neutral_steps << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   target structure id: " << p.struct_id << endl;
	s << "#" << endl;
	return s;
}

Parameters getParams( int ac, char **av )
{
	if ( ac != 7 )
	{
		cout << "Start program like this:" << endl;
		cout << "  " << av[0] << " <prot length> <free_energy_cutoff> <repetitions> <neutral steps> <random seed> [<struct id>|-1]" << endl;
		exit (-1);
	}

	Parameters p;
	int i = 1;
	p.protein_length = atoi( av[i++] );
	p.free_energy_cutoff = atof( av[i++] );
	p.repetitions = atoi( av[i++] );
	p.neutral_steps = atoi( av[i++] );
	p.random_seed = atoi( av[i++] );
	p.struct_id = atoi( av[i++] );

	return p;
}


// do neutral evolution and find next sequence fitting the specified criteria
Gene getNextSequence( CodingDNA s, Folder &b, const Parameters &p, ostream &out, int nfolded )
{
	double mutation_rate = 1.0/(3.*p.protein_length);
	SimpleMutator mut(mutation_rate);
	int count = 0;
	auto_ptr<FoldInfo> fdata;
	Protein prot( p.protein_length );
	do {
		CodingDNA s2 = s;
		bool changed = false;
		do {
			changed = mut.mutate(s2);
		} while (!changed);

		if ( !s2.encodesFullLength() )
			continue;

		prot = s2.translate();
		fdata = auto_ptr<FoldInfo>( b.fold(prot) );
			
		if ( p.struct_id >= 0 && p.struct_id != fdata->getStructure() )
			continue;
		
		if ( fdata->getDeltaG() >= p.free_energy_cutoff )
			continue;

		// we found a valid sequence
		s = s2;
		count += 1;
	}
	while ( count <= p.neutral_steps );

	out << s << " " << fdata->getDeltaG() << " " << fdata->getStructure() << " " << FolderUtil::calcNeutrality( b, prot, p.free_energy_cutoff )
	  << " " << (b.getNumFolded()-nfolded) << endl;
	return s;
}


int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );

	// set random seed
	Random::seed( p.random_seed );

	int side_length = (int)(sqrt(float(p.protein_length)));
	// initialize the protein folder
	CompactLatticeFolder b(side_length);

	cout << p;
	cout << "# <sequence> <free energy> <structure id> <neutrality> <nfolded>" << endl;

	CodingDNA s;
	int nfolded = 0;
	// get starting sequence
	if ( p.struct_id < 0 )
		s = FolderUtil::getSequence(b, 3*p.protein_length, p.free_energy_cutoff);
	else
		s = FolderUtil::getSequenceForStructure(b, 3*p.protein_length, p.free_energy_cutoff, p.struct_id);

	
	for ( int i=0; i<p.repetitions; i++ )
	{
		s = getNextSequence( s, b, p, cout, nfolded );
		nfolded = b.getNumFolded();
	}
}




