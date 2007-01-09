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


#ifndef FOLDER_UTIL_HH
#define FOLDER_UTIL_HH

#include "coding-sequence.hh"
#include "folder.hh"
#include "mutator.hh"

#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

class FolderUtil
{
public:
	/**
	 * Calculates the neutrality of the given protein. Cutoff is the free energy cutoff
	 * below which the protein folds.
	 **/
	static double calcNeutrality( const Folder &b, Protein p, double cutoff )
	{
		auto_ptr<FoldInfo> fold_data( b.fold(p) );
		if ( fold_data->getDeltaG() > cutoff )
			return 0;
		int structure_id = fold_data->getStructure();

		int count = 0;
		// go through all positions in the protein
		for ( unsigned int i=0; i<p.length(); i++ )	{
			// go through all possible point mutations
			char oldaa = p[i];
			for (int j=0; j<20; j++) {
				char newaa = GeneticCodeUtil::AMINO_ACIDS[j];
				if (newaa == oldaa)
					continue;
				p[i] = newaa;
				// sequence folds into correct structure with low free energy?
				fold_data = auto_ptr<FoldInfo>( b.fold(p) );
				if (fold_data->getStructure() == structure_id && fold_data->getDeltaG() < cutoff) {
					count += 1;
				}
			}
			p[i] = oldaa;
		}
		return count / (19.0*p.length());
	}

	/**
	 * Finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
	 */
	static CodingDNA getSequenceForStructure( const Folder &b, unsigned int length, double deltag_cutoff, const int struct_id )
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		CodingDNA g(length);
		auto_ptr<FoldInfo> fdata;
		double mutation_rate = 1.0/length;
		SimpleMutator mut(mutation_rate);
		bool found = false;
		double min_free_energy_for_starting = max(300.0, deltag_cutoff);

		//int nfolded = b.getNumFolded();
		// find sequence that encodes the desired structure
		//int q=0;
		do {
			g = CodingDNA::createRandomNoStops( length );
			//cout << g << endl;
			Protein p = g.translate();
			fdata = auto_ptr<FoldInfo>( b.fold(p) );
			found = (fdata->getStructure() == struct_id && fdata->getDeltaG() <= min_free_energy_for_starting);
			//cout << fdata->getStructure() << "\t" << fdata->getDeltaG() << "\t" << g << endl;
			//cout << q++ << " " << g << " " << p << " " << fdata->getStructure() << " " << fdata->getDeltaG() << endl;
		} while ( !found );
		//cout << "nf: " << (b.getNumFolded() - nfolded) << endl;

		//cout << "hey" << endl;

		int fail_count = 0;
		int total_fail_count = 0;
		G = fdata->getDeltaG();

		// optimize the sequence for stability
		do {
			CodingDNA g2 = g;
			bool changed = false;
			do {
				changed = mut.mutate(g2);
			} while (!changed);

			if (!g2.encodesFullLength()) {
				fail_count++;
				total_fail_count++;
			}
			else {
				Protein p = g2.translate();
				fdata = auto_ptr<FoldInfo>( b.fold(p) );
				//cout << fdata->getStructure() << "\t" << fdata->getDeltaG() << "\t" << g << endl;
				if (fdata->getStructure() == struct_id && fdata->getDeltaG() <= G-0.001) {
					// we found an improved sequence. grab it, and reset failure count
					g = g2;
					G = fdata->getDeltaG();
					fail_count = 0;
					//cout << G << endl;
				}
			}

			// start again with random genotype if search is not successful after 50000 failures
			if ( fail_count > 50000 || total_fail_count > 1e6 )	{
				found = false;
				do {
					g = CodingDNA::createRandomNoStops( length );
					Protein p = g.translate();
					fdata = auto_ptr<FoldInfo>( b.fold(p) );
					found = (fdata->getStructure() == struct_id && fdata->getDeltaG() <= min_free_energy_for_starting);
				} while ( !found );
				G = fdata->getDeltaG();
				fail_count = 0;
				total_fail_count = 0;
			}
		} while( G > deltag_cutoff );

		return g;
	}

	/**
	 * Finds a random sequence with folding energy smaller than cutoff
	 */
	static CodingDNA getSequence( const Folder &b, unsigned int length, double free_energy_cutoff)
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		CodingDNA g(length);
		auto_ptr<FoldInfo> fdata;
		double mutation_rate = 1.0/length;
		SimpleMutator mut(mutation_rate);
		bool found = false;
		double min_free_energy_for_starting = max(300.0, free_energy_cutoff);

		// find sequence that encodes the desired structure
		do {
			g = CodingDNA::createRandomNoStops( length );
			//cout << g << endl;
			Protein p = g.translate();
			//cout << p << endl;
			fdata = auto_ptr<FoldInfo>( b.fold(p) );
			found = (fdata->getDeltaG() <= min_free_energy_for_starting);
			//cout << fdata.first << "\t" << fdata.second << "\t" << g << endl;
		} while ( !found );

		int fail_count = 0;
		int total_fail_count = 0;
		G = fdata->getDeltaG();

		// optimize the sequence for stability
		do {
			CodingDNA g2 = g;
			bool changed = false;
			do {
				changed = mut.mutate(g2);
			} while (!changed);

			if (!g2.encodesFullLength()) {
				fail_count++;
				total_fail_count++;
			}
			else {
				Protein p = g2.translate();
				fdata = auto_ptr<FoldInfo>( b.fold(p) );
				//cout << fdata.first << "\t" << fdata.second << "\t" << G << endl;
				if (fdata->getDeltaG() <= G-0.001) {
					// we found an improved sequence. grab it, and reset failure count
					g = g2;
					G = fdata->getDeltaG();
					fail_count = 0;
					//cout << G << endl;
				}
			}

			// start again with random genotype if search is not successful after 50000 failures
			if ( fail_count > 50000 || total_fail_count > 1e6 )	{
				found = false;
				do {
					g = CodingDNA::createRandomNoStops( length );
					Protein p = g.translate();
					fdata = auto_ptr<FoldInfo>( b.fold(p) );
					found = (fdata->getDeltaG() <= min_free_energy_for_starting);
				} while ( !found );
				G = fdata->getDeltaG();
				fail_count = 0;
				total_fail_count = 0;
			}
		} while( G > free_energy_cutoff );

		return g;
	}

};


#endif //FOLDER_UTIL_HH
