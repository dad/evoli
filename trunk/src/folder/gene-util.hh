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


#ifndef GENE_UTIL_HH
#define GENE_UTIL_HH

#include "random.hh"
#include "codon.hh"
#include "genetic-code.hh"
#include "translator.hh"
#include "folder.hh"
#include "protein.hh"

#include <algorithm>
#include <iostream>
#include <memory>

using namespace std;

class GeneUtil
{
public:

	/**
	 * Calculates the total number of sites in the gene.
	 **/
	static double calcTotalSites( const Gene& g ) {
		return (double)g.length();
	}

	/**
	 * Calculates the number of synonymous sites in the genotype.
	 **/
	static double calcSynonymousSites( const Gene & g )
	{
		Gene::const_iterator it = g.begin(), e = g.end();
		double s = 0;
		for ( ; it != e; it++ )
			s += GeneticCodeUtil::calcSynonymousSites( *it );
		return s;
	}

	/**
	 * Calculates the number of synonymous and nonsynonymous sites in surface and core.
	 **/
	static void calcSNSitesSurfaceCore( double &NSurf, double &NCore, double &SSurf, double &SCore,
					    const Gene &g, const vector<int> &surface )
	{
		unsigned int e = g.codonLength();
		assert( e == surface.size() );
		NSurf = NCore = SSurf = SCore = 0;
		double tmp_N, tmp_S;

		for ( unsigned int i=0; i<e; i++ )
		{
			tmp_S = GeneticCodeUtil::calcSynonymousSites( g[i] );
			tmp_N = 3 - tmp_S;

			if ( surface[i] )
			{
				NSurf += tmp_N;
				SSurf += tmp_S;
			}
			else
			{
				NCore += tmp_N;
				SCore += tmp_S;
			}
		}
	}

	/**
	 * Calculates the number of synonymous and nonsynonymous substitutions
	 * from genotype 1 to genotype 2 (these numbers are not! normalized by
	 * the number of synonymous/nonsynonymous sites).
	 **/
	static void calcDnDs( double &dn, double &ds, const Gene &g1, const Gene &g2 )
	{
		assert( g1.length() == g2.length() );
		int e = g1.codonLength();
		dn = ds =0;
		double tmp_dn, tmp_ds;

		for ( int i=0; i<e; i++ )
		{
			GeneticCodeUtil::calcDnDs( tmp_dn, tmp_ds, g1[i], g2[i] );
			//			 CodonUtil::printCodon( cout, g1[i] );
			//			 cout << " ";
			//			 CodonUtil::printCodon( cout, g2[i] );
			//			 cout << " " << tmp_dn << " " << tmp_ds << endl;
			dn += tmp_dn;
			ds += tmp_ds;
		}
	}

	/**
	 * Same as calcDnDs, but distinguishes between surface/core residues.
	 * The vector "surface" is expected to be of the same length as the protein,
	 * and must contain a 1 at surface positions and a 0 at core positions.
	 **/
	static void calcDnDsSurfaceCore( double &dnSurf, double &dnCore, double &dsSurf, double &dsCore,
					 const Gene &g1, const Gene &g2, const vector<int> &surface )
	{
		assert( g1.codonLength() == g2.codonLength() );
		assert( g1.codonLength() == surface.size() );
		int e = g1.codonLength();
		dnSurf = dnCore = dsSurf = dsCore = 0;
		double tmp_dn, tmp_ds;

		for ( int i=0; i<e; i++ )
		{
			GeneticCodeUtil::calcDnDs( tmp_dn, tmp_ds, g1[i], g2[i] );
			//			 CodonUtil::printCodon( cout, g1[i] );
			//			 cout << " ";
			//			 CodonUtil::printCodon( cout, g2[i] );
			//			 cout << " " << tmp_dn << " " << tmp_ds << endl;
			if ( surface[i] )
			{
				dnSurf += tmp_dn;
				dsSurf += tmp_ds;
			}
			else
			{
				dnCore += tmp_dn;
				dsCore += tmp_ds;
			}
		}
	}

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
		for ( unsigned int i=0; i<p.length(); i++ )
		{

			// go through all possible point mutations
			// (avoid operator %, which can be very slow)
			const int res = p[i];
			int tempres = res + 1;
			while ( tempres < 20 )
			{
				p[i] = tempres;
				// sequence folds into correct structure with low free energy?
				fold_data = auto_ptr<FoldInfo>( b.fold(p) );
				if (fold_data->getStructure() == structure_id && fold_data->getDeltaG() < cutoff) {
					count += 1;
				}
				tempres++;
			}

			tempres = 0;
			while ( tempres < res )
			{
				p[i] = tempres;
				// sequence folds into correct structure with low free energy?
				fold_data = auto_ptr<FoldInfo>( b.fold(p) );
				if (fold_data->getStructure() == structure_id && fold_data->getDeltaG() < cutoff) {
					count += 1;
				}
				tempres++;
			}
			p[i] = res;
		}
		return count / (19.0*p.length());
	}


	/**
	 * Calculates the fraction of optimal codons from the set of optimal codons.
	 **/
	static double calcFop( const Gene &g, const vector<bool>& is_optimal ) {
		int count = 0, tot=0;
		Gene::const_iterator it = g.begin(), e = g.end();
		int methionine = 14;
		int tryptophan = 58;
		int codon = -1;

		for ( ; it != e; it++ )	{
			codon = *it;
			if ( !(codon == methionine || codon == tryptophan)) {
				tot += 1;
				if (is_optimal[codon]) {
					count += 1;
				}
			}
		}
		return count / (double)tot;
	}

	/**
	 * As calcFop, but considers surface and core regions separately
	 **/
	static void calcFopSurfaceCore( double &Fop_surf, double &Fop_core, const Gene &g, const double *codon_costs, const vector<int> &surface )
	{
		assert( g.codonLength() == surface.size() );

		Fop_surf = Fop_core = 0;

		if ( codon_costs )
		{
			int sopt = 0, copt = 0, scount = 0, ccount = 0;

			for ( uint i=0; i<g.codonLength(); i++ )
			{
				if ( surface[i] )
				{
					scount += 1;
					if ( codon_costs[g[i]] == 0 )
						sopt += 1;
				}
				else
				{
					ccount += 1;
					if ( codon_costs[g[i]] == 0 )
						copt += 1;
				}
			}
			Fop_surf = (double) sopt / (double) scount;
			Fop_core = (double) copt / (double) ccount;
		}

	}

	/**
	 * Produces a sequence encoding the same amino acid sequence as gene g, but with randomly chosen codons.
	 * Returns the original sequence if it doesn't translate properly.
	 **/
	static Gene randomizeCodons( const Gene &g ) {
		Gene result(g);
		if ( g.encodesFullLength() ) {
			Protein p = g.translate();
			for ( unsigned int i=0; i<p.length(); i++ ) {
				int num_alts = GeneticCodeUtil::residueToAllCodonsTable[p[i]][0];
				// Now choose from alternatives at random
				int randbin = Random::rint( num_alts );
				result[i] = GeneticCodeUtil::residueToAllCodonsTable[p[i]][randbin+1];
			}
		}
		return result;
	}

	static Gene reverseTranslate(const Protein& p) {
		Gene g(p.length()*3);
		Protein::const_iterator it = p.begin();
		
		for (uint i=0; it != p.end(); it++, i++) {
			int aa = *it;
			int codon = GeneticCodeUtil::residueToCodonTable[aa+1];
			g[i] = codon;
		}
		
		return g;
	}

	/**
	 * Finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
	 */
	static Gene getSequenceForStructure( const Folder &b, unsigned int length, double free_energy_cutoff, const int struct_id )
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		Gene g(length);
		auto_ptr<FoldInfo> fdata;
		double mutation_rate = 1.0/length;
		bool found = false;
		double min_free_energy_for_starting = max(300.0, free_energy_cutoff);

		// find sequence that encodes the desired structure
		do {
			g = Gene::createRandomNoStops( length );
			//cout << g << endl;
			Protein p = g.translate();
			//cout << p << endl;
			fdata = auto_ptr<FoldInfo>( b.fold(p) );
			found = (fdata->getStructure() == struct_id && fdata->getDeltaG() <= min_free_energy_for_starting);
			//cout << fdata.first << "\t" << fdata.second << "\t" << g << endl;
		} while ( !found );

		int fail_count = 0;
		int total_fail_count = 0;
		G = fdata->getDeltaG();

		// optimize the sequence for stability
		do {
			Gene g2 = g;
			bool changed = false;
			do {
				changed = g2.mutate( mutation_rate );
			} while (!changed);

			if (!g2.encodesFullLength()) {
				fail_count++;
				total_fail_count++;
			}
			else {
				Protein p = g2.translate();
				fdata = auto_ptr<FoldInfo>( b.fold(p) );
				//cout << fdata.first << "\t" << fdata.second << "\t" << G << endl;
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
					g = Gene::createRandomNoStops( length );
					Protein p = g.translate();
					fdata = auto_ptr<FoldInfo>( b.fold(p) );
					found = (fdata->getStructure() == struct_id && fdata->getDeltaG() <= min_free_energy_for_starting);
				} while ( !found );
				G = fdata->getDeltaG();
				fail_count = 0;
				total_fail_count = 0;
			}
		} while( G > free_energy_cutoff );

		return g;
	}

	/**
	 * Finds a random sequence with folding energy smaller than cutoff
	 */
	static Gene getSequence( const Folder &b, unsigned int length, double free_energy_cutoff)
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		Gene g(length);
		auto_ptr<FoldInfo> fdata;
		double mutation_rate = 1.0/length;
		bool found = false;
		double min_free_energy_for_starting = max(300.0, free_energy_cutoff);

		// find sequence that encodes the desired structure
		do {
			g = Gene::createRandomNoStops( length );
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
			Gene g2 = g;
			bool changed = false;
			do {
				changed = g2.mutate( mutation_rate );
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
					g = Gene::createRandomNoStops( length );
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


#endif //GENE_UTIL_HH
