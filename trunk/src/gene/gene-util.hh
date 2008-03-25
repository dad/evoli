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


#ifndef GENE_UTIL_HH
#define GENE_UTIL_HH

#include "random.hh"
#include "codon.hh"
#include "genetic-code.hh"
#include "translator.hh"
#include "protein.hh"
#include "mutator.hh"

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
	static double calcTotalSites( const CodingDNA& g ) {
		return (double)g.length();
	}

	/**
	 * Calculates the number of synonymous sites in the genotype.
	 **/
	static double calcSynonymousSites( const CodingDNA & g )
	{
		CodingRNA rna = g.transcribe();
		double s = 0;
		for (unsigned int i=0; i<rna.codonLength(); i++) {
			s += GeneticCodeUtil::calcSynonymousSites( rna.getCodon(i) );
		}
		return s;
	}

	/**
	 * Calculates the number of synonymous and nonsynonymous sites in surface and core.
	 **/
	static void calcSNSitesSurfaceCore( double &NSurf, double &NCore, double &SSurf, double &SCore,
					    const CodingDNA &g, const vector<int> &surface )
	{
		unsigned int e = g.codonLength();
		assert( e == surface.size() );
		NSurf = NCore = SSurf = SCore = 0;
		double tmp_N, tmp_S;

		for ( unsigned int i=0; i<e; i++ )
		{
			tmp_S = GeneticCodeUtil::calcSynonymousSites( g.getCodon(i) );
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
	static pair<double,double> calcDnDs( const CodingDNA &g1, const CodingDNA &g2 )
	{
		assert( g1.length() == g2.length() );
		int e = g1.codonLength();
		double dn =0, ds =0;
		for ( int i=0; i<e; i++ )
		{
			pair<double,double> dnds = GeneticCodeUtil::calcDnDs( g1.getCodon(i), g2.getCodon(i) );
			//			 CodonUtil::printCodon( cout, g1[i] );
			//			 cout << " ";
			//			 CodonUtil::printCodon( cout, g2[i] );
			//			 cout << " " << tmp_dn << " " << tmp_ds << endl;
			dn += dnds.first;
			ds += dnds.second;
		}
		return pair<double,double>(dn,ds);
	}

	/**
	 * Same as calcDnDs, but distinguishes between surface/core residues.
	 * The vector "surface" is expected to be of the same length as the protein,
	 * and must contain a 1 at surface positions and a 0 at core positions.
	 **/
	static void calcDnDsSurfaceCore( double &dnSurf, double &dnCore, double &dsSurf, double &dsCore,
					 const CodingDNA &g1, const CodingDNA &g2, const vector<int> &surface )
	{
		assert( g1.codonLength() == g2.codonLength() );
		assert( g1.codonLength() == surface.size() );
		int e = g1.codonLength();
		dnSurf = dnCore = dsSurf = dsCore = 0;

		for ( int i=0; i<e; i++ )
		{
			pair<double,double> dnds = GeneticCodeUtil::calcDnDs( g1.getCodon(i), g2.getCodon(i) );
			if ( surface[i] )
			{
				dnSurf += dnds.first;
				dsSurf += dnds.second;
			}
			else
			{
				dnCore += dnds.first;
				dsCore += dnds.second;
			}
		}
	}


	/**
	 * Calculates the fraction of optimal codons from the set of optimal codons.
	 **/
	static double calcFop( const CodingDNA &g, const vector<bool>& is_optimal ) {
		int count = 0, tot=0;
		char methionine = 'M';
		char tryptophan = 'W';

		for (unsigned int i=0; i<g.codonLength(); i++) {
			Codon c = g.getCodon(i);
			unsigned int index = GeneticCodeUtil::codonToIndex(c);
			assert( index >= 0 && index < is_optimal.size() );
			//assert(is_optimal.size() == 64);
			char aa = GeneticCodeUtil::geneticCode(c);
			if ( !(aa == methionine || aa == tryptophan)) {
				tot += 1;
				if (is_optimal[index]) {
					count += 1;
				}
			}
		}
		return count / (double)tot;
	}

	/**
	 * As calcFop, but considers surface and core regions separately
	 **/
	static void calcFopSurfaceCore( double &Fop_surf, double &Fop_core, const CodingDNA &g, const double *codon_costs, const vector<int> &surface )
	{
		assert( g.codonLength() == surface.size() );

		Fop_surf = Fop_core = 0;

		if ( codon_costs )
		{
			int sopt = 0, copt = 0, scount = 0, ccount = 0;

			for ( uint i=0; i<g.codonLength(); i++ ) {
				Codon ci = g.getCodon(i);
				int index = GeneticCodeUtil::codonToIndex(ci);
				if ( surface[i] )
				{
					scount += 1;
					if ( codon_costs[index] == 0 )
						sopt += 1;
				}
				else
				{
					ccount += 1;
					if ( codon_costs[index] == 0 )
						copt += 1;
				}
			}
			Fop_surf = (double) sopt / (double) scount;
			Fop_core = (double) copt / (double) ccount;
		}

	}

	/**
	 * Produces a sequence encoding the same amino acid sequence as gene g, but with randomly chosen codons.
	 * @param g A gene sequence.
	 * @return A gene sequence encoding the same protein sequence as g, with randomly selected codons.
	 * Returns the original sequence if it doesn't translate properly.
	 **/
	static CodingDNA randomizeCodons( const CodingDNA &g ) {
		CodingDNA result(g);
		if ( g.encodesFullLength() ) {
			Protein p = g.translate();
			result = reverseTranslate(p);
		}
		return result;
	}

	/**
	 * Turn a protein sequence into a DNA sequence.
	 *
	 * @param prot The protein to reverse-translate.
	 * @return A CodingDNA sequence which, when transcribed and translated, will be equivalent to prot.
	 */
	static CodingDNA reverseTranslate(const Protein& prot) {
		// DAD: strange...size of the GeneticCodeUtil hash_map RNACodonToAA is 139, not 64!  Figure this out.
		hash_map<const Codon, char, hash_codon> gc(GeneticCodeUtil::codonAAPairs + 64, GeneticCodeUtil::codonAAPairs+128);
		// DAD: probably only want to do this once!
		hash_map<char, vector<Codon>, hash<char> > AAToDNACodons;
		hash_map<const Codon, char, hash_codon>::iterator map_it = gc.begin();
		//int i = 0;
		for (; map_it != gc.end(); map_it++) {
			pair<const Codon, char> p = *map_it;
			char key = p.second;
			Codon codon = p.first;
			//cout << ++i << " " << key << " " << codon << endl;

			hash_map<char, vector<Codon>, hash<char> >::iterator found_it = AAToDNACodons.find(key);
			if (found_it == AAToDNACodons.end()) {
				// Insert a new vector
				AAToDNACodons[key] = vector<Codon>();
				AAToDNACodons[key].push_back(codon);
			}
			else {
				// Append next codon
				AAToDNACodons[key].push_back(codon);
				//cout << key << " " << codon << " " << AAToDNACodons[key].size() << endl;
			}
		}
		CodingDNA g(prot.length()*3);
		//cout << prot << endl;
		Protein::const_iterator it = prot.begin();
		int which_codon = 0;
		for (; it != prot.end(); it++, which_codon++) {
			char aa = *it;
			//cout << which_codon << " " << aa << endl;
			//continue;
			vector<Codon>& v = AAToDNACodons[aa];
			int rint = Random::rint(v.size());
			//cout << which_codon << tab << aa << tab << rint << tab << v.size() << tab << flush;
			Codon codon = v[rint];
			//cout << codon << tab << gc[codon] << endl;
			g.setCodon(which_codon, codon);
		}

		return g;
	}

};


#endif //GENE_UTIL_HH
