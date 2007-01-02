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


#include "translator.hh"

#include "random.hh"
#include "tools.hh"
#include "genetic-code.hh"
#include "protein.hh"


Translator::Translator( double mutation_prob)
		: m_mutation_prob( mutation_prob )
{}

Translator::Translator()
		: m_mutation_prob( 0 )
{}

bool Translator::translateErrorFree( const CodingRNA &g, Protein& residue_sequence ) const {
	bool no_stop = true;

	int i = 0;
	for (int i=0; i<g.codonLength(); i++) {
		Codon ci = g.getCodon(i);
		char residue = GeneticCodeUtil::geneticCode(ci);
		no_stop = (residue != GeneticCodeUtil::STOP );
		residue_sequence[i] = residue;
		//cout << ci << " " << i << " " << residue << endl;
	}
	return no_stop;
}

bool Translator::translate( const CodingRNA &g, Protein& residue_sequence ) const {
	if ( m_mutation_prob == 0 )
		return translateErrorFree( g, residue_sequence );

	bool no_stop = true;

	//int i = 0;
	//for ( CodingRNA::const_iterator it = g.begin(); it != g.end() && no_stop; it++ )	{
	for ( unsigned int i=0; i<g.codonLength(); i++) {
		Codon ci = g.getCodon(i);
		char residue = GeneticCodeUtil::geneticCode(ci);
		if ( Random::runif() < m_mutation_prob )
			residue = ( residue + Random::rint( 20 ) ) % 20;
		no_stop = (residue != GeneticCodeUtil::STOP );
		residue_sequence[i++] = residue;
	}
	return no_stop;
}

int Translator::translateWeighted( const CodingRNA &g, Protein& residue_sequence, const vector<vector<pair<double, char> > >& weights,
									       const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated)
{
	double mut_weight_total = 0.0;
	for ( CodingRNA::const_iterator it = g.begin(); it != g.end(); it++ )	{
		mut_weight_total += (1.0 + prefCodons[*it]*(nonPrefCodonPenalty-1));
	}

	truncated = false;
	int numErrors = 0;
	for ( unsigned int i=0; i<g.codonLength() && !truncated; i++) {
		Codon ci = g.getCodon(i);
		int codon_index = GeneticCodeUtil::codonToIndex(ci);
		char residue = GeneticCodeUtil::geneticCode(ci);
		residue_sequence[i] = residue;

		if ( residue == GeneticCodeUtil::STOP ) {
			truncated = true;
		}
		else {
			double threshold = m_mutation_prob*(1.0 + prefCodons[codon_index]*(nonPrefCodonPenalty-1))/(mut_weight_total/g.codonLength());
			double rand = Random::runif();
			if ( rand < threshold ) {
				// Weight the outcomes of a missense substitution.
				rand = Random::runif();
				double targ = 0.0;
				// The first member of the pair is a cumulative probability; the second is
				// the residue resulting from the error.
				for (unsigned int j=0; j<weights[codon_index].size() && (rand > targ); j++) {
					pair<double, char> p = weights[codon_index][j];
					residue_sequence[i] = p.second;
					targ = p.first;
				}
				if (residue_sequence[i] == GeneticCodeUtil::STOP) { // truncation error
					truncated = true;
				}
				if (residue_sequence[i] != residue) {
					numErrors++;
				}
			}
		}
	}
	return numErrors;
}

int Translator::translateRelativeWeighted( const CodingRNA &g, Protein& residue_sequence, const double error_weight,
										   const vector<vector<pair<double, char> > >& weights, const double* prefCodons, 
										   const double nonPrefCodonPenalty, bool& truncated)
{
	truncated = false;
	int numErrors = 0;
	// For an average gene encoding a folded protein, the sum of the
	// site_weights over all codons should be equal to error_weight,
	// and thus the per-codon probability of a translation error
	// (possibly synonymous) will be given by m_mutation_prob.  CodingRNAs
	// with higher site_weight sums are more likely to be
	// mistranslated.
	//
	// When codons are translated with equal accuracy (no codon
	// preference), error_weight should simply be the length of the
	// gene, each codon's site_weight should be 1.0, and the per-codon
	// probability of error is exactly given by m_mutation_prob.
	double avg_error_per_site_weight = error_weight/g.codonLength();
	for ( int i=0; i<g.codonLength() && !truncated; i++) {
		Codon ci = g.getCodon(i);
		int codon_index = GeneticCodeUtil::codonToIndex(ci);
		char residue = GeneticCodeUtil::geneticCode(ci);
		
		residue_sequence[i] = residue;

		// Compute the site weight, accounting for codon preference.
		// Each site weight is proportional to the probability that an
		// error occurs at this codon relative to other codons in the gene.
		double site_weight = (1.0 + prefCodons[codon_index]*(nonPrefCodonPenalty-1));

		// With probability threshold_prob, make a translation error (possibly synonymous).
		double threshold_prob = m_mutation_prob * site_weight / avg_error_per_site_weight;
		double rand = Random::runif();
		if ( rand < threshold_prob ) {
			// We've made an error.  Now determine what it is.
			rand = Random::runif();
			double targ = 0.0;
			// The first member of the pair is a cumulative probability; the second is
			// the residue resulting from the error.
			// Find the event corresponding to rand.
			for (unsigned int j=0; j<weights[codon_index].size() && (rand > targ); j++) {
				pair<double, char> p = weights[codon_index][j];
				residue_sequence[i] = p.second;
				targ = p.first;
			}
			// Tabulate the results
			// Only count an error if the polymerized amino acid differs from
			// the natively encoded residue.
			if (residue_sequence[i] != residue) {
				numErrors++;
			}
		}
		// Check for truncation.  May not be an error!
		if ( residue_sequence[i] == GeneticCodeUtil::STOP ) {
			truncated = true;
		}
	}
	return numErrors;
}
