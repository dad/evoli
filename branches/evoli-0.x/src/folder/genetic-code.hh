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


#ifndef GENETIC_CODE_HH
#define GENETIC_CODE_HH

#include <iostream>
#include <map>

using namespace std;

class GeneticCodeUtil
{
private:
	// All this is needed for the function calc DnDs.
	static bool m_setup;
	static double m_dnTable[64][64];
	static double m_dsTable[64][64];
	static void calcDnDsPrivate( double &dn, double &ds, int codon1, int codon2 );
	static pair<double, double> calcDnDsWeightedPrivate( int codon1, int codon2, double rho );
public:
	/**
	 * The mapping from integer to residue
	 **/
	static const char* residues[20];
	/**
	 * The mapping from integer to residue single-letter codes
	 **/
	static const char* residueLetters[21];

	/**
	 * The genetic code
	 * we represent each base by two binary digits:
	 *   A: 00, C: 01, G: 10, U: 11
	 * Then, each codon is an integer between 0 and 63. For example,
	 * GUG = 101110 = 46. At position 46 in the array, we therefore need a 5
	 * (because GUG codes for the residue VAL, which has index 5 according to the
	 * mapping defined above.
	 **/
	static const int geneticCode[64];

	/**
	 * The mapping from characters to residues
	 **/
    static map<const char, int, less<const char> > letter_to_residue_map;

	/**
	 * Reverse lookup for genetic code. Contains all codons for a given residue. The first number
	 * in the array for each residue is the number of codons.
	 **/
	static const int residueToAllCodonsTable[20][7];

	/**
	 * Table to convert from residues to a single (arbitrarily chosen) codon.
	 */
	static const int residueToCodonTable[21];

	/**
	 * Table to figure out whether a particular amino acid can occur
	 * during mistranslation of a codon if only singe-base mismatches are allowed.
	 **/
	static const int singleSubstsTranslErrors[64][20];

	/**
	 * Translates the given codon into a residue, and outputs the result
	 * to the given ostream.
	 **/
	static void printResidue( ostream &s, int codon );

	/**
	 * Translates given codon into a residue, and returns the residue one-character code.
	 **/
	static char residueLetter( int codon );

	/**
	 * Prints a table of the genetic code.
	 **/
	static void printGeneticCode( ostream &s );


	/**
	 * Calculates the number of synonymous sites in the codon. The variable
	 * sites can be used in order to do this analysis on select bases:
	 *  sites = 7 (default) implies all bases are considered
	 *  sites = 1 implies only the last base is considered
	 *  sites = 2 implies only the middle base is considered
	 *  sites = 4 implies only the first base is considered
	 *  etc.
	 **/
	static double calcSynonymousSites( int codon, int sites = 7 );

	/**
	Calculates the number of synonymous sites according to the mutational
	opportunity definition, i.e., corrected for a potential transition to
	transversion bias.
	\param codon The codon to analyze.
	\param rho Ratio of transitions to transversions. Note that no
	transition to transversion bias corresponds to rho=.5. 
	**/
	static double calcSynMutationOpportunity( int codon, double rho );

	/**
	Calculates the number of nonsynonymous sites according to the
	mutational opportunity definition, i.e., corrected for a potential
	transition to transversion bias. Mutations to stop codons are not counted.
	\param codon The codon to analyze.
	\param rho Ratio of transitions to transversions. Note that no
	transition to transversion bias corresponds to rho=.5. 
	**/
	static double calcNonsynMutationOpportunity( int codon, double rho );


	/**
	 * Calculates the number of changes in nonsynonymous (dn) and synonymous
	 * (ds) substitutions between codon1 and codon2. If the two codons differ
	 * by more than one base, all possible paths are weighted equally.
	 * Stop codons are not excluded from the analysis.
	 **/
	static void calcDnDs( double &dn, double &ds, int codon1, int codon2 );

	/**
	Calculates the number of changes in nonsynonymous (dn) and synonymous
	(ds) substitutions per mutational opportunity between codon1 and codon2.
	If the two codons differ by more than one base, all possible paths are weighted
	according to the transition/transversion ratio rho. Stop codons are excluded
	from the analysis. Note that mutations are assumed to lead from codon1 to
	codon2. If you interchange the two, the results will in general change.
	\warning This function has not been fully debugged, and may contain mistakes (though none are known). Use at your own peril.
	\param codon1 The originating codon.
	\param codon2 The destination codon.
	\param rho The transition/transversion ratio (a ratio of .5 means both are
	equally frequent).
	\return A pair containing dn and ds (in this order).
	 **/
	static pair<double, double> calcDnDsWeighted( int codon1, int codon2, double rho );

};

#endif






