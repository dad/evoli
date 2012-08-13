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


#ifndef GENETIC_CODE_HH
#define GENETIC_CODE_HH

#include "protein.hh"
#include "codon.hh"
#include <iostream>
// uncomment next line for older versions of gcc
#include <unordered_map>
// comment out next two lines for older versions of gcc
//#include <ext/unordered_map>
//using namespace __gnu_cxx;

using namespace std;

// Helper class to hash Codons
struct hash_codon {
	size_t operator()(const Codon c) const {
		hash<const char*> h;
		return h(c.c_str());
	}
};

class GeneticCodeUtil
{
public:
	typedef unordered_map<const Codon, char, hash_codon > CodonMap;
	typedef unordered_map<const char*, double, hash<const char*> > StringDoubleMap;
private:
	// All this is needed for the function calc DnDs.
	static double m_dnTable[][64];
	static double m_dsTable[][64];
	static void calcDnDsPrivate( double &dn, double &ds, Codon codon1, Codon codon2 );
	static pair<double, double> calcDnDsWeightedPrivate( Codon codon1, Codon codon2, double rho );

	/**
	 * Pairing of amino-acid letters with integer indices.
	 **/
	static const pair<const char, int> aaLetterIndices[21];

	/**
	 * Mapping from amino-acid letters to integer indices
	 **/
	static const unordered_map<const char, int, less<const char> > aminoAcidLetterToIndexMap;

	/**
	 * Mapping from proper RNA codons to amino acids.
	 **/
	static const CodonMap RNACodonToAA;

	/**
	 * Mapping from DNA pseudo-codons to amino acids.
	 **/
	static const CodonMap DNACodonToAA;

	/**
	 * Mapping from all codons (RNA and DNA) to amino acids.
	 **/
	static const CodonMap codonToAA;

public:
	static const char STOP;
	static const char INVALID_AA;
	static const char INVALID_NT;
	static const int INVALID_INDEX;
	static const char* AMINO_ACIDS;
	static const char* RNA_NUCLEOTIDES;
	static const char* DNA_NUCLEOTIDES;
	static const char* ALL_NUCLEOTIDES;

	static const pair<const Codon,char> codonAAPairs[128];

	/**
	 * Retrieve the amino acid corresponding to this codon.
	 **/
	static char geneticCode(const Codon);

	/**
	 * Retrieve an index for this amino-acid letter code
	 **/
	static int aminoAcidLetterToIndex(char aa);

	/**
	 * Retrieve an amino-acid letter code for this index
	 **/
	static char indexToAminoAcidLetter(int index);

	/**
	 * Retrieve an index for this codon.
	 **/
	static int codonToIndex(Codon codon);

	/**
	 * Retrieve the codon for this index.
	 **/
	static Codon indexToCodon(int index);

  /**
   * Is this nucleotide change a transition?
   * @param nt1 The first nucleotide
   * @param nt2 The second nucleotide
   * @return Whether the change from nt1 to nt2 is a transition.
   **/
  static bool isTransition( char nt1, char nt2 );

	/**
	 * Mapping from amino acids to codons
	 **/
	//static unordered_map<char, vector<const Codon>, hash<char> > AAToRNACodon;

	/**
	 * Reverse lookup for genetic code. Contains all codons for a given residue. The first number
	 * in the array for each residue is the number of codons.
	 **/
	static const int residueToAllCodonsTable[20][7];

	/**
	 * Table to figure out whether a particular amino acid can occur
	 * during mistranslation of a codon if only singe-base mismatches are allowed.
	 **/
	static const int singleSubstsTranslErrors[64][20];

	/**
	 * Translates the given codon into a residue, and outputs the result
	 * to the given ostream.
	 **/
	static void printResidue( ostream &s, Codon codon );

	/**
	 * Translates given codon into a residue, and returns the residue one-character code.
	 **/
	static char residueLetter( Codon codon );

	/**
	 * Prints a table of the genetic code.
	 **/
	//static void printGeneticCode( ostream &s );

	/**
	 * Calculates the number of synonymous sites in the codon. The variable
	 * sites can be used in order to do this analysis on select bases:
	 *  sites = 7 (default) implies all bases are considered
	 *  sites = 1 implies only the last base is considered
	 *  sites = 2 implies only the middle base is considered
	 *  sites = 4 implies only the first base is considered
	 *  etc.
	 **/
	static double calcSynonymousSites( Codon codon, int sites = 7 );

	/**
	Calculates the number of synonymous sites according to the mutational
	opportunity definition, i.e., corrected for a potential transition to
	transversion bias.
	\param codon The codon to analyze.
	\param rho Ratio of transitions to transversions. Note that no
	transition to transversion bias corresponds to rho=.5.
	**/
	static double calcSynMutationOpportunity( Codon codon, double rho );

	/**
	Calculates the number of nonsynonymous sites according to the
	mutational opportunity definition, i.e., corrected for a potential
	transition to transversion bias. Mutations to stop codons are not counted.
	\param codon The codon to analyze.
	\param rho Ratio of transitions to transversions. Note that no
	transition to transversion bias corresponds to rho=.5.
	**/
	static double calcNonsynMutationOpportunity( Codon codon, double rho );


	/**
	 * Calculates the number of changes in nonsynonymous (dn) and synonymous
	 * (ds) substitutions between codon1 and codon2. If the two codons differ
	 * by more than one base, all possible paths are weighted equally.
	 * Stop codons are not excluded from the analysis.
	 \warning This function has not been fully debugged, and may contain mistakes (though none are known). Use at your own peril.
	 \param codon1 The originating codon.
	 \param codon2 The destination codon.
	 \return A pair containing dn and ds (in this order).
	 **/
	static pair<double, double> calcDnDs( Codon codon1, Codon codon2 );

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
	static pair<double, double> calcDnDsWeighted( Codon codon1, Codon codon2, double rho );

};

#endif






