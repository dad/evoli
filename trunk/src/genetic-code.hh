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
	 * Calculates the number of changes in nonsynonymous (dn) and synonymous
	 * (ds) substitutions between codon1 and codon2. If the two codons differ
	 * by more than one base, all possible paths are weighted equally.
	 * Stop codons are not excluded from the analysis.
	 **/
	static void calcDnDs( double &dn, double &ds, int codon1, int codon2 );


};

#endif






