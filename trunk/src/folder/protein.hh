#ifndef PROTEIN_HH
#define PROTEIN_HH

#include <vector>
#include <string>
#include "sequence.hh"
#include "genetic-code.hh"

using namespace std;

typedef pair<int, int> Contact;

/** \brief A class holding a protein sequence.
*/
class Protein : public Sequence {
private:
	Protein();
protected:
public:
	Protein(const int length);
	Protein(const Protein& p);
	Protein(const string& s);
	~Protein() {}

	/**
	Calculates the Hamming distance to another protein. Assumes that 
	the two proteins are of equal lengths.
	@return The hamming distance to the protein given as argument.
	*/
	int distance(const Protein& p) const;

	/**
	Converts the protein sequence into a string useful for output etc.
	@return A string object containing the protein sequence, using the single-letter amino-acid alphabet.
	*/
	string toString() const;
};

class Translator;

/** \brief A class holding a protein-coding DNA sequence.
Attention: The current implementation is based on codons, not individual
nucleotides. As a consequence, the iterators inherited from \ref Sequence iterate over codons, not nucleotides. In the future, we should add a nucleotide_iterator or something similar to iterate over individual nucleotides.
*/
class Gene : public Sequence {
private:
protected:
public:
	Gene();
	Gene(const int length);
	Gene(const Gene& g);
	Gene(const string& s);
	~Gene() {}

	/**
	@return The length of the DNA sequence in nucleotides.
	*/
	virtual uint length(void) const { return 3*m_sequence.size(); }

	/**
	@return The length of the DNA sequence in codons.
	*/
	virtual uint codonLength(void) const { return m_sequence.size(); }

	/**
	Translates the DNA sequence into a protein sequence using the given \ref Translator object. Useful for non-standard genetic codes, error-prone translation, etc.
	@param t \ref Translator object specifying the details of translation.
	@return The protein sequence.
	*/
	Protein translate(const Translator& t) const;

	/**
	Translates the DNA sequence into a protein sequence using the standard genetic code.
	@return The protein sequence.
	*/
	Protein translate(void) const;

	/**
	Creates a random gene. The resulting gene may contain stop codons.
	@param length Length of the desired gene, in nucleotides.
	@return The random gene.
	*/
	static Gene createRandom(const int length);

	/**
	As @ref createRandom, but now stop codons are avoided.
	@param length Length of the desired gene, in nucleotides.
	@return The random gene.
	*/
	static Gene createRandomNoStops(const int length);
	bool mutate(const double prob);

	/**
	Tests whether the DNA sequence contains any stop codons.
	@return True if there are no stop codons (the entire gene is coding sequence), False otherwise.
	*/
	bool encodesFullLength(void) const;

	/**
	Accessor function, returns the nucleotide at a given position.
	@param index The position in the DNA sequence; counting starts from zero.
	@return The resulting nucleotide, i.e., 'A', 'C', 'G', or 'U'.
	*/
	char getBase(const uint index) const;

	/**
	Converts the DNA sequence into a string useful for output etc.
	@return A string object containing the DNA sequence, using the single-letter nucleotide abbreviations.
	*/
	string toString() const;

};


/**
 * Printing a Gene.
 **/
inline ostream & operator<<( ostream &s, const Gene &g ) {
	s << g.toString();
	return s;
}

/**
 * Printing a Protein.
 **/
inline ostream & operator<<( ostream &s, const Protein &p ) {
	for ( Protein::const_iterator it=p.begin(); it != p.end(); it++) {
		s << GeneticCodeUtil::residueLetters[*it+1]; // << "  ";
	}
	return s;
}

inline istream & operator>>( istream &s, Gene & g ) {
	string str;
	s >> str;
	g = Gene(str);
	return s;
}


#endif //PROTEIN_HH
