/*
This file is part of the E.voli project.  Copyright (C) 2004, 2005,
2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>, Allan Drummond
<drummond@alumni.princeton.edu>

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


#ifndef PROTEIN_HH
#define PROTEIN_HH

#include <vector>
#include <string>
#include <cassert>
#include "sequence.hh"
#include "codon.hh"

using namespace std;

typedef pair<int, int> Contact;

/** \brief A class holding a protein sequence.
*/
class Protein : public Sequence {
protected:
public:
	Protein(unsigned int length);
	Protein(unsigned int length, char val);
	Protein(const string& s) : Sequence(s) {}
	~Protein() {}

	/**
	Converts the protein sequence into a string useful for output etc.
	@return A string object containing the protein sequence, using the single-letter amino-acid alphabet.
	*/
	string toString() const { return *this; }

	/**
	Creates a random protein.
	@param length Length of the desired protein, in residues
	@return The random protein.
	*/
	static Protein createRandom(unsigned int length);

};

class Translator;
class CodingDNA;
typedef CodingDNA CodingRNA;

/** \brief A class holding a protein-coding DNA sequence.
*/
class CodingDNA : public Sequence {
private:
protected:
public:
	CodingDNA();
	CodingDNA(const CodingDNA& g);
	CodingDNA(const string& s) : Sequence(s) {}
	CodingDNA(unsigned int length);
	CodingDNA(unsigned int length, char val);
	virtual ~CodingDNA() {}

	/**
	@return The length of the DNA sequence in codons.
	*/
	virtual uint codonLength(void) const { return size()/3; }

	/**
	Translates the DNA sequence into a protein sequence using the
	given \ref Translator object. Useful for non-standard genetic
	codes, error-prone translation, etc.  
	
	\warning The function doesn't test whether the protein is actually
	translatable, and will return an invalid protein if translation
	fails. Always check with \ref encodesFullLength() first whether
	translation is safe.

	@param t \ref Translator object specifying the details of translation.
	@return The protein sequence.
	*/
	Protein translate(const Translator& t) const;

	/**
	Translates the DNA sequence into a protein sequence using the
	standard genetic code.  

	\warning The function doesn't test whether the protein is actually
	translatable, and will return an invalid protein if translation
	fails. Always check with \ref encodesFullLength() first whether
	translation is safe.

	@return	The protein sequence.
	*/
	Protein translate(void) const;

	/**
	Creates a random coding sequence. The resulting gene may contain stop codons.
	@param length Length of the desired gene, in nucleotides.
	@return The random gene.
	*/
	static CodingDNA createRandom(unsigned int length);

	/**
	As @ref createRandom, but now stop codons are avoided.
	@param length Length of the desired gene, in nucleotides.
	@return The random gene.
	*/
	static CodingDNA createRandomNoStops(unsigned int length);

	/**
	Tests whether the DNA sequence contains any stop codons.
	@return True if there are no stop codons (the entire gene is coding sequence), False otherwise.
	*/
	bool encodesFullLength(void) const;

	/**
	 * Retrieves the specified codon.
	 * @return A sequence
	 **/
	Codon getCodon(unsigned int codon_index) const {
		assert(codon_index*3 < length()-2);
		return Codon(*this, 3*codon_index);
	}

	/**
	 * Replaces the specified codon.
	 * @return A sequence
	 **/
	void setCodon(unsigned int codon_index, const Codon& codon);

	/**
	 * Transcribes DNA into RNA.
	 * @return A CodingRNA, with instances of T in the original CodingDNA replaced by U.
	 **/
	CodingRNA transcribe() const;
};

// DAD: Just for now
typedef CodingDNA Gene;


#endif //PROTEIN_HH
