/*
This file is part of the E.voli project.
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


#include "protein.hh"
#include "translator.hh"
#include "codon.hh"
#include "genetic-code.hh"
#include <sstream>


Codon Codon::transcribe() const {
	Codon res = *this;
	for (int i=0; i<3; i++) {
		if (res[i] == 'T')
			res[i] = 'U';
	}
	return res;
}

Protein::Protein(unsigned int length) : Sequence(length, 'A') {
}

Protein::Protein(unsigned int length, char val) : Sequence(length, val) {
}

Protein Protein::createRandom(unsigned int length) {
	Protein p(length, 'A');
	Protein::iterator pit = p.begin();
	for (; pit != p.end(); pit++) {
		char aa = GeneticCodeUtil::indexToAminoAcidLetter(Random::rint(20)+1);
		*pit = aa;
	}
	return p;
}


CodingDNA::CodingDNA() : Sequence(0, 'A') {
}

CodingDNA::CodingDNA(const CodingDNA& dna) : Sequence(dna) {
}

CodingDNA::CodingDNA(unsigned int length) : Sequence(length, 'A') {
}

CodingDNA::CodingDNA(unsigned int length, char val) : Sequence(length, val) {
}

CodingDNA CodingDNA::createRandom(unsigned int length ) {
	CodingDNA g( length, 'A' );
	string nts("ATGC");

	for (unsigned int j=0; j<length; j++) {
		char nt = nts[Random::rint( 4 )];
		g[j] = nt;
	}
	return g;
}

CodingDNA CodingDNA::createRandomNoStops(unsigned int length ) {
	assert( length % 3 == 0 );
	CodingDNA g( length, 'A' );
	CodingDNA::iterator it = g.begin();
	string nts("ATGC");

	for (unsigned int j=0; j<length/3; j++) {
		do {
			for (unsigned int k=0; k<3; k++) {
				char nt = nts[Random::rint( 4 )];
				g[3*j+k] = nt;
			}
		} while (GeneticCodeUtil::geneticCode(g.transcribe().getCodon(j)) == GeneticCodeUtil::STOP);
	}
	return g;
}

bool CodingDNA::encodesFullLength(void) const {
	bool full_length = (length() % 3)==0;
	CodingRNA rna = transcribe();
	//cout << rna.codonLength() << endl;
	for (int i=0; i<rna.codonLength() && full_length; i++) {
		full_length = (GeneticCodeUtil::geneticCode(rna.getCodon(i)) != GeneticCodeUtil::STOP);
		//cout << i << " " << full_length << endl;
	}
	return full_length;
}

Protein CodingDNA::translate(const Translator& t) const {
	Protein prot(codonLength());
	CodingRNA rna = transcribe();
	//cout << rna << endl;
	bool translation_successful = t.translate(rna, prot);
	//cout << prot << endl;
	assert( translation_successful ); // This function should never be used on genes that don't translate correctly
	return prot;
}

Protein CodingDNA::translate(void) const {
	int len = codonLength();
	Translator t;
	Protein prot(len);
	CodingRNA rna = transcribe();
	//cout << rna << endl;
	bool translation_successful = t.translateErrorFree(rna, prot);
	//cout << prot << endl;
	assert( translation_successful ); // This function should never be used on genes that don't translate correctly
	//cout << this << tab << "translation" << endl;
	return prot;
}

CodingRNA CodingDNA::transcribe() const {
	CodingRNA rna(*this);
	string T("T");
	string U("U");
	int pos = 0;
	pos = rna.find(T, pos);
	while (pos>=0) {
		rna.replace(pos, 1, U);
		pos = rna.find(T, pos);
	}
	return rna;
}

void CodingDNA::setCodon(unsigned int codon_index, const Codon& codon) {
	//cout << codon_index << " " << codon << endl;
	assert(codon_index*3 < length()-2);
	replace(codon_index*3, 3, codon);
}


