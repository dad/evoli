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


#include "genetic-code.hh"
#include "codon.hh"
#include <vector>
#include <ostream>

const char GeneticCodeUtil::STOP = '*';
const char GeneticCodeUtil::INVALID_AA = '@';
const char GeneticCodeUtil::INVALID_NT = '@';
const int GeneticCodeUtil::INVALID_INDEX = -2;
const char* GeneticCodeUtil::AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY";
const char* GeneticCodeUtil::RNA_NUCLEOTIDES = "ACGU";
const char* GeneticCodeUtil::DNA_NUCLEOTIDES = "ACGT";
const char* GeneticCodeUtil::ALL_NUCLEOTIDES = "ACGTU";

// The genetic code
const pair<const Codon,char> GeneticCodeUtil::codonAAPairs[128] =  {
	pair<const Codon,char>(Codon("AAA"), 'K'),
	pair<const Codon,char>(Codon("AAC"), 'N'),
	pair<const Codon,char>(Codon("AAG"), 'K'),
	pair<const Codon,char>(Codon("AAU"), 'N'),
	pair<const Codon,char>(Codon("ACA"), 'T'),
	pair<const Codon,char>(Codon("ACC"), 'T'),
	pair<const Codon,char>(Codon("ACG"), 'T'),
	pair<const Codon,char>(Codon("ACU"), 'T'),
	pair<const Codon,char>(Codon("AGA"), 'R'),
	pair<const Codon,char>(Codon("AGC"), 'S'),
	pair<const Codon,char>(Codon("AGG"), 'R'),
	pair<const Codon,char>(Codon("AGU"), 'S'),
	pair<const Codon,char>(Codon("AUA"), 'I'),
	pair<const Codon,char>(Codon("AUC"), 'I'),
	pair<const Codon,char>(Codon("AUG"), 'M'),
	pair<const Codon,char>(Codon("AUU"), 'I'),
	pair<const Codon,char>(Codon("CAA"), 'Q'),
	pair<const Codon,char>(Codon("CAC"), 'H'),
	pair<const Codon,char>(Codon("CAG"), 'Q'),
	pair<const Codon,char>(Codon("CAU"), 'H'),
	pair<const Codon,char>(Codon("CCA"), 'P'),
	pair<const Codon,char>(Codon("CCC"), 'P'),
	pair<const Codon,char>(Codon("CCG"), 'P'),
	pair<const Codon,char>(Codon("CCU"), 'P'),
	pair<const Codon,char>(Codon("CGA"), 'R'),
	pair<const Codon,char>(Codon("CGC"), 'R'),
	pair<const Codon,char>(Codon("CGG"), 'R'),
	pair<const Codon,char>(Codon("CGU"), 'R'),
	pair<const Codon,char>(Codon("CUA"), 'L'),
	pair<const Codon,char>(Codon("CUC"), 'L'),
	pair<const Codon,char>(Codon("CUG"), 'L'),
	pair<const Codon,char>(Codon("CUU"), 'L'),
	pair<const Codon,char>(Codon("GAA"), 'E'),
	pair<const Codon,char>(Codon("GAC"), 'D'),
	pair<const Codon,char>(Codon("GAG"), 'E'),
	pair<const Codon,char>(Codon("GAU"), 'D'),
	pair<const Codon,char>(Codon("GCA"), 'A'),
	pair<const Codon,char>(Codon("GCC"), 'A'),
	pair<const Codon,char>(Codon("GCG"), 'A'),
	pair<const Codon,char>(Codon("GCU"), 'A'),
	pair<const Codon,char>(Codon("GGA"), 'G'),
	pair<const Codon,char>(Codon("GGC"), 'G'),
	pair<const Codon,char>(Codon("GGG"), 'G'),
	pair<const Codon,char>(Codon("GGU"), 'G'),
	pair<const Codon,char>(Codon("GUA"), 'V'),
	pair<const Codon,char>(Codon("GUC"), 'V'),
	pair<const Codon,char>(Codon("GUG"), 'V'),
	pair<const Codon,char>(Codon("GUU"), 'V'),
	pair<const Codon,char>(Codon("UAA"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("UAC"), 'Y'),
	pair<const Codon,char>(Codon("UAG"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("UAU"), 'Y'),
	pair<const Codon,char>(Codon("UCA"), 'S'),
	pair<const Codon,char>(Codon("UCC"), 'S'),
	pair<const Codon,char>(Codon("UCG"), 'S'),
	pair<const Codon,char>(Codon("UCU"), 'S'),
	pair<const Codon,char>(Codon("UGA"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("UGC"), 'C'),
	pair<const Codon,char>(Codon("UGG"), 'W'),
	pair<const Codon,char>(Codon("UGU"), 'C'),
	pair<const Codon,char>(Codon("UUA"), 'L'),
	pair<const Codon,char>(Codon("UUC"), 'F'),
	pair<const Codon,char>(Codon("UUG"), 'L'),
	pair<const Codon,char>(Codon("UUU"), 'F'),
	// DNA pseudo-codons
	pair<const Codon,char>(Codon("AAA"), 'K'),
	pair<const Codon,char>(Codon("AAC"), 'N'),
	pair<const Codon,char>(Codon("AAG"), 'K'),
	pair<const Codon,char>(Codon("AAT"), 'N'),
	pair<const Codon,char>(Codon("ACA"), 'T'),
	pair<const Codon,char>(Codon("ACC"), 'T'),
	pair<const Codon,char>(Codon("ACG"), 'T'),
	pair<const Codon,char>(Codon("ACT"), 'T'),
	pair<const Codon,char>(Codon("AGA"), 'R'),
	pair<const Codon,char>(Codon("AGC"), 'S'),
	pair<const Codon,char>(Codon("AGG"), 'R'),
	pair<const Codon,char>(Codon("AGT"), 'S'),
	pair<const Codon,char>(Codon("ATA"), 'I'),
	pair<const Codon,char>(Codon("ATC"), 'I'),
	pair<const Codon,char>(Codon("ATG"), 'M'),
	pair<const Codon,char>(Codon("ATT"), 'I'),
	pair<const Codon,char>(Codon("CAA"), 'Q'),
	pair<const Codon,char>(Codon("CAC"), 'H'),
	pair<const Codon,char>(Codon("CAG"), 'Q'),
	pair<const Codon,char>(Codon("CAT"), 'H'),
	pair<const Codon,char>(Codon("CCA"), 'P'),
	pair<const Codon,char>(Codon("CCC"), 'P'),
	pair<const Codon,char>(Codon("CCG"), 'P'),
	pair<const Codon,char>(Codon("CCT"), 'P'),
	pair<const Codon,char>(Codon("CGA"), 'R'),
	pair<const Codon,char>(Codon("CGC"), 'R'),
	pair<const Codon,char>(Codon("CGG"), 'R'),
	pair<const Codon,char>(Codon("CGT"), 'R'),
	pair<const Codon,char>(Codon("CTA"), 'L'),
	pair<const Codon,char>(Codon("CTC"), 'L'),
	pair<const Codon,char>(Codon("CTG"), 'L'),
	pair<const Codon,char>(Codon("CTT"), 'L'),
	pair<const Codon,char>(Codon("GAA"), 'E'),
	pair<const Codon,char>(Codon("GAC"), 'D'),
	pair<const Codon,char>(Codon("GAG"), 'E'),
	pair<const Codon,char>(Codon("GAT"), 'D'),
	pair<const Codon,char>(Codon("GCA"), 'A'),
	pair<const Codon,char>(Codon("GCC"), 'A'),
	pair<const Codon,char>(Codon("GCG"), 'A'),
	pair<const Codon,char>(Codon("GCT"), 'A'),
	pair<const Codon,char>(Codon("GGA"), 'G'),
	pair<const Codon,char>(Codon("GGC"), 'G'),
	pair<const Codon,char>(Codon("GGG"), 'G'),
	pair<const Codon,char>(Codon("GGT"), 'G'),
	pair<const Codon,char>(Codon("GTA"), 'V'),
	pair<const Codon,char>(Codon("GTC"), 'V'),
	pair<const Codon,char>(Codon("GTG"), 'V'),
	pair<const Codon,char>(Codon("GTT"), 'V'),
	pair<const Codon,char>(Codon("TAA"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("TAC"), 'Y'),
	pair<const Codon,char>(Codon("TAG"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("TAT"), 'Y'),
	pair<const Codon,char>(Codon("TCA"), 'S'),
	pair<const Codon,char>(Codon("TCC"), 'S'),
	pair<const Codon,char>(Codon("TCG"), 'S'),
	pair<const Codon,char>(Codon("TCT"), 'S'),
	pair<const Codon,char>(Codon("TGA"), GeneticCodeUtil::STOP),
	pair<const Codon,char>(Codon("TGC"), 'C'),
	pair<const Codon,char>(Codon("TGG"), 'W'),
	pair<const Codon,char>(Codon("TGT"), 'C'),
	pair<const Codon,char>(Codon("TTA"), 'L'),
	pair<const Codon,char>(Codon("TTC"), 'F'),
	pair<const Codon,char>(Codon("TTG"), 'L'),
	pair<const Codon,char>(Codon("TTT"), 'F')
};

const GeneticCodeUtil::CodonMap GeneticCodeUtil::RNACodonToAA(GeneticCodeUtil::codonAAPairs, GeneticCodeUtil::codonAAPairs+64);
const GeneticCodeUtil::CodonMap GeneticCodeUtil::DNACodonToAA(GeneticCodeUtil::codonAAPairs+64, GeneticCodeUtil::codonAAPairs+128);
const GeneticCodeUtil::CodonMap GeneticCodeUtil::codonToAA(GeneticCodeUtil::codonAAPairs, GeneticCodeUtil::codonAAPairs+128);

// Reverse mapping, from AA to vector<Codon>, would be nice, but for some reason doesn't work properly.

const pair<const char, int> GeneticCodeUtil::aaLetterIndices[21] =
{
	pair<char, int>('*', -1),
	pair<char, int>('C', 0),
	pair<char, int>('M', 1),
	pair<char, int>('F', 2),
	pair<char, int>('I', 3),
	pair<char, int>('L', 4),
	pair<char, int>('V', 5),
	pair<char, int>('W', 6),
	pair<char, int>('Y', 7),
	pair<char, int>('A', 8),
	pair<char, int>('G', 9),
	pair<char, int>('T', 10),
	pair<char, int>('S', 11),
	pair<char, int>('Q', 12),
	pair<char, int>('N', 13),
	pair<char, int>('E', 14),
	pair<char, int>('D', 15),
	pair<char, int>('H', 16),
	pair<char, int>('R', 17),
	pair<char, int>('K', 18),
	pair<char, int>('P', 19)

};

const map<const char, int, less<const char> > GeneticCodeUtil::aminoAcidLetterToIndexMap(aaLetterIndices, aaLetterIndices+sizeof(aaLetterIndices)/sizeof(aaLetterIndices[0]));

int GeneticCodeUtil::aminoAcidLetterToIndex(char aa) {
	map<const char, int, less<const char> >::const_iterator it = aminoAcidLetterToIndexMap.find(aa);
	if (it != aminoAcidLetterToIndexMap.end()) {
		return (*it).second;
	}
	else {
		return INVALID_INDEX;
	}
}

char GeneticCodeUtil::indexToAminoAcidLetter(int index) {
	assert( index>=-1 && index <=19 );
	return aaLetterIndices[index+1].first;
}

int GeneticCodeUtil::codonToIndex(Codon codon) {
	// DAD: currently (potentially) PAINFULLY slow.
	for (int i=0; i<128; i++) {
		Codon c = codonAAPairs[i].first;
		if (c == codon) {
			return i;
		}
	}
	// Codon was not found.  Error.
	return -1;
}

Codon GeneticCodeUtil::indexToCodon(int index) {
	assert (index >= 0 && index < 64);
	return codonAAPairs[index].first;
}


const int GeneticCodeUtil::residueToAllCodonsTable[20][7] =
{ // first number is number of codons, then followed by actual codons, and finished off with -1.
	//  0 Cysteine       C
	{ 2, 57, 59, -1, -1, -1, -1 },
	//  1 Methionine     M
	{ 1, 14, -1, -1, -1, -1, -1 },
	//  2 Phenylalanine  F
	{ 2, 61, 63, -1, -1, -1, -1 },
	//  3 Isoleucine     I
	{ 3, 12, 13, 15, -1, -1, -1 },
	//  4 Leucine	L
	{ 6, 28, 29, 30, 31, 60, 62 },
	//  5 Valine	 V
	{ 4, 44, 45, 46, 47, -1, -1 },
	//  6 Tryptophan     W
	{ 1, 58, -1, -1, -1, -1, -1 },
	//  7 Tyrosine       Y
	{ 2, 49, 51, -1, -1, -1, -1 },
	//  8 Alanine	A
	{ 4, 36, 37, 38, 39, -1, -1 },
	//  9 Glycine	G
	{ 4, 40, 41, 42, 43, -1, -1 },
	// 10 Threonine      T
	{ 4, 4, 5, 6, 7, -1, -1 },
	// 11 Serine	 S
	{ 6, 9, 11, 52, 53, 54, 55 },
	// 12 Glutamine      Q
	{ 2, 16, 18, -1, -1, -1, -1 },
	// 13 Asparagine     N
	{ 2, 1, 3, -1, -1, -1, -1 },
	// 14 Glutamic Acid  E
	{ 2, 32, 34, -1, -1, -1, -1 },
	// 15 Aspartic Acid  D
	{ 2, 33, 35, -1, -1, -1, -1 },
	// 16 Histidine      H
	{ 2, 17, 19, -1, -1, -1, -1 },
	// 17 Arginine       R
	{ 6, 8, 10, 24, 25, 26, 27 },
	// 18 Lysine	 K
	{ 2, 0, 2, -1, -1, -1, -1 },
	// 19 Proline	P
	{ 4, 20, 21, 22, 23, -1, -1 }
};


const int GeneticCodeUtil::singleSubstsTranslErrors[64][20] =
{
//		 AAA
	{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0 },
//		 AAC
	{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0 },
//		 AAG
	{ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0 },
//		 AAU
	{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0 },
//		 ACA
	{ 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1 },
//		 ACC
	{ 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 },
//		 ACG
	{ 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1 },
//		 ACU
	{ 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1 },
//		 AGA
	{ 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0 },
//		 AGC
	{ 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0 },
//		 AGG
	{ 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0 },
//		 AGU
	{ 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0 },
//		 AUA
	{ 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0 },
//		 AUC
	{ 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0 },
//		 AUG
	{ 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0 },
//		 AUU
	{ 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0 },
//		 CAA
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1 },
//		 CAC
	{ 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1 },
//		 CAG
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1 },
//		 CAU
	{ 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1 },
//		 CCA
	{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0 },
//		 CCC
	{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0 },
//		 CCG
	{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0 },
//		 CCU
	{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0 },
//		 CGA
	{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
//		 CGC
	{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1 },
//		 CGG
	{ 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 },
//		 CGU
	{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1 },
//		 CUA
	{ 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1 },
//		 CUC
	{ 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1 },
//		 CUG
	{ 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1 },
//		 CUU
	{ 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1 },
//		 GAA
	{ 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
//		 GAC
	{ 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0 },
//		 GAG
	{ 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 },
//		 GAU
	{ 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0 },
//		 GCA
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1 },
//		 GCC
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1 },
//		 GCG
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1 },
//		 GCU
	{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1 },
//		 GGA
	{ 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
//		 GGC
	{ 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0 },
//		 GGG
	{ 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0 },
//		 GGU
	{ 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0 },
//		 GUA
	{ 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
//		 GUC
	{ 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
//		 GUG
	{ 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 },
//		 GUU
	{ 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 },
//		 UAA
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UAC
	{ 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0 },
//		 UAG
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UAU
	{ 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0 },
//		 UCA
	{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
//		 UCC
	{ 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
//		 UCG
	{ 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
//		 UCU
	{ 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1 },
//		 UGA
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UGC
	{ 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0 },
//		 UGG
	{ 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0 },
//		 UGU
	{ 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0 },
//		 UUA
	{ 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UUC
	{ 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UUG
	{ 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
//		 UUU
	{ 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 }
};

char GeneticCodeUtil::geneticCode(const Codon codon) {
	CodonMap::const_iterator it = codonToAA.find(codon);
	if (it != codonToAA.end())
		return (*it).second;
	else
		return INVALID_AA;
}

void GeneticCodeUtil::printResidue( ostream &s, Codon codon )
{
	s << geneticCode(codon);
}

char GeneticCodeUtil::residueLetter( Codon codon )
{
	return geneticCode(codon);
}

/*
void GeneticCodeUtil::printGeneticCode( ostream &s )
{
	s << "Genetic code, sorted according to frequency:" << endl;
	vector<int> counts;
	counts.resize(20);

	for ( int i=0; i<64; i++ ) {
		int p = geneticCode[i];
		if ( p>=0 )
			counts[p]+=1;
	}
	for ( int k=7; k>0; k-- ) {
		for ( int p=0; p<20; p++ ) {
			if ( counts[p] == k ) {
				s << residues[p] << ": ";
				for ( int i=0; i<64; i++ ) {
					int q = geneticCode[i];
					if ( q == p ) {
						s << " ";
						CodonUtil::printCodon( s, i );
					}
				}
				s << endl;
			}
		}
	}
	s << "Stop:";
	for ( int i=0; i<64; i++ )
		if ( geneticCode[i] < 0 ) {
			s << " ";
			CodonUtil::printCodon( s, i );
		}
	s << endl;
}
*/

int synHelper(Codon& codon, Codon& varcodon, unsigned int ntpos) {
	const char* nts = GeneticCodeUtil::RNA_NUCLEOTIDES;
	int numSyn = 0;
	for ( int i=0; i<4; i++ ) {
		// Don't consider non-mutations
		if (codon[ntpos] == nts[i])
			continue;
		varcodon[ntpos] = nts[i];
		if ( GeneticCodeUtil::geneticCode(codon) == GeneticCodeUtil::geneticCode(varcodon))
			numSyn += 1;
		// restore nucleotide
		varcodon[ntpos] = codon[ntpos];
	}
	return numSyn;
}

double GeneticCodeUtil::calcSynonymousSites( Codon codon, int sites )
{
	int numSyn = 0;
	Codon varcodon(codon);

	if ( ( sites & 4 ) != 0 )
		numSyn += synHelper(codon, varcodon, 0);
				
	if ( ( sites & 2 ) != 0 )
		numSyn += synHelper(codon, varcodon, 1);
				
	if ( ( sites & 1 ) != 0 )
		numSyn += synHelper(codon, varcodon, 2);

	return (double) numSyn / 3.;
}


double GeneticCodeUtil::calcSynMutationOpportunity( Codon codon, double rho )
{
	double S = 0;
	double wti = 2*rho; // weight for transitions, have to multiply by two because there are two times as many transversions
	double wtv = 1; // weight for transversions
	const char* nts = GeneticCodeUtil::RNA_NUCLEOTIDES;

	Codon varcodon(codon);

	char res = geneticCode(codon);
	if (res == INVALID_AA) {
		// Return an invalid number.  Need to specify semantics of this function.
		return -1.0;
	}

	// For each codon position...
	for ( unsigned int ntpos=0; ntpos<3; ntpos++) {
		// For each alternative nucleotide...
		for ( int i=0; i<4; i++ ) {
			char new_nt = nts[i];
			// Don't consider "mutations" to the current nucleotide
			if (new_nt == codon[ntpos])
				continue;
			// Replace the nucleotide
			varcodon[ntpos] = new_nt;
			// If synonymous...
			if ( res == geneticCode(varcodon)) {
				// Accumulate weights depending on whether mutation is transition or transversion
				if ( CodonUtil::isTransition( codon[ntpos], varcodon[ntpos] ) )
					S += wti;
				else
					S += wtv;
			}
			// Restore the old nucleotide before continuing
			varcodon[ntpos] = codon[ntpos];
		}
	}
	return S/(2.+2*rho);
}

double GeneticCodeUtil::calcNonsynMutationOpportunity( Codon codon, double rho )
{
	double N = 0;
	double wti = 2*rho; // weight for transitions, have to multiply by two because there are two times as many transversions
	double wtv = 1; // weight for transversions
	string nts("AUGC"); // Alternative nucleotides

	Codon varcodon(codon);

	char res = geneticCode(codon);
	if (res == INVALID_AA) {
		// Return an invalid number.  Need to specify semantics of this function.
		return -1.0;
	}

	// For each codon position...
	for ( unsigned int ntpos=0; ntpos<3; ntpos++) {
		// For each alternative nucleotide...
		for ( int i=0; i<4; i++ ) {
			char new_nt = nts[i];
			// Don't consider "mutations" to the current nucleotide
			if (new_nt == codon[ntpos])
				continue;
			// Replace the nucleotide
			varcodon[ntpos] = new_nt;
			// If nonsynonymous...
			if ( res != geneticCode(varcodon)) {
				// Accumulate weights depending on whether mutation is transition or transversion
				if ( CodonUtil::isTransition( codon[ntpos], varcodon[ntpos] ) )
					N += wti;
				else
					N += wtv;
			}
			// Restore the old nucleotide before continuing
			varcodon[ntpos] = codon[ntpos];
		}
	}
	return N/(2.+2*rho);
}

bool GeneticCodeUtil::m_setup=false;
hash_map<const char*, double, hash<const char*> > GeneticCodeUtil::m_dnLookup;
hash_map<const char*, double, hash<const char*> > GeneticCodeUtil::m_dsLookup;
double GeneticCodeUtil::m_dnTable[64][64];
double GeneticCodeUtil::m_dsTable[64][64];

void GeneticCodeUtil::calcDnDs( double &dn, double &ds, Codon codon1, Codon codon2 ) {
	// This implementation is not thread-safe.
	if ( !m_setup )	{
		m_setup = true;
		for ( int i=0; i<128; i++) {
			for ( int j=0; j<128; j++ )	{
				Codon c1 = codonAAPairs[i].first;
				Codon c2 = codonAAPairs[j].first;
				calcDnDsPrivate( dn, ds, c1, c2 );
				const char* key = c1.append(c2).c_str();
				m_dnLookup[key] = dn;
				m_dsLookup[key] = ds;
				//m_dnTable[i][j] = dn;
				//m_dsTable[i][j] = ds;
//				 CodonUtil::printCodon( cout, i );
//				 cout << " ";
//				 CodonUtil::printCodon( cout, j );
//				 cout << " " << dn << " " << ds << endl;
			}
		}
	}

	//dn = m_dnTable[codon1][codon2];
	//ds = m_dsTable[codon1][codon2];
	const char* key = codon1.append(codon2).c_str();
	hash_map<const char*, double, hash<const char*> >::const_iterator it = m_dnLookup.find(key);
	if (it != m_dnLookup.end()) {
		dn = (*it).second;
	}
	it = m_dsLookup.find(key);
	if (it != m_dsLookup.end()) {
		ds = (*it).second;
	}
}


void GeneticCodeUtil::calcDnDsPrivate( double &dn, double &ds, Codon codon1, Codon codon2 ) {
	dn = 0; // the number of nonsynonymous substitutions between the two codons
	ds = 0; // the number of synonymous substitutions between the two codons

	//Codon c1 = codon1, c2 = codon2; // holds the two codons as arrays of single letters
	int subst_positions[3]; // holds the positions at which the two codons differ
	int differences = 0; // number of positions at which the two codons differ
	Codon tmpc; // temporary codon (for paths)
	char res1, res2, tmpres1, tmpres2;

	res1 = geneticCode(codon1);
	res2 = geneticCode(codon2);

	// first, count total number of differences. Record them in subst_positions
	if ( codon1[0] != codon2[0] )
	{
		subst_positions[differences] = 0;
		differences += 1;
	}
	if ( codon1[1] != codon2[1] )
	{
		subst_positions[differences] = 1;
		differences += 1;
	}
	if ( codon1[2] != codon2[2] )
	{
		subst_positions[differences] = 2;
		differences += 1;
	}

	if ( differences == 0 ) // no differences, we are done
		return;

	if ( differences == 1 )
	{ // this case is simple
		if ( res1 == res2 )
			ds = 1;
		else
			dn = 1;
		return;
	}

	if ( differences == 2 )
	{ // there are only two possible cases here, still simple
		// path 1
		tmpc = codon1; // copy the initial codon
		// subst 1
		tmpc[subst_positions[0]] = codon2[subst_positions[0]];
		tmpres1 = geneticCode(tmpc); //[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
		if ( res1 == tmpres1 )
			ds += 1;
		else
			dn += 1;

		// subst 2 leads to codon 2
		if ( tmpres1 == res2 )
			ds += 1;
		else
			dn += 1;

		// path 2
		tmpc = codon1; // copy the initial codon

		// subst 1
		tmpc[subst_positions[1]] = codon2[subst_positions[1]];
		tmpres1 = geneticCode(tmpc); 
		if ( res1 == tmpres1 )
			ds += 1;
		else
			dn += 1;
		// subst 2 leads to codon 2
		if ( tmpres1 == res2 )
			ds += 1;
		else
			dn += 1;
		ds /= 2.;
		dn /= 2.;

		return;
	}

	// finally, if all letters are different, we have 6 different cases
	// these are the 6 orderings of 3 mutations (3!).
	// now, we don't need the array subst_positions any more, because
	// we have to do all possible substitutions anyway

	// Set up all possible paths for first two mutations
	vector<pair<int,int> > mutation_paths(6);
	mutation_paths.push_back(pair<int,int>(0,1));
	mutation_paths.push_back(pair<int,int>(1,0));
	mutation_paths.push_back(pair<int,int>(1,2));
	mutation_paths.push_back(pair<int,int>(2,1));
	mutation_paths.push_back(pair<int,int>(0,2));
	mutation_paths.push_back(pair<int,int>(2,0));
	
	vector<pair<int,int> >::const_iterator mutit = mutation_paths.begin();
	for (; mutit != mutation_paths.end(); mutit++) {
		pair<int,int> mut = *mutit;
		tmpc = codon1; // copy the initial codon
		// subst 1
		tmpc[mut.first] = codon2[mut.first];
		tmpres1 = geneticCode(tmpc);
		if ( res1 == tmpres1 )
			ds += 1;
		else
			dn += 1;
		// subst 2
		tmpc[mut.second] = codon2[mut.second];
		tmpres2 = geneticCode(tmpc);
		if ( tmpres2 == tmpres1 )
			ds += 1;
		else
			dn += 1;
		// subst 3 always leads to codon 2
		if ( tmpres2 == res2 )
			ds += 1;
		else
			dn += 1;
	}

	ds /= 6.;
	dn /= 6.;

	return;
}


pair<double, double> GeneticCodeUtil::calcDnDsWeightedPrivate( Codon codon1, Codon codon2, double rho )
{
	#ifndef NDEBUG
	// check that function is called correctly
	//int cala, calb, calc, cbla, cblb, cblc;
	//CodonUtil::codonToLetters( cala, calb, calc, codon1 );
	//CodonUtil::codonToLetters( cbla, cblb, cblc, codon2 );
	// we expect exactly one base difference
	assert( ( codon1[0] != codon2[0] && codon1[1] == codon2[1] && codon1[2] == codon2[2] ) ||
			( codon1[0] == codon2[0] && codon1[1] != codon2[1] && codon1[2] == codon2[2] ) ||
			( codon1[0] == codon2[0] && codon1[1] == codon2[1] && codon1[2] != codon2[2] ) );
	#endif

	double S = calcSynMutationOpportunity( codon1, rho );
	double N = calcNonsynMutationOpportunity( codon1, rho );
	double dn = 0., ds = 0.;

	char res1 = geneticCode( codon1 );
	char res2 = geneticCode( codon2 );
	// we set distances involving stop codons to zero, in order to not count them
	if ( res1 == GeneticCodeUtil::STOP || res2 == GeneticCodeUtil::STOP )
		return pair<double, double>( 0., 0. );

	// is the mutation synonymous or not?
	if ( res1 == res2 )
		ds = 1./S;
	else
		dn = 1./N;
	return pair<double, double>( dn, ds );
}


pair<double, double> GeneticCodeUtil::calcDnDsWeighted( Codon codon1, Codon codon2, double rho )
{
	double dn = 0; // the number of nonsynonymous substitutions between the two codons
	double ds = 0; // the number of synonymous substitutions between the two codons

	//int c1[3], c2[3]; // holds the two codons as arrays of single letters
	int subst_positions[3]; // holds the positions at which the two codons differ
	int differences = 0; // number of positions at which the two codons differ
	Codon tmpc, tmpc2; // temporary codon (for paths)
	int res1, res2, tmpres1, tmpres2;

	res1 = geneticCode(codon1);
	res2 = geneticCode(codon2);

	//CodonUtil::codonToLetters( c1[0], c1[1], c1[2], codon1 );
	//CodonUtil::codonToLetters( c2[0], c2[1], c2[2], codon2 );

	// first, count total number of differences. Record them in subst_positions
	if ( codon1[0] != codon2[0] )
	{
		subst_positions[differences] = 0;
		differences += 1;
	}
	if ( codon1[1] != codon2[1] )
	{
		subst_positions[differences] = 1;
		differences += 1;
	}
	if ( codon1[2] != codon2[2] )
	{
		subst_positions[differences] = 2;
		differences += 1;
	}

	if ( differences == 0 ) // no differences, we are done
		return pair<double, double>( 0., 0. );

	if ( differences == 1 )
	{ // this case is simple
		return calcDnDsWeightedPrivate( codon1, codon2, rho );
	}


	pair<double, double> p;
	if ( differences == 2 )
	{ // there are only two possible cases here, still simple
		// path 1
		tmpc = codon1; // initialize codon.
		// subst 1
		tmpc[subst_positions[0]] = codon2[subst_positions[0]];
		p = calcDnDsWeightedPrivate( codon1, tmpc, rho );
		dn = p.first;
		ds = p.second;

		// subst 2 leads to codon 2
		p = calcDnDsWeightedPrivate( tmpc, codon2, rho );
		dn += p.first;
		ds += p.second;

		// path 2
		tmpc = codon1; // copy the initial codon
		// subst 1
		tmpc[subst_positions[1]] = codon2[subst_positions[1]];
		p = calcDnDsWeightedPrivate( codon1, tmpc, rho );
		dn += p.first;
		ds += p.second;

		// subst 2 leads to codon 2
		p = calcDnDsWeightedPrivate( tmpc, codon2, rho );
		dn += p.first;
		ds += p.second;

		ds /= 2.;
		dn /= 2.;

		return pair<double, double>( dn, ds );
	}

	// finally, if all letters are different, we have 6 different cases
	// these are the 6 orderings of 3 mutations (3!).
	// now, we don't need the array subst_positions any more, because
	// we have to do all possible substitutions anyway

	// Set up all possible paths for first two mutations
	vector<pair<int,int> > mutation_paths(6);
	mutation_paths.push_back(pair<int,int>(0,1));
	mutation_paths.push_back(pair<int,int>(1,0));
	mutation_paths.push_back(pair<int,int>(1,2));
	mutation_paths.push_back(pair<int,int>(2,1));
	mutation_paths.push_back(pair<int,int>(0,2));
	mutation_paths.push_back(pair<int,int>(2,0));
	
	vector<pair<int,int> >::const_iterator mutit = mutation_paths.begin();
	for (; mutit != mutation_paths.end(); mutit++) {
		pair<int,int> mut = *mutit;
		tmpc = codon1; // copy the initial codon
		// subst 1
		tmpc[mut.first] = codon2[mut.first];
		p = calcDnDsWeightedPrivate( codon1, tmpc, rho );
		dn = p.first;
		ds = p.second;
		// subst 2
		tmpc2 = tmpc;
		tmpc2[mut.second] = codon2[mut.second];
		p = calcDnDsWeightedPrivate( tmpc, tmpc2, rho );
		dn += p.first;
		ds += p.second;
		// subst 3 always leads to codon 2
		p = calcDnDsWeightedPrivate( tmpc2, codon2, rho );
		dn += p.first;
		ds += p.second;
	}

	ds /= 6.;
	dn /= 6.;

	return pair<double, double>( dn, ds );
}
