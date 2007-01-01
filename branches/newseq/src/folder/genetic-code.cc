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

// The genetic code
const pair<const Codon,char> GeneticCodeUtil::codon_aa_pairs[64] = 
{
	pair<const Codon,char>(Codon("AAA"), 'K'),
	pair<const Codon,char>(Codon("AAC"), 'N'),
	pair<const Codon,char>(Codon("AAG"), 'K'),
	pair<const Codon,char>(Codon("AAU"), 'N'),
	pair<const Codon,char>(Codon("ACA"), 'T'),
	pair<const Codon,char>(Codon("ACC"), 'T'),
	pair<const Codon,char>(Codon("ACG"), 'T'),
	pair<const Codon,char>(Codon("ACU"), 'T'),
	pair<const Codon,char>(Codon("AGA"), 'R'),
	pair<const Codon,char>(Codon("AGC"), 'R'),
	pair<const Codon,char>(Codon("AGG"), 'R'),
	pair<const Codon,char>(Codon("AGU"), 'R'),
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
	pair<const Codon,char>(Codon("CGA"), 'S'),
	pair<const Codon,char>(Codon("CGC"), 'R'),
	pair<const Codon,char>(Codon("CGG"), 'S'),
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
	pair<const Codon,char>(Codon("UUU"), 'F')
};

hash_map<const Codon, char, hash_codon > GeneticCodeUtil::RNACodonToAA(GeneticCodeUtil::codon_aa_pairs, GeneticCodeUtil::codon_aa_pairs+64);

/*typedef hash_map<char, vector<const Codon>, hash<char> > aa_codon_map;
hash_map<char, vector<const Codon>, hash<char> > GeneticCodeUtil::AAToRNACodon();
static {
	for (int i=0; i<64; i++) {
		pair<const Codon, char> p = codon_aa_pairs[i];
		char key = p.first;
		aa_codon_map::iterator it = GeneticCodeUtil::AAToRNACodon.find(key);
		if (it == GeneticCodeUtil::AAToRNACodon.end()) {
			// Insert a new vector
			GeneticCodeUtil::AAToRNACodon[key] = vector<const Codon>(1,p.second);
		}
		else {
			// Append next codon
			GeneticCodeUtil::AAToRNACodon[key].push_back(p.second);
		}
	}
	};*/

pair<const char, int> letterResidues[21] =
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

map<const char, int, less<const char> > GeneticCodeUtil::letter_to_residue_map(letterResidues, letterResidues+sizeof(letterResidues)/sizeof(letterResidues[0]));

// this is the mapping from integer to residue that we use
const char *GeneticCodeUtil::residues[20] =
{
	"CYS", //  0 Cysteine       C
	"MET", //  1 Methionine     M
	"PHE", //  2 Phenylalanine  F
	"ILE", //  3 Isoleucine     I
	"LEU", //  4 Leucine	L
	"VAL", //  5 Valine	 V
	"TRP", //  6 Tryptophan     W
	"TYR", //  7 Tyrosine       Y
	"ALA", //  8 Alanine	A
	"GLY", //  9 Glycine	G
	"THR", // 10 Threonine      T
	"SER", // 11 Serine	 S
	"GLN", // 12 Glutamine      Q
	"ASN", // 13 Asparagine     N
	"GLU", // 14 Glutamic Acid  E
	"ASP", // 15 Aspartic Acid  D
	"HIS", // 16 Histidine      H
	"ARG", // 17 Arginine       R
	"LYS", // 18 Lysine	 K
	"PRO"  // 19 Proline	P
};

const char *GeneticCodeUtil::residueLetters[21] =
{
	"*", // -1 STOP			*
	"C", //  0 Cysteine       C
	"M", //  1 Methionine     M
	"F", //  2 Phenylalanine  F
	"I", //  3 Isoleucine     I
	"L", //  4 Leucine	L
	"V", //  5 Valine	 V
	"W", //  6 Tryptophan     W
	"Y", //  7 Tyrosine       Y
	"A", //  8 Alanine	A
	"G", //  9 Glycine	G
	"T", // 10 Threonine      T
	"S", // 11 Serine	 S
	"Q", // 12 Glutamine      Q
	"N", // 13 Asparagine     N
	"E", // 14 Glutamic Acid  E
	"D", // 15 Aspartic Acid  D
	"H", // 16 Histidine      H
	"R", // 17 Arginine       R
	"K", // 18 Lysine	 K
	"P"  // 19 Proline	P
};

// The genetic code
const int GeneticCodeUtil::geneticCode[64] =
{
	18, // AAA -> LYS 0
	13, // AAC -> ASN 1
	18, // AAG -> LYS 2
	13, // AAU -> ASN 3
	10, // ACA -> THR 4
	10, // ACC -> THR 5
	10, // ACG -> THR 6
	10, // ACU -> THR 7
	17, // AGA -> ARG 8
	11, // AGC -> SER 9
	17, // AGG -> ARG 10
	11, // AGU -> SER 11
	3,  // AUA -> ILE 12
	3,  // AUC -> ILE 13
	1,  // AUG -> MET 14
	3,  // AUU -> ILE 15
	12, // CAA -> GLN 16
	16, // CAC -> HIS 17
	12, // CAG -> GLN 18
	16, // CAU -> HIS 19
	19, // CCA -> PRO 20
	19, // CCC -> PRO 21
	19, // CCG -> PRO 22
	19, // CCU -> PRO 23
	17, // CGA -> ARG 24
	17, // CGC -> ARG 25
	17, // CGG -> ARG 26
	17, // CGU -> ARG 27
	4,  // CUA -> LEU 28
	4,  // CUC -> LEU 29
	4,  // CUG -> LEU 30
	4,  // CUU -> LEU 31
	14, // GAA -> GLU 32
	15, // GAC -> ASP 33
	14, // GAG -> GLU 34
	15, // GAU -> ASP 35
	8,  // GCA -> ALA 36
	8,  // GCC -> ALA 37
	8,  // GCG -> ALA 38
	8,  // GCU -> ALA 39
	9,  // GGA -> GLY 40
	9,  // GGC -> GLY 41
	9,  // GGG -> GLY 42
	9,  // GGU -> GLY 43
	5,  // GUA -> VAL 44
	5,  // GUC -> VAL 45
	5,  // GUG -> VAL 46
	5,  // GUU -> VAL 47
	-1, // UAA -> STOP48
	7,  // UAC -> TYR 49
	-1, // UAG -> STOP50
	7,  // UAU -> TYR 51
	11, // UCA -> SER 52
	11, // UCC -> SER 53
	11, // UCG -> SER 54
	11, // UCU -> SER 55
	-1, // UGA -> STOP56
	0,  // UGC -> CYS 57
	6,  // UGG -> TRP 58
	0,  // UGU -> CYS 59
	4,  // UUA -> LEU 60
	2,  // UUC -> PHE 61
	4,  // UUG -> LEU 62
	2   // UUU -> PHE 63
};



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


const int GeneticCodeUtil::residueToCodonTable[21] =
{
	48, // UAA -> STOP 49 * 0
	57,  // UGC -> CYS 58 C 1
	14,  // AUG -> MET 15 M 2
	61,  // UUC -> PHE 62 F 3
	12,  // AUA -> ILE 13 I 4
	31,  // CUU -> LEU 32 L 5
	47,  // GUU -> VAL 48 V 6
	58,  // UGG -> TRP 59 W 7
	49,  // UAC -> TYR 50 Y 8
	39,  // GCU -> ALA 40 A 9
	40,  // GGA -> GLY 41 G 10
	4,   // ACA -> THR 5  T 11
	9,   // AGC -> SER 10 S 12
	16,  // CAA -> GLN 17 Q 13
	1,   // AAC -> ASN 2  N 14
	32,  // GAA -> GLU 33 E 15
	33,  // GAC -> ASP 34 D 16
	17,  // CAC -> HIS 18 H 17
	8,   // AGA -> ARG 9  R 18
	0,   // AAA -> LYS 1  K 19
	23,  // CCU -> PRO 24 P 20
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

void GeneticCodeUtil::printResidue( ostream &s, int codon )
{
	assert( codon < 64 );

	int r = geneticCode[codon];
	if ( r < 0 )
		s << "STOP";
	else
		s << residues[r];
}

char GeneticCodeUtil::residueLetter( int codon )
{
	assert( codon < 64 );

	int r = geneticCode[codon];
	return *residueLetters[r+1];
}


void GeneticCodeUtil::printGeneticCode( ostream &s )
{
	s << "Genetic code, sorted according to frequency:" << endl;
	vector<int> counts;
	counts.resize(20);

	for ( int i=0; i<64; i++ )
	{
		int p = geneticCode[i];
		if ( p>=0 )
			counts[p]+=1;
	}
		for ( int k=7; k>0; k-- )
	{
		for ( int p=0; p<20; p++ )
		{
			if ( counts[p] == k )
			{
				s << residues[p] << ": ";
				for ( int i=0; i<64; i++ )
				{
					int q = geneticCode[i];
					if ( q == p )
					{
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
		if ( geneticCode[i] < 0 )
		{
			s << " ";
			CodonUtil::printCodon( s, i );
		}
	s << endl;
}


double GeneticCodeUtil::calcSynonymousSites( int codon, int sites )
{
	int l1, l2, l3, ltmp;
	int numSyn = 0;

	CodonUtil::codonToLetters( l1, l2, l3, codon );
	if ( ( sites & 4 ) != 0 )
	{
		for ( int i=1; i<4; i++ )
		{
			ltmp = (l1+i) & 3;
			if ( geneticCode[codon] == geneticCode[CodonUtil::lettersToCodon( ltmp, l2, l3 )] )
				numSyn += 1;
		}
	}

	if ( ( sites & 2 ) != 0 )
	{
		for ( int i=1; i<4; i++ )
		{
			ltmp = (l2+i) & 3;
			if ( geneticCode[codon] == geneticCode[CodonUtil::lettersToCodon( l1, ltmp, l3 )] )
				numSyn += 1;
		}
	}
	if ( ( sites & 1 ) != 0 )
	{
		for ( int i=1; i<4; i++ )
		{
			ltmp = (l3+i) & 3;
			if ( geneticCode[codon] == geneticCode[CodonUtil::lettersToCodon( l1, l2, ltmp )] )
				numSyn += 1;
		}
	}

	return (double) numSyn / 3.;
}


double GeneticCodeUtil::calcSynMutationOpportunity( int codon, double rho )
{
	int l1, l2, l3, ltmp;
	double S = 0;
	double wti = 2*rho; // weight for transitions, have to multiply by two because there are two times as many transversions
	double wtv = 1; // weight for transversions

	CodonUtil::codonToLetters( l1, l2, l3, codon );
	int res = geneticCode[codon];

	for ( int i=1; i<4; i++ )
	{
		ltmp = (l1+i) & 3;
		if ( res == geneticCode[CodonUtil::lettersToCodon( ltmp, l2, l3 )] )
		{
			if ( CodonUtil::isTransition( l1, ltmp ) )
				S += wti;
			else
				S += wtv;
		}
	}
	for ( int i=1; i<4; i++ )
	{
		ltmp = (l2+i) & 3;
		if ( res == geneticCode[CodonUtil::lettersToCodon( l1, ltmp, l3 )] )
		{
			if ( CodonUtil::isTransition( l2, ltmp ) )
				S += wti;
			else
				S += wtv;
		}
	}
	for ( int i=1; i<4; i++ )
	{
		ltmp = (l3+i) & 3;
		if ( res == geneticCode[CodonUtil::lettersToCodon( l1, l2, ltmp )] )
		{
			if ( CodonUtil::isTransition( l3, ltmp ) )
				S += wti;
			else
				S += wtv;
		}
	}
	return S/(2.+2*rho);
}

double GeneticCodeUtil::calcNonsynMutationOpportunity( int codon, double rho )
{
	int l1, l2, l3, ltmp;
	double N = 0;
	double wti = 2*rho; // weight for transitions, have to multiply by two because there are two times as many transversions
	double wtv = 1.; // weight for transversions

	CodonUtil::codonToLetters( l1, l2, l3, codon );
	int res = geneticCode[codon];
	int tmpres;

	for ( int i=1; i<4; i++ )
	{
		ltmp = (l1+i) & 3;
		tmpres = geneticCode[CodonUtil::lettersToCodon( ltmp, l2, l3 )];
		if ( tmpres >= 0 && res != tmpres )
		{
			if ( CodonUtil::isTransition( l1, ltmp ) )
				N += wti; // transitions occur rho times more frequent
			else
				N += wtv;
		}
	}
	for ( int i=1; i<4; i++ )
	{
		ltmp = (l2+i) & 3;
		tmpres = geneticCode[CodonUtil::lettersToCodon( l1, ltmp, l3 )];
		if ( tmpres >= 0 && res != tmpres )
		{
			if ( CodonUtil::isTransition( l2, ltmp ) )
				N += wti; // transitions occur rho times more frequent
			else
				N += wtv;
		}
	}
	for ( int i=1; i<4; i++ )
	{
		ltmp = (l3+i) & 3;
		tmpres = geneticCode[CodonUtil::lettersToCodon( l1, l2, ltmp )];
		if ( tmpres >= 0 && res != tmpres )
		{
			if ( CodonUtil::isTransition( l3, ltmp ) )
				N += wti; // transitions occur rho times more frequent
			else
				N += wtv;
		}
	}
	return N/(2.+2*rho);
}


bool GeneticCodeUtil::m_setup=false;
double GeneticCodeUtil::m_dnTable[64][64];
double GeneticCodeUtil::m_dsTable[64][64];


void GeneticCodeUtil::calcDnDs( double &dn, double &ds, int codon1, int codon2 )
{
	if ( !m_setup )
	{
		m_setup = true;
		for ( int i=0; i<64; i++ )
			for ( int j=0; j<64; j++ )
			{
				calcDnDsPrivate( dn, ds, i, j );
				m_dnTable[i][j] = dn;
				m_dsTable[i][j] = ds;
//				 CodonUtil::printCodon( cout, i );
//				 cout << " ";
//				 CodonUtil::printCodon( cout, j );
//				 cout << " " << dn << " " << ds << endl;
			}
	}

	dn = m_dnTable[codon1][codon2];
	ds = m_dsTable[codon1][codon2];
}


void GeneticCodeUtil::calcDnDsPrivate( double &dn, double &ds, int codon1, int codon2 )
{
	dn = 0; // the number of nonsynonymous substitutions between the two codons
	ds = 0; // the number of synonymous substitutions between the two codons

	int c1[3], c2[3]; // holds the two codons as arrays of single letters
	int subst_positions[3]; // holds the positions at which the two codons differ
	int differences = 0; // number of positions at which the two codons differ
	int tmpc[3]; // temporary codon (for paths)
	int res1, res2, tmpres1, tmpres2;

	res1 = geneticCode[codon1];
	res2 = geneticCode[codon2];

	CodonUtil::codonToLetters( c1[0], c1[1], c1[2], codon1 );
	CodonUtil::codonToLetters( c2[0], c2[1], c2[2], codon2 );

	// first, count total number of differences. Record them in subst_positions
	if ( c1[0] != c2[0] )
	{
		subst_positions[differences] = 0;
		differences += 1;
	}
	if ( c1[1] != c2[1] )
	{
		subst_positions[differences] = 1;
		differences += 1;
	}
	if ( c1[2] != c2[2] )
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
		else dn = 1;
		return;
	}

	if ( differences == 2 )
	{ // there are only two possible cases here, still simple
		// path 1
		for ( int i=0; i<3; i++ ) // copy the initial codon
			tmpc[i] = c1[i];
		// subst 1
		tmpc[subst_positions[0]] = c2[subst_positions[0]];
		tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
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
		for ( int i=0; i<3; i++ ) // copy the initial codon
			tmpc[i] = c1[i];

		// subst 1
		tmpc[subst_positions[1]] = c2[subst_positions[1]];
		tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
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
	// now, we don't need the array subst_positions any more, because
	// we have to do all possible substitutions anyway
	// path 1
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[0] = c2[0];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[1] = c2[1];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	// path 2
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[0] = c2[0];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[2] = c2[2];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	// path 3
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[1] = c2[1];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[0] = c2[0];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	// path 4
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[1] = c2[1];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[2] = c2[2];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	// path 5
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[2] = c2[2];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[0] = c2[0];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	// path 6
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[2] = c2[2];
	tmpres1 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( res1 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 2
	tmpc[1] = c2[1];
	tmpres2 = geneticCode[ CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] ) ];
	if ( tmpres2 == tmpres1 )
		ds += 1;
	else
		dn += 1;
	// subst 3 leads to codon 2
	if ( tmpres2 == res2 )
		ds += 1;
	else
		dn += 1;

	ds /= 6.;
	dn /= 6.;

	return;
}


pair<double, double> GeneticCodeUtil::calcDnDsWeightedPrivate( int codon1, int codon2, double rho )
{
	#ifndef NDEBUG
	// check that function is called correctly
	int cala, calb, calc, cbla, cblb, cblc;
	CodonUtil::codonToLetters( cala, calb, calc, codon1 );
	CodonUtil::codonToLetters( cbla, cblb, cblc, codon2 );
	// we expect exatly one base difference
	assert( ( cala != cbla && calb == cblb && calc == cblc ) ||
			( cala == cbla && calb != cblb && calc == cblc ) ||
			( cala == cbla && calb == cblb && calc != cblc ) );
	#endif

	double S = calcSynMutationOpportunity( codon1, rho );
	double N = calcNonsynMutationOpportunity( codon1, rho );
	double dn = 0., ds = 0.;

	int res1 = geneticCode[ codon1 ];
	int res2 = geneticCode[ codon2 ];
	// we set distances involving stop codons to zero, in order to not count them
	if ( res1<0 || res2<0 )
		return pair<double, double>( 0., 0. );

	// is the mutation synonymous or not?
	if ( res1 == res2 )
		ds = 1./S;
	else
		dn = 1./N;
	return pair<double, double>( dn, ds );
}


pair<double, double> GeneticCodeUtil::calcDnDsWeighted( int codon1, int codon2, double rho )
{
	double dn = 0; // the number of nonsynonymous substitutions between the two codons
	double ds = 0; // the number of synonymous substitutions between the two codons

	int c1[3], c2[3]; // holds the two codons as arrays of single letters
	int subst_positions[3]; // holds the positions at which the two codons differ
	int differences = 0; // number of positions at which the two codons differ
	int tmpc[3]; // temporary codon (for paths)
	int tmpcodon1, tmpcodon2; // temporary codon
	int res1, res2, tmpres1, tmpres2;

	res1 = geneticCode[codon1];
	res2 = geneticCode[codon2];

	CodonUtil::codonToLetters( c1[0], c1[1], c1[2], codon1 );
	CodonUtil::codonToLetters( c2[0], c2[1], c2[2], codon2 );

	// first, count total number of differences. Record them in subst_positions
	if ( c1[0] != c2[0] )
	{
		subst_positions[differences] = 0;
		differences += 1;
	}
	if ( c1[1] != c2[1] )
	{
		subst_positions[differences] = 1;
		differences += 1;
	}
	if ( c1[2] != c2[2] )
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
		for ( int i=0; i<3; i++ ) // copy the initial codon
			tmpc[i] = c1[i];
		// subst 1
		tmpc[subst_positions[0]] = c2[subst_positions[0]];
		tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
		p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
		dn = p.first;
		ds = p.second;

		// subst 2 leads to codon 2
		p = calcDnDsWeightedPrivate( tmpcodon1, codon2, rho );
		dn += p.first;
		ds += p.second;

		// path 2
		for ( int i=0; i<3; i++ ) // copy the initial codon
			tmpc[i] = c1[i];

		// subst 1
		tmpc[subst_positions[1]] = c2[subst_positions[1]];
		tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
		p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
		dn += p.first;
		ds += p.second;

		// subst 2 leads to codon 2
		p = calcDnDsWeightedPrivate( tmpcodon1, codon2, rho );
		dn += p.first;
		ds += p.second;

		ds /= 2.;
		dn /= 2.;

		return pair<double, double>( dn, ds );
	}

	// finally, if all letters are different, we have 6 different cases
	// now, we don't need the array subst_positions any more, because
	// we have to do all possible substitutions anyway
	// path 1
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[0] = c2[0];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn = p.first;
	ds = p.second;

	// subst 2
	tmpc[1] = c2[1];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	// path 2
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[0] = c2[0];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn += p.first;
	ds += p.second;

	// subst 2
	tmpc[2] = c2[2];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	// path 3
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[1] = c2[1];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn += p.first;
	ds += p.second;

	// subst 2
	tmpc[2] = c2[2];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	// path 4
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[1] = c2[1];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn += p.first;
	ds += p.second;

	// subst 2
	tmpc[0] = c2[0];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	// path 5
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[2] = c2[2];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn += p.first;
	ds += p.second;

	// subst 2
	tmpc[0] = c2[0];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	// path 6
	for ( int i=0; i<3; i++ ) // copy the initial codon
		tmpc[i] = c1[i];
	// subst 1
	tmpc[2] = c2[2];
	tmpcodon1 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( codon1, tmpcodon1, rho );
	dn += p.first;
	ds += p.second;

	// subst 2
	tmpc[1] = c2[1];
	tmpcodon2 = CodonUtil::lettersToCodon( tmpc[0], tmpc[1], tmpc[2] );
	p = calcDnDsWeightedPrivate( tmpcodon1, tmpcodon2, rho );
	dn += p.first;
	ds += p.second;

	// subst 3 leads to codon 2
	p = calcDnDsWeightedPrivate( tmpcodon2, codon2, rho );
	dn += p.first;
	ds += p.second;

	ds /= 6.;
	dn /= 6.;

	return pair<double, double>( dn, ds );
}
