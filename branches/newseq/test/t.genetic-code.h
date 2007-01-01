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


#ifndef _T_GENETIC_CODE_H__
#define _T_GENETIC_CODE_H__
#include "cutee.h"
#include "protein.hh"
#include "genetic-code.hh"
#include <iostream>
using namespace std;

struct TEST_CLASS( genetic_code_basic ) {
	void TEST_FUNCTION( aa_from_codon ) {
		Codon c("AUG");
		char aa = GeneticCodeUtil::RNACodonToAA[c];
		TEST_ASSERT( aa == 'M' );
		// Systematically test all codons
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AGA")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAA")] == '*' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAG")] == '*' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AAA")] == 'K' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AAG")] == 'K' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AAU")] == 'N' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AAC")] == 'N' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("ACA")] == 'T' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("ACG")] == 'T' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("ACU")] == 'T' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("ACC")] == 'T' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AGA")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AGG")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AGU")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AGC")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AUA")] == 'I' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AUG")] == 'M' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AUU")] == 'I' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("AUC")] == 'I' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CAA")] == 'Q' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CAG")] == 'Q' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CAU")] == 'H' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CAC")] == 'H' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CCA")] == 'P' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CCG")] == 'P' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CCU")] == 'P' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CCC")] == 'P' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CGA")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CGG")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CGU")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CGC")] == 'R' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CUA")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CUG")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CUU")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("CUC")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GAA")] == 'E' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GAG")] == 'E' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GAU")] == 'D' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GAC")] == 'D' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GCA")] == 'A' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GCG")] == 'A' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GCU")] == 'A' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GCC")] == 'A' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GGA")] == 'G' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GGG")] == 'G' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GGU")] == 'G' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GGC")] == 'G' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GUA")] == 'V' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GUG")] == 'V' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GUU")] == 'V' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("GUC")] == 'V' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAA")] == '*' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAG")] == '*' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAU")] == 'Y' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UAC")] == 'Y' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UCA")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UCG")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UCU")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UCC")] == 'S' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UGA")] == '*' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UGG")] == 'W' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UGU")] == 'C' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UGC")] == 'C' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UUA")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UUG")] == 'L' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UUU")] == 'F' );
		TEST_ASSERT( GeneticCodeUtil::RNACodonToAA[Codon("UUC")] == 'F' );

	}
};

#endif //_T_GENETIC_CODE_H__
