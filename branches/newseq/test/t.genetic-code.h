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
		char aa = GeneticCodeUtil::geneticCode(c);
		TEST_ASSERT( aa == 'M' );
		// Systematically test all codons
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AGA")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAA")) == '*' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAG")) == '*' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AAA")) == 'K' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AAG")) == 'K' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AAU")) == 'N' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AAC")) == 'N' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("ACA")) == 'T' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("ACG")) == 'T' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("ACU")) == 'T' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("ACC")) == 'T' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AGA")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AGG")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AGU")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AGC")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AUA")) == 'I' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AUG")) == 'M' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AUU")) == 'I' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("AUC")) == 'I' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CAA")) == 'Q' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CAG")) == 'Q' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CAU")) == 'H' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CAC")) == 'H' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CCA")) == 'P' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CCG")) == 'P' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CCU")) == 'P' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CCC")) == 'P' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CGA")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CGG")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CGU")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CGC")) == 'R' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CUA")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CUG")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CUU")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("CUC")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GAA")) == 'E' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GAG")) == 'E' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GAU")) == 'D' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GAC")) == 'D' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GCA")) == 'A' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GCG")) == 'A' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GCU")) == 'A' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GCC")) == 'A' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GGA")) == 'G' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GGG")) == 'G' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GGU")) == 'G' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GGC")) == 'G' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GUA")) == 'V' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GUG")) == 'V' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GUU")) == 'V' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("GUC")) == 'V' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAA")) == '*' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAG")) == '*' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAU")) == 'Y' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UAC")) == 'Y' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UCA")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UCG")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UCU")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UCC")) == 'S' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UGA")) == '*' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UGG")) == 'W' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UGU")) == 'C' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UGC")) == 'C' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UUA")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UUG")) == 'L' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UUU")) == 'F' );
		TEST_ASSERT( GeneticCodeUtil::geneticCode(Codon("UUC")) == 'F' );

	}
	void TEST_FUNCTION( aa_index_conversion ) {
		const char *aas = "ACDEFGHIKLMNPQRSTVWY";
		for (int i=0; i<20; i++) {
			char aa = aas[i];
			int ind = GeneticCodeUtil::aminoAcidLetterToIndex(aa);
			char aa2 = GeneticCodeUtil::indexToAminoAcidLetter(ind);
			TEST_ASSERT(aa == aa2);
		}
	}
	void TEST_FUNCTION( codon_index_conversion ) {
		const char *aas = "ACDEFGHIKLMNPQRSTVWY";
		for (int i=0; i<64; i++) {
			Codon c1 = GeneticCodeUtil::indexToCodon(i);
			int c1_index = GeneticCodeUtil::codonToIndex(c1);
			Codon c2 = GeneticCodeUtil::indexToCodon(c1_index);
			TEST_ASSERT(i == c1_index);
			TEST_ASSERT(c1 == c2);
		}
	}
	void TEST_FUNCTION( syn_sites_zero ) {
		Codon c1("AUG");
		double syn_sites = 0;
		syn_sites = GeneticCodeUtil::calcSynonymousSites(c1, 7);
		TEST_ASSERT(syn_sites == 0.0);
	}

	void TEST_FUNCTION( syn_sites_arginine_2fold ) {
		Codon c1("AGA");
		double syn_sites = 0;
		syn_sites = GeneticCodeUtil::calcSynonymousSites(c1, 7);
		//cout << endl << c1 << " " << syn_sites << endl;
		TEST_ASSERT(syn_sites == 2.0/3.0);
	}
	void TEST_FUNCTION( syn_sites_arginine_4fold_dna ) {
		Codon c1("CGT");
		double syn_sites = 0;
		syn_sites = GeneticCodeUtil::calcSynonymousSites(c1, 7);
		//cout << endl << c1 << " " << syn_sites << endl;
		TEST_ASSERT(syn_sites == 3.0/3.0);
	}
	void TEST_FUNCTION( syn_sites_arginine_4fold_rna ) {
		Codon c1("CGU");
		double syn_sites = 0;
		syn_sites = GeneticCodeUtil::calcSynonymousSites(c1, 7);
		//cout << endl << c1 << " " << syn_sites << endl;
		TEST_ASSERT(syn_sites == 3.0/3.0);
	}


};

#endif //_T_GENETIC_CODE_H__
