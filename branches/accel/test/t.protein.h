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


#ifndef _T_PROTEIN_H__
#define _T_PROTEIN_H__
#include "cutee.h"
#include "protein.hh"
#include "gene-util.hh"
#include "compact-lattice-folder.hh"
#include "genetic-code.hh"
#include <iostream>
using namespace std;

struct TEST_CLASS( protein_gene_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	int countAminoAcids(char aa, const Protein& p){
	  int count = 0;
	  for (int i=0; i<p.size(); i++) {
		if (p[i] == aa) {
		  count++;
		}
	  }
	  return count;
	}


	
	void TEST_FUNCTION( sequence_from_string ) {
		string str("ACTGCT");
		Sequence s(str);
		TEST_ASSERT( str == s );
	}

	void TEST_FUNCTION( transcribe ) {
		string str("ACTGCT");
		CodingDNA dna(str);
		CodingRNA r = dna.transcribe();
		TEST_ASSERT( r.toString() == string("ACUGCU") );
	}
	void TEST_FUNCTION( length_test )
	{
		CodingDNA g = CodingDNA::createRandom(gene_length);
		TEST_ASSERT( g.length()==gene_length );
		return;
	}
	void TEST_FUNCTION( reverse_translate )
	{
		CodingDNA g = CodingDNA::createRandomNoStops(gene_length);
		//cout << g << endl;
		Protein p = g.translate();
		//cout << p << endl;
		Protein p2 = GeneUtil::reverseTranslate(p).translate();
		//cout << p2 << endl;
		TEST_ASSERT( p == p2 );
		return;
	}
	void TEST_FUNCTION( randomize_codons )
	{
		CodingDNA g = CodingDNA::createRandomNoStops(gene_length);
		for (int i=0; i<100; i++) {
			//cout << g << endl;
			Protein p = g.translate();
			//cout << p << endl;
			Protein p2 = GeneUtil::randomizeCodons(g).translate();
			//cout << p2 << endl;
			TEST_ASSERT_M( p == p2, p + "\n" + p2 );
		}
		return;
	}

	void TEST_FUNCTION( gene_string_gene )
	{
		CodingDNA g = CodingDNA::createRandom(gene_length);
		string str(g.toString());
		CodingDNA g2(str);
		TEST_ASSERT( g == g2 );
		return;
	}
	void TEST_FUNCTION( translate_known )
	{
		CodingDNA dna("ATGTGGGGG");
		CodingRNA rna = dna.transcribe();
		Translator t(0);
		Protein p(rna.codonLength());
		t.translate(rna, p);
		string str(p.toString());
		TEST_ASSERT( str == "MWG" );
		return;
	}
	void TEST_FUNCTION( full_length_true )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		TEST_ASSERT( g.encodesFullLength() );
		return;
	}
	void TEST_FUNCTION( full_length_false )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		Codon c("TAG");
		g.setCodon(3, c);
		TEST_ASSERT( g.getCodon(3) == c );
		TEST_ASSERT( !g.transcribe().encodesFullLength() );
		return;
	}
	void TEST_FUNCTION( full_length_protein )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		//cout << g << endl;
		Protein p = g.translate();
		//cout << p << endl;
		TEST_ASSERT( g.codonLength() == p.length() );
		return;
	}
	void TEST_FUNCTION( protein_from_string )
	{
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p = g.translate();
		string s = p.toString();
		Protein p2(s);
		//cout << p << endl;
		//cout << p2 << endl;
		TEST_ASSERT( p2 == p );
		return;
	}
	void TEST_FUNCTION( protein_from_strings_equality ) {
		Gene g = Gene::createRandomNoStops(gene_length);
		Protein p1 = g.translate();
		Protein p2(p1.toString());
		TEST_ASSERT( p1 == p2 );
	}
	void TEST_FUNCTION( random_protein ) {
	  string letters = "ACDEFGHIKLMNPQRSTVWY";
	  Protein p = Protein::createRandom(1000);
	  for (int i=0; i<letters.size(); i++) {
		char aa = letters[i];
		int count = countAminoAcids(aa, p);
		// Let this be random; the probability of failure is 
		// extraordinarily low for a sequence of this length.
		TEST_ASSERT(count > 0);
	  }

	}
	// Test getDifferences()
};


#endif //_T_PROTEIN_H__
