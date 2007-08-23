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


#ifndef _T_MUTATOR_H__
#define _T_MUTATOR_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "protein.hh"
#include "gene-util.hh"
#include "tools.hh"
#include <fstream>
#include <cmath>

struct TEST_CLASS( mutator_basic )
{
	void TEST_FUNCTION( is_mutated ) {
		CodingDNA g = CodingDNA::createRandomNoStops(75);
		SimpleMutator mut(0.01);
		for (int i=0; i<100; i++) {
			CodingDNA g2 = g;
			bool changed = mut.mutate(g2);
			TEST_ASSERT( changed == (g2 != g) );
		}
	}

	void TEST_FUNCTION( mutator_push ) {
		vector<double> AtoCGT;
		vector<double> CtoGTA;
		vector<double> GtoTAC;
		vector<double> TtoACG;
		AtoCGT.push_back( 1. );	// A -> C
		AtoCGT.push_back( 0. );	// A -> G
		AtoCGT.push_back( 0. );	// A -> T
		CtoGTA.push_back( 1. );	// C -> G
		CtoGTA.push_back( 0. );	// C -> T
		CtoGTA.push_back( 0. );	// C -> A
		GtoTAC.push_back( 1. );	// G -> T
		GtoTAC.push_back( 0. );	// G -> A
		GtoTAC.push_back( 0. );	// G -> C
		TtoACG.push_back( 0. );	// T -> A
		TtoACG.push_back( 0. );	// T -> C
		TtoACG.push_back( 0. );	// T -> G

		Polymerase p2( 1, AtoCGT, CtoGTA, GtoTAC, TtoACG );
		// Every step leads to a certain mutation.  After one step, all A's should be C's;
		// after two steps, G's, and after three steps, T's.  Similarly all
		vector<char> steps;
		steps.push_back('A');
		steps.push_back('C');
		steps.push_back('G');
		steps.push_back('T');
		steps.push_back('T');

		int L = 10;
		int reps = steps.size()-1;
		NucleotideSequence s( L, 'A' );
		for ( int i=0; i<reps; i++ ) {
			bool mutated = p2.mutate(s);
			TEST_ASSERT( reps>2 || mutated ); // Probability of a mutation is 1!
			for (int j=0; j<L; j++) {
				//cout << s[j] << " " << steps[i+1] << endl;
				TEST_ASSERT( s[j] == steps[i+1] );
			}
		}
	}
	void TEST_FUNCTION( is_poisson_accurate ) {
	  int size = 1000;
	  double sum_j = 0;
	  double sum_sq_j = 0;
	  int MAX = 1000;
	  for (int i=0; i<MAX; i++) {
	    CodingDNA g = CodingDNA::createRandomNoStops(size);
	    SimpleMutator mut(0.01);
	    double j = 0;
	    CodingDNA g2 = g;
	    mut.mutate(g2);
	    for (int k=0; k<size; k++) {
	      if (g[k]!=g2[k]){ 
		j++;
	      }
	    }
	    cout << "There were: " << j << " mutations in the " << i+1 << "th sequence. " <<  endl; 
	    sum_j += j;
	    sum_sq_j += sum_j*sum_j;
	  }
	  cout << " The sum of all the mutations is: " << sum_j << endl;
	  cout << " The sum_squared of all the mutations is: " << sum_sq_j << endl;
	  Accumulator ac(sum_j, sum_sq_j, MAX);
	  cout << " The mean of the mutations is: " << ac.mean() << endl;
	  cout << " The standard deviation of the mutations is: " << ac.stdev() << endl;
	  cout << " The standard error of the mutations is: " << ac.stderror() << endl;
	}
	
};

#endif // _T_MUTATOR_H__
