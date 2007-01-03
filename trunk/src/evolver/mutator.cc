#include "mutator.hh"
#include "codon.hh"
#include "genetic-code.hh"
#include "random.hh"

#include <iostream>
using namespace std;

////////////////////
// SimpleMutator methods
////////////////////

SimpleMutator::SimpleMutator(double mutation_rate) {
	m_mutation_rate = mutation_rate;
}

SimpleMutator::~SimpleMutator() {
}

bool SimpleMutator::mutate(CodingDNA& seq) {
	bool changed = false;
	const char* mutA = "CGT";
	const char* mutC = "GTA";
	const char* mutG = "TAC";
	const char* mutT = "ACG";
	for (int i=0; i<seq.length(); i++) {
		char old_nt = seq[i];
		if (Random::runif() < m_mutation_rate) {
			changed = true;
			int j = Random::rint( 3 );
			switch( seq[i] ){
			case 'A':
				seq[i] = mutA[j]; break;
			case 'C':
				seq[i] = mutC[j]; break;
			case 'G':
				seq[i] = mutG[j]; break;
			case 'T':
				seq[i] = mutT[j]; break;
			default:
				assert( false ); // should never get here
			}
		}
	}
	return changed;
}



////////////////////
// Polymerase methods
////////////////////

Polymerase::Polymerase(double mutation_rate) {
	m_mutation_rate = mutation_rate;
	// DAD: implement
	// m_mutation_matrix = 
}

Polymerase::Polymerase(double mutation_rate, vector<vector<double> >& mutation_matrix) {
	m_mutation_rate = mutation_rate;
	m_mutation_matrix = mutation_matrix;
}

Polymerase::~Polymerase() {
}

bool Polymerase::mutate(CodingDNA& seq) {
	// DAD: implement mutation matrix
	bool changed = false;
	for (int i=0; i<seq.length(); i++) {
		char old_nt = seq[i];
		if (Random::runif() < m_mutation_rate) {
			do {
				seq[i] = GeneticCodeUtil::DNA_NUCLEOTIDES[Random::rint( 4 )];
			} while (seq[i] == old_nt);
			changed = true;
		}
	}
	return changed;
}
