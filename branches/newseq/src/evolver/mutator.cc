#include "mutator.hh"
#include "codon.hh"


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
	string nts("ACGT");
	CodingDNA::iterator it = seq.begin();
	for (int i=0; i<seq.length(); i++) {
		if (Random::runif() < m_mutation_rate) {
			char nt = nts[Random::rint( 4 )];
			seq[i] = nt;
			changed = true;
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
	string nts("ACGT");
	CodingDNA::iterator it = seq.begin();
	for (int i=0; i<seq.length(); i++) {
		if (Random::runif() < m_mutation_rate) {
			char nt = nts[Random::rint( 4 )];
			seq[i] = nt;
			changed = true;
		}
	}
	return changed;
}
