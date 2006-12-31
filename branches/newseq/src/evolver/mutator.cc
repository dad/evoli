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

bool SimpleMutator::mutate(Gene& seq) {
	bool changed = false;
	Gene::iterator it = seq.begin();
	for ( ; it != seq.end(); it++) {
		int codon = *it;
		(*it) = CodonUtil::mutateCodon( m_mutation_rate, codon );
		if ( !changed && (*it) != codon ) {
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

bool Polymerase::mutate(Gene& seq) {
	// DAD: implement mutation matrix
	bool changed = false;
	Gene::iterator it = seq.begin();
	for ( ; it != seq.end(); it++) {
		int codon = *it;
		(*it) = CodonUtil::mutateCodon( m_mutation_rate, codon );
		if ( !changed && (*it) != codon ) {
			changed = true;
		}
	}
	return changed;
	
}
