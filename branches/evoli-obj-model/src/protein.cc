#include "protein.hh"
#include "protein-folder.hh"
#include "translator.hh"
#include "genotype.hh"
#include <sstream>


Sequence::Sequence(const vector<int>& v) {
	m_sequence = v;
}

int& Sequence::operator[](const int index) {
	m_modified = true;
	return m_sequence[index];
}

int Sequence::operator[](const int index) const {
	return m_sequence[index];
}



Protein::Protein(const int length) : Sequence(length) {
	m_free_energy = 1000000;
	m_structure_id = -1;
	m_modified = true;
}

Protein::Protein(const Protein& p) : Sequence(p.m_sequence) {
	m_free_energy = p.m_free_energy;
	m_structure_id = p.m_structure_id;
	m_modified = p.m_modified;
}

int Protein::distance(const Protein& p) const {
	int diffs = 0;
	Protein::const_iterator qit = this->begin();
	Protein::const_iterator pit = p.begin();

	for (; qit != this->end() && pit != p.end(); qit++, pit++) {
		if (*pit != *qit) {
			diffs++;
		}
	}
	return diffs;
}

pair<int, double> Protein::fold(const ProteinFolder& folder) {
	if (m_structure_id < 0 || m_modified) {
		m_free_energy = folder.foldProtein(*this);
		m_structure_id = folder.getLastFoldedProteinStructureID();
		m_modified = false;
	}
	return pair<int, double>(m_structure_id, m_free_energy);
}

vector<Contact> Protein::getContacts(const ProteinFolder& folder) {
	fold(folder);
	return folder.getStructure(m_structure_id)->getInteractingPairs();
}

Gene Protein::reverseTranslate() const {
	Gene g(length()*3);
	Protein::const_iterator it = begin();

	for (uint16 i=0; it != end(); it++, i++) {
		int aa = *it;
		int codon = GeneticCodeUtil::residueToCodonTable[aa+1];
		g[i] = codon;
	}

	return g;
}

bool Protein::operator==(const Protein& p) const {
	bool identical = (length() == p.length());
	Protein::const_iterator qit = begin();
	Protein::const_iterator pit = p.begin();

	for (; qit != end() && pit != p.end() && identical; qit++, pit++) {
		identical = (*pit == *qit);
	}
	return identical;
}

string Protein::toString() const {
	string res;
	for ( Protein::const_iterator it=begin(); it != end(); it++) {
		res += GeneticCodeUtil::residueLetters[*it+1];
	}
	return res;
}



Gene::Gene() : Sequence(0) {
	m_modified = true;
}

Gene::Gene(const int length) : Sequence(length/3) {
	m_modified = true;
}

Gene::Gene(const Gene& g) : Sequence(g.m_sequence) {
	m_modified = true;
}

Gene::Gene(const string& seq_string) : Sequence(0) {
	stringstream s(seq_string);
	char c1, c2, c3;

	// read leading whitespace
	do {
		s.get( c1 );
	}
	while ( s && ( ( c1 == ' ' ) || ( c1 == '\n' ) || ( c1 == '\t' ) ) );

	// read until next whitespace
	while ( s && ( c1 != ' ' ) && ( c1 != '\n' ) && ( c1 != '\t' ) ) {
		s.get( c2 );
		s.get( c3 );
		m_sequence.push_back( CodonUtil::lettersToCodon( c1, c2, c3 ) );
		s.get( c1 );
	}
}

char Gene::getBase(const uint16 index) const {
	char res = 'X';
	int triplet_pos = index % 3;
	int codon = (*this)[(index-triplet_pos)/3];
	//cout << codon << " ";
	//CodonUtil::printCodon(cout, codon);
	//cout << " ";
	int x = (codon >> (2*(2-triplet_pos))) & 3;
	switch( x )	{
		case 0:
			res = 'A';
			break;
		case 1:
			res = 'C';
			break;
		case 2:
			res = 'G';
			break;
		case 3:
			res = 'U';
			break;
	}
	//cout << res << endl;
	return res;
}

string Gene::toString() const {
	string res;
	for (uint16 i=0; i<m_sequence.size()*3; i++) {
		res += getBase(i);
	}
	return res;
}

Gene Gene::createRandom(const int length ) {
	Gene g( length );
	//vector<int>::iterator it = g.begin();
	Gene::iterator it = g.begin();

	for ( ; it != g.end(); it++) {
		*it = static_cast<int>( 64*drand48() );
	}
	return g;
}

Gene Gene::createRandomNoStops(const int length ) {
	Gene g( length );
	//vector<int>::iterator it = g.begin();
	Gene::iterator it = g.begin();

	for ( ; it != g.end(); it++) {
		do {
			*it = static_cast<int>( 64*myRand() );
		} while (GeneticCodeUtil::geneticCode[*it] < 0);
	}
	return g;
}

bool Gene::encodesFullLength(void) const {
	bool no_stop = true;
	Gene::const_iterator it = begin();
	while (it != end() && no_stop) {
		no_stop = (GeneticCodeUtil::geneticCode[*it] >= 0);
		it++;
	}
	return no_stop;
}

bool Gene::mutate(const double prob) {
	bool changed = false;
	Gene::iterator it = begin();
	for ( ; it != end(); it++) {
		int codon = *it;
		(*it) = CodonUtil::mutateCodon( prob, codon );
		if ( !changed && (*it) != codon ) {
			changed = true;
			m_modified = true;
		}
	}
	return changed;
}

Gene::operator const Genotype&(void) const {
	return m_sequence;
}

Protein Gene::translate(const Translator& t) const {
	Protein prot(length()/3);
	t.translate(*this, prot);
	return prot;
}

Protein Gene::translate(void) const {
	int len = length()/3;
	Translator t;
	Protein prot(len);
	t.translate(*this, prot);
	return prot;
}

Gene Gene::getSequenceForStructure( ProteinFolder &b, double free_energy_cutoff, const int struct_id )
{
	// find a random sequence with folding energy smaller than cutoff
	double G;
	int length = 3*b.getProteinLength();
	Gene g(length);
	pair<int, double> fdata;
	bool found = false;
	double min_free_energy_for_starting = max(0.0, free_energy_cutoff);

	// find sequence that encodes the desired structure
	do {
		g = Gene::createRandomNoStops( length );
		//cout << g << endl;
		Protein p = g.translate();
		//cout << p << endl;
		fdata = p.fold(b);
		found = (fdata.first == struct_id && fdata.second <= min_free_energy_for_starting);
		//cout << fdata.first << "\t" << fdata.second << "\t" << g << endl;
	} while ( !found );

	int fail_count = 0;
	int total_fail_count = 0;
	G = fdata.second;

	// optimize the sequence for stability
	do {
		Gene g2 = g;
		bool changed = false;
		do {
			changed = g2.mutate( 0.02 );
		} while (!changed);

		if (!g2.encodesFullLength()) {
			fail_count++;
			total_fail_count++;
		}
		else {
			Protein p = g2.translate();
			fdata = p.fold(b);
			//cout << fdata.first << "\t" << fdata.second << "\t" << G << endl;
			if (fdata.first == struct_id && fdata.second <= G-0.001) {
				// we found an improved sequence. grab it, and reset failure count
				g = g2;
				G = fdata.second;
				fail_count = 0;
				//cout << G << endl;
			}
		}

		// start again with random genotype if search is not successful after 50000 failures
		if ( fail_count > 50000 || total_fail_count > 1e6 )	{
			found = false;
			do {
				g = Gene::createRandomNoStops( length );
				Protein p = g.translate();
				fdata = p.fold(b);
				found = (fdata.first == struct_id && fdata.second <= min_free_energy_for_starting);
			} while ( !found );
			G = fdata.first;
			fail_count = 0;
			total_fail_count = 0;
		}
	} while( G > free_energy_cutoff );

	return g;
}

Gene Gene::getSequence( ProteinFolder &b, double free_energy_cutoff)
{
	// find a random sequence with folding energy smaller than cutoff
	double G;
	int length = 3*b.getProteinLength();
	Gene g(length);
	pair<int, double> fdata;
	bool found = false;
	double min_free_energy_for_starting = max(0.0, free_energy_cutoff);

	// find sequence that encodes the desired structure
	do {
		g = Gene::createRandomNoStops( length );
		//cout << g << endl;
		Protein p = g.translate();
		//cout << p << endl;
		fdata = p.fold(b);
		found = (fdata.second <= min_free_energy_for_starting);
		//cout << fdata.first << "\t" << fdata.second << "\t" << g << endl;
	} while ( !found );

	int fail_count = 0;
	int total_fail_count = 0;
	G = fdata.second;

	// optimize the sequence for stability
	do {
		Gene g2 = g;
		bool changed = false;
		do {
			changed = g2.mutate( 0.02 );
		} while (!changed);

		if (!g2.encodesFullLength()) {
			fail_count++;
			total_fail_count++;
		}
		else {
			Protein p = g2.translate();
			fdata = p.fold(b);
			//cout << fdata.first << "\t" << fdata.second << "\t" << G << endl;
			if (fdata.second <= G-0.001) {
				// we found an improved sequence. grab it, and reset failure count
				g = g2;
				G = fdata.second;
				fail_count = 0;
				//cout << G << endl;
			}
		}

		// start again with random genotype if search is not successful after 50000 failures
		if ( fail_count > 50000 || total_fail_count > 1e6 )	{
			found = false;
			do {
				g = Gene::createRandomNoStops( length );
				Protein p = g.translate();
				fdata = p.fold(b);
				found = (fdata.second <= min_free_energy_for_starting);
			} while ( !found );
			G = fdata.first;
			fail_count = 0;
			total_fail_count = 0;
		}
	} while( G > free_energy_cutoff );

	return g;
}

