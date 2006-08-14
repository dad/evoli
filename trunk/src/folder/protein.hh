#ifndef PROTEIN_HH
#define PROTEIN_HH

#include <vector>
#include <utility>
#include <iostream>
#include "genotype.hh"
using namespace std;

class ProteinFolder;
class Translator;
class Gene;

typedef pair<int, int> Contact;

class Sequence {
private:
	Sequence();
protected:
	vector<int> m_sequence;
	bool m_modified;

	virtual ~Sequence(void) {}

public:
	Sequence(const vector<int>& v);
	Sequence(const int length) : m_sequence(length) {}

	virtual uint16 length() const { return m_sequence.size(); }
	bool modified() const { return m_modified; }
	int& operator[](const int index);
	int operator[](const int index) const;
	void clear() {
		m_sequence.clear();
		m_modified = true;
	}

	typedef vector<int>::iterator iterator;
	typedef vector<int>::const_iterator const_iterator;

	iterator begin() { return m_sequence.begin(); }
	iterator end() { return m_sequence.end(); }

	const_iterator begin() const { return m_sequence.begin(); }
	const_iterator end() const { return m_sequence.end(); }
};

class Protein : public Sequence {
private:
	int m_structure_id;
	double m_free_energy;

	Protein();
protected:
public:
	Protein(const int length);
	Protein(const Protein& p);
	Protein(const string& s);
	~Protein() {}

	pair<int, double> fold(const ProteinFolder& folder);
	int distance(const Protein& p) const;
	vector<Contact> getContacts(const ProteinFolder& folder);
	string toString() const;

	Gene reverseTranslate() const;
	bool operator==(const Protein& p) const;
};



class Gene : public Sequence {
private:
protected:
public:
	Gene();
	Gene(const int length);
	Gene(const Gene& g);
	Gene(const string& s);
	~Gene() {}

	virtual uint16 length(void) const { return 3*m_sequence.size(); }
	virtual uint16 codonLength(void) const { return m_sequence.size(); }
	Protein translate(const Translator& t) const;
	Protein translate(void) const;

	static Gene createRandom(const int length);
	static Gene createRandomNoStops(const int length);
	bool mutate(const double prob);
	bool encodesFullLength(void) const;
	char getBase(const uint16 index) const;
	string toString() const;

	operator const Genotype&(void) const;

	/**
	 * Finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
	 */
	static Gene getSequenceForStructure( ProteinFolder &b, double free_energy_cutoff, const int struct_id );
	/**
	 * Finds a random sequence with folding energy smaller than cutoff
	 */
	static Gene getSequence( ProteinFolder &b, double free_energy_cutoff);
};

/*
ostream & operator<<( ostream &s, const Gene& g );
ostream & operator<<( ostream &s, const Protein& p );
*/
/**
 * Printing a Gene.
 **/
inline ostream & operator<<( ostream &s, const Gene &g ) {
	s << g.toString();
	return s;
}

/**
 * Printing a Protein.
 **/
inline ostream & operator<<( ostream &s, const Protein &p ) {
	for ( Protein::const_iterator it=p.begin(); it != p.end(); it++) {
		s << GeneticCodeUtil::residueLetters[*it+1]; // << "  ";
	}
	return s;
}

inline istream & operator>>( istream &s, Gene & g )
{
        string str;
		s >> str;
		g = Gene(str);
		return s;
}


#endif //PROTEIN_HH
