#ifndef PROTEIN_HH
#define PROTEIN_HH

#include <vector>
#include <string>
#include "tools.hh"
#include "sequence.hh"
#include "genetic-code.hh"

using namespace std;

typedef pair<int, int> Contact;

class Protein : public Sequence {
private:
	Protein();
protected:
public:
	Protein(const int length);
	Protein(const Protein& p);
	Protein(const string& s);
	~Protein() {}

	int distance(const Protein& p) const;
	string toString() const;
};

class Translator;

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

inline istream & operator>>( istream &s, Gene & g ) {
	string str;
	s >> str;
	g = Gene(str);
	return s;
}


#endif //PROTEIN_HH
