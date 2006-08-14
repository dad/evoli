#ifndef GENOTYPE_HH
#define GENOTYPE_HH

#include "codon.hh"
#include "genetic-code.hh"
#include "tools.hh"

#include <vector>

using namespace std;

typedef vector<int> Genotype;
typedef vector<int>::const_iterator GenotypeCIterator;
typedef vector<int>::iterator GenotypeIterator;


inline ostream & operator<<( ostream &s, const Genotype &g )
{
        GenotypeCIterator it = g.begin(), e = g.end();

        for ( ; it != e; it++ )
        {
                CodonUtil::printCodon( s, (*it) );
                //     s << " ";
        }
        return s;
}

inline istream & operator>>( istream &s, Genotype & g )
{
        char c1, c2, c3;

        g.clear();

        // read leading whitespace
        do
        {
                s.get( c1 );
        }
        while ( s && ( ( c1 == ' ' ) || ( c1 == '\n' ) || ( c1 == '\t' ) ) );

        // read until next whitespace
        while ( s && ( c1 != ' ' ) && ( c1 != '\n' ) && ( c1 != '\t' ) )
        {
                s.get( c2 );
                s.get( c3 );
                g.push_back( CodonUtil::lettersToCodon( c1, c2, c3 ) );
                s.get( c1 );
        }
        return s;
}

inline bool operator==(const Genotype& g1, const Genotype& g2) {
	bool res = (g1.size() == g2.size());
	for (unsigned int i=0; i<g1.size() && res; i++) {
		res = (g1[i] == g2[i]);
	}
	return res;
}


#endif
