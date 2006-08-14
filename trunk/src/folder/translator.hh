#ifndef TRANSLATOR_HH
#define TRANSLATOR_HH

#include <vector>
#include "protein.hh"

class Translator
{
private:
	double m_mutation_prob;

	Translator( const Translator & );
	Translator& operator=( const Translator & );
public:
	Translator( double mutation_prob);
	Translator();

	bool translateErrorFree( const Gene &g, Protein& residue_sequence ) const;
	bool translate( const Gene &g, Protein& residue_sequence ) const;

	int translateWeighted( const Gene &g, Protein& residue_sequence, const vector<vector<pair<double, int> > >& weights,
						   const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);
	int translateRelativeWeighted( const Gene &g, Protein& residue_sequence, const double relative_gene_weight,
								   const vector<vector<pair<double, int> > >& weights, const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);
};



#endif
