#ifndef TRANSLATOR_HH
#define TRANSLATOR_HH

#include "genotype.hh"


class Translator
{
private:
	double m_mutation_prob;
	int m_L;

	Translator( const Translator & );
	Translator& operator=( const Translator & );
public:
	Translator( double mutation_prob, int length );
	Translator( int length );

	bool translateErrorFree( const Genotype &, int *residue_sequence );
	bool translate( const Genotype &, int *residue_sequence );
	int translateWeighted( const Genotype &g, int *residue_sequence, const vector<vector<pair<double, int> > >& weights,
		const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);

	int translateRelativeWeighted( const Genotype &g, int *residue_sequence, const double relative_gene_weight,
		const vector<vector<pair<double, int> > >& weights,	const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated);

};













#endif
