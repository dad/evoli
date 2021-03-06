#include "translator.hh"


#include "tools.hh"
#include "genetic-code.hh"


Translator::Translator( double mutation_prob, int length )
		: m_mutation_prob( mutation_prob ), m_L( length )
{}

Translator::Translator( int length )
		: m_mutation_prob( 0 ), m_L( length )
{}

bool Translator::translateErrorFree( const Genotype &g, int *residue_sequence )
{
	bool result = true;

	for ( int i=0; i<m_L; i++ )
	{
		int residue = GeneticCodeUtil::geneticCode[g[i]];
		if ( residue < 0 )
		{
			result = false;
			break;
		}
		else
		{
			residue_sequence[i] = residue;
		}
	}

	return result;

}

bool Translator::translate( const Genotype &g, int *residue_sequence ) {
	if ( m_mutation_prob == 0 )
		return translateErrorFree( g, residue_sequence );

	bool result = true;

	for ( int i=0; i<m_L; i++ ) {
		int residue = GeneticCodeUtil::geneticCode[g[i]];
		if ( residue < 0 ) {
			result = false;
			break;
		}
		else {
			if ( myRand() < m_mutation_prob )
				residue = ( residue
					    + (int) (20*myRand() )) % 20;
			residue_sequence[i] = residue;
		}
	}

	return result;


}

int Translator::translateWeighted( const Genotype &g, int *residue_sequence, const vector<vector<pair<double, int> > >& weights,
	const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated)
{
	double mut_weight_total = 0.0;
	for ( unsigned int i=0; i<g.size(); i++ ) {
		mut_weight_total += (1.0 + prefCodons[g[i]]*(nonPrefCodonPenalty-1));
	}

	truncated = false;
	int numErrors = 0;
	for ( int i=0; i<m_L && !truncated; i++ ) {
		int residue = GeneticCodeUtil::geneticCode[g[i]];
		residue_sequence[i] = residue;

		if ( residue < 0 ) {
			truncated = true;
		}
		else {
			double threshold = m_mutation_prob*(1.0 + prefCodons[g[i]]*(nonPrefCodonPenalty-1))/(mut_weight_total/g.size());
			double rand = myRand();
			if ( rand < threshold ) {
				// Weight the outcomes of a missense substitution.
				rand = myRand();
				double targ = 0.0;
				// The first member of the pair is a cumulative probability; the second is
				// the residue resulting from the error.
				for (unsigned int j=0; j<weights[g[i]].size() && (rand > targ); j++) {
					pair<double, int> p = weights[g[i]][j];
					residue_sequence[i] = p.second;
					targ = p.first;
				}
				if (residue_sequence[i] < 0) { // truncation error
					truncated = true;
				}
				if (residue_sequence[i] != residue) {
					numErrors++;
				}
			}
		}
	}
	return numErrors;
}

int Translator::translateRelativeWeighted( const Genotype &g, int *residue_sequence, const double relative_gene_weight,
	const vector<vector<pair<double, int> > >& weights,	const double* prefCodons, const double nonPrefCodonPenalty, bool& truncated)
{
	truncated = false;
	int numErrors = 0;
	// Hopefully, relative_gene_weight is tabulated in a reasonable way.  If there are no accuracy
	// differences between codons, relative_gene_weight should be equal to the number of codons g.size().
	double relative_site_weight = g.size()/relative_gene_weight;
	for ( int i=0; i<m_L && !truncated; i++ ) {
		int residue = GeneticCodeUtil::geneticCode[g[i]];
		residue_sequence[i] = residue;

		// With probability threshold_prob, make a translation error (possibly synonymous).
		double site_weight = (1.0 + prefCodons[g[i]]*(nonPrefCodonPenalty-1));
		double threshold_prob = m_mutation_prob * site_weight * relative_site_weight;
		double rand = myRand();
		if ( rand < threshold_prob ) {
			// We've made an error.  Now determine what it is.
			rand = myRand();
			double targ = 0.0;
			// The first member of the pair is a cumulative probability; the second is
			// the residue resulting from the error.
			// Find the event corresponding to rand.
			for (unsigned int j=0; j<weights[g[i]].size() && (rand > targ); j++) {
				pair<double, int> p = weights[g[i]][j];
				residue_sequence[i] = p.second;
				targ = p.first;
			}
			// Tabulate the results
			if (residue_sequence[i] < 0) { // truncation error
				truncated = true;
			}
			if (residue_sequence[i] != residue) {
				numErrors++;
			}
		}
	}
	return numErrors;
}
