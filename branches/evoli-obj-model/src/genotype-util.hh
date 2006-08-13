#ifndef GENOTYPE_UTIL_HH
#define GENOTYPE_UTIL_HH

#include "codon.hh"
#include "genetic-code.hh"
#include "tools.hh"
#include "genotype.hh"
#include "translator.hh"
#include "protein-folder.hh"

#include <algorithm>
#include <iostream>

using namespace std;

class GenotypeUtil
{
public:
	/**
	 * Outputs the given genotype to the given stream.
	 *
	 * This function is deprecated. Better use << operator.
	 **/
	static void printGenotype( ostream & s, const Genotype & g )
	{
		s << g << endl;
   	}

	static void printProtein( ostream& s, const Genotype &g ) {
		for ( unsigned int i=0; i<g.size(); i++ ) {
			s << GeneticCodeUtil::residueLetters[ GeneticCodeUtil::geneticCode[g[i]]+1];
		}
		s << endl;
	}

	/**
	 * Reads genotype from the given stream.
	 * This function is deprecated. Better use >> operator.
	 **/
	static void readGenotype( istream & s, Genotype & g )
	{
		s >> g;
	}

	/**
	 * Creates a random genotype of the given length. Length is here
	 * measured in codons, not in bases!
	 **/
	static Genotype createRandomGenotype( int length )
	{
		Genotype g;
		g.resize( length );
		GenotypeIterator it = g.begin(), e = g.end();
		for ( ; it != e; it++ )
			(*it) =  static_cast<int>( 64*drand48() );
		return g;
	}

	/**
	 * Creates a random genotype of the given length. Length is here
	 * measured in codons, not in bases!
	 **/
	static Genotype createRandomGenotypeNoStops( int length )
	{
		Genotype g;
		g.resize( length );
		GenotypeIterator it = g.begin(), e = g.end();
		int codon = 0;
		for ( ; it != e; it++ ) {
			do {
				codon = static_cast<int>( 64*drand48() );
			} while (GeneticCodeUtil::geneticCode[codon] < 0); // eliminate stop codons
			*it = codon;
		}
		return g;
	}


	/**
	 * Mutates the given genotype with the per-site mutation rate U.
	 * A base is equally likely to be substituted by all other bases.
	 *
	 * Returns true if a mutation did occur, otherwise returns false.
	 **/
	static bool mutateGenotype( Genotype & g, double U )
	{
		bool changed = false;
		GenotypeIterator it = g.begin(), e = g.end();
		for ( ; it != e; it++ )
		{
			int c = (*it);
			(*it) = CodonUtil::mutateCodon( U, c );
			if ( !changed && (*it) != c )
				changed = true;
		}
		return changed;
	}

	/**
	 * Finds a random sequence with folding energy smaller than cutoff.
	 */
	static Genotype getSequence( ProteinFolder &b, double free_energy_cutoff )
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		Genotype g, g2;
		pair<double, int> fdata;
		bool found = false;
		Translator t(0, b.getProteinLength());
		int *seq = new int[b.getProteinLength()];
		double min_free_energy_for_starting = max(0.0, free_energy_cutoff);

		// find sequence with at least a reasonable starting stability
		do {
			g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
			fdata = GenotypeUtil::translateAndFold(b, g);
			found = (fdata.first <= min_free_energy_for_starting);
		} while ( !found );

		int fail_count = 0;
		int total_fail_count = 0;
		G = fdata.first;

		// optimize the sequence for stability
		do {
			g2 = g;
			bool changed = false;
			do {
				changed = GenotypeUtil::mutateGenotype( g2, 0.02 );
			} while (!changed);

			if (t.translateErrorFree(g2, seq)) {
				fdata = GenotypeUtil::translateAndFold(b, g2);
				if (fdata.first < G) {
					// we found an improved sequence. grab it, and reset failure count
					g = g2;
					G = fdata.first;
					fail_count = 0;
				}
			}
			else {
				fail_count++;
				total_fail_count++;
			}

			// start again with random genotype if search is not successful after 50000 failures
			if ( fail_count > 20000 || total_fail_count > 1e5 )	{
				found = false;
				do  {
					g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
					fdata = GenotypeUtil::translateAndFold(b, g);
					found = (fdata.first <= min_free_energy_for_starting);
				} while ( !found );
				G = fdata.first;
				// Reset counts
				fail_count = 0;
				total_fail_count = 0;
			}
		}
		while( G > free_energy_cutoff );

		delete [] seq;
		return g;
	}


	/**
	 * Finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
	 */
	static Genotype getSequenceForStructure( ProteinFolder &b, double free_energy_cutoff, const int struct_id )
	{
		// find a random sequence with folding energy smaller than cutoff
		double G;
		Genotype g, g2;
		pair<double, int> fdata;
		bool found = false;
		Translator t(0, b.getProteinLength());
		int *seq = new int[b.getProteinLength()];
		double min_free_energy_for_starting = max(0.0, free_energy_cutoff);

		// find sequence that encodes the desired structure
		do {
			g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
			found = t.translateErrorFree(g, seq) && b.isFoldedBelowThreshold(seq, struct_id, min_free_energy_for_starting);
		} while ( !found );

		int fail_count = 0;
		int total_fail_count = 0;
		G = fdata.first;

		// optimize the sequence for stability
		do {
			g2 = g;
			bool changed = false;
			do {
				changed = GenotypeUtil::mutateGenotype( g2, 0.02 );
			} while (!changed);

			if (t.translateErrorFree(g2, seq) && b.isFoldedBelowThreshold(seq, struct_id, G-0.001)) {
				// we found an improved sequence. grab it, and reset failure count
				g = g2;
				G = b.foldProtein(seq);
				fail_count = 0;
			}
			else {
				fail_count++;
				total_fail_count++;
			}

			// start again with random genotype if search is not successful after 50000 failures
			if ( fail_count > 50000 || total_fail_count > 1e6 )	{
				found = false;
				do  {
					g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
					found = t.translateErrorFree(g, seq) && b.isFoldedBelowThreshold(seq, struct_id, min_free_energy_for_starting);
				} while ( !found );
				G = fdata.first;
				fail_count = 0;
				total_fail_count = 0;
			}
		}
		while( G > free_energy_cutoff );

		delete [] seq;
		return g;
	}

	/**
	* Randomly replaces the codon at a random position with another
	* (not identical) codon, such that the amino acid has changed. Stop
	* codons are excluded.  All amino acids are equally likely to
	* occur.
	*
	* Returns the position at which the mutation was done.
	**/
	static int nonsynPointMutation( Genotype & g )
	{
		int pos = (int) (g.size()*myRand());
		int offset = (int) (19*myRand()) + 1;
		int aa = GeneticCodeUtil::geneticCode[g[pos]];
		// now mutate
		aa = ((aa + offset) % 20);
		g[pos]=GeneticCodeUtil::residueToCodonTable[aa+1];
		return pos;
	}


	/**
	 * Calculates the number of synonymous sites in the genotype.
	 **/
	static double calcSynonymousSites( const Genotype & g )
	{
		GenotypeCIterator it = g.begin(), e = g.end();
		double s = 0;
		for ( ; it != e; it++ )
			s += GeneticCodeUtil::calcSynonymousSites( *it );
		return s;
	}

	/**
	 * Calculates the total number of sites in the genotype.
	 **/
	static int calcTotalSites( const Genotype & g )
	{
		return 3*g.size();
	}

	/**
	 * Calculates the number of synonymous and nonsynonymous sites in surface and core.
	 **/
	static void calcSNSitesSurfaceCore( double &NSurf, double &NCore, double &SSurf, double &SCore,
					    const Genotype &g, const vector<int> &surface )
	{
		assert( g.size() == surface.size() );
		int e = g.size();
		NSurf = NCore = SSurf = SCore = 0;
		double tmp_N, tmp_S;

		for ( int i=0; i<e; i++ )
		{
			tmp_S = GeneticCodeUtil::calcSynonymousSites( g[i] );
			tmp_N = 3 - tmp_S;

			if ( surface[i] )
			{
				NSurf += tmp_N;
				SSurf += tmp_S;
			}
			else
			{
				NCore += tmp_N;
				SCore += tmp_S;
			}
		}
	}

	/**
	 * Calculates the number of synonymous and nonsynonymous substitutions
	 * from genotype 1 to genotype 2 (these numbers are not! normalized by
	 * the number of synonymous/nonsynonymous sites).
	 **/
	static void calcDnDs( double &dn, double &ds, const Genotype &g1, const Genotype &g2 )
	{
		assert( g1.size() == g2.size() );
		int e = g1.size();
		dn = ds =0;
		double tmp_dn, tmp_ds;

		for ( int i=0; i<e; i++ )
		{
			GeneticCodeUtil::calcDnDs( tmp_dn, tmp_ds, g1[i], g2[i] );
			//			 CodonUtil::printCodon( cout, g1[i] );
			//			 cout << " ";
			//			 CodonUtil::printCodon( cout, g2[i] );
			//			 cout << " " << tmp_dn << " " << tmp_ds << endl;
			dn += tmp_dn;
			ds += tmp_ds;
		}
	}

	/**
	 * Same as calcDnDs, but distinguishes between surface/core residues.
	 * The vector "surface" is expected to be of the same length as the protein,
	 * and must contain a 1 at surface positions and a 0 at core positions.
	 **/
	static void calcDnDsSurfaceCore( double &dnSurf, double &dnCore, double &dsSurf, double &dsCore,
					 const Genotype &g1, const Genotype &g2, const vector<int> &surface )
	{
		assert( g1.size() == g2.size() );
		assert( g1.size() == surface.size() );
		int e = g1.size();
		dnSurf = dnCore = dsSurf = dsCore = 0;
		double tmp_dn, tmp_ds;

		for ( int i=0; i<e; i++ )
		{
			GeneticCodeUtil::calcDnDs( tmp_dn, tmp_ds, g1[i], g2[i] );
			//			 CodonUtil::printCodon( cout, g1[i] );
			//			 cout << " ";
			//			 CodonUtil::printCodon( cout, g2[i] );
			//			 cout << " " << tmp_dn << " " << tmp_ds << endl;
			if ( surface[i] )
			{
				dnSurf += tmp_dn;
				dsSurf += tmp_ds;
			}
			else
			{
				dnCore += tmp_dn;
				dsCore += tmp_ds;
			}
		}
	}

	/**
	 * Calculates the neutrality of the given genotype. Cutoff is the free energy cutoff
	 * below which the protein folds.
	 **/
	static double calcNeutrality( ProteinFolder &b, const Genotype &g, double cutoff )
	{
		int l = b.getProteinLength();
		Translator t( 0, l );
		int *seq = new int[l];
		double G;
		int id;

		if ( t.translateErrorFree( g, seq ) )
		{
			G = b.foldProtein( seq );
			id = b.getLastFoldedProteinStructureID();
		}
		else
			return 0;

		if ( G > cutoff )
			return 0;

		int count = 0;
		// go through all positions in the genome
		for ( int i=0; i<l; i++ )
		{

			// go through all possible point mutations
			// (avoid operator %, which can be very slow)
			const int res = seq[i];
			int tempres = res + 1;
			while ( tempres < 20 )
			{
				seq[i] = tempres;
				// sequence folds into correct structure with low free energy?
				if ( b.foldProtein( seq ) < cutoff )
					if ( b.getLastFoldedProteinStructureID() == id )
						count += 1;
				tempres++;
			}

			tempres = 0;
			while ( tempres < res )
			{
				seq[i] = tempres;
				if ( b.foldProtein( seq ) < cutoff )
					if ( b.getLastFoldedProteinStructureID() == id )
						count += 1;
				tempres++;
			}
			seq[i] = res;
		}

		delete [] seq;

		return (double) count / (double) (19*l);
	}

	/**
	 * Calculates the number (not fraction) of possible neutral substitutions at the given site.
	 * Cutoff is the free energy cutoff below which the protein folds.
	 **/
	static pair<int, int> calcSiteNeutrality( ProteinFolder &b, const Genotype &g, double cutoff, int site )
	{
		int l = b.getProteinLength();
		Translator t( 0, l );
		int *seq = new int[l];
		double G;
		int id;

		if ( t.translateErrorFree( g, seq ) )
		{
			G = b.foldProtein( seq );
			id = b.getLastFoldedProteinStructureID();
		}
		else
			return pair<int, int>( 0, 0 );

		if ( G > cutoff )
			return pair<int, int>( 0, 0 );

		int count = 0, count_res = 0;

		// go through all possible point mutations
		// (avoid operator %, which can be very slow)
		const int res = seq[site];
		int tempres = res + 1;
		while ( tempres < 20 )
		{
			seq[site] = tempres;
			// sequence folds into correct structure with low free energy?
			if ( b.foldProtein( seq ) < cutoff )
				if ( b.getLastFoldedProteinStructureID() == id )
				{
					count += 1;
					if ( GeneticCodeUtil::singleSubstsTranslErrors[g[site]][tempres] == 1 )
						count_res += 1;
				}
			tempres++;
		}

		tempres = 0;
		while ( tempres < res )
		{
			seq[site] = tempres;
			if ( b.foldProtein( seq ) < cutoff )
				if ( b.getLastFoldedProteinStructureID() == id )
				{
					count += 1;
					if ( GeneticCodeUtil::singleSubstsTranslErrors[g[site]][tempres] == 1 )
						count_res += 1;
				}
			tempres++;
		}
		seq[site] = res;

		delete [] seq;
		return pair<int, int>( count, count_res );
	}

	/**
	 * Translates the DNA sequence into a residue sequence, and tries to fold it.
	 * The pair is returned contains the free energy of folding (pair.first) and the
	 * id of the protein structure (pair.second). If the sequence doesn't translate (because
	 * it contains a stop codon), then pair.first=10000000 and pair.second=-1.
	 **/
	static pair<double, int> translateAndFold( ProteinFolder &b, const Genotype &g )
	{
		int l = b.getProteinLength();
		Translator t( 0, l );
		int *seq = new int[l];
		int id;
		double G;

		if ( t.translateErrorFree( g, seq ) )
		{
			G = b.foldProtein( seq );
			id = b.getLastFoldedProteinStructureID();
		}
		else
		{
			G = 10000000;
			id = -1;
		}

		delete [] seq;
		return pair<double, int>( G, id );
	}

	/**
	 * Calculates the fraction of optimal codons from the set of optimal codons.
	 **/
	static double calcFop( const Genotype &g, const vector<bool>& is_optimal ) {
		int count = 0, tot=0;
		GenotypeCIterator it = g.begin(), e = g.end();
		int methionine = 14;
		int tryptophan = 58;
		int codon = -1;

		for ( ; it != e; it++ )	{
			codon = *it;
			if ( !(codon == methionine || codon == tryptophan)) {
				tot += 1;
				if (is_optimal[codon]) {
					count += 1;
				}
			}
		}
		return count / (double)tot;
	}

	/**
	 * Finds the lowest-cost codon for each amino acid and designates it optimal.
	 **/
	static vector<bool> getOptimalCodons(const double *codon_costs, const vector<vector<pair<double, int> > >& weights, const double unopt_codon_cost) {
		// Compute fraction accurately translated (facc) for each codon.  Optimal codons
		// have facc significantly higher than average for the family.
		int max_reps = 1000;
		Translator t(0.5,1); // Translator makes errors half the time.
		int res[1];
		Genotype g(1);
		bool truncated = false;

		vector<pair<int,int> > family_accuracies(20, pair<int,int>(0,0));  // Overall (accurate,total) pairs for all amino-acid families
		vector<pair<int,int> > codon_accuracies(64, pair<int,int>(0,0));   // (accurate,total) pairs for all codons
		vector<int> family_degeneracies(20,0);

		for (int codon=0; codon<64; codon++) {
			int aa = GeneticCodeUtil::geneticCode[codon];
			if (aa >= 0) {  // if not stop codon
				family_degeneracies[aa]++;
				for (int reps=0; reps<max_reps; reps++) {
					g[0] = codon;
					t.translateWeighted(g, res, weights, codon_costs, unopt_codon_cost, truncated);
					// Increment totals
					family_accuracies[aa].second++;
					codon_accuracies[codon].second++;
					if (res[0] == aa) {  // Properly translated
						// Increment accuracy counts
						family_accuracies[aa].first++;
						codon_accuracies[codon].first++;
					}
				}
			}
		}

		// Now we've tabulated overall accuracies and individual accuracies.
		// First designate as optimal any codon with the highest accuracy in its class, ensuring
		// that at least one optimal codon is found.
		vector<int> max_facc_codons(20, -1);  // The best codons so far.
		vector<double> aa_max_faccs(20, -1.0);  // The best codons' faccs.
		vector<bool> is_optimal(64, false);

		for (int codon=0; codon<64; codon++) {
			int aa = GeneticCodeUtil::geneticCode[codon];
			if (aa >= 0) {  // if not stop codon
				pair<int,int> acc = codon_accuracies[codon];
				double facc = acc.first/(double)acc.second;
				if (facc > aa_max_faccs[aa]) {
					aa_max_faccs[aa] = facc;
					max_facc_codons[aa] = codon;
				}
			}
		}

		for (int aa=0; aa<20; aa++) {
			if (family_degeneracies[aa]>1) { // Don't bother with nondegenerate codons
				is_optimal[max_facc_codons[aa]] = true;
			}
		}

		for (int codon=0; codon<64; codon++) {
			int aa = GeneticCodeUtil::geneticCode[codon];
			if (is_optimal[codon]) {
				int tot = family_accuracies[aa].second - codon_accuracies[codon].second;
				int acc = family_accuracies[aa].first - codon_accuracies[codon].first;
				int myacc = codon_accuracies[codon].first;
				cout << "# Opt for " << GeneticCodeUtil::residueLetters[aa+1] << " = ";
				CodonUtil::printCodon(cout, codon);
				cout << ", acc = " << myacc/(double)max_reps << " vs. family acc = " << acc/(double)tot << endl;
			}
		}

		// Now designate optimal any codons with accuracies significantly higher than the family average
		// By significantly, we will mean that the probability of obtaining the given
		// number of accurate translations in max_reps tries is < 0.01/family_degeneracy given the
		// observed proportion of accurate translations for the family and a binomial distribution.

		// Disallow more than 2 optimal codons per family
		vector<int> aa_optimal_codons(20,0);

		max_facc_codons.clear();
		max_facc_codons.resize(20,-1);
		aa_max_faccs.clear();
		aa_max_faccs.resize(20,-1.0);

		for (int codon=0; codon<64; codon++) {
			int aa = GeneticCodeUtil::geneticCode[codon];
			// Don't both reconsidering optimal codons or nondegenerate codons
			if (is_optimal[codon] || family_degeneracies[aa] < 2) {
				continue;
			}
			if (aa >= 0) {  // if not stop codon
				pair<int,int> acc = codon_accuracies[codon];
				pair<int,int> tacc = family_accuracies[aa];
				double facc = tacc.first/(double)tacc.second;

				double p_threshold = 0.01/family_degeneracies[aa];  // Bonferroni correction for multiple tests
				double p = 0.0; // probability that at least acc.first accurate translations would occur by chance given family accuracy
				for (int succ=acc.first; succ<=max_reps; succ++) {
					p += dbinom(succ, max_reps, facc);
				}
				if (p < p_threshold) {
					// Store the second best
					if (facc > aa_max_faccs[aa]) {
						aa_max_faccs[aa] = facc;
						max_facc_codons[aa] = codon;
					}
				}
			}
		}

		// Having identified the second-best, in cases where accuracy was significantly better, designate them optimal, too.
		for (int aa=0; aa<20; aa++) {
			if (family_degeneracies[aa]>1 && max_facc_codons[aa]!=-1) { // Don't bother with nondegenerate families or families with no significantly better codons
				is_optimal[max_facc_codons[aa]] = true;
			}
		}


		// debugging
		/*
		int num_optimal_codons = 0;
		for (int j=0; j<64; j++) {
			CodonUtil::printCodon(cout, j);
			cout << "\t" << GeneticCodeUtil::residueLetters[GeneticCodeUtil::geneticCode[j]+1] << "\t" << is_optimal[j] << "\t" << unopt_codon_cost*codon_costs[j] << endl;
			if (is_optimal[j]) {
				num_optimal_codons++;
			}
		}
		cout << "Total of " << num_optimal_codons << " optimal codons." << endl;
		*/

		return is_optimal;
	}

	/**
	 * As calcFop, but considers surface and core regions separately
	 **/
	static void calcFopSurfaceCore( double &Fop_surf, double &Fop_core, const Genotype &g, const double *codon_costs, const vector<int> &surface )
	{
		assert( g.size() == surface.size() );

		Fop_surf = Fop_core = 0;

		if ( codon_costs )
		{
			int sopt = 0, copt = 0, scount = 0, ccount = 0;

			for ( uint i=0; i<g.size(); i++ )
			{
				if ( surface[i] )
				{
					scount += 1;
					if ( codon_costs[g[i]] == 0 )
						sopt += 1;
				}
				else
				{
					ccount += 1;
					if ( codon_costs[g[i]] == 0 )
						copt += 1;
				}
			}
			Fop_surf = (double) sopt / (double) scount;
			Fop_core = (double) copt / (double) ccount;
		}

	}

	/**
	 * Produces a sequence encoding the same amino acid sequence as gene g, but with randomly chosen codons.
	 * Returns the original sequence if it doesn't translate properly.
	 **/
	static Genotype randomizeCodons( const Genotype &g ) {
		int l = g.size();
		Translator t( 0, l );
		int *seq = new int[l];
		Genotype result = g;

		if ( t.translateErrorFree( g, seq ) ) {
			for ( int i=0; i<l; i++ ) {
				int num_alts = GeneticCodeUtil::residueToAllCodonsTable[seq[i]][0];
				// Now choose from alternatives at random
				int randbin = (int)(myRand()*num_alts);
				result[i] = GeneticCodeUtil::residueToAllCodonsTable[seq[i]][randbin+1];
			}
		}
		delete [] seq;
		return result;
	}



};



#endif
