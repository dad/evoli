#include "fitness-evaluator.hh"

#include "genetic-code.hh"
#include "protein-folder.hh"
#include "translator.hh"
#include "gene-util.hh"

#include <cmath>
#include <iomanip>

FitnessEvaluator::FitnessEvaluator()
{}

FitnessEvaluator::~FitnessEvaluator()
{}

ProteinFreeEnergyFitness::ProteinFreeEnergyFitness( ProteinFolder *protein_folder )
	: m_protein_folder( protein_folder )
{
}

ProteinFreeEnergyFitness::~ProteinFreeEnergyFitness() {
}

double ProteinFreeEnergyFitness::getFitness( const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		return getFitness(p);
	}
	else
		return 0;
}

double ProteinFreeEnergyFitness::getFitness( Protein &p ) {
	pair<int, double> pf = p.fold(*m_protein_folder);
	return exp(-pf.second);
}



ProteinStructureFitness::ProteinStructureFitness( ProteinFolder *protein_folder, int protein_structure_ID, double max_free_energy )
	: m_protein_folder( protein_folder ),
	  m_protein_structure_ID( protein_structure_ID ),
	  m_max_free_energy( max_free_energy )
{
}

ProteinStructureFitness::~ProteinStructureFitness() {
}

double ProteinStructureFitness::getFitness( const Gene &g ) {
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		return getFitness(p);
	}
	else
		return 0;
}

double ProteinStructureFitness::getFitness( Protein &p ) {
	pair<int, double> pf = p.fold(*m_protein_folder);
	if ( pf.second > m_max_free_energy )
		return 0;

	if ( pf.first != m_protein_structure_ID )
		return 0;
	return 1;
}

struct greater_pair_first
{
	bool operator()(const pair<double, int>& p1, const pair<double, int>& p2)
	{
		return p1.first > p2.first;
	}
};

ErrorproneTranslation::ErrorproneTranslation() {
	m_protein_folder = NULL;
	m_max_free_energy = 0.0;
	m_tr_cost = 0.0;
	m_ca_cost = 0.0;
	m_error_rate = 0.0;
	m_protein_structure_ID = -1;
	m_last_struct_id = 0;
	m_last_free_energy = 0;
	m_last_sensitivity = 0;
	m_last_sensitivity_no_stop = 0;
	m_accuracy_weight = 1;
	m_error_weight = 1;
}

void ErrorproneTranslation::init( ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight )
{
	m_protein_folder = protein_folder;
	m_protein_length = m_protein_folder->getProteinLength();
	m_residue_sequence = new int[m_protein_length];

	m_max_free_energy = max_free_energy;
	m_tr_cost = tr_cost;
	m_ca_cost = ca_cost;
	m_error_rate = error_rate;
	m_accuracy_weight = accuracy_weight;
	m_error_weight = error_weight;
	m_protein_structure_ID = protein_structure_ID;

	// This weight matrix contains the probability that a given codon will be mistranslated as a given amino acid,
	// including a stop (at pos 20).
	m_weight_matrix.resize( 64 );

	//This weight matrix contains the sorted cumulative probability that
	//a given codon will be mistranslated as a given amino acid, including a stop.
	m_cum_weight_matrix.resize(64);

	for (int c=0; c<64; c++)
	{
		m_weight_matrix[c].resize( 21 ); // 20 aa's plus stop
		m_cum_weight_matrix[c].resize(21, pair<double,int>(0.0,-2));
	}
	vector<vector<double> > codon_to_codon_probabilities(64);
	for ( int c = 0; c<64; c++ )
	{
		codon_to_codon_probabilities[c].resize(64);
		double total_weight = 0.0;
		// now loop over all target codons
		// compute normalization
		for ( int tc = 0; tc<64; tc++ )
		{
			total_weight += calcWeight(c, tc);
		}
		// compute probabilities of transitions to each codon
		for ( int tc = 0; tc<64; tc++ )
		{
			double w = calcWeight(c, tc);
			codon_to_codon_probabilities[c][tc] = w/total_weight;
		}
	}
	// Now compute the codon-to-aa probabilities
	for ( int c=0; c<64; c++ ) {
		for (int aa=0; aa<21; aa++)	{
			if (aa<20) {
				m_cum_weight_matrix[c][aa].second = aa;
			}
			else {
				m_cum_weight_matrix[c][aa].second = -1;
			}
		}
		for (int tc=0; tc<64; tc++)	{
			int aa = GeneticCodeUtil::geneticCode[tc];
			if (aa < 0)	{
				aa = 20;
			}

			m_weight_matrix[c][aa] += codon_to_codon_probabilities[c][tc];
			m_cum_weight_matrix[c][aa].first += codon_to_codon_probabilities[c][tc];
		}
		// Sort the weights in preparation for computing ordered cumulative probabilities
		vector<pair<double,int> >& v = m_cum_weight_matrix[c];
		sort(v.begin(), v.end(), greater_pair_first());
		// Calculate the cumulative probabilities.
		double tot = 0.0;
		for (unsigned int i=0; i<m_cum_weight_matrix[c].size(); i++) {
			pair<double, int>&p = m_cum_weight_matrix[c][i];
			p.first += tot;
			tot = p.first;
		}
	}


	//    print the result for good measure
//	 for ( int c=0; c<64; c++ )
//	 {
//		 CodonUtil::printCodon( cout, c );
//		 for ( int a=0; a<=20; a++ )
//			 cout << " " << m_weight_matrix[c][a];
//		 cout << endl;
//	 }
}

ErrorproneTranslation::~ErrorproneTranslation()
{
	delete [] m_residue_sequence;
}


bool ErrorproneTranslation::sequenceFolds(Protein& p)
{
	// test if residue sequence folds into correct structure and has correct free energy
	pair<int, double> fold_data = p.fold(*m_protein_folder);
	
	m_last_free_energy = fold_data.second;
	m_last_struct_id = fold_data.first;

	if ( m_last_free_energy > m_max_free_energy )
		return false;  // free energy above cutoff

	if ( m_last_struct_id != m_protein_structure_ID )
		return false; // sequence folds into the wrong structure

	return true;
}


double ErrorproneTranslation::calcWeight( int co, int ct )
{
	if ( co == ct )
		return 0;

	double stop_factor = 1;
	if ( GeneticCodeUtil::geneticCode[ct] < 0 ) // does target codon code for stop?
		stop_factor = 0.33; // then probability is reduced to approx. one third.

	int lo1, lo2, lo3;
	int lt1, lt2, lt3;

	CodonUtil::codonToLetters( lo1, lo2, lo3, co );
	CodonUtil::codonToLetters( lt1, lt2, lt3, ct );

	if ( lo1 != lt1 && lo2 == lt2 && lo3 == lt3 )
	{
		if ( CodonUtil::isTransition( lo1, lt1 ) )
			return 1.*stop_factor;
		else
			return 0.5*stop_factor;
	}

	if ( lo2 != lt2 && lo1 == lt1 && lo3 == lt3 )
	{
		if ( CodonUtil::isTransition( lo2, lt2 ) )
			return 0.5*stop_factor;
		else
			return 0.1*stop_factor;
	}

	if ( lo3 != lt3 && lo1 == lt1 && lo2 == lt2 )
	{
		if ( CodonUtil::isTransition( lo3, lt3 ) )
			return 1.*stop_factor;
		else
			return 1.*stop_factor;
	}

	return 0;
}

/**
 * Assays protein encoded by the gene g for folding.  Calls sequenceFolds(), which may
 * be overridden in interesting ways.
 */
bool ErrorproneTranslation::getFolded( const Gene &g ) {
	// first test if protein translates at all, and if so, if it folds
	bool res = g.encodesFullLength();
	if ( res ) {
		Protein p = g.translate();
		res = sequenceFolds(p);
	}
	return res;
}

double ErrorproneTranslation::getFitness( const Gene &g ) {
	//cout << &g << tab << "begin getFitness" << endl;
	double fitness = 1.0;
	if ( m_tr_cost > 0 ) {
		double ffold, frob, facc, ftrunc;
		calcOutcomes(g, facc, frob, ftrunc, ffold);
		if (ffold > 0.0 ) {
			fitness = exp( - m_tr_cost * (1-ffold) / ffold );
		}
		else {
			fitness = 0.0;
		}
	}
	//cout << &g << tab << "end getFitness" << endl;
	return fitness;
}

void ErrorproneTranslation::setTargetAccuracyOfRandomGenes(const Gene& seed_genotype, const double facc, const int num_equil, const int num_rand) {
	// sets m_error_weight and m_accuracy_weight
	getWeightsForTargetAccuracy(seed_genotype, facc, m_error_rate, m_accuracy_weight, m_error_weight, num_equil, num_rand);
}

double ErrorproneTranslation::estimateErrorRateFromAccuracy(const double facc, const double accuracy_weight, const double error_weight) const {
	double error_rate_estimate = (1.0 - pow(facc, 1.0/m_protein_length))*(error_weight/accuracy_weight);
	return error_rate_estimate;
}

void ErrorproneTranslation::getWeightsForTargetAccuracy(const Gene& seed_genotype, const double target_accuracy, double& error_rate, 
														double& accuracy_weight, double& error_weight, const int num_equil, const int num_rand) {
	// Use for this is: we will specify an error rate for translation to match biological observations.  However, that error rate
	// cannot always hold when accuracy can change.  Therefore we would like to set up the system so that
	// a random folded sequence (e.g. w/ codons and aa usage constrained by folding only, no translational selection)
	// has the desired error rate.  Thus we will enumerate a large set of genes enconding folded sequences and determine what the
	// average sequence weight (m_error_weight) and sequence weight controlling for accuracy (m_accuracy_weight) of those genes are.
	Accumulator random_accuracy_weight; // running average
	Accumulator random_error_weight; // running average

	Gene g(seed_genotype);
	if (!getFolded(g)) {
		cerr << "# ERROR: getRandomWeights() requires a folded seed_genotype." << endl;
		return;
	}

	// Neutrally evolve for tot_equil steps, then for tot_rand steps, recording mutation weights
	int nrand=0, nequil=0;
	while ( nrand < num_rand ) {
		double rand = myRand();
		int randpos = (int)(rand * m_protein_length+1)-1;
		// go through all possible point mutations
		int from_codon = g[randpos];
		int to_codon = from_codon;
		do {
			rand = myRand();
			to_codon = (int)(64*rand+1)-1;
		} while (to_codon == from_codon);

		g[randpos] = to_codon;
		if (getFolded(g)) {
			nequil++;
			if (nequil > num_equil) {
				// Compute the weight
				double acc_weight_total = 0.0;
				double err_weight_total = 0.0;
				for ( unsigned int i=0; i<g.codonLength(); i++ ) {
					double last_event_prob = 0.0;

					// Compute probability of synonymous error
					double p_synonymous = 0.0;
					for (unsigned int outcome=0; outcome<m_cum_weight_matrix[g[i]].size(); outcome++) {
						pair<double, int> event = m_cum_weight_matrix[g[i]][outcome];
						// events record cumulative probabilities, so get the probability.
						double event_prob = event.first - last_event_prob;
						if (event.second == GeneticCodeUtil::geneticCode[g[i]]) {
							// Synonymous error
							p_synonymous += event_prob;
						}
						last_event_prob = event_prob;
					}
					double site_weight = (1.0 + m_codon_cost[g[i]]*(m_ca_cost-1));
					acc_weight_total += site_weight * (1-p_synonymous);
					err_weight_total += site_weight;
				}
				random_accuracy_weight += acc_weight_total;
				random_error_weight += err_weight_total;
				nrand++;
			}
		}
		else {
			g[randpos] = from_codon;
		}
	}

	accuracy_weight = random_accuracy_weight.value();
	error_weight = random_error_weight.value();
	error_rate = estimateErrorRateFromAccuracy(target_accuracy, accuracy_weight, error_weight);
}


/*
 * Computes the expected value of various translational outcomes from the gene g.
 */
double ErrorproneTranslation::calcOutcomes( const Gene &g, double &frac_accurate, double &frac_robust, double &frac_truncated, double &frac_folded ) {
	Protein prot = g.translate();
	// Assess folding of native sequence.
	bool native_seq_folds = sequenceFolds(prot);

	if (!native_seq_folds) {
		frac_accurate = pow((1-m_error_rate*m_accuracy_weight/m_error_weight), (double)m_protein_length);
		frac_robust = 0.0;
		frac_folded = 0.0;
		frac_truncated = 0.0;  // Not true...is there a fast way to estimate this (will rarely be used for anything)?
		return 0.0;
	}

	// Initialize our target properties
	double p_acc = 1.0;
	double p_fold = 1.0;
	double p_notrunc = 1.0;
	double inv_site_error_weight = 1.0/(m_error_weight/g.codonLength());
	int codon = -2;
	// Iterate over all sites
	for ( int i=0; i<prot.length(); i++ ) {
		const int old_res = prot[i];
		codon = g[i];
		// Probability of an event at this site
		double site_weight = 1.0 + m_codon_cost[codon]*(m_ca_cost-1);
		double p_codon_error = m_error_rate * site_weight * inv_site_error_weight;
		// Now go through possible translation outcomes at this codon.
		double last_prob = 0.0;
		double p_error_folds = 0.0;
		double p_synonymous = 0.0;
		double p_trunc_site = 0.0;
		for (unsigned int noutcome=0; noutcome<21; noutcome++) {
			// The first member of the pair is a cumulative probability; the second is
			// the residue resulting from the error.
			pair<double, int> p = m_cum_weight_matrix[codon][noutcome];
			prot[i] = p.second;
			// Probability of this particular event given an error at this site
			double p_outcome = (p.first-last_prob);
			if (p_outcome < 1e-6) { // Ignore really low-probability events
				break;
			}
			bool trunc = (prot[i] < 0); // truncation error.
			if (trunc) {
				p_trunc_site += p_outcome;
			}
			if (prot[i] != old_res) { // Change from wildtype sequence
				if (!trunc && sequenceFolds(prot)) { // Truncated means misfolded
					p_error_folds += p_outcome;
				}
			}
			else { // No change from wildtype sequence
				p_synonymous += p_outcome;
			}
			last_prob += p_outcome;
		}
		// to be accurate, either have no error, or have a synonymous error
		// p_acc = prod_i (1 - p_codon_error(i)) + p_codon_error(i) * p_synonymous(i)
		p_acc  *= (1.0 - p_codon_error * (1 - p_synonymous));
		// to be folded, either have no error, a synonymous error, or a nonsynonymous error that doesn't disrupt folding
		// p_fold = prod_i (1 - p_codon_error(i)) + p_codon_error(i) * p_synonymous(i) + p_codon_error(i) * (1 - p_synonymous(i)) * p_error_folds(i))
		// this is an approximation, and encodes the assumption that tolerated errors are tolerated together.
		p_fold *= (1.0 - p_codon_error * (1 - p_synonymous - p_error_folds));
		// to be untruncated, have no truncation
		// p_notrunc = prod_i (1- p_codon_error(i) * p_trunc(i))
		p_notrunc *= (1.0 - p_codon_error * p_trunc_site);
		prot[i] = old_res;
	}

	// Compute fractions of interest.
	frac_accurate = p_acc;   // No error at protein level
	frac_folded   = p_fold;  // Folded
	frac_robust   = 0.0;     // Folded despite mistranslation
	if ((1-frac_accurate)>0) {
		frac_robust = (frac_folded - frac_accurate)/(1-frac_accurate);
	}
	frac_truncated = 1.0 - p_notrunc;

	// Return fitness
	return exp( - m_tr_cost * (1-frac_folded) / frac_folded );
}

/**
 * Compute the actual fraction functional by translating num_to_fold proteins using
 * the error spectrum of this translator.
 * This function assumes the native protein folds.
 */
double ErrorproneTranslation::countOutcomes( const Gene &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded) {
	if (num_to_fold <= 0) {
		num_truncated = 0;
		num_accurate = 0;
		num_folded = 0;
		num_robust = 0;
		return (getFolded(g) ? 1.0 : 0.0);
	}
	bool native_seq_folds = getFolded(g);

	Translator t(m_error_rate);
	Protein p = g.translate();
	int numFolded = 0;
	int numTrunc = 0;
	int numMistranslated = 0;
	for (int i=0; i<num_to_fold; i++) {
		bool truncated = false;
		int numErrors = t.translateRelativeWeighted(g, p, m_error_weight, m_cum_weight_matrix, m_codon_cost, m_ca_cost, truncated);
		if (truncated) {
			numTrunc++;
		}
		if (numErrors > 0) {
			numMistranslated++;
			if (!truncated) {
				if (sequenceFolds(p)) {
					numFolded++;
				}
			}
		}
		else {
			if (native_seq_folds) {
				numFolded++;
			}
		}
	}
	//cout << "ff: " << numNotTrunc << " not trunc, " << numFolded << " folded of " << num_to_fold << endl;
	num_truncated = numTrunc;
	num_accurate = num_to_fold-numMistranslated;
	num_folded = numFolded;
	num_robust = numFolded - (native_seq_folds ? 1 : 0)*num_accurate;
	double frac_folded = numFolded/(double)num_to_fold;

	// return fitness
	return exp( - m_tr_cost * (1-frac_folded) / frac_folded );
}


/**
 * Compute the actual fraction functional by translating num_to_fold proteins using
 * the error spectrum of this translator.
 * This function assumes the native protein folds.
 */
void ErrorproneTranslation::stabilityOutcomes( const Gene &g, const int num_to_fold, vector<double>& ddgs ) {
	getFolded(g);

	Translator t(m_error_rate);
	Protein p = g.translate();
	for (int i=0; i<num_to_fold;) {
		bool truncated = false;
		int numErrors = t.translateRelativeWeighted(g, p, m_error_weight, m_cum_weight_matrix, m_codon_cost, m_ca_cost, truncated);
		// Only record mistranslations that are folded into the correct structure
		if (numErrors>0 && !truncated) {
			pair<int,double> fold_data = p.fold(*m_protein_folder);
			if (fold_data.first == m_protein_structure_ID) {
				ddgs.push_back(fold_data.second);
				i++;
			}
		}
	}
}


void ErrorproneTranslation::setCodonCosts() {
	// This algorithm weights each codon by the alphabetical order of its
	// one-letter code.
	// Codon cost matrix
	double codon_costs[64] = {
		1., // AAA -> LYS 0
		0., // AAC -> ASN 1
		0., // AAG -> LYS 2
		1., // AAU -> ASN 3
		1., // ACA -> THR 4
		0., // ACC -> THR 5
		1., // ACG -> THR 6
		0., // ACU -> THR 7
		0., // AGA -> ARG 8
		1., // AGC -> SER 9
		1., // AGG -> ARG 10
		1., // AGU -> SER 11
		1., // AUA -> ILE 12
		0., // AUC -> ILE 13
		1., // AUG -> MET 14
		0., // AUU -> ILE 15
		0., // CAA -> GLN 16
		0., // CAC -> HIS 17
		1., // CAG -> GLN 18
		1., // CAU -> HIS 19
		0., // CCA -> PRO 20
		1., // CCC -> PRO 21
		1., // CCG -> PRO 22
		1., // CCU -> PRO 23
		1., // CGA -> ARG 24
		1., // CGC -> ARG 25
		1., // CGG -> ARG 26
		1., // CGU -> ARG 27
		1., // CUA -> LEU 28
		1., // CUC -> LEU 29
		1., // CUG -> LEU 30
		1., // CUU -> LEU 31
		0., // GAA -> GLU 32
		0., // GAC -> ASP 33
		1., // GAG -> GLU 34
		1., // GAU -> ASP 35
		1., // GCA -> ALA 36
		1., // GCC -> ALA 37
		1., // GCG -> ALA 38
		0., // GCU -> ALA 39
		1., // GGA -> GLY 40
		1., // GGC -> GLY 41
		1., // GGG -> GLY 42
		0., // GGU -> GLY 43
		1., // GUA -> VAL 44
		0., // GUC -> VAL 45
		1., // GUG -> VAL 46
		0., // GUU -> VAL 47
		1., // UAA -> STOP48
		0., // UAC -> TYR 49
		1., // UAG -> STOP50
		1., // UAU -> TYR 51
		1., // UCA -> SER 52
		0., // UCC -> SER 53
		1., // UCG -> SER 54
		0., // UCU -> SER 55
		1., // UGA -> STOP56
		1., // UGC -> CYS 57
		1., // UGG -> TRP 58
		0., // UGU -> CYS 59
		1., // UUA -> LEU 60
		0., // UUC -> PHE 61
		0., // UUG -> LEU 62
		1.  // UUU -> PHE 63
	};
	// Reset codon costs
	for (int codon=0; codon<64; codon++) {
		m_codon_cost[codon] = codon_costs[codon];
	}
	for (int codon=0; codon<64; codon++) {
		int aa = GeneticCodeUtil::geneticCode[codon]+1;
		char c = *(GeneticCodeUtil::residueLetters[aa]);
		if (c == '*') {
			c = (char)100;
		}
		float weight = (float)((int)c-(int)'A'+1.0)/20.0;  // scaled from 0 to 1.
		m_codon_cost[codon] += weight;
		/*if (m_codon_cost[codon] == 0) {
			m_codon_cost[codon] = m_codon_cost[codon] + weight - 0.5;
		}
		else {
			m_codon_cost[codon] = m_codon_cost[codon] + weight - 0.5;
		}*/
	}
}

void ErrorproneTranslation::printCodonCosts(ostream& os) {
	for (int codon=0; codon<64; codon++) {
		os << "# ";
		CodonUtil::printCodon(cout, codon);
		os << "\t" << GeneticCodeUtil::residueLetters[GeneticCodeUtil::geneticCode[codon]+1] << "\t" << m_codon_cost[codon] << endl;
	}

}

double* ErrorproneTranslation::getCodonCosts() const { return ErrorproneTranslation::m_codon_cost; }

vector<vector<pair<double, int> > > ErrorproneTranslation::getTranslationWeights() const { return m_cum_weight_matrix; }

vector<bool> ErrorproneTranslation::getOptimalCodons(bool print_report) const {
	// Compute fraction accurately translated (facc) for each codon.  Optimal codons
	// have facc significantly higher than average for the family.
	int max_reps = 1000;
	Translator t(0.5); // Translator makes errors half the time.
	Gene g(3);
	Protein p(1);
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
				int n_errors = t.translateRelativeWeighted(g, p, m_error_weight/m_protein_length, m_cum_weight_matrix, m_codon_cost, m_ca_cost, truncated);
				// Increment totals
				family_accuracies[aa].second++;
				codon_accuracies[codon].second++;
				if (n_errors==0) {  // Properly translated
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


	// Now designate optimal any codons with accuracies significantly higher than the family average
	// By significantly, we will mean that the fold-improvement in accuracy over the rest of the
	// family is within 3-fold of the m_ca_cost factor.


	if (m_ca_cost > 1) {
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
				double facc = (tacc.first-acc.first)/(double)(tacc.second-acc.second);
				double myfacc = acc.first/(double)acc.second;
				if ((1-facc)/(1-myfacc) > m_ca_cost/3) { // Arbitrary fold cutoff.
					is_optimal[codon] = true;
				}
			}
		}
	}
	else {
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
				double facc = (tacc.first-acc.first)/(double)(tacc.second-acc.second);

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
	}


	if (print_report) {
		int num_optimal_codons = 0;
		for (int codon=0; codon<64; codon++) {
			int aa = GeneticCodeUtil::geneticCode[codon];
			if (is_optimal[codon]) {
				int tot = family_accuracies[aa].second - codon_accuracies[codon].second;
				int acc = family_accuracies[aa].first - codon_accuracies[codon].first;
				int myacc = codon_accuracies[codon].first;
				cout << "# Opt for " << GeneticCodeUtil::residueLetters[aa+1] << " = ";
				CodonUtil::printCodon(cout, codon);
				double myfacc = myacc/(double)max_reps;
				double facc = acc/(double)tot;
				cout << setprecision(3) << ", err = " << (1-myfacc) << " vs. family acc = " << facc << "\t(" << setprecision(4) << ((1-facc)/(1-myfacc)) << "-fold\tfewer errors)" << endl;
				num_optimal_codons++;
			}
		}
		cout << "# Found " << num_optimal_codons << " optimal codons." << endl;
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





// Codon cost matrix
double ErrorproneTranslation::m_codon_cost[64] = {
		1., // AAA -> LYS 0
		0., // AAC -> ASN 1
		0., // AAG -> LYS 2
		1., // AAU -> ASN 3
		1., // ACA -> THR 4
		0., // ACC -> THR 5
		1., // ACG -> THR 6
		0., // ACU -> THR 7
		0., // AGA -> ARG 8
		1., // AGC -> SER 9
		1., // AGG -> ARG 10
		1., // AGU -> SER 11
		1., // AUA -> ILE 12
		0., // AUC -> ILE 13
		1., // AUG -> MET 14
		0., // AUU -> ILE 15
		0., // CAA -> GLN 16
		0., // CAC -> HIS 17
		1., // CAG -> GLN 18
		1., // CAU -> HIS 19
		0., // CCA -> PRO 20
		1., // CCC -> PRO 21
		1., // CCG -> PRO 22
		1., // CCU -> PRO 23
		1., // CGA -> ARG 24
		1., // CGC -> ARG 25
		1., // CGG -> ARG 26
		1., // CGU -> ARG 27
		1., // CUA -> LEU 28
		1., // CUC -> LEU 29
		1., // CUG -> LEU 30
		1., // CUU -> LEU 31
		0., // GAA -> GLU 32
		0., // GAC -> ASP 33
		1., // GAG -> GLU 34
		1., // GAU -> ASP 35
		1., // GCA -> ALA 36
		1., // GCC -> ALA 37
		1., // GCG -> ALA 38
		0., // GCU -> ALA 39
		1., // GGA -> GLY 40
		1., // GGC -> GLY 41
		1., // GGG -> GLY 42
		0., // GGU -> GLY 43
		1., // GUA -> VAL 44
		0., // GUC -> VAL 45
		1., // GUG -> VAL 46
		0., // GUU -> VAL 47
		1., // UAA -> STOP48
		0., // UAC -> TYR 49
		1., // UAG -> STOP50
		1., // UAU -> TYR 51
		1., // UCA -> SER 52
		0., // UCC -> SER 53
		1., // UCG -> SER 54
		0., // UCU -> SER 55
		1., // UGA -> STOP56
		1., // UGC -> CYS 57
		1., // UGG -> TRP 58
		0., // UGU -> CYS 59
		1., // UUA -> LEU 60
		0., // UUC -> PHE 61
		0., // UUG -> LEU 62
		1.  // UUU -> PHE 63
};




AccuracyOnlyTranslation::AccuracyOnlyTranslation() : m_target_sequence(0) {
}

AccuracyOnlyTranslation::~AccuracyOnlyTranslation() {
}

void AccuracyOnlyTranslation::init( ProteinFolder *protein_folder, const int structure_id, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight ) {
	m_target_sequence = Protein(protein_folder->getProteinLength());
	//t.translateErrorFree( target_gene_sequence, m_target_sequence );
	//m_last_free_energy = protein_folder->foldProtein( m_target_sequence );
	//m_last_struct_id = protein_folder->getLastFoldedProteinStructureID();

	// Superclass init for the rest.
	ErrorproneTranslation::init( protein_folder, structure_id, max_free_energy, tr_cost, ca_cost, error_rate, accuracy_weight, error_weight );
}


bool AccuracyOnlyTranslation::sequenceFolds(Protein& p) {
	// test if residue sequence is the same as the initialized sequence.
	//cout << "AOT::sF" << endl;
	return m_target_sequence == p;
}

bool AccuracyOnlyTranslation::getFolded(const Gene &g) {
	bool res = g.encodesFullLength();
	// first test if protein translates at all, and if so, if it folds
	if ( res ) {
		Protein p = g.translate();
		res = ErrorproneTranslation::sequenceFolds(p);
	}
	return res;
}

double AccuracyOnlyTranslation::getFitness( const Gene &g ) {
	// test if protein translates at all, and if so, if it folds
	if (!getFolded(g)) {
		return 0.0;
	}

	//cout << "seq folded" << endl;

	if ( m_tr_cost > 0 )
	{
		double ffold, frob, facc, ftrunc;
		calcOutcomes(g, facc, frob, ftrunc, ffold);
		return exp( - m_tr_cost * (1-ffold)/ ffold );
	}
	else return 1.;
}

RobustnessOnlyTranslation::RobustnessOnlyTranslation() {
}

RobustnessOnlyTranslation::~RobustnessOnlyTranslation() {
}

void RobustnessOnlyTranslation::init(ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate ) {

	double length = (double)protein_folder->getProteinLength();
	ErrorproneTranslation::init( protein_folder, protein_structure_ID, max_free_energy, tr_cost, ca_cost, error_rate, length, length );
	// Fixed fraction mistranslated based on error rate
	m_fraction_accurate = 1.0 - pow((1.0 - error_rate), length);
}

double RobustnessOnlyTranslation::getFitness( const Gene &g )
{
	double fitness = 0.0;
	// test if protein translates at all, and if so, if it folds
	if ( g.encodesFullLength() ) {
		Protein p = g.translate();
		if ( ErrorproneTranslation::sequenceFolds(p) ) {
			if ( m_tr_cost > 0 ) {
				// Actual fraction folded will be (1-m_fraction_mistranslated) [all fold] + m_fraction_mistranslated*nu [neutral point mutations]
				double nu = GeneUtil::calcNeutrality(*m_protein_folder, p, m_max_free_energy);
				double ffold = m_fraction_accurate + nu*(1 - m_fraction_accurate);
				double fitness = exp( - m_tr_cost * (1.0 - ffold) / ffold );
			}
			else {
				fitness = 1.0;
			}
		}
	}
	return fitness;
}

/*
 * Computes the expected value of various translational outcomes from the gene g.
 */
double RobustnessOnlyTranslation::calcOutcomes( const Gene &g, double &frac_accurate, double &frac_robust, double &frac_truncated, double &frac_folded ) {
	Protein prot = g.translate();
	// Assess folding of native sequence.
	bool native_seq_folds = sequenceFolds(prot);

	if (!native_seq_folds) {
		frac_accurate = m_fraction_accurate;
		frac_robust = 0.0;
		frac_folded = 0.0;
		frac_truncated = 0.0;  // Not true...is there a fast way to estimate this (will rarely be used for anything)?
		return 0.0;
	}

	double nu = GeneUtil::calcNeutrality(*m_protein_folder, prot, m_max_free_energy);
	double ffold = m_fraction_accurate + nu*(1 - m_fraction_accurate);
	double fitness = exp( - m_tr_cost * (1.0 - ffold) / ffold );

	// Compute fractions of interest.
	frac_accurate = m_fraction_accurate;   // No error at protein level
	frac_folded   = ffold;  // Folded
	frac_robust   = 0.0;     // Folded despite mistranslation
	if ((1-frac_accurate)>0) {
		frac_robust = (frac_folded - frac_accurate)/(1-frac_accurate);
	}
	frac_truncated = 0.0;

	// Return fitness
	return fitness;
}

/**
 * Compute the actual fraction functional by translating num_to_fold proteins using
 * the error spectrum of this translator.
 * This is stubbed out in RobustnessOnlyTranslation.
 */
double RobustnessOnlyTranslation::countOutcomes( const Gene &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded) {
	num_truncated = 0;
	num_accurate = 0;
	num_folded = 0;
	num_robust = 0;
	return 0.0;
}

NeutralFitness::NeutralFitness()
{}


NeutralFitness::~NeutralFitness()
{}


double fixation_probability(int N, double s) {
	double res = 0.0;
	double eps = 1e-10;
	if (s > -eps && s < eps) {
		// completely neutral evolution
		res = 1.0/N;
	}
	else if (s <= -1.0) {
		// mutant has zero fitness.
		res = 0.0;
	}
	else {
		// Sella and Hirsh probability of fixation
		double delta = log(1.0+s);
		double num = 1.0 - exp(-2.0 * delta);
		double denom = 1.0 - exp(-2.0 * N * delta);
		if (denom != 0.0) {
			res = num/denom;
		}
	}
	return res;
}

