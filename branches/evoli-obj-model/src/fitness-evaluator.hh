#ifndef FITNESS_EVALUATOR_HH
#define FITNESS_EVALUATOR_HH

#include "protein.hh"

class ProteinFolder;

class FitnessEvaluator {
private:
	FitnessEvaluator( const FitnessEvaluator & );
	FitnessEvaluator& operator=( const FitnessEvaluator & );
public:
	FitnessEvaluator();
	virtual ~FitnessEvaluator();

	virtual double getFitness( const Gene & ) = 0;
	virtual double getFitness( Protein & ) = 0;
};


class ProteinFreeEnergyFitness : public FitnessEvaluator {
private:
	ProteinFolder *m_protein_folder;
public:
	ProteinFreeEnergyFitness( ProteinFolder *protein_folder );
	~ProteinFreeEnergyFitness();

	double getFitness( const Gene &g );
	double getFitness( Protein &p );
};


class ProteinStructureFitness : public FitnessEvaluator {
private:
	ProteinFolder *m_protein_folder;
	int m_protein_structure_ID;
	double m_max_free_energy;

public:
	ProteinStructureFitness( ProteinFolder *protein_folder, int protein_structure_ID, double max_free_energy );
	~ProteinStructureFitness();

	double getFitness( const Gene &g );
	double getFitness( Protein &p );
};


class NeutralFitness : public FitnessEvaluator {
public:
	NeutralFitness();
	~NeutralFitness();

	double getFitness( const Gene & ) {
		return 1;
	}
	double getFitness( Protein & ) {
		return 1;
	}
};

class Accumulator {
private:
	double m_running_sum;
	int m_count;
public:
	Accumulator() { m_running_sum=0.0; m_count=0; }
	Accumulator(double sum, int count) { m_running_sum = sum; m_count = count; }

	double value() const { return m_running_sum/m_count; }
	void add(double val) { m_running_sum += val; m_count++; }
	void reset() { m_running_sum = 0.0; m_count = 0; }

	int count() const { return m_count; }
	double sum() const { return m_running_sum; }

	operator double() { return value(); }
	void operator+=(double x) { add(x); }
};


class ErrorproneTranslation : public FitnessEvaluator
{
protected:
	ProteinFolder *m_protein_folder;
	double m_max_free_energy;
	double m_tr_cost;
	double m_ca_cost;
	double m_error_rate;

	int m_protein_structure_ID;
	int m_protein_length;
	int *m_residue_sequence;

	int m_last_struct_id;
	double m_last_free_energy;
	double m_last_sensitivity;
	double m_last_sensitivity_no_stop;

	double m_accuracy_weight;
	double m_error_weight;

	vector<vector<double> > m_weight_matrix;
	vector<vector<pair<double, int> > > m_cum_weight_matrix;

	double calcWeight( int co, int ct );
	virtual bool sequenceFolds(Protein& p);

	/**
	 * Compute the translational accuracy-related gene weights of a large set of random genotypes encoding folded proteins.
	 */
	void setRandomWeights(const Gene& seed_genotype, const int num_equil=5000, const int num_rand=1000);

public:
	ErrorproneTranslation();
	virtual ~ErrorproneTranslation();

	void init(ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );

	void changeStructure( const int structureID ) {
		m_protein_structure_ID = structureID;
		m_last_free_energy = 999;
	}

	double getFitness( const Gene &g );
	double getFitness( Protein& p ) { return getFitness( p.reverseTranslate() ); }
	virtual bool getFolded( const Gene &g );

	/**
	 * Set costs for codons.  For example, may want to have variation in codon costs
	 * between amino acids.
	 */
	virtual void setCodonCosts();
	void printCodonCosts(ostream& os);
	double* getCodonCosts() const;
	vector<vector<pair<double, int> > > getTranslationWeights() const;
	vector<bool> getOptimalCodons(bool print_report = true) const;

	/**
	 * Compute the per-site error rate estimated to achieve the desired fraction accurate over a large set of random (unselected) genotypes encoding folded proteins.
	 */
	void getWeightsForTargetAccuracy(const Gene& seed_genotype, const double target_accuracy, double& error_rate, 
									 double& accuracy_weight, double& error_weight, const int num_equil, const int num_rand);
	void setTargetAccuracyOfRandomGenes(const Gene& seed_genotype, const double facc, const int num_equil, const int num_rand);
	void setErrorRate(const double error_rate) { m_error_rate = error_rate; }
	void setErrorWeights(const double error_rate, const double accuracy_weight, const double error_weight) {
		m_error_rate = error_rate;
		m_accuracy_weight = accuracy_weight;
		m_error_weight = error_weight;
	}
	/**
	 * Compute the per-site error rate estimated to achieve the desired fraction accurate.
	 * Calling setRandomWeights(...) first is recommended.
	 */
	double estimateErrorRateFromAccuracy(const double base_fraction_accurate, const double accuracy_weight, const double error_weight) const;

	/**
	 * Compute the estimated fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator.
	 * Returns the estimated fitness.
	 */
	virtual double calcOutcomes( const Gene &g, double &frac_accurate, double &frac_robust, double& frac_truncated, double &frac_folded );

	/**
	 * Record the actual fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator, by translating num_to_fold proteins.
	 * Returns the estimated fitness that would result.
	 */
	virtual double countOutcomes(const Gene &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);

	/**
	 * Record stabilities of mistranslated proteins.
	 */
	virtual void stabilityOutcomes( const Gene &g, const int num_to_fold, vector<double>& ddgs );

	/**
	 * Returns the translational error probability per codon.
	 */
	double getErrorRate() const {
		return m_error_rate;
	}

	// these functions must be called after getFitness
	int getLastStructId() const	{
			return m_last_struct_id;
	}

	/*
	* Returns the free energy of the sequence that was most-recently evaluated.
	*/
	double getLastFreeEnergy() const {
			return m_last_free_energy;
	}

	/*
	* Returns the sensitivity of the sequence that was most-recently evaluated.
	*/
	double getLastSensitivity() const
	{
			return m_last_sensitivity;
	}

	/*
	* Returns the sensitivity, excluding contributions from premature stop codons,
	* of the sequence that was most-recently evaluated.
	*/
	double getLastSensitivityNoStop() const	{
			return m_last_sensitivity_no_stop;
	}

	ProteinFolder* getFolder() {
		return m_protein_folder;
	}

	/**
	* Codon cost matrix.
	* We represent each base by two binary digits:
	*   A: 00, C: 01, G: 10, U: 11
	* Then, each codon is an integer between 0 and 63. For example,
	* GUG = 101110 = 46.
	*
	* The value at the individual positions gives the cost associated with the given codon.
	* A value of zero means no cost, while a value of one means maximal cost.
	**/
	static double m_codon_cost[64];
};

/**
 * AccuracyOnlyTranslation implements a case where translation is error-prone and any
 * mistranslation to an incorrect amino acid leads to a consequent fitness cost.
 */

class AccuracyOnlyTranslation : public ErrorproneTranslation {
protected:
	Protein m_target_sequence;
	bool sequenceFolds(Protein& p);

public:
	AccuracyOnlyTranslation();
	virtual ~AccuracyOnlyTranslation();

	void init( ProteinFolder *protein_folder, const int target_structure_id, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	double getFitness( const Gene &g );
	bool getFolded( const Gene &g );
};

/**
 * RobustnessOnlyTranslation implements a case where translation always yields a
 * fixed percentage of mistranslated polypeptides, which may then fold properly.
 */

class RobustnessOnlyTranslation : public ErrorproneTranslation {
protected:
	double m_fraction_accurate;

public:
	RobustnessOnlyTranslation();
	virtual ~RobustnessOnlyTranslation();

	void init( ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate );
	double getFitness( const Gene &g );
	/**
	 * Compute the estimated fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator.
	 * Returns the estimated fitness.
	 */
	virtual double calcOutcomes( const Gene &g, double &frac_accurate, double &frac_robust, double& frac_truncated, double &frac_folded );

	/**
	 * Record the actual fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator, by translating num_to_fold proteins.
	 * Returns the estimated fitness that would result.
	 */
	virtual double countOutcomes(const Gene &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);
};

double fixation_probability(int N, double s);

#endif
