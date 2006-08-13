#ifndef FITNESS_EVALUATOR_HH
#define FITNESS_EVALUATOR_HH

#include "genotype.hh"

class ProteinFolder;

class FitnessEvaluator
{
private:

        FitnessEvaluator( const FitnessEvaluator & );
        FitnessEvaluator& operator=( const FitnessEvaluator & );
public:
        FitnessEvaluator();
        virtual ~FitnessEvaluator();

        virtual double getFitness( const Genotype & ) = 0;

};



class ProteinFreeEnergyFitness : public FitnessEvaluator
{
private:
        ProteinFolder *m_protein_folder;
        int m_protein_length;
        int *m_residue_sequence;

public:
        ProteinFreeEnergyFitness( ProteinFolder *protein_folder );
        ~ProteinFreeEnergyFitness();

        double getFitness( const Genotype &g );


};


class ProteinStructureFitness : public FitnessEvaluator
{
private:
        ProteinFolder *m_protein_folder;
        int m_protein_structure_ID;
        double m_max_free_energy;

        int m_protein_length;
        int *m_residue_sequence;

public:
        ProteinStructureFitness( ProteinFolder *protein_folder, int protein_structure_ID, double max_free_energy );
        ~ProteinStructureFitness();

        double getFitness( const Genotype &g );
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
	double calcSensitivity( const Genotype &g );
	virtual bool sequenceFolds();

	/**
	 * Compute the translational accuracy-related gene weights of a large set of random genotypes encoding folded proteins.
	 */
	void setRandomWeights(const Genotype& seed_genotype, const int num_equil=5000, const int num_rand=1000);

public:
	ErrorproneTranslation();
	virtual ~ErrorproneTranslation();

	void init(ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );

	void changeStructure( const int structureID ) {
		m_protein_structure_ID = structureID;
		m_last_free_energy = 999;
	}

	double getFitnessSensitivity( const Genotype &g );

	double getFitness( const Genotype &g );
	virtual bool getFolded( const Genotype &g );

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
	void getWeightsForTargetAccuracy(const Genotype& seed_genotype, const double target_accuracy, double& error_rate, double& accuracy_weight, double& error_weight, const int num_equil, const int num_rand);
	void setTargetAccuracyOfRandomGenes(const Genotype& seed_genotype, const double facc, const int num_equil, const int num_rand);
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
	virtual double calcOutcomes( const Genotype &g, double &frac_accurate, double &frac_robust, double& frac_truncated, double &frac_folded );

	/**
	 * Record the actual fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator, by translating num_to_fold proteins.
	 * Returns the estimated fitness that would result.
	 */
	virtual double countOutcomes(const Genotype &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);

	/**
	 * Record stabilities of mistranslated proteins.
	 */
	virtual void stabilityOutcomes( const Genotype &g, const int num_to_fold, vector<double>& ddgs );

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
	int* m_target_sequence;
	bool sequenceFolds();

public:
	AccuracyOnlyTranslation();
	virtual ~AccuracyOnlyTranslation();

	void init( ProteinFolder *protein_folder, const int target_structure_id, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	double getFitness( const Genotype &g );
	bool getFolded( const Genotype &g );
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
	double getFitness( const Genotype &g );
	/**
	 * Compute the estimated fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator.
	 * Returns the estimated fitness.
	 */
	virtual double calcOutcomes( const Genotype &g, double &frac_accurate, double &frac_robust, double& frac_truncated, double &frac_folded );

	/**
	 * Record the actual fractions accurately translated, folded despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator, by translating num_to_fold proteins.
	 * Returns the estimated fitness that would result.
	 */
	virtual double countOutcomes(const Genotype &g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);
};

/**
 * StabilityConstraint implements a case where translation is error-prone and any
 * mistranslation leading to a protein with stability above or below certain cutoffs incurs a cost.
 */

class StabilityConstraint : public ErrorproneTranslation {
protected:
	double m_min_free_energy;
	bool sequenceFolds();

public:
	StabilityConstraint();
	virtual ~StabilityConstraint();

	void init( ProteinFolder *protein_folder, const int protein_structure_ID, const double max_free_energy,
		const double min_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
};


class NeutralFitness : public FitnessEvaluator {
public:
	NeutralFitness();
	~NeutralFitness();

	double getFitness( const Genotype & ) {
			return 1;
	}
};

double fixation_probability(int N, double s);


class FitnessDensityEvaluator {
public:
	FitnessDensityEvaluator() {}
	~FitnessDensityEvaluator() {}
	double getFitnessDensity(const Genotype& g, FitnessEvaluator& fe, unsigned int pop_size) const;
	//double getFitnessDensityNonsyn(const Genotype& g, FitnessEvaluator& fe, unsigned int pop_size) const;
	//double getFitnessDensitySyn(const Genotype& g, FitnessEvaluator& fe, unsigned int pop_size) const;
};


#endif
