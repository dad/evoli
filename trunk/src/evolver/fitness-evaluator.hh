/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


#ifndef FITNESS_EVALUATOR_HH
#define FITNESS_EVALUATOR_HH

#include "protein.hh"
#include "folder.hh"

class FitnessEvaluator {
private:
	FitnessEvaluator( const FitnessEvaluator& );
	FitnessEvaluator& operator=( const FitnessEvaluator& );

public:
	FitnessEvaluator();
	virtual ~FitnessEvaluator();

	virtual double getFitness( const Gene& ) = 0;
	virtual double getFitness( const Protein& ) = 0;
};


class ProteinFreeEnergyFitness : public FitnessEvaluator {
private:
	ProteinFreeEnergyFitness();

	Folder *m_protein_folder;
public:
	ProteinFreeEnergyFitness( Folder *protein_folder );
	virtual ~ProteinFreeEnergyFitness();

	double getFitness( const Gene& g );
	double getFitness( const Protein& p );
};


class ProteinStructureFitness : public FitnessEvaluator {
private:
	ProteinStructureFitness();

	Folder *m_protein_folder;
	int m_protein_structure_ID;
	double m_max_free_energy;

public:
	ProteinStructureFitness( Folder *protein_folder, int protein_structure_ID, double max_free_energy );
	virtual ~ProteinStructureFitness();

	void setFreeEnergyCutoff(double cutoff) { m_max_free_energy = cutoff; }
	double getFreeEnergyCutoff() const { return m_max_free_energy; }

	double getFitness( const Gene& g );
	double getFitness( const Protein& p );
};


class NeutralFitness : public FitnessEvaluator {
public:
	NeutralFitness();
	virtual ~NeutralFitness();

	double getFitness( const Gene& ) {
		return 1;
	}
	double getFitness( const Protein& ) {
		return 1;
	}
};


/** \brief A \ref FitnessEvaluator in which the fitness of a gene sequence is determined from the amount of misfolded protein generated under error-prone mistranslation.
*/
class ErrorproneTranslation : public FitnessEvaluator
{
protected:
	Folder *m_protein_folder;
	double m_max_free_energy;
	double m_tr_cost;
	double m_ca_cost;
	double m_error_rate;

	int m_protein_structure_ID;
	int m_protein_length;

	double m_accuracy_weight;
	double m_error_weight;

	/**
	This weight matrix contains the probability that a given codon will
	be mistranslated as a given amino acid, including a stop (at pos 20).
	*/
	vector<vector<double> > m_weight_matrix;

	/**
	This weight matrix contains the sorted cumulative probability that
	a given codon will be mistranslated as a given amino acid, including
	a stop.
	*/
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

	void init(Folder *protein_folder, const int length, const int protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );

	void changeStructure( const int structureID ) {
		m_protein_structure_ID = structureID;
	}

	double getFitness( const Gene& g );
	double getFitness( const Protein& p );

	/**
	Tests whether the protein encoded by the \ref Gene g folds
	correctly. Calls \ref sequenceFolds(), which may be overridden
	in interesting ways.
	@param g The \ref Gene to test.
	@return True if the protein folds correctly, False otherwise.
	*/
	virtual bool getFolded( const Gene& g );

	/**
	Set fitness costs for codons. The implemented cost scheme weights
	each codon by the alphabetical order of its one-letter code.
	*/
	virtual void setCodonCosts();

	/**
	Writes the currently defined codon costs to a stream.
	*/
	void printCodonCosts(ostream& os);

	/**
	Returns a pointer to an array holding the codon costs.
	\warning Usage of this pointer is potentially unsafe, because
	the pointer will become invalid upon destruction of the \ref ErrorproneTranslation class.
	*/
	double* getCodonCosts() const;
	vector<vector<pair<double, int> > > getTranslationWeights() const;

	/**
	@return A vector of bools indicating whether a given codon
	is optimal or not (True means codon is optimal). 
	*/
	vector<bool> getOptimalCodons(bool print_report = true) const;

	/**
	Compute the per-site error rate estimated to achieve the desired fraction accurate over a large set of random (unselected) genotypes encoding folded proteins.

	The use for this function is the following: We will specify an error
	rate for translation to match biological observations. However, that
	error rate cannot always hold when translation accuracy can change.
	Therefore, we would like to set up the system so that a randomly
	chosen sequence (e.g. w/ codons and aa usage constrained by
	folding only, no translational selection) has the desired error rate.
	Thus, we will enumerate a large set of genes encoding folded
	sequences and determine what the average sequence weight (returned in variable \ref error_weight) and sequence weight controlling for accuracy (returned in variable \ref accuracy_weight) of those genes are.
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
	Compute the per-site error rate estimated to achieve the desired fraction accurate. It is recommended to call the function \ref setRandomWeights(...) first.
	 */
	double estimateErrorRateFromAccuracy(const double base_fraction_accurate, const double accuracy_weight, const double error_weight) const;

	/**
	Compute the estimated fractions accurately translated, folded
	despite mistranslation, truncated and folded, using the error
	spectrum of this FitnessEvaluator.
	@return The estimated fitness.
	*/
	virtual double calcOutcomes( const Gene& g, double& frac_accurate, double& frac_robust, double& frac_truncated, double& frac_folded );

	/**
	Record the actual fractions accurately translated, folded despite
	mistranslation, truncated and folded, using the error spectrum of this
	FitnessEvaluator, by translating a specified number of proteins.
	@param num_to_fold The number of proteins that should be translated to determine the error spectrum.
	@return The fitness that would result from the generated set of proteins.
	 */
	virtual double countOutcomes(const Gene& g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);

	/**
	 * Record stabilities of mistranslated proteins.
	 */
	virtual void stabilityOutcomes( const Gene& g, const int num_to_fold, vector<double>& ddgs );

	/**
	 * Returns the translational error probability per codon.
	 */
	double getErrorRate() const {
		return m_error_rate;
	}


	/**
	@return A pointer to the currently active protein folder.
	*/
	Folder* getFolder() {
		return m_protein_folder;
	}

	/**
	* Codon cost matrix.
	* We represent each base by two binary digits:
	*   A: 00, C: 01, G: 10, U: 11.
	* Then, each codon is an integer between 0 and 63. For example,
	* GUG = 101110 = 46.
	*
	* The value at the individual positions gives the cost associated with the given codon.
	* A value of zero means no cost, while a value of one means maximal cost.
	**/
	static double m_codon_cost[64];
};

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which all fitness cost is due to mistranslation only.

This class implements a case of \ref ErrorproneTranslation  where translation is error-prone and any mistranslation to an incorrect amino acid leads to a consequent fitness cost, regardless of whether the translated protein would fold correctly or not.
 */

class AccuracyOnlyTranslation : public ErrorproneTranslation {
protected:
	Protein m_target_sequence;
	bool sequenceFolds(Protein& p);

public:
	AccuracyOnlyTranslation();
	virtual ~AccuracyOnlyTranslation();

	void init( Folder* protein_folder, const int length, const int target_structure_id, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	double getFitness( const Gene& g );
	bool getFolded( const Gene& g );
};

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which all fitness differences are due to differing protein robustness to translation errors.

This class implements a case of \ref ErrorproneTranslation where translation always yields a fixed percentage of mistranslated polypeptides, which may then fold properly.
*/

class RobustnessOnlyTranslation : public ErrorproneTranslation {
protected:
	double m_fraction_accurate;

public:
	RobustnessOnlyTranslation();
	virtual ~RobustnessOnlyTranslation();

	void init( Folder* protein_folder, const int length, const int protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate );
	double getFitness( const Gene& g );

	/**
	Compute the estimated fractions accurately translated, folded
	despite mistranslation, truncated and folded, using the error
	spectrum of this FitnessEvaluator.
	@return The estimated fitness.
	 */
	virtual double calcOutcomes( const Gene& g, double& frac_accurate, double& frac_robust, double& frac_truncated, double& frac_folded );

	/**
	\warning This function is currently not implemented and sets all values to zero!

	Record the actual fractions accurately translated, folded despite
	mistranslation, truncated and folded, using the error spectrum of this
	FitnessEvaluator, by translating a specified number of proteins.
	@param num_to_fold The number of proteins that should be translated to determine the error spectrum.
	@return The fitness that would result from the generated set of proteins.
	 */
	virtual double countOutcomes(const Gene& g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);
};

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which a minimum number of misfolded proteins are required before any fitness cost is incurred.
 */

class CutoffErrorproneTranslation : public ErrorproneTranslation {
protected:
	/**
	 * The multiplier to convert translational robustness cost into number of misfolded proteins.
	 **/
	double m_cost_constant;
	/**
	 * The minimum number of misfolded proteins required to get a toxic effect.
	 **/
	int m_toxicity_cutoff;
public:
	/**
	@param cost_constant The multiplier to convert translational robustness cost into number of misfolded proteins.
	@param toxicity_cutoff The minimum number of misfolded proteins required to get a toxic effect.
	*/
	CutoffErrorproneTranslation(double cost_constant, int toxicity_cutoff);
	
	virtual ~CutoffErrorproneTranslation();
	double getFitness( const Gene& g );
};

/**
 * Computes the probability of fixation of a mutant with fitness advantage s in a population of size N.
 *
 * @param s The fitness advantage relative to the wild type.
 * @param N The population size.
 **/
double fixation_probability(int N, double s);

#endif
