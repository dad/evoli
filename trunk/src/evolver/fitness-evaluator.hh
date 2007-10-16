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
#include "gene-util.hh"

class FitnessEvaluator {
private:
	FitnessEvaluator( const FitnessEvaluator& );
	FitnessEvaluator& operator=( const FitnessEvaluator& );

public:
	FitnessEvaluator();
	virtual ~FitnessEvaluator();

	virtual double getFitness( const CodingDNA& ) = 0;
	virtual double getFitness( const Protein& ) = 0;
};

class ProteinFreeEnergyFitness : public FitnessEvaluator {
private:
	ProteinFreeEnergyFitness();

	Folder *m_protein_folder;
public:
	ProteinFreeEnergyFitness( Folder *protein_folder );
	virtual ~ProteinFreeEnergyFitness();

	double getFitness( const CodingDNA& g );
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
	double getDeltaGCutoff() const { return m_max_free_energy; }

	double getFitness( const CodingDNA& g );
	double getFitness( const Protein& p );
};

/**
 * \brief A \ref FitnessEvaluator that assigns equal fitness to all sequences.
 */
class NeutralFitness : public FitnessEvaluator {
public:
	NeutralFitness();
	virtual ~NeutralFitness();

	double getFitness( const CodingDNA& ) {
		return 1;
	}
	double getFitness( const Protein& ) {
		return 1;
	}
};


/** \brief A \ref FitnessEvaluator in which the fitness of a gene sequence is determined from the amount of misfolded protein generated under error-prone mistranslation.

**/

class ErrorproneTranslation : public FitnessEvaluator
{
private:
	ErrorproneTranslation();

protected:
	Folder *m_protein_folder;
	double m_max_free_energy;
	double m_tr_cost;
	double m_ca_cost;
	double m_error_rate;

	StructureID m_protein_structure_ID;
	int m_protein_length;

	double m_accuracy_weight;
	double m_error_weight;

	/**
	This weight matrix contains the probability that a given codon will
	be mistranslated as a given amino acid, including a stop (at pos 20).
	*/
	//vector<vector<double> > m_weight_matrix;

	/**
	This weight matrix contains the sorted cumulative probability that
	a given codon will be mistranslated as a given amino acid, including
	a stop.
	*/
	vector<vector<pair<double, char> > > m_cum_weight_matrix;

	double calcWeight( int co, int ct );
	virtual bool sequenceFolds(Protein& p);

	/**
	 * Compute the translational accuracy-related gene weights of a large set of random genotypes encoding folded proteins.
	 */
	void setRandomWeights(const CodingDNA& seed_genotype, const int num_equil=5000, const int num_rand=1000);

	/**
	 * Build the weight matrices that will be used in translating genes.
	 *
	 * Note: must be called before getWeightsForTargetAccuracy is called.
	 **/
	void buildWeightMatrix();

public:

	/**
	 * \brief Create new ErrorproneTranslation object with pre-determined sequence weights.
	 *
	 * @param protein_folder An initialized \ref Folder object.
	 * @param protein_length Protein sequence length in amino acids.
	 * @param protein_structure_ID Structure identifier for the native structure all folded proteins must attain.
	 * @param max_free_energy Maximum free energy of folding for a folded protein.
	 * @param tr_cost Cost factor used to convert fraction of misfolded proteins f into fitness cost via fitness = exp{- tr_cost*f/(1-f).
	 * @param ca_cost Codon adapation cost.  This cost represents the average fold-decrease in codon accuracy for non-optimal codons relative to optimal synonymous codons.
	 * @param error_rate The per-codon base translational error rate.  Non-optimal codons will have a higher error rate, determined by ca_cost.
	 * @param accuracy_weight The reference average sequence weight for a large ensemble of folded proteins, excluding the possibility of synonymous errors.  @see getWeightsForTargetAccuracy.
	 * @param error_weight The reference average sequence weight for a large ensemble of folded proteins, accounting for both synonymous and nonsynonymous errors.  @see getWeightsForTargetAccuracy.
	 **/
	ErrorproneTranslation(Folder *protein_folder, const int protein_length, const StructureID protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );

	/**
	 * \brief Create new ErrorproneTranslation object which will self-determine sequence weights to achieve a target fraction accurate for randomly chosen genes.
	 *
	 * @param protein_folder An initialized \ref Folder object.
	 * @param protein_length Protein sequence length in amino acids.
	 * @param protein_structure_ID Structure identifier for the native structure all folded proteins must attain.
	 * @param max_free_energy Maximum free energy of folding for a folded protein.
	 * @param tr_cost Cost factor for
	 * @param ca_cost Codon adapation cost.  This cost represents the average fold-decrease in codon accuracy for non-optimal codons relative to optimal synonymous codons.
	 * @param target_fraction_accurate Desired probability that an average folded protein will be translated with out errors.
	 **/
	ErrorproneTranslation(Folder *protein_folder, const int protein_length, const StructureID protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double target_fraction_accurate );

	/**
	 * \brief Create new ErrorproneTranslation object from a preexisting object.
	 **/
	ErrorproneTranslation(const ErrorproneTranslation& ept) {
		m_protein_folder = ept.m_protein_folder;
		m_protein_length = ept.m_protein_length;
		m_protein_structure_ID = ept.m_protein_structure_ID;
		m_max_free_energy = ept.m_max_free_energy;
		m_ca_cost = ept.m_ca_cost;
		m_error_rate = ept.m_error_rate;
		m_accuracy_weight = ept.m_accuracy_weight;
		m_error_weight = ept.m_error_weight;
		buildWeightMatrix();
	}

	/**
	 * \brief Destroy this ErrorproneTranslation object.
	 **/
	virtual ~ErrorproneTranslation();

	/**
	 * \brief Change the target structure ID that defines a folded protein.
	 **/
	void changeStructure( const StructureID structure_ID ) {
		m_protein_structure_ID = structure_ID;
	}

	double getFitness( const CodingDNA& g );
	double getFitness( const Protein& p );

	/**
	Tests whether the protein encoded by the \ref CodingDNA g folds
	correctly. Calls \ref sequenceFolds(), which may be overridden
	in interesting ways.
	@param g The \ref CodingDNA to test.
	@return True if the protein folds correctly, False otherwise.
	*/
	virtual bool getFolded( const CodingDNA& g );

	/**
	Set fitness costs for codons. The implemented cost scheme weights
	each codon by the alphabetical order of its one-letter code.
	*/
	//virtual void setCodonCosts();

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
	vector<vector<pair<double, char> > > getTranslationWeights() const;

	/**
	@return A vector of bools indicating whether a given codon
	is optimal or not (True means codon is optimal).
	*/
	vector<bool> getOptimalCodons(bool print_report = true) const;

	/**
	Compute the per-site error rate estimated to achieve the desired fraction accurate over a large set of random (unselected) genotypes encoding folded proteins.

	Our goal is to specify an error rate for translation to match
	biological estimates, such as that 81% of proteins have at
	least one error on average (estimate from Drummond et al., PNAS
	102:14338-14343 (2005)).

	When codons have different accuracies, it is difficult to set
	an error rate for translation per site that achieves this
	overall target accuracy.  Some genes will be more or less
	likely to be mistranslated than others based on their codon
	composition, and particular sites within each gene will be
	similarly more or less likely to be mistranslated.

	In addition, some errors may be synonymous -- that is, a codon
	may be translated as if it were another codon, but the result
	may be that the same amino acid is used.  Such synonymous
	errors will generally not be reflected in biologically observed
	error rates.

	Our strategy is to characterize the propensity of a large set
	of random genes that encode folded proteins to be mistranslated
	based on their codon composition.  (Folded proteins are
	required because amino acid composition will influence codon
	composition.)  Each gene will be assigned a weight derived from
	its codon composition.  The average weight over this large,
	random gene set will then be used as a reference weight for
	subsequent translation of arbitrary (non-random) genes.

	The basic weight unit is the fold decrease in accuracy at a site
	relative to a site bearing a preferred codon,
	[1 + unpreferred<1|0>*(unpreferred_penalty-1)].  For example, if a
	site has a preferred codon, the site weight is 1; if a site has an
	unpreferred codon and the accuracy difference between preferred
	and unpreferred codons is 6, then the site weight is 6.

	The "error weight" for a gene is the sum of its site
	weights. CodingDNAs with higher error weights are more likely to be
	mistranslated, though in some cases the translation errors may be
	synonymous.

	The average error weight for a large set of random genes
	encoding folded proteins sets a reference point for translation
	of arbitrary genes, and is stored in error_weight.

	Correspondingly, the "accuracy weight" for a gene is the sum of
	each site weight multiplied by the probability of a
	nonsynonymous change given an error at that site.  The average
	accuracy weight again sets a reference point for translation of
	arbitrary genes, and is stored in accuracy_weight.

	Finally, given an accuracy weight and an error weight, we can
	set an error rate for translation that will, on average over a
	large random sample, yield the desired fraction of genes which
	are accurately translated, whether through avoidance of error
	or synonymous errors.  This error rate is returned in error_rate.

	@param num_equil Number of steps to evolve, preserving stable folding only, before recording weights.
	@param num_rand Number of steps (genes) over which average weights will be accumulated.
	@param target_fraction_accurate Desired fraction of random genes encoding folded proteins that will be translated without amino-acid errors.
	@param seed_genotype Starting gene sequence, which must encode a folded protein consistent with initialization values of this ErrorproneTranslation object.
	@param error_rate Stores computed error rate (see above).
	@param error_weight Stores computed error weight (see above).
	@param accuracy_weight Stores computed accuracy weight (see above).
	 */
	void getWeightsForTargetAccuracy(const CodingDNA& seed_genotype, const double target_fraction_accurate, double& error_rate,
									 double& accuracy_weight, double& error_weight, const int num_equil, const int num_rand);

	/**
	 * Given a [[DAD: is this used?]]
	 **/
	void setTargetAccuracyOfRandomGenes(const CodingDNA& seed_genotype, const double target_fraction_accurate, const int num_equil, const int num_rand);

	/**
	 * Sets the translational error rate per codon directly.
	 **/
	void setErrorRate(const double error_rate) { m_error_rate = error_rate; }

	/**
	 * Sets the translational cost directly.
	 **/
	void setMisfoldingCost(const double cost) { m_tr_cost = cost; }

	/**
	 * Sets the translational error rate, accuracy weight and error weight directly.
	 * @see getWeightsForTargetAccuracy for detailed description of each quantity.
	 **/
	void setErrorWeights(const double error_rate, const double accuracy_weight, const double error_weight) {
		m_error_rate = error_rate;
		m_accuracy_weight = accuracy_weight;
		m_error_weight = error_weight;
	}

	/**
	 * Compute the per-codon error rate estimated to achieve the
	 * desired fraction accurate, given the weights provided.
	 */
	double estimateErrorRateFromAccuracy(const double target_fraction_accurate, const double accuracy_weight, const double error_weight) const;

	/**
	 * Compute the expected fraction accurately translated,
	 * given the error rate and weights provided.
	 */
	double estimateAccuracyFromErrorRate(const double error_rate, const double accuracy_weight, const double error_weight) const;

	/**
	 * Compute the estimated fractions accurately translated, folded
	 * despite mistranslation, truncated and folded, using the error
	 * spectrum of this FitnessEvaluator.  If the natively encoded
	 * sequence does not fold, these values are invalid.  @return The
	 * estimated fitness.
	 *
	 * @param g The gene to assay.
	 * @param frac_accurate The estimated fraction of accurately translated (proper amino acid sequence) proteins from this gene.
	 * @param frac_robust The estimated fraction of proteins that fold properly from this gene.
	 * @param frac_truncated The estimated fraction of truncated proteins from this gene.
	 * @param frac_folded The estimated fraction of folded proteins from this gene.
	 * @return The fitness that would result from these estimated translation outcomes.
	 **/
	virtual double calcOutcomes( const CodingDNA& g, double& frac_accurate, double& frac_robust, double& frac_truncated, double& frac_folded );

	/**
	 *	Record the actual fractions accurately translated, folded despite
	 *	mistranslation, truncated and folded, using the error spectrum of this
	 *	FitnessEvaluator, by translating a specified number of proteins.
	 *
	 * @param g The gene to assay.
	 * @param num_to_fold The number of proteins that should be translated to determine the error spectrum.
	 * @param num_accurate The observed number of accurately translated (proper amino acid sequence) proteins from this gene.
	 * @param num_robust The observed number of proteins that fold properly from this gene.
	 * @param num_truncated The observed number of truncated proteins from this gene.
	 * @param num_folded The observed number of folded proteins from this gene.
	 * @return The fitness that would result from the generated set of proteins.
	 **/
	virtual double countOutcomes(const CodingDNA& g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);

	/**
	 * Record stabilities of mistranslated proteins.
	 */
	virtual void stabilityOutcomes( const CodingDNA& g, const int num_to_fold, vector<double>& ddgs );

	/**
	 * Record proteins coming off the simulated ribosome.
	 */
	virtual void translationOutcomes( const CodingDNA& g, const int num_to_translate, bool only_mistranslated, vector<Protein>& proteins );

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

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which fitness is 1 if the native protein folds stably, and 0 otherwise.

This class implements a case of \ref ErrorproneTranslation where translation is error-prone but mistranslation has no cost.
 */
class FoldingOnlyFitness : public ErrorproneTranslation {
private:
	FoldingOnlyFitness();

protected:

public:
	FoldingOnlyFitness( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy,
						const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	FoldingOnlyFitness( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy,
						const double tr_cost, const double ca_cost, const double target_fraction_accurate );

  virtual ~FoldingOnlyFitness();

  double getFitness( const CodingDNA& g );
  double getFitness( const Protein& p );
  bool getFolded( const CodingDNA& g );
};

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which fitness differences are due to different effective translational accuracies influenced by the gene sequence.

This class implements a case of \ref ErrorproneTranslation where translation is error-prone and all mistranslated proteins are costly, independent of whether they fold.
 */
class AccuracyOnlyTranslation : public ErrorproneTranslation {
private:
	AccuracyOnlyTranslation();

protected:
	Protein m_target_sequence;
	bool sequenceFolds(Protein& p);

public:
	AccuracyOnlyTranslation( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	virtual ~AccuracyOnlyTranslation();
	//void setTargetSequence( const Protein& p);

	double getFitness( const CodingDNA& g );
    double getFitness( const Protein& p ) { return getFitness( GeneUtil::reverseTranslate(p) ); }
	bool getFolded( const CodingDNA& g );
};

/** \brief A \ref FitnessEvaluator, derived from \ref ErrorproneTranslation, in which all fitness differences are due to differing protein robustness to translation errors.

This class implements a case of \ref ErrorproneTranslation where translation always yields a fixed percentage of mistranslated polypeptides, which may then fold properly.
*/

class RobustnessOnlyTranslation : public ErrorproneTranslation {
private:
	RobustnessOnlyTranslation();
protected:
	double m_fraction_accurate;

public:

	RobustnessOnlyTranslation( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate );
	RobustnessOnlyTranslation( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy,
		const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight );
	virtual ~RobustnessOnlyTranslation();

	double getFitness( const CodingDNA& g );
    double getFitness( const Protein& p ) { return getFitness( GeneUtil::reverseTranslate(p) ); }

	/**
	Compute the estimated fractions accurately translated, folded
	despite mistranslation, truncated and folded, using the error
	spectrum of this FitnessEvaluator.
	@return The estimated fitness.
	 */
	virtual double calcOutcomes( const CodingDNA& g, double& frac_accurate, double& frac_robust, double& frac_truncated, double& frac_folded );

	/**
	 * \warning This function is currently not implemented and sets all values to zero!
	 *
	 *	Record the actual fractions accurately translated, folded despite
	 *	mistranslation, truncated and folded, using the error spectrum of this
	 *	FitnessEvaluator, by translating a specified number of proteins.
	 *
	 * @param g The gene to assay.
	 * @param num_to_fold The number of proteins that should be translated to determine the error spectrum.
	 * @param num_accurate The observed fraction of accurately translated (proper amino acid sequence) proteins from this gene.
	 * @param num_robust The observed fraction of proteins that fold properly from this gene.
	 * @param num_truncated The observed fraction of truncated proteins from this gene.
	 * @param num_folded The observed fraction of folded proteins from this gene.
	 * @return The fitness that would result from the generated set of proteins.
	 **/
	virtual double countOutcomes(const CodingDNA& g, const int num_to_fold, int& num_accurate, int& num_robust, int& num_truncated, int& num_folded);
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
	 * @param protein_folder An initialized \ref Folder object.
	 * @param length Protein sequence length in amino acids.
	 * @param protein_structure_ID Structure identifier for the native structure all folded proteins must attain.
	 * @param max_free_energy Maximum free energy of folding for a folded protein.
	 * @param tr_cost Cost factor used to convert fraction of misfolded proteins f into fitness cost via fitness = exp{- tr_cost*f/(1-f).
	 * @param ca_cost Codon adapation cost.  This cost represents the average fold-decrease in codon accuracy for non-optimal codons relative to optimal synonymous codons.
	 * @param error_rate The per-codon base translational error rate.  Non-optimal codons will have a higher error rate, determined by ca_cost.
	 * @param accuracy_weight The reference average sequence weight for a large ensemble of folded proteins, excluding the possibility of synonymous errors.  @see getWeightsForTargetAccuracy.
	 * @param error_weight The reference average sequence weight for a large ensemble of folded proteins, accounting for both synonymous and nonsynonymous errors.  @see ErrorproneTranslation::getWeightsForTargetAccuracy.
	 * @param cost_constant The multiplier to convert translational robustness cost into number of misfolded proteins.
	 * @param toxicity_cutoff The minimum number of misfolded proteins required to get a toxic effect.
	 **/
	CutoffErrorproneTranslation( Folder* protein_folder, const int length, const StructureID protein_structure_ID, const double max_free_energy, const double tr_cost, const double ca_cost, const double error_rate, const double accuracy_weight, const double error_weight, double cost_constant, int toxicity_cutoff);

	virtual ~CutoffErrorproneTranslation();
    double getFitness( const CodingDNA& g );
    double getFitness( const Protein& p ) { return getFitness( GeneUtil::reverseTranslate(p) ); }
};

/**
 * Computes the probability of fixation of a mutant with fitness advantage s in a population of size N.
 *
 * @param s The fitness advantage relative to the wild type.
 * @param N The population size.
 **/
double fixation_probability(int N, double s);

#endif
