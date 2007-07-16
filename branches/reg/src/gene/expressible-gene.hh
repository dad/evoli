#ifndef _EXPRESSIBLE_GENE_HH__
#define _EXPRESSIBLE_GENE_HH__

#include "coding-sequence.hh"
#include "upstream.hh"

class ExpressibleGene {
protected:
	Promoter m_promoter;
	CodingDNA m_coding_seq;
	double m_expression_level;
public:
	ExpressibleGene(const Promoter& prom, const CodingDNA& dna) : m_promoter(prom), m_coding_seq(dna) {
		m_expression_level = 0.0;
	}

	ExpressibleGene(const CodingDNA& gene) : m_promoter(22), m_coding_seq(gene) {
		m_expression_level = 0.0;
	}

	CodingDNA& getCodingDNA() { return m_coding_seq; }

	const CodingDNA& getCodingDNA() const { return m_coding_seq; }

	/**
	 * @return The expression level of the gene when its promoter is bound by the given transcription factor.
	 */
	double express(const TranscriptionFactor& tf, const BindingInteraction& bi, double max_binding_energy) const {
		double expr = 0.0;
		double max_expr = 1000; // Arbitrary
		//cout << m_promoter << endl;
		if (bi.bindsBelowCutoff(m_promoter, tf, max_binding_energy)) {
			expr = max_expr;
			//double time_bound = bi.getTimeBound(m_promoter, tf);
			//double time_polymerase_bound = bi.getPolymeraseAffinity();
		}
		return expr;
	}

	/**
	 * @return This gene's promoter.
	 */
	const Promoter& getPromoter() const {
		return m_promoter;
	}

	/**
	 * @return This gene's promoter.
	 */
	Promoter& getPromoter() {
		return m_promoter;
	}

	/**
	 * \brief Print this gene.
	 **/
	void print(ostream& os) const {
		os << m_promoter << " " << m_coding_seq;
	}
};

ostream& operator<<(ostream& os, const ExpressibleGene& g);

#endif // _EXPRESSIBLE_GENE_HH__

