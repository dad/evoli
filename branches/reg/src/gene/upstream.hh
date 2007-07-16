#ifndef _UPSTREAM_HH__
#define _UPSTREAM_HH__

#include "protein.hh"
#include "codon.hh"
#include <fstream>
#include <cmath>

using namespace std;

class Promoter : public NucleotideSequence {
public:
	Promoter(unsigned int length) : NucleotideSequence(length) {}
	Promoter(const string& s) : NucleotideSequence(s) {}
	const Sequence getBindingSequence() const {
		return *this;
	}
};


class TranscriptionFactor : public Protein {
public:
	TranscriptionFactor(const string& seq) : Protein(seq) {
	}

	/**
	 *
	 */
	const Sequence getBindingSurface() const {
		return *this;
	}

	/**
	 * Return the affinity of this TF for RNA polymerase.
	 *
	 * TFs recruit polymerase to promoters.
	 * @return the expected proportion of time polymerase occupies the promoter when this TF is bound.
	 */
	double getPolymeraseAffinity() {
		return 1.0;
	}
};

class BindingInteraction {
private:
	vector<vector<double> > m_energy_matrix;

	/**
	 * @return Proportion of time bound, given binding energy
	 */
	double getTimeBound(double energy) const {
		double kT = 0.6;
		return 1.0 - exp(-energy/kT);
	}

public:
	BindingInteraction(const string& matrix_filename) {
		ifstream fin(matrix_filename.c_str());
		if (fin.good()) {
			readEnergyMatrix(fin);
		}
		else {
			cerr << "Bad energy matrix file " << matrix_filename << endl;
		}
		fin.close();
	}

	BindingInteraction(istream& fin) {
		if (fin.good()) {
			readEnergyMatrix(fin);
		}
	}

	/**
	 * @return Energy of the promoter-TF interaction.
	 */
	double getBindingEnergy(const Promoter& promoter, const TranscriptionFactor& tfactor) const {
		double energy = 0.0;
		// Get binding surface from tfactor
		//Sequence binding_surface = tfactor.getBindingSequence();
		// Get binding sequence from promoter
		const NucleotideSequence& bound_sequence = promoter.getBindingSequence();
		// Presently, binding is independent of tfactor sequence.
		// Compute additive energy using energy matrix
		for (unsigned int i=0; i<bound_sequence.length(); i++) {
			int bsi = Codon::baseToInt(bound_sequence[i]);
			assert(bsi>=0 && bsi<4);
			//cout << i << " " << bound_sequence[i] << " " << bsi << " " << energy << " " << m_energy_matrix[i][bsi] << endl;
			energy += m_energy_matrix[i][bsi];
		}
		return energy;
	}

	/**
	 * @return Whether the TF binds the promoter with an energy below max_binding_energy.
	 */
	bool bindsBelowCutoff(const Promoter& promoter, const TranscriptionFactor& tfactor, double max_binding_energy) const {
		double energy = getBindingEnergy(promoter, tfactor);
		return energy <= max_binding_energy;
	}

	/**
	 * @return Proportion of time bound
	 */
	double getTimeBound(const Promoter& promoter, const TranscriptionFactor& tfactor) const {
		// In reality, we should compute a free energy which accounts for all other
		// possible binding sites for this tfactor.  For now, just use energy as free energy.

		// The energy for the best (consensus) sequence is 0 in the Brown/Callan matrix.
		// Shift down so that a good sequence (E=7) is bound most of the time.
		double energy_shift = -9;
		double energy = getBindingEnergy(promoter, tfactor);
		
		return getTimeBound(energy + energy_shift);
	}

	/**
	 * Reads Px4 matrix of energies, where P is the number of sites.
	 */
	void readEnergyMatrix(istream& fin) {
		unsigned int num_rows;
		fin >> num_rows;
		m_energy_matrix.resize(num_rows);
		for (unsigned int i=0; i<num_rows; i++) {
			vector<double> rowvec(4);
			// Order in Brown/Callan's matrix: ACGT
			for (int j=0; j<4; j++) {
				fin >> rowvec[j];
			}
			// Order in E.voli parlance: ACGT(U)
			m_energy_matrix[i] = rowvec;
		}
		assert(m_energy_matrix.size()==num_rows);
	}

	void print(ostream& os) const {
		// rows
		for (unsigned int i=0; i<m_energy_matrix.size(); i++) {
			// columns
			for (unsigned int j=0; j<m_energy_matrix[i].size(); j++) {
				os << m_energy_matrix[i][j] << "\t";
			}
			os << endl;
		}
	}
};

ostream& operator<<(ostream& os, const BindingInteraction& bi);

#endif // _UPSTREAM_HH__
