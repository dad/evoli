/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@alumni.princeton.edu>

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


#ifndef MUTATOR_HH
#define MUTATOR_HH

#include "tools.hh"
#include "coding-sequence.hh"

using namespace std;

/**
 * \brief Implements a simple mutation model in which all point mutations are equally likely.
 **/
class SimpleMutator {
private:
	double m_mutation_rate;
public:
	SimpleMutator(double mutation_rate);
	virtual ~SimpleMutator();

	/**
	 * Returns the mutation rate.
	 * @return The per-nucleotide mutation rate.
	 */
	double getMutationRate() const { return m_mutation_rate; }

	/**
	 * Mutates a nucleotide sequence with a pre-specified probability.
	 * @return Whether any mutations occurred.
	 **/
	virtual bool mutate(NucleotideSequence& dna) const;
};

/**
 * \brief Implements a more intricate mutation model in which probabilities of mutations from and to each nucleotide may be specified.
 **/
class Polymerase {
private:
	double m_mutation_rate;  //!< overall mutation frequency
	vector<double> m_mut_rate_ACGT; //!< mutation rates for each of the four bases
	vector<double> m_AtoCGT; //!< mutation distr. from base A
	vector<double> m_CtoGTA; //!< mutation distr. from base C
	vector<double> m_GtoTAC; //!< mutation distr. from base G
	vector<double> m_TtoACG; //!< mutation distr. from base T

protected:
	/**
	Turns the vector of relative mutation frequencies into a properly normalized (cummulative)
	distribution.
	*/
	void normalize( vector<double> & p );

	/**
	Prepares mutation matrix (normalizes, turns into cumulative).
	*/
	void prepareMutationMatrix();

	/**
	Chooses a mutation at random, given the cumulative distribution vector p.
	(Assumes specifically that there are only three possible states a mutation can lead to.)
	*/
	int chooseMutation( const vector<double> & p ) const;
public:
	/**
	The default constructor uses uniform relative mutation frequencies.
	\param mutation_rate Overall mutation rate.
	*/
	Polymerase(double mutation_rate);
	/**
	The mutation frequencies from A to C, G, or T etc. need not be normalized. However, normalization is done over all four mutation frequencies. I.e., if AtoCGT is ( 2, 1, 1 ) and
	CtoGTA is ( 4, 2, 2 ), then overall mutations from C are twice as frequent as mutations from A.
	\param mutation_rate Overall mutation rate.
	\param AtoCGT Relative mutation frequencies from A to C, G, or T.
	\param CtoGTA As AtoCGT.
	\param GtoTAC As AtoCGT.
	\param TtoACG As AtoCGT.
	*/
	Polymerase(double mutation_rate, const vector<double>& AtoCGT, const vector<double>& CtoGTA,  const vector<double>& GtoTAC,  const vector<double>& TtoACG );

	/**
	The mutation frequencies from A to C, G, or T etc. need not be normalized. However, normalization is done over all four mutation frequencies. I.e., if AtoCGT is ( 2, 1, 1 ) and
	CtoGTA is ( 4, 2, 2 ), then overall mutations from C are twice as frequent as mutations from A.
	\param mutation_rate Overall mutation rate.
	\param GCtoAT Relative mutation frequencies from G->A and C->T
	\param ATtoGC As GCtoAT.
	\param GCtoTA As GCtoAT.
	\param GCtoCG As GCtoAT.
	\param ATtoCG As GCtoAT.
	\param ATtoTA As GCtoAT.
	*/
	Polymerase(double mutation_rate, double GCtoAT, double ATtoGC, double GCtoTA, double GCtoCG, double ATtoCG, double ATtoTA );

	virtual ~Polymerase();

	/**
	 * Returns the overall mutation rate.
	 * @return The average per-nucleotide mutation rate.
	 */
	double getMutationRate() const { return m_mutation_rate; }

	/**
	 * Mutates a nucleotide sequence with the pre-specified probabilities.
	 * @return Whether any mutations occurred.
	 **/
	virtual bool mutate(NucleotideSequence& dna) const;
};

#endif // MUTATOR_HH
