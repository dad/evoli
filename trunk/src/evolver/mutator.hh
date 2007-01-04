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
	double getMutationRate() { return m_mutation_rate; }
	
	/**
	 * Mutates a gene with the pre-specified probability.
	 * @return Whether any mutations occurred.
	 **/
	virtual bool mutate(CodingDNA& gene);
};

/**
 * \brief Implements a more intricate mutation model in which probabilities of mutations from and to each nucleotide may be specified.
 **/
class Polymerase {
private:
	double m_mutation_rate;
	vector<vector<double> > m_mutation_matrix;
public:
	Polymerase(double mutation_rate);
	Polymerase(double mutation_rate, vector<vector<double> >& mutation_matrix);
	virtual ~Polymerase();

	/**
	 * Returns the mutation rate.
	 * @return The per-nucleotide mutation rate.
	 */
	double getMutationRate() { return m_mutation_rate; }

	/**
	 * Mutates a gene with the pre-specified probabilities.
	 * @return Whether any mutations occurred.
	 **/
	virtual bool mutate(CodingDNA& gene);
};

#endif // MUTATOR_HH
