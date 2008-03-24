/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>,
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


#ifndef POPULATION_HH
#define POPULATION_HH

#include "genebank.hh"
#include "random.hh"

#include <fstream>
#include <cassert>

using namespace std;



/** \brief The class that handles the population of evolving organisms.
 */

template <typename Organism, typename FitnessEvaluator, typename Mutator>
class Population {
public:
	typedef GenebankEntry<Organism> const *const * const_iterator;
	typedef GenebankEntry<Organism> *const * iterator;

private:
	GenebankEntry<Organism> **m_pop[2]; // the actual population
	double *m_selectionBins;
	int m_buffer; // the currently active population

	int m_maxN; // maximal population size
	int m_N;  // current population size
	int m_time; // the time in generations
	FitnessEvaluator *m_fitness_evaluator; // needed to evaluate organism fitnesses
	Mutator *m_mutator; // introduces mutations into organisms

	Genebank<Organism> m_genebank;

	Population();
	Population( const Population & );
	const Population & operator=( const Population & );
public:
	Population( int maxN );
	virtual ~Population();

	// manipulators

	/**
	 * Initializes the population with the given genotype, stores the FitnessEvaluator and the Mutator.
	 **/
	void init( const Organism &g, FitnessEvaluator *fe, Mutator* mut );

	/**
	 * Create a new offspring based on a given one, with the correct
	 * mutations etc.
	 **/
	GenebankEntry<Organism>* createOffspring( GenebankEntry<Organism>* e );

	/**
	 * Sets the population size. N must be smaller than maxN given to the
	 * constructor. The population size is changed by resampling the population,
	 * that is, after the change the new population is a random sample of the
	 * old one.
	 **/
	void setPopulationSize( int N );

	/**
	 * Retrieves the population size.
	 * @return The population size.
	 **/
	int getPopulationSize() { return m_N; }

	/**
	 * Retrieves the fitness evaluator.
	 * @return The fitness evaluator.
	 **/
	FitnessEvaluator* getFitnessEvaluator() { return m_fitness_evaluator; }

	/**
	 * Retrieves the mutation rate.
	 * @return The mutation rate.
	 **/
	double getMutationRate() { return m_mutator.getMutationRate(); }

	/**
	 * Advances the population one time step.
	 **/
	void evolve();


	/**
	 * Calculates the average fitness of the population.
	 **/
	double getAveFitness() const;

	/**
	 * Returns this population's Genebank.
	 **/
	Genebank<Organism>* getGenebank() { return &m_genebank; }

	/**
	 * Returns the number of generations that this population has experienced.
	 **/
	int getNumGenerations() { return m_time; }

	/**
	 * Outputs the contents of the genebank (which contains the complete ancestry of the population)
	 * into the stream given.
	 **/
	void printGenebank( ostream &s ) const
	{
		m_genebank.print( s );
	}

	/**
	 * Outputs the contents of the current population to stdout.
	 **/
	void printPopulation() const
	{
		//iterator over the population
		for ( const_iterator g = begin(); g != end(); g++)
			cout << (*g)->getId() << endl;
	}

	/**
	 * An iterator to the first genotype in the population. The iterator
	 * is only valid until the next call of evolve().
	 **/
	const_iterator begin() const
	{
		return m_pop[m_buffer];
	}

	iterator begin()
	{
		return m_pop[m_buffer];
	}

	/**
	 * An iterator to one after the last genotype in the population. See begin().
	 **/
	const_iterator end() const
	{
		return m_pop[m_buffer] + m_N;
	}

	iterator end()
	{
		return m_pop[m_buffer] + m_N;
	}
	const Organism& operator[](int i) const {
		return (*(m_pop[m_buffer] + i))->getOrganism();
	}

	/**
	Calculates the last point in time before the any descendants of the most-recent common
	ancestor begin to branch off. That is, if for example all organisms at time t are direct
	descendants of the most-recent common ancestor, then this function returns t-1.
	*/
	int calcCoalescenceTime()
	{
		// first, make sure that all coalescent genotypes are marked as such
		m_genebank.checkCoalescence( *(begin()) );

		// store current time
		int min_time = m_time;
		const_iterator it = begin();
		for (; it != end(); it++) {
			//for ( int i=0; i<m_population.size(); i++ ) {
			int t = m_time-1;
			const GenebankEntry<Organism> *e = *it;

			while( !e->isCoalescent() ) {
				t = e->getBirthTime()-1;
				e = e->getParent();
				assert( e!= 0 );
			}

			if ( min_time > t )
				min_time = t;
		}
		return min_time;
	}


protected:
	/**
	 * This function deletes all the memory that the class has allocated.
	 **/
	virtual void clearMemory();
};


///////////////////
// Population methods
///////////////////

template <typename Organism, typename FitnessEvaluator, typename Mutator>
Population<Organism, FitnessEvaluator, Mutator>::Population( int maxN )
	: m_maxN( maxN ), m_N( maxN ), m_time( 0 )
{
	m_pop[0] = new GenebankEntry<Organism>*[m_maxN];
	m_pop[1] = new GenebankEntry<Organism>*[m_maxN];
	for ( int i=0; i<m_maxN; i++ )
        {
			m_pop[0][i] = 0;
			m_pop[1][i] = 0;
        }
	m_selectionBins = new double[m_maxN];
	m_buffer = 0;
}

template <typename Organism, typename FitnessEvaluator, typename Mutator>
Population<Organism, FitnessEvaluator, Mutator>::~Population()
{
	clearMemory();
}

template <typename Organism, typename FitnessEvaluator, typename Mutator>
void Population<Organism, FitnessEvaluator, Mutator>::init( const Organism &g, FitnessEvaluator *fe, Mutator *mut )
{
	m_mutator = mut;
	m_time = 0;
	assert( fe != 0 );
	m_fitness_evaluator = fe;

	// empty the genebank
	m_genebank.clear();

	// get the parent fitness
	double f = m_fitness_evaluator->getFitness( g );
	assert( f > 0.0 );

	// create the parent of all organisms
	GenebankEntry<Organism> *parent = m_genebank.createOrganism( g, f, 0, 0 );
	parent->setCoalescent();

	// initialize population with given organism
	for ( int i=0; i<m_N; i++ )
		m_pop[m_buffer][i] = m_genebank.createOrganism( g, f, parent, 0 );

	// the parent is not part of the population anymore
	m_genebank.removeOrganism( parent );
}

template <typename Organism, typename FitnessEvaluator, typename Mutator>
GenebankEntry<Organism>* Population<Organism, FitnessEvaluator, Mutator>::createOffspring( GenebankEntry<Organism>* e )
{
	Organism g = e->getOrganism();

	bool mutated = m_mutator->mutate(g);
	if ( !mutated ) {
		m_genebank.addOrganism( e );
		return e;
	}

	assert( m_fitness_evaluator != 0 );
	double f =  m_fitness_evaluator->getFitness( g );
	assert( f > 0.0 );
	return m_genebank.createOrganism( g, f, e, m_time );
}

template <typename Organism, typename FitnessEvaluator, typename Mutator>
void Population<Organism, FitnessEvaluator, Mutator>::setPopulationSize( int N )
{
	if ( N > m_maxN )
		return;

	int outBuffer, index;

	outBuffer = (m_buffer ^ 1);

	// resample the population
	for ( int i=0; i<N; i++ ) {
		index = Random::rint( m_N );
			m_pop[outBuffer][i] = m_pop[m_buffer][index];
			m_genebank.addOrganism( m_pop[m_buffer][index] );
	}

	// un-register the organisms of the old generation
	for (int i=0; i<m_N; i++) {
		m_genebank.removeOrganism( m_pop[m_buffer][i] );
		m_pop[m_buffer][i] = 0; // make sure we don't have dangling pointers
	}

	m_N = N;

	m_buffer = outBuffer;
}

template <typename Organism, typename FitnessEvaluator, typename Mutator>
void Population<Organism, FitnessEvaluator, Mutator>::evolve()
{
	double fitnessSum;
	int outBuffer, index;

	outBuffer = (m_buffer ^ 1);
	fitnessSum=0;

	// advance generation time by one.
	m_time += 1;

	// create the selection bins
	iterator g = begin(); //iterator over the population
	iterator e = end();
	double *bin = m_selectionBins;
	for ( ; g != e; g++) {
		fitnessSum += (*g)->getFitness();
		(*bin) = fitnessSum; bin++;
	}
	assert(fitnessSum > 0.0);
	if (fitnessSum <= 0.0) {
		cout << m_time << ": got <= 0 fitness sum " << fitnessSum << endl;
	}

	// now normalize the bins
	for (int i=0; i<m_N; i++) {
		m_selectionBins[i] /= fitnessSum;
	}

	// now do the selection/mutation step
	for ( int i=0; i<m_N; i++ ) {
		index = Random::randintFromDistr( m_selectionBins, m_N );
		double f = m_pop[m_buffer][index]->getFitness();
		assert( f > 0.0);
		if (f <= 0.0) {
			// Need to decide what to do -- besides printing -- when this happens.
			//cout << "reproducing with fitness <= 0.0" << endl;
		}

		m_pop[outBuffer][i] = createOffspring( m_pop[m_buffer][index] );
	}

	// un-register the organisms of the old generation
	for (int i=0; i<m_N; i++) {
		m_genebank.removeOrganism( m_pop[m_buffer][i] );
		m_pop[m_buffer][i] = 0; // make sure we don't have dangling pointers
	}

	m_buffer = outBuffer;
}


template <typename Organism, typename FitnessEvaluator, typename Mutator>
double Population<Organism, FitnessEvaluator, Mutator>::getAveFitness() const
{
	const_iterator bi = begin();
	const_iterator ei = end();
	double aveF = 0;


	for ( ; bi!=ei; bi++ )
		aveF += (*bi)->getFitness();

	return aveF/m_N;
}


template <typename Organism, typename FitnessEvaluator, typename Mutator>
void Population<Organism, FitnessEvaluator, Mutator>::clearMemory()
{
	delete [] m_pop[0];
	m_pop[0] = 0;
	delete [] m_pop[1];
	m_pop[1] = 0;
	delete [] m_selectionBins;
	m_selectionBins = 0;

}


#endif
