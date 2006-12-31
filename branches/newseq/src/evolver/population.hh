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


#ifndef POPULATION_HH
#define POPULATION_HH

#include "tools.hh"
#include "protein.hh"
#include "random.hh"
#include "gene-util.hh"
#include "fitness-evaluator.hh"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>

#include <map>

using namespace std;

class FitnessEvaluator;

/** \brief Helper class used by Genebank
 */

template <typename Organism>
class GenebankEntry
{
private:
	const Organism m_organism; // the stored organism
	const double m_fitness; // the fitness of this organism
	GenebankEntry<Organism> *const m_parent; // the parent organism
	const int m_id; // organism id
	const int m_birthTime; // the time at which the organism first arose
	int m_count; // number of references to this entry currently in the system
	bool m_coalescent;
	mutable bool m_tagged;

	GenebankEntry();
	GenebankEntry( const GenebankEntry<Organism> & g );
	const GenebankEntry<Organism> & operator=( const GenebankEntry<Organism> &g );
public:
	GenebankEntry( const Organism &organism, double fitness, GenebankEntry<Organism>* parent, int id, int birthTime ) :
		m_organism( organism ), m_fitness( fitness ),
		m_parent( parent ),
		m_id( id ), m_birthTime( birthTime ), m_count( 1 ),
        m_coalescent( false ), m_tagged( false ) {}

	void incrementCount()
	{
		m_count += 1;
	}
	bool decrementCount()
	{
		m_count -= 1;
		return m_count == 0;
	}

	void setCoalescent( bool c = true )
	{
		m_coalescent = c;
	}

	void setTagged( bool c = true ) const
	{
		m_tagged = c;
	}

	int getId() const
	{
		return m_id;
	}
	GenebankEntry* getParent() const
	{
		return m_parent;
	}
	int getBirthTime() const
	{
		return m_birthTime;
	}
	const Organism & getOrganism() const
	{
		return m_organism;
	}

	double getFitness() const
	{
		return m_fitness;
	}
	int getCount() const
	{
		return m_count;
	}
	bool isCoalescent() const
	{
		return m_coalescent;
	}
	bool isTagged() const
	{
		return m_tagged;
	}
};


/** \brief Class that keeps track of the different organisms and their ancestry.
 */


template <typename Organism>
class Genebank
{
private:
	int m_maxId;
	map<int, GenebankEntry<Organism>*> m_organismMap;


	Genebank( const Genebank & );
	const Genebank & operator=( const Genebank & );
public:
	Genebank();
	virtual ~Genebank();

	/**
	 * Resets the genebank to an empty state.
	 **/
	void clear();

	void addOrganism( GenebankEntry<Organism>* g );
	void removeOrganism( GenebankEntry<Organism>* g );

	/**
	 * Creates a new organism with the given characteristics. The organism
	 * does not need to be added.
	 **/
	GenebankEntry<Organism>* createOrganism( const Organism &g, double fitness, GenebankEntry<Organism>* parent, int birthTime );

	/**
	 * Reliably finds the last coalescent organism if the given organism is from
	 * the current generation.
	 **/
	bool checkCoalescence( GenebankEntry<Organism>*g );

	/**
	 * Print the state of the genebank into the given stream.
	 **/
	void print( ostream & s ) const;
};


/**
 * \brief Analyzes various results from an evolving population and the history stored in a Genebank object.
 *
 * \ref GenebankAnalyzer carries out evolutionary calculations on a population.  Typical usage is as follows:
 * \verbatim
 Population<Gene, ProteinStructureFitness, SimpleMutator> pop;
 //...initialize population, etc.
 GenebankAnalyzer<Gene> analyzer( pop.getGenebank() );
 //...evolve population
 analyzer.prepareCoalescenceCalcs(pop.begin(), pop.end(), pop.getNumGenerations());
 analyzer.analyzeDnDs(...);
 int coal_time = analyzer.calcCoalescenceTime();
 \endverbatim
 **/
template <typename Organism>
class GenebankAnalyzer {
public:
	typedef GenebankEntry<Organism> const *const * const_iterator;
	typedef GenebankEntry<Organism> *const * iterator;

private:
	Genebank<Organism>* m_genebank;
	iterator m_population_begin; // beginning of the population
	iterator m_population_end;  // end of the population
	int m_time; // number of generations of evolution.
	
public:
	/**
	 * Make a new GenebankAnalyzer object.
	 */
	GenebankAnalyzer(Genebank<Organism>* genebank);
	~GenebankAnalyzer() {}

	/**
	 * Count the number of synonymous and non-synonymous substitutions along the line of descent, from the
	 * most-recent common ancestor backwards for 'window_size' generations. The function outputs an error (and
	 * returns false) if there are not at least 'window_size' generations available from the most-recent
	 * common ancestor to the beginning of the simulation. The function also calculates the average fitness
	 * and the average fraction of optimal codons along the way.
	 *
	 * The vector is_optimal indicates codon optimality, allowing the fraction of optimal codons to be computed.
	 *
	 * Return values are given in the variables 'ave_dn', 'ave_ds', 'ave_N', 'ave_S', 'ave_f', 'ave_fop'.
	 **/
	bool analyzeDnDs( int window_size, double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f, double &ave_fop, const vector<bool>& is_optimal );

	/**
	 * Needs to be called before any of the functions concerning
	 * coalescence distance or Hamming distance can be used.
	 * @param begin An iterator pointing to the beginning of the evolving population
	 * @param end An iterator pointing to the end of the evolving population
	 * @param time The number of generations for which the population has evolved.
	 **/
	void prepareCoalescenceCalcs(iterator begin, iterator end, int time);

	/**
	 * Calculates the last point in time before the any descendants of the most-recent common ancestor
	 * begin to branch off. That is, if for example all organisms at time t are direct descendants of
	 * the most-recent common ancestor, then this function returns t-1.
	 *
	 * The function prepareCoalescenceCalcs must have been called before this one is called.
	 **/
	int calcCoalescenceTime();
};

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


protected:
	/**
	 * This function deletes all the memory that the class has allocated.
	 **/
	virtual void clearMemory();
};


template <typename Organism>
Genebank<Organism>::Genebank()
{
	clear();
}

template <typename Organism>
Genebank<Organism>::~Genebank()
{
	clear();
}

template <typename Organism>
void Genebank<Organism>::clear()
{
	m_maxId = 0;
	// DAD: Why does this fail to compile?
	/*map<int, GenebankEntry<Organism> * >::iterator it = m_organismMap.begin();
	for ( ; it!=m_organismMap.end(); it++ )
		delete (*it).second;
		m_organismMap.clear();*/
}

template <typename Organism>
void Genebank<Organism>::addOrganism( GenebankEntry<Organism>* g )
{
	g->incrementCount();
}

template <typename Organism>
void Genebank<Organism>::removeOrganism( GenebankEntry<Organism> *g )
{
	if ( g==0 ) return;

	if ( g->decrementCount() )
        {
			// this organism has lost all its references. We have to remove it.
			m_organismMap.erase( g->getId() );

			// recursively remove all reference counts from parents
			removeOrganism( g->getParent() );
			delete g;
        }
}

template <typename Organism>
GenebankEntry<Organism>* Genebank<Organism>::createOrganism( const Organism &g, double fitness, GenebankEntry<Organism>* parent, int birthTime )
{
	GenebankEntry<Organism> *e = new GenebankEntry<Organism>( g, fitness, parent, m_maxId, birthTime );
	m_organismMap[m_maxId] = e;
	m_maxId += 1;

	// the parent gets one more reference
	if ( parent )
		parent->incrementCount();

	return e;
}

template <typename Organism>
bool Genebank<Organism>::checkCoalescence( GenebankEntry<Organism> *g ) {
	if ( g->isCoalescent() ) {
		//    cout << g->getId() << " is already coalescent." << endl;
		return true;
	}
	//else
	//  cout << g->getId() << " is not coalescent. Checking parent..." << endl;
	
	
	GenebankEntry<Organism> *parent = g->getParent();
	
	assert( g != 0 );
	
	if ( checkCoalescence( parent ) ) {
		// if parent is coalescent and has only a count of 1, then we are
		// coalescent as well
		
		//  cout << "parent " << parent->getId() << " is coalescent, checking count..." << endl;
		if ( parent->getCount() == 1 ) {
			// cout << "count is 1. " << g->getId() << " is coalescent." << endl;
			g->setCoalescent();
			return true;
		}
		else {
			//cout << "count is > 1. " << g->getId() << " is not coalescent." << endl;
			return false;
		}
	}
	else {
		//cout << "parent " << parent->getId() << " is not coalescent. We are done." << endl;
		return false;
	}
}


template <typename Organism>
void Genebank<Organism>::print( ostream &s ) const
{
	s << "#--- Genebank ---\n#<Id> <parentId> <Count> <birth time> <coalescent> <fitness> <organism>" << endl;

	const char* tab = "\t";
	// DAD: Why does this fail to compile?
	/*map<int, GenebankEntry<Organism>*>::const_iterator it = m_organismMap.begin();
	for ( ; it!=m_organismMap.end(); it++ )
        {
			int parentId;
			GenebankEntry<Organism> *g = (*it).second;
			if ( g->getParent() == 0 )
				parentId = -1;
			else
				parentId = g->getParent()->getId();

			char c = 'n';
			if ( g->isCoalescent() )
				c = 'y';

			s << g->getId() << tab << parentId << tab << g->getCount() << tab << g->getBirthTime() << tab << c << tab << g->getFitness() << tab << g->getOrganism() << endl;
			}*/
	s << endl;
}

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
	for ( ; g != e; g++)
        {
			fitnessSum += (*g)->getFitness();
			(*bin) = fitnessSum; bin++;
        }
	// now normalize the bins
	for (int i=0; i<m_N; i++)
		m_selectionBins[i] /= fitnessSum;


	// now do the selection/mutation step
	for ( int i=0; i<m_N; i++ )
        {
			index = Random::randintFromDistr( m_selectionBins, m_N );
			m_pop[outBuffer][i] =  createOffspring( m_pop[m_buffer][index] );
        }

	// un-register the organisms of the old generation
	for (int i=0; i<m_N; i++)
        {
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

//////////////////
// GenebankAnalyzer methods
//////////////////

template <typename Organism>
GenebankAnalyzer<Organism>::GenebankAnalyzer(Genebank<Organism>* genebank) {
	m_genebank = genebank;
}

template <typename Organism>
void GenebankAnalyzer<Organism>::prepareCoalescenceCalcs(iterator begin, iterator end, int time)
{
	m_population_begin = begin;
	m_population_end = end;
	m_time = time;
	m_genebank->checkCoalescence( *m_population_begin );
}

template <typename Organism>
int GenebankAnalyzer<Organism>::calcCoalescenceTime() {
	int min_time = m_time;
	const_iterator it = m_population_begin;
	for (; it != m_population_end; it++) {
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


template <typename Organism>
bool GenebankAnalyzer<Organism>::analyzeDnDs( int window_size, double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, 
											  double &ave_f, double &ave_fop, const vector<bool>& is_optimal )
{
	vector<double> dn_vect;
	vector<double> ds_vect;
	vector<double> nsyn_sites_vect;
	vector<double> syn_sites_vect;
	vector<double> fitness_vect;
	vector<double> fop_vect;

	// first, make sure that all coalescent genotypes are marked as such
	prepareCoalescenceCalcs(m_population_begin, m_population_end, m_time);

	const GenebankEntry<Organism> *e = *m_population_begin, *e2;

	// find the first coalescent genotype in the population
	while( !e->isCoalescent() ) {
		e = e->getParent();
		assert( e != 0 );
	}
	
	// we end our analysis with this genotype, so get its birth time
	int end_time = calcCoalescenceTime();

	if ( end_time < window_size ) {
		cerr << "Window size exceeds time from start of simulation to most-recent common ancestor" << endl;
		return false;
	}

	// this defines the number of entries we need in the dn_vect and ds_vect variables
	dn_vect.resize( end_time + 1 );
	ds_vect.resize( end_time + 1 );
	nsyn_sites_vect.resize( end_time + 1 );
	syn_sites_vect.resize( end_time + 1 );
	fitness_vect.resize( end_time + 1 );
	fop_vect.resize( end_time + 1 );

	// record number of syn. and nonsyn. sites for last genotype
	int b_time = e->getBirthTime();
	double f = e->getFitness();
	double fop = GeneUtil::calcFop( e->getOrganism(), is_optimal );
	double sites = GeneUtil::calcTotalSites( e->getOrganism() );
	double syn_sites = GeneUtil::calcSynonymousSites( e->getOrganism() );
	double nsyn_sites = sites - syn_sites;
	for ( int i=b_time+1; i<end_time + 1; i++ ) {
		dn_vect[i]=0; // nothing happened here
		ds_vect[i]=0; // nothing happened here
		nsyn_sites_vect[i] = nsyn_sites;
		syn_sites_vect[i] = syn_sites;
		fitness_vect[i] = f;
		fop_vect[i] = fop;
	}
	nsyn_sites_vect[b_time] = nsyn_sites; // arrays for these quantities need to be shifted by 1
	syn_sites_vect[b_time] = syn_sites;
	fitness_vect[b_time] = f;
	fop_vect[b_time] = fop;

	// now go backwards and record all the changes
	int last_b_time = b_time;
	while ( 1 ) {
		e2 = e->getParent();
		if ( e2 == 0 ) // we have reached the final ancestor, so we are done
			break;
		b_time = e2->getBirthTime();
		f = e2->getFitness();
		fop = GeneUtil::calcFop( e->getOrganism(), is_optimal );
		sites = GeneUtil::calcTotalSites( e2->getOrganism() );
		syn_sites = GeneUtil::calcSynonymousSites( e2->getOrganism() );
		nsyn_sites = sites - syn_sites;
		double dn, ds;
		GeneUtil::calcDnDs( dn, ds, e2->getOrganism(), e->getOrganism() );
		
		// and record
		dn_vect[last_b_time] = dn;
		ds_vect[last_b_time]= ds;
		nsyn_sites_vect[b_time] = nsyn_sites;
		syn_sites_vect[b_time] = syn_sites;
		fitness_vect[b_time] = f;
		fop_vect[b_time] = fop;
		
		for ( int i=b_time+1; i<last_b_time; i++ ) {
			dn_vect[i]=0; // nothing happened here
			ds_vect[i]=0; // nothing happened here
			nsyn_sites_vect[i] = nsyn_sites;
			syn_sites_vect[i] = syn_sites;
			fitness_vect[i] = f;
			fop_vect[i] = fop;
		}
		e = e2;
		last_b_time = b_time;
	}

	// now that we have all the changes, analyze

	ave_dn = 0;
	ave_ds = 0;
	ave_N = 0;
	ave_S = 0;
	ave_f = 0;
	ave_fop = 0;

	for ( int i=end_time-window_size+1; i<=end_time; i++ ) {
		//                cout << i << " " << dn_vect[i] << " " << ds_vect[i] << " " << nsyn_sites_vect[i] << " " << syn_sites_vect[i] << " " << fitness_vect[i] << endl;
		ave_dn += dn_vect[i];
		ave_ds += ds_vect[i];
		ave_N += nsyn_sites_vect[i];
		ave_S += syn_sites_vect[i];
		ave_f += fitness_vect[i];
		ave_fop += fop_vect[i];
	}
	
	ave_N /= (double) window_size;
	ave_S /= (double) window_size;
	ave_f /= (double) window_size;
	ave_fop /= (double) window_size;

	return true;
}

#endif
