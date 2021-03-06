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


#ifndef GENEBANK_HH
#define GENEBANK_HH

#include "tools.hh"
#include "gene-util.hh" // for GenebankAnalyzer

#include <iostream>
#include <fstream>
#include <cassert>
#include <map>

using namespace std;

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
public:
	typedef map<int, GenebankEntry<Organism> * > EntryMap;
	//typedef EntryMap::iterator EntryMapIterator;
	//typedef EntryMap::const_iterator const_EntryMapIterator;
private:
	int m_maxId;
	EntryMap m_organismMap;


	Genebank( const Genebank & );
	const Genebank & operator=( const Genebank & );
public:
	Genebank() {
		clear();
	}

	virtual ~Genebank(){
		clear();
	}

	/**
	 * Resets the genebank to an empty state.
	 **/
	void clear()
	{
		m_maxId = 0;
		typename EntryMap::iterator it = m_organismMap.begin();
		for ( ; it!=m_organismMap.end(); it++ )
			delete (*it).second;
			m_organismMap.clear();
	}

	void addOrganism( GenebankEntry<Organism>* g ) {
		g->incrementCount();
	}

	void removeOrganism( GenebankEntry<Organism>* g )
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

	/**
	 * Creates a new organism with the given characteristics. The organism
	 * does not need to be added.
	 **/
	GenebankEntry<Organism>* createOrganism( const Organism &g, double fitness, GenebankEntry<Organism>* parent, int birthTime )
	{
		GenebankEntry<Organism> *e = new GenebankEntry<Organism>( g, fitness, parent, m_maxId, birthTime );
		m_organismMap[m_maxId] = e;
		m_maxId += 1;

		// the parent gets one more reference
		if ( parent )
			parent->incrementCount();

		return e;
	}



	/**
	 * Reliably finds the last coalescent organism if the given organism is from
	 * the current generation.
	 **/
	bool checkCoalescence( GenebankEntry<Organism>*g )
	{
		if ( g->isCoalescent() ) {
			return true;
		}
	
		GenebankEntry<Organism> *parent = g->getParent();
		assert( g != 0 );
	
		if ( checkCoalescence( parent ) ) {
		// if parent is coalescent and has only a count of 1, then we are
		// coalescent as well	
			if ( parent->getCount() == 1 ) {
				g->setCoalescent();
				return true;
			}
			else {
				return false;
			}
		}
		else {
			return false;
		}
	}


	/**
	 * Print the state of the genebank into the given stream.
	 **/
	void print( ostream & s ) const
	{
		s << "#--- Genebank ---\n#<Id> <parentId> <Count> <birth time> <coalescent> <fitness> <organism>" << endl;

		const char* tab = "\t";
		typename EntryMap::const_iterator it = m_organismMap.begin();
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
			}
		s << endl;
	}
};


/**
 * \brief Analyzes various results from an evolving population and the history stored in a Genebank object.
 *
 * \ref GenebankAnalyzer carries out evolutionary calculations on a population.  Typical usage is as follows:
 * \verbatim
 Population<Gene, ProteinStructureFitness, SimpleMutator> pop;
 //...initialize population, etc.
 GenebankAnalyzer<Gene> analyzer( pop.getGenebank(), pop.begin(), pop.end(), pop.getNumGenerations() );
 //...evolve population
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
	iterator m_population_begin; ///< beginning of the population
	iterator m_population_end;  ///< end of the population
	int m_time; ///< number of generations of evolution.
	int m_coalescence_time; ///< birth time of the last-recent common ancestor
	
public:
	/**
	 * Make a new GenebankAnalyzer object.
	 * @param begin An iterator pointing to the beginning of the evolving population
	 * @param end An iterator pointing to the end of the evolving population
	 * @param end_time The number of generations for which the population has evolved.
	 * @param coalescence_time The birth time of the last-recent common ancestor
	 */
	GenebankAnalyzer(Genebank<Organism>* genebank, iterator begin, iterator end, int end_time, int coalescence_time )
		: m_genebank( genebank), m_population_begin( begin ), m_population_end( end ), m_time( end_time ), m_coalescence_time( coalescence_time )
	{
	}
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


};


//////////////////
// GenebankAnalyzer methods
//////////////////



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

	// This shouldn't be necessary, but just to be on the save side
	m_genebank->checkCoalescence( *m_population_begin );

	const GenebankEntry<Organism> *e = *m_population_begin, *e2;

	// find the first coalescent genotype in the population
	while( !e->isCoalescent() ) {
		e = e->getParent();
		assert( e != 0 );
	}
	
	// we end our analysis with this genotype, so get its birth time
	int end_time = m_coalescence_time;

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
	while ( true ) {
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
		
		pair<double,double> dnds = GeneUtil::calcDnDs( e2->getOrganism(), e->getOrganism() );
		
		// and record
		dn_vect[last_b_time] = dnds.first;
		ds_vect[last_b_time]= dnds.second;
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
