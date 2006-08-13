#include "population.hh"
#include "genotype-util.hh"
#include "fitness-evaluator.hh"

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>


Genebank::Genebank()
{
        clear();
}

Genebank::~Genebank()
{
        clear();
}

void Genebank::clear()
{
        m_maxId = 0;
        map<int, GenebankEntry*>::iterator it = m_genotypeMap.begin();
        for ( ; it!=m_genotypeMap.end(); it++ )
                delete (*it).second;
        m_genotypeMap.clear();
}

void Genebank::addGenotype( GenebankEntry* g )
{
        g->incrementCount();
}

void Genebank::removeGenotype( GenebankEntry *g )
{
        if ( g==0 ) return;

        if ( g->decrementCount() )
        {
                // this genotyp has lost all its references. We have to remove it.
                m_genotypeMap.erase( g->getId() );

                // recursively remove all reference counts from parents
                removeGenotype( g->getParent() );
                delete g;
        }
}

GenebankEntry* Genebank::createGenotype( const Genotype &g, double fitness, GenebankEntry* parent, int birthTime )
{
        GenebankEntry *e = new GenebankEntry( g, fitness, parent, m_maxId, birthTime );
        m_genotypeMap[m_maxId] = e;
        m_maxId += 1;

        // the parent gets one more reference
        if ( parent )
                parent->incrementCount();

        return e;
}

bool Genebank::checkCoalescence( GenebankEntry *g )
{
        if ( g->isCoalescent() )
        {
                //    cout << g->getId() << " is already coalescent." << endl;
                return true;
        }
        //else
        //  cout << g->getId() << " is not coalescent. Checking parent..." << endl;


        GenebankEntry *parent = g->getParent();

        assert( g != 0 );

        if ( checkCoalescence( parent ) )
        {
                // if parent is coalescent and has only a count of 1, then we are
                // coalescent as well

                //  cout << "parent " << parent->getId() << " is coalescent, checking count..." << endl;
                if ( parent->getCount() == 1 )
                {
                        // cout << "count is 1. " << g->getId() << " is coalescent." << endl;
                        g->setCoalescent();
                        return true;
                }
                else
                {
                        //cout << "count is > 1. " << g->getId() << " is not coalescent." << endl;
                        return false;
                }
        }
        else
        {
                //cout << "parent " << parent->getId() << " is not coalescent. We are done." << endl;
                return false;
        }
}


void Genebank::print( ostream &s ) const
{
        s << "#--- Genebank ---\n#<Id> <parentId> <Count> <birth time> <coalescent> <fitness> <genotype>" << endl;

		const char* tab = "\t";
        map<int, GenebankEntry*>::const_iterator it = m_genotypeMap.begin();
        for ( ; it!=m_genotypeMap.end(); it++ )
        {
                int parentId;
                GenebankEntry *g = (*it).second;
                if ( g->getParent() == 0 )
                        parentId = -1;
                else
                        parentId = g->getParent()->getId();

                char c = 'n';
                if ( g->isCoalescent() )
                        c = 'y';

                s << g->getId() << tab << parentId << tab << g->getCount() << tab << g->getBirthTime() << tab << c << tab << g->getFitness() << tab << g->getGenotype() << endl;
        }
        s << endl;
}



Population::Population( int maxN )
                : m_maxN( maxN ), m_N( maxN ), m_U( 0 ), m_time( 0 )
{
        m_pop[0] = new GenebankEntry*[m_maxN];
        m_pop[1] = new GenebankEntry*[m_maxN];
        for ( int i=0; i<m_maxN; i++ )
        {
                m_pop[0][i] = 0;
                m_pop[1][i] = 0;
        }
        m_selectionBins = new double[m_maxN];
        m_buffer = 0;
}


Population::~Population()
{
        clearMemory();
}


void Population::init( const Genotype &g, FitnessEvaluator *fe, double U )
{
        m_U = U;
        m_time = 0;
        assert( fe != 0 );
        m_fitness_evaluator = fe;

        // empty the genebank
        m_genebank.clear();

        // get the parent fitness
        double f = m_fitness_evaluator->getFitness( g );

        // create the parent of all genotypes
        GenebankEntry *parent = m_genebank.createGenotype( g, f, 0, 0 );
        parent->setCoalescent();

        // initialize population with given genotype
        for ( int i=0; i<m_N; i++ )
                m_pop[m_buffer][i] = m_genebank.createGenotype( g, f, parent, 0 );

        // the parent is not part of the population anymore
        m_genebank.removeGenotype( parent );
}


GenebankEntry* Population::createOffspring( GenebankEntry* e )
{
        Genotype g = e->getGenotype();

        bool mutated = GenotypeUtil::mutateGenotype( g, m_U );
        if ( !mutated )
        {
                m_genebank.addGenotype( e );
                return e;
        }

        assert( m_fitness_evaluator != 0 );
        double f =  m_fitness_evaluator->getFitness( g );
        return m_genebank.createGenotype( g, f, e, m_time );
}


void Population::setPopulationSize( int N )
{
        if ( N > m_maxN )
                return;

        int outBuffer, index;

        outBuffer = (m_buffer ^ 1);

        // resample the population
        for ( int i=0; i<N; i++ )
        {
                index = (int) (m_N*myRand());
                m_pop[outBuffer][i] = m_pop[m_buffer][index];
                m_genebank.addGenotype( m_pop[m_buffer][index] );
        }

        // un-register the genotypes of the old generation
        for (int i=0; i<m_N; i++)
        {
                m_genebank.removeGenotype( m_pop[m_buffer][i] );
                m_pop[m_buffer][i] = 0; // make sure we don't have dangling pointers
        }

        m_N = N;

        m_buffer = outBuffer;
}

void Population::evolve()
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
                index = chooseAtRandom( m_selectionBins, m_N );
                m_pop[outBuffer][i] =  createOffspring( m_pop[m_buffer][index] );
        }

        // un-register the genotypes of the old generation
        for (int i=0; i<m_N; i++)
        {
                m_genebank.removeGenotype( m_pop[m_buffer][i] );
                m_pop[m_buffer][i] = 0; // make sure we don't have dangling pointers
        }

        m_buffer = outBuffer;
}

void Population::prepareCoalescenceCalcs()
{
        m_genebank.checkCoalescence( m_pop[m_buffer][0] );
}

int Population::calcCoalescenceTime()
{
        int min_time = m_time;

        for ( int i=0; i<m_N; i++ )
        {
                int t = m_time-1;
                const GenebankEntry *e = m_pop[m_buffer][i];


                while( !e->isCoalescent() )
                {
                        t = e->getBirthTime()-1;
                        e = e->getParent();
                        assert( e!= 0 );
                }

                if ( min_time > t )
                        min_time = t;
        }

        return min_time;
}

double Population::getAveFitness() const
{
        const_iterator bi = begin();
        const_iterator ei = end();
        double aveF = 0;


        for ( ; bi!=ei; bi++ )
                aveF += (*bi)->getFitness();

        return aveF/m_N;
}


void Population::clearMemory()
{
        delete [] m_pop[0];
        m_pop[0] = 0;
        delete [] m_pop[1];
        m_pop[1] = 0;

        delete [] m_selectionBins;
        m_selectionBins = 0;

}




bool Population::analyzeDnDs( int window_size, double &ave_dn, double &ave_ds, double &ave_N, double &ave_S, double &ave_f, double &ave_fop, const vector<bool>& is_optimal )
{
        vector<double> dn_vect;
        vector<double> ds_vect;
        vector<double> nsyn_sites_vect;
        vector<double> syn_sites_vect;
        vector<double> fitness_vect;
        vector<double> fop_vect;

        // first, make sure that all coalescent genotypes are marked as such
        prepareCoalescenceCalcs();

        const GenebankEntry *e = m_pop[m_buffer][0], *e2;

        // find the first coalescent genotype in the population
        while( !e->isCoalescent() )
        {
                e = e->getParent();
                assert( e!= 0 );
        }

        // we end our analysis with this genotype, so get its birth time
        int end_time = calcCoalescenceTime();

        if ( end_time < window_size )
        {
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
        //double fop = GenotypeUtil::calcFop( e->getGenotype(), codon_costs );
        double fop = GenotypeUtil::calcFop( e->getGenotype(), is_optimal );
        double sites = GenotypeUtil::calcTotalSites( e->getGenotype() );
        double syn_sites = GenotypeUtil::calcSynonymousSites( e->getGenotype() );
        double nsyn_sites = sites - syn_sites;
        for ( int i=b_time+1; i<end_time + 1; i++ )
        {
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
        while ( 1 )
        {
                e2 = e->getParent();
                if ( e2 == 0 ) // we have reached the final ancestor, so we are done
                        break;
                b_time = e2->getBirthTime();
                f = e2->getFitness();
                fop = GenotypeUtil::calcFop( e->getGenotype(), is_optimal );
                sites = GenotypeUtil::calcTotalSites( e2->getGenotype() );
                syn_sites = GenotypeUtil::calcSynonymousSites( e2->getGenotype() );
                nsyn_sites = sites - syn_sites;
                double dn, ds;
                GenotypeUtil::calcDnDs( dn, ds, e2->getGenotype(), e->getGenotype() );

                // and record
                dn_vect[last_b_time] = dn;
                ds_vect[last_b_time]= ds;
                nsyn_sites_vect[b_time] = nsyn_sites;
                syn_sites_vect[b_time] = syn_sites;
                fitness_vect[b_time] = f;
                fop_vect[b_time] = fop;

                for ( int i=b_time+1; i<last_b_time; i++ )
                {
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

        for ( int i=end_time-window_size+1; i<=end_time; i++ )
        {
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



