#ifndef POPULATION_HH
#define POPULATION_HH

#include "tools.hh"
#include "protein.hh"

#include <map>

using namespace std;

class FitnessEvaluator;

/** \brief Helper class used by Genebank
*/

class GenebankEntry
{
private:
        const Gene m_genotype; // the stored genotype
        const double m_fitness; // the fitness of this genotype
        GenebankEntry *const m_parent; // the parent genotype
        const int m_id; // genotype id
        const int m_birthTime; // the time at which the genotype first arose
        int m_count; // number of references to this entry currently in the system
        bool m_coalescent;
        mutable bool m_tagged;

        GenebankEntry();
        GenebankEntry( const GenebankEntry & g );
        const GenebankEntry & operator=( const GenebankEntry &g );
public:
        GenebankEntry( const Gene &genotype, double fitness, GenebankEntry* parent, int id, int birthTime ) :
                        m_genotype( genotype ), m_fitness( fitness ),
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
        const Gene & getGenotype() const
        {
                return m_genotype;
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


/** \brief Class that keeps track of the different genotypes and their ancestry.
*/


class Genebank
{
private:
        int m_maxId;
        map<int, GenebankEntry*> m_genotypeMap;


        Genebank( const Genebank & );
        const Genebank & operator=( const Genebank & );
public:
        Genebank();
        virtual ~Genebank();

        /**
         * Resets the genebank to an empty state.
         **/
        void clear();

        void addGenotype( GenebankEntry* g );
        void removeGenotype( GenebankEntry* g );

        /**
         * Creates a new genotype with the given characteristics. The genotype
         * does not need to be added.
         **/
        GenebankEntry* createGenotype( const Gene &g, double fitness, GenebankEntry* parent, int birthTime );

        /**
         * Reliably finds the last coalescent genotype if the given genotype is from
         * the current generation.
         **/
        bool checkCoalescence( GenebankEntry *g );

        /**
         * Print the state of the genebank into the given stream.
         **/
        void print( ostream & s ) const;
};

/** \brief The class that handles the population of evolving proteins.
*/


class Population
{
public:
        typedef GenebankEntry const *const * const_iterator;
        typedef GenebankEntry *const * iterator;

private:
        GenebankEntry **m_pop[2]; // the actual population
        double *m_selectionBins;
        int m_buffer; // the currently active population

        int m_maxN; // maximal population size
        int m_N;  // current population size
        double m_U; // the (genomic) mutation rate
        int m_time; // the time in generations
        FitnessEvaluator *m_fitness_evaluator; // needed to translate genotypes into
        // fitnesses

        Genebank m_genebank;

        Population();
        Population( const Population & );
        const Population & operator=( const Population & );
public:
        Population( int maxN );
        virtual ~Population();

        // manipulators

        /**
         * Initializes the population with the given genotype and sets
         * the mutation rate to U.
         **/
        void init( const Gene &g, FitnessEvaluator *fe, double U );

        /**
         * Create a new offspring based on a given one, with the correct
         * mutations etc.
         **/
        GenebankEntry* createOffspring( GenebankEntry* e );

        /**
         * Sets the population size. N must be smaller than maxN given to the
         * constructor. The population size is changed by resampling the population,
         * that is, after the change the new population is a random sample of the
         * old one.
         **/
        void setPopulationSize( int N );

        /**
         * Advances the population one time step.
         **/
        void evolve();

        /**
         * Needs to be called before any of the functions concerning
         * coalescence distance or Hamming distance can be used.
         **/
        void prepareCoalescenceCalcs();

        // accessors
        /**
         * Calculates the last point in time before the any descendants of the most-recent common ancestor
         * begin to branch off. That is, if for example all organisms at time t are direct descendants of
         * the most-recent common ancestor, then this function returns t-1.
         *
         * The function prepareCoalescenceCalcs must have been called before this one is called.
         **/
        int calcCoalescenceTime();

        /**
         * Calculates the average fitness of the population.
         **/
        double getAveFitness() const;

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
	const Gene& operator[](int i) const {
		return (*(m_pop[m_buffer] + i))->getGenotype();
	}


protected:
        /**
         * This function deletes all the memory that the class has allocated.
         **/
        virtual void clearMemory();
};


#endif
