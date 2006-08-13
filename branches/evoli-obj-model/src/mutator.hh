#ifndef MUTATOR_HH
#define MUTATOR_HH

#include "tools.hh"
#include "protein-folder.hh"
#include "translator.hh"

#include <vector>
#include <iostream>

/**
*  This class is used to create mutated amino acid sequences and test them for viabilities. Note that
* it does not work on DNA sequences!
*/
class Mutator
{
private:
        int m_L;
        mutable int *m_sequence;
        mutable ProteinFolder &m_folder;
        double m_cutoff;

        /**
        * Copies the given sequence into the private variable m_sequence.
        *
        */
        void copySequence( const int *sequence ) const
        {
                // copy sequence
                const int *src = sequence;
                int *tar = m_sequence;
                for ( int i=0; i<m_L; i++ )
                {
                        *tar = *src;
                        src++; tar++;
                }
        }

        /**
        * Takes the given aa sequence, and stores a mutated copy in the variable m_sequence. The mutation
        * positions are taken from the vector mut_list.
        */
        void mutateSequence( const int *sequence, const vector<int> &mut_list ) const
        {
                // first, copy the sequence
                copySequence( sequence );

                // now make mutations
                vector<int>::const_iterator it = mut_list.begin(), e = mut_list.end();

                for ( ; it!= e; it++ )
                {
                        m_sequence[*it] += (int) (19 * myRand() ) + 1;
                        if ( m_sequence[*it] >= 20 )
                                m_sequence[*it] -= 20;

                }
        }


        Mutator( const Mutator & );
        Mutator& operator=( const Mutator & );
public:
        Mutator( ProteinFolder &folder, double cutoff );
        ~Mutator();

        /**
        * Creates a list of n distinct mutation positions that lie
        * between 0 and length-1.
        */
        vector<int> createMutationList( int n ) const
        {
                vector<int> v;
                vector<int>::iterator it, e;

                // this algorithm is somewhat inefficient but foolproof
                int i = 0;
                while ( i < n )
                {
                        int pos = int( ( m_L ) * myRand() );

                        it = v.begin();
                        e = v.end();
                        bool exists = false;
                        for ( ; it!= e; it++ )
                        {
                                if ( *it == pos )
                                {
                                        exists = true;
                                        break;
                                }
                        }

                        if ( !exists )
                        {
                                v.push_back( pos );
                                i += 1;
                        }
                }
                return v;
        }

        /**
        * Makes all one mutants to the given sequence, and returns in a pair the number of correctly folding
        * one mutants (pair.first) and the total number of mutants (pair.second).
        */

        pair<int, int> calcNeutralityRaw( const int *sequence ) const
        {
                // copy the sequence
                copySequence( sequence );

                double G = m_folder.foldProtein( sequence );
                int id = m_folder.getLastFoldedProteinStructureID();

                if ( G > m_cutoff )
                        return pair<int, int>( 0, 1 );

                int count = 0;
                // go through all positions in the genome
                for ( int i=0; i<m_L; i++ )
                {

                        // go through all possible point mutations
                        // (avoid operator %, which can be very slow)
                        const int res = m_sequence[i];
                        int tempres = res + 1;
                        while ( tempres < 20 )
                        {
                                m_sequence[i] = tempres;
                                // sequence folds into correct structure with low free energy?
                                if ( m_folder.foldProtein( m_sequence ) < m_cutoff )
                                        if ( m_folder.getLastFoldedProteinStructureID() == id )
                                                count += 1;
                                tempres++;
                        }

                        tempres = 0;
                        while ( tempres < res )
                        {
                                m_sequence[i] = tempres;
                                if ( m_folder.foldProtein( m_sequence ) < m_cutoff )
                                        if ( m_folder.getLastFoldedProteinStructureID() == id )
                                                count += 1;
                                tempres++;
                        }
                        m_sequence[i] = res;
                }

                return pair<int, int> ( count, 19*m_L );
        }

        /**
        * Same as calcNeutralityRaw, but returns the fraction functional rather than a pair with count and sample size.
        */
        double calcNeutrality( const int *sequence ) const
        {
                pair<int, int> p = calcNeutralityRaw( sequence );
                return (double) p.first / (double) p.second;
        }

        vector<double> calcDDGDistribution( const Genotype &g )
        {
                int l = m_folder.getProteinLength();
                Translator t( 0, l );
                int *seq = new int[l];

                vector<double> result;
                if ( t.translateErrorFree( g, seq ) )
                {
                        result = calcDDGDistribution( seq );
                }

                return result;
        }

        vector<double> calcDDGDistribution( const int *sequence ) const
        {
                vector<double> result;
                // copy the sequence
                copySequence( sequence );

                double G = m_folder.foldProtein( sequence );

                if ( G > m_cutoff )
                        return result;

                // go through all positions in the genome
                for ( int i=0; i<m_L; i++ )
                {

                        // go through all possible point mutations
                        // (avoid operator %, which can be very slow)
                        const int res = m_sequence[i];
                        int tempres = res + 1;
                        while ( tempres < 20 )
                        {
                                m_sequence[i] = tempres;
                                // sequence folds into correct structure with low free energy?
                                result.push_back( m_folder.foldProtein( m_sequence ) - G );
                                tempres++;
                        }

                        tempres = 0;
                        while ( tempres < res )
                        {
                                m_sequence[i] = tempres;
                                result.push_back( m_folder.foldProtein( m_sequence ) - G );
                                tempres++;
                        }
                        m_sequence[i] = res;
                }

                sort( result.begin(), result.end() );
                return result;
        }

        /**
        * Makes all n mutants to the given sequence, and returns in a pair the number of correctly folding
        * n mutants (pair.first) and the total number of mutants tested (pair.second). The second quantity is either
        * sample size, or the total number of all possible n mutants if n is small.
        */

        pair<int, int> sampleMutantsRaw( const int *sequence, int n, int sample_size ) const
        {
                if ( n == 1 )
                        return calcNeutralityRaw( sequence );

                double G = m_folder.foldProtein( sequence );
                int id = m_folder.getLastFoldedProteinStructureID();

                if ( G > m_cutoff )
                        return pair<int, int> ( 0, 1 );

                int count = 0;
                for ( int i=0; i<sample_size; i++ )
                {
                        mutateSequence( sequence, createMutationList( n ) );

                        if ( m_folder.foldProtein( m_sequence ) < m_cutoff )
                                if ( m_folder.getLastFoldedProteinStructureID() == id )
                                        count += 1;
                }
                return pair<int, int> ( count, sample_size );
        }

        /**
        * Same as sampleMutantsRaw, but returns the fraction functional rather than a pair with count and sample size.
        */
        double sampleMutants( const int *sequence, int n, int sample_size ) const
        {
                pair<int, int> p =  sampleMutantsRaw( sequence, n, sample_size );
                return (double) p.first / (double) p.second;
        }
};









#endif
