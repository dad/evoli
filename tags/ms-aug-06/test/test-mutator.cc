#include "mutator.hh"
#include "fitness-evaluator.hh"
#include "genotype-util.hh"

#include <iostream>
#include <iterator>
#include <algorithm>


void testCreateMutationList( ProteinFolder &b, double cutoff )
{
        int length = b.getProteinLength();

        Mutator m( b, cutoff );

        for ( int i=0; i<100; i++ )
        {
                vector<int> v = m.createMutationList( length - 1 );
                sort( v.begin(), v.end() );
                copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
                cout << endl;
        }

        for ( int i=0; i<100; i++ )
        {
                vector<int> v = m.createMutationList( 3 );
                sort( v.begin(), v.end() );
                copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
                cout << endl;
        }
}

Genotype findStartingGenotype( ProteinFolder &b, double cutoff )
{
        double G, G2;
        Genotype g, g2;
        ProteinFreeEnergyFitness fetmp( &b );

        int length = b.getProteinLength();
        g = GenotypeUtil::createRandomGenotype( length );
        G = -log( fetmp.getFitness( g ) );
        int count = 1;
        do
        {
                g2 = g;
                GenotypeUtil::mutateGenotype( g2, 0.02 );
                G2 = -log( fetmp.getFitness( g2 ) );
                if ( G2 < G )
                {
                        g = g2;
                        G = G2;
                }
                count += 1;
                if ( count > 50000 )
                { // start again with random genotype if search is not successful
                        g = GenotypeUtil::createRandomGenotype( length );
                        G = -log( fetmp.getFitness( g ) );
                        count = 1;
                }
                //    cout << "# " << count << " " << G << endl;
        }
        while( G > cutoff );

        return g;
}

double testRandomTwoMutants( ProteinFolder &b, double cutoff, int structureID, Genotype g, int sample_size )
{
        ProteinStructureFitness fe( &b, structureID, cutoff );
        int count = 0;
        for ( int i=0; i<sample_size; i++ )
        {
                Genotype g2 = g;
                int pos1 = GenotypeUtil::nonsymPointMutation( g2 );
                int pos2 = GenotypeUtil::nonsymPointMutation( g2 );

                if ( pos1 == pos2 )
                        i--;
                else if ( fe.getFitness( g2 ) > 0 )
                        count += 1;
        }
        return (double) count / (double) sample_size;
}

void testSampleMutants( ProteinFolder &b, double cutoff )
{
        Genotype g = findStartingGenotype( b, cutoff - 1 );
        int structureID = b.getLastFoldedProteinStructureID();
        int length = b.getProteinLength();
        Translator t( length );
        int *seq = new int[length];

        if ( !t.translateErrorFree( g, seq ) )
        {
                cerr << "found sequence doesn't translate. something is wrong!" << endl;
                exit(-1);
        }

        Mutator m( b, cutoff );
        vector<double> v = m.calcDDGDistribution( seq );
        copy( v.begin(), v.end(), ostream_iterator<double>( cout, " " ) );
        cout << endl;

        cout << "Neutrality (genotype util): " << endl;
        cout << GenotypeUtil::calcNeutrality( b, g, cutoff ) << endl;
        cout << "Random two mutants (genotype util): " << endl;
        cout << testRandomTwoMutants( b, cutoff, structureID, g, 10000 ) << endl;
        cout << "Neutrality and higher order terms (mutator): " << endl;
        cout << 1 << " " << m.sampleMutants( seq, 1, 0 ) << endl;

        cout << 2 << " " << m.sampleMutants( seq, 2, 10000 ) << endl;
        cout << 3 << " " << m.sampleMutants( seq, 3, 10000 ) << endl;
        cout << 4 << " " << m.sampleMutants( seq, 4, 10000 ) << endl;
        cout << 5 << " " << m.sampleMutants( seq, 5, 10000 ) << endl;
        cout << 6 << " " << m.sampleMutants( seq, 6, 10000 ) << endl;

        // this takes approx. 3 min on 2.8 GHz Pentium 4:
        //cout << 7 << " " << m.sampleMutants( seq, 7, 1000000 ) << endl;

        delete [] seq;
}


int main()
{
        const int size = 5;
        double maxFreeEnergy = -5.0;

        srand48(3);

        ProteinFolder b(size);
        b.enumerateStructures();
        //testCreateMutationList( b, maxFreeEnergy );
        testSampleMutants( b, maxFreeEnergy );

}




