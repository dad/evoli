#include "mutator.hh"
#include "fitness-evaluator.hh"
#include "genotype-util.hh"

#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <string>

struct Parameters
{
        double G_cutoff;
        double G_extra;
        int repetitions;
        int random_seed;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   free energy cutoff: " << p.G_cutoff << endl;
        s << "#   min. extra free energy: " << p.G_extra << endl;
        s << "#   repetitions: " << p.repetitions << endl;
        s << "#   random seed: " << p.random_seed << endl;
        s << "#" << endl;
        return s;
}

Parameters getParams( int ac, char **av )
{
        if ( ac != 5 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <free energy cutoff> <min. extra free energy> <repetitions> <random seed>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.G_cutoff = atof( av[i++] );
        p.G_extra = atof( av[i++] );
        p.repetitions = atoi( av[i++] );
        p.random_seed = atoi( av[i++] );

        return p;
}


pair<Genotype, double> findStartingGenotype( ProteinFolder &b, double cutoff )
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

        return pair<Genotype, double> ( g, G );
}


void doMutations( const Mutator &m, const int *sequence, int order, int sample_size, ostream &s )
{
        pair<int, int> p = m.sampleMutantsRaw( sequence, order, sample_size );
        s << order << " ";
        s << (double) p.first / (double) p.second << " ";
        s << p.first << " " << p.second << endl;
}

void analyzeSequence( ProteinFolder &b, Parameters p )
{
        // seed the random number generator
        srand48( p.random_seed );

        stringstream base_name;
        base_name << "id" << p.random_seed;
        string fname_ddG = base_name.str() + "-ddg.dat";
        string fname_ffunct = base_name.str() + "-fract-func.dat";

        ofstream ddG_file( fname_ddG.c_str(), ios::out );
        ofstream ffunct_file( fname_ffunct.c_str(), ios::out );

        // find the initial starting sequence
        pair<Genotype, double> gG = findStartingGenotype( b, p.G_cutoff - p.G_extra );
        Genotype g = gG.first;
        int structureID = b.getLastFoldedProteinStructureID();
        int length = b.getProteinLength();
        Translator t( length );
        int *seq = new int[length];

        if ( !t.translateErrorFree( g, seq ) )
        {
                cerr << "found sequence doesn't translate. something is wrong!" << endl;
                return;
        }

        Mutator m( b, p.G_cutoff );
        stringstream header;
        header << "# Sequence: " << g << "\n# Free energy: " << gG.second << endl;
        header << "# Cutoff: " << p.G_cutoff << endl;
        header << "# Structure ID: " << structureID << endl;
        ddG_file << header.str();
        ffunct_file << header.str();
        
        vector<double> v = m.calcDDGDistribution( seq );
        copy( v.begin(), v.end(), ostream_iterator<double>( ddG_file, "\n" ) );
        ddG_file << endl;
                /*
        ffunct_file << "0 1" << endl;
        doMutations( m, seq, 1, 0, ffunct_file );
        doMutations( m, seq, 2, 10000, ffunct_file );
        doMutations( m, seq, 3, 10000, ffunct_file );
        doMutations( m, seq, 4, 10000, ffunct_file );
        doMutations( m, seq, 5, 10000, ffunct_file );
        doMutations( m, seq, 6, 10000, ffunct_file );
                */
        // A sample of 1,000,000 takes approximately 3 min.
        // Hence, we can easily go out to 10,000,000, and do the whole thing in
        // maybe an hour.
        ///*
        ffunct_file << "0 1" << endl;
        /*
        doMutations( m, seq, 1, 0, ffunct_file );
        doMutations( m, seq, 2, 10000, ffunct_file );
        doMutations( m, seq, 3, 50000, ffunct_file );
        doMutations( m, seq, 4, 100000, ffunct_file );
        doMutations( m, seq, 5, 500000, ffunct_file );
        doMutations( m, seq, 6, 1000000, ffunct_file );
        doMutations( m, seq, 7, 5000000, ffunct_file );
        doMutations( m, seq, 8, 10000000, ffunct_file );
        */
        //*/

        delete [] seq;
}


int main( int ac, char **av )
{
        Parameters p = getParams( ac, av );

        ProteinFolder b( 5 );
        b.enumerateStructures();

        for ( int i=0; i<p.repetitions; i++ )
        {
                cout << "[" << i+1 << flush;
                analyzeSequence( b, p );
                cout << "/" << p.repetitions << "] " << flush;
                p.random_seed += 1;
        }
        cout << endl;
}




