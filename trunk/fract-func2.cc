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
        int intermediate_steps;
        int repetitions;
        int random_seed;
        Genotype g;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   free energy cutoff: " << p.G_cutoff << endl;
        s << "#   starting sequence: " << p.g << endl;
        s << "#   intermediate steps: " << p.intermediate_steps << endl;
        s << "#   repetitions: " << p.repetitions << endl;
        s << "#   random seed: " << p.random_seed << endl;
        s << "#" << endl;
        return s;
}

Parameters getParams( int ac, char **av )
{
        if ( ac != 6 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <free energy cutoff> <sequence file> <intermediate steps> <repetitions> <random seed>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.G_cutoff = atof( av[i++] );
        // read initial sequence
        ifstream in( av[i++] );
        in >> p.g;
        p.intermediate_steps = atoi( av[i++] );
        p.repetitions = atoi( av[i++] );
        p.random_seed = atoi( av[i++] );

        return p;
}


void doMutations( const Mutator &m, const int *sequence, int order, int sample_size, ostream &s )
{
        pair<int, int> p = m.sampleMutantsRaw( sequence, order, sample_size );
        s << order << " ";
        s << (double) p.first / (double) p.second << " ";
        s << p.first << " " << p.second << endl;
}

void analyzeSequence( ProteinFolder &b, Parameters p, double freeEnergy, int structureID )
{
        // seed the random number generator
        srand48( p.random_seed );

        stringstream base_name;
        base_name << "id" << p.random_seed;
        string fname_ddG = base_name.str() + "-ddg.dat";
        string fname_ffunct = base_name.str() + "-fract-func.dat";

        ofstream ddG_file( fname_ddG.c_str(), ios::out );
        ofstream ffunct_file( fname_ffunct.c_str(), ios::out );

        Genotype g = p.g;
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
        header << "# Sequence: " << g << "\n# Free energy: " << freeEnergy << endl;
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
    
        doMutations( m, seq, 1, 0, ffunct_file );
        doMutations( m, seq, 2, 10000, ffunct_file );
        doMutations( m, seq, 3, 50000, ffunct_file );
        doMutations( m, seq, 4, 100000, ffunct_file );
        doMutations( m, seq, 5, 500000, ffunct_file );
        doMutations( m, seq, 6, 1000000, ffunct_file );
        doMutations( m, seq, 7, 5000000, ffunct_file );
        doMutations( m, seq, 8, 10000000, ffunct_file );
        //*/

        delete [] seq;
}


int main( int ac, char **av )
{
        Parameters p = getParams( ac, av );

        ProteinFolder b( 5 );
        b.enumerateStructures();
        pair<double, int> fp = GenotypeUtil::translateAndFold( b, p.g );
        if ( fp.first > p.G_cutoff || fp.second < 0 )
        {
                cerr << "Input sequence does not translate or is not stable with given free energy cutoff!" << endl;
                cerr << "Exiting. No data written." << endl;
                exit( - 1 );
        }
        int structureID = fp.second;

        for ( int i=0; i<p.repetitions; i++ )
        {
                cout << "[" << i+1 << flush;
                analyzeSequence( b, p, fp.first, fp.second );
                // do intermediate steps
                int count = 0;
                Genotype g;
                while( count < p.intermediate_steps )
                {
                        g = p.g;
                        GenotypeUtil::nonsymPointMutation( g );
                        fp = GenotypeUtil::translateAndFold( b, g );
                        if ( fp.first > p.G_cutoff || fp.second != structureID )
                                continue;

                        p.g = g;
                        count += 1;
                }
                cout << "/" << p.repetitions << "] " << flush;
                p.random_seed += 1;
        }
        cout << endl;
}




