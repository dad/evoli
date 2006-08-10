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
        int samplings;
        int intermediate_steps;
        int random_seed;
        Genotype g;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Sequence: " << p.g << endl;
        s << "# Cutoff: " << p.G_cutoff << endl;
        s << "# Samplings: " << p.samplings << endl;
        s << "# Intermediate steps: " << p.intermediate_steps << endl;
        s << "# Random seed: " << p.random_seed << endl;
        s << "#" << endl;
        return s;
}

Parameters getParams( int ac, char **av )
{
        if ( ac != 6 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <free energy cutoff> <samplings> <intermediate steps> <random seed> <sequence file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.G_cutoff = atof( av[i++] );
        p.samplings = atoi( av[i++] );
        p.intermediate_steps = atoi( av[i++] );
        p.random_seed = atoi( av[i++] );

        // read initial sequence
        ifstream in( av[i] );
        in >> p.g;

        return p;
}



void analyzeSequence( ProteinFolder &b, Parameters p )
{
        // seed the random number generator
        srand48( p.random_seed );

        // find the initial starting sequence
        pair<double, int> fp = GenotypeUtil::translateAndFold( b, p.g );

        if ( fp.first > p.G_cutoff || fp.second < 0 )
        {
                cerr << "Input sequence does not translate or is not stable with given free energy cutoff!" << endl;
                cerr << "Exiting. No data written." << endl;
                exit( - 1 );
        }
        int structureID = fp.second;

        Mutator m( b, p.G_cutoff );

        Genotype g = p.g, g2;
        cout << fp.first << " " << fp.second << " " << g << endl; // always start with the initial sequence
        vector<double> result = m.calcDDGDistribution( g );
        int count = 1;
        while( count < p.samplings )
        {
                // do intermediate steps
                int count2 = 0;
                while( count2 < p.intermediate_steps )
                {
                        g2 = g;
                        GenotypeUtil::nonsymPointMutation( g2 );
                        fp = GenotypeUtil::translateAndFold( b, g2 );
                        if ( fp.first > p.G_cutoff || fp.second != structureID )
                                continue;

                        g = g2;
                        count2 += 1;
                }

                cout << fp.first << " " << fp.second << " " << g << endl;
                vector<double> v1 = result;
                vector<double> v2 = m.calcDDGDistribution( g );
                result.resize( v1.size() + v2.size() );
                merge( v1.begin(), v1.end(), v2.begin(), v2.end(), result.begin() );

                count += 1;
        }

        // output the result
        stringstream name;
        name << "id" << p.random_seed << "-ddg" << p.samplings << ".dat";
        ofstream out( name.str().c_str(), ios::out );
        out << p;
        copy( result.begin(), result.end(), ostream_iterator<double>( out, "\n" ) );
}


int main( int ac, char **av )
{
        Parameters p = getParams( ac, av );

        ProteinFolder b( 5 );
        b.enumerateStructures();

        analyzeSequence( b, p );
}




