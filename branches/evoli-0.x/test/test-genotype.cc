#include "genotype.hh"
#include "genotype-util.hh"

#include <iostream>
#include <sstream>
#include <string>

int main()
{
        int reps = 10;
        double changes = 0;
        double U = 0.01;

        for ( int i=0; i<reps; i++ )
        {
                Genotype g = GenotypeUtil::createRandomGenotype( 20 );
                Genotype g2 = g;

                bool change = GenotypeUtil::mutateGenotype( g2, U );

                GenotypeUtil::printGenotype( cout, g );
                GenotypeUtil::printGenotype( cout, g2 );

                cout << "Synonymous sites g1: " << GenotypeUtil::calcSynonymousSites( g ) << endl;
                cout << "Synonymous sites g2: " << GenotypeUtil::calcSynonymousSites( g2 ) << endl;

                double dn, ds;
                GenotypeUtil::calcDnDs( dn, ds, g, g2 );

                cout << dn << " " << ds;
                if ( change )
                        cout << " Mutation" << endl;
                else
                        cout << endl;

                changes += dn + ds;
        }

        cout << "Average number of changes: " << changes/ (double) reps << endl;

        cout << "\ntest for readGenotype:" << endl;

        Genotype g = GenotypeUtil::createRandomGenotype( 20 );
        cout << " random genotype:\n    " << g << endl;
        cout << " transfer to string, add some random stuff:\n    ";
        stringstream ss;
        ss << g << " 1 2 3 \n";
        cout << ss.str();
          g = GenotypeUtil::createRandomGenotype( 20 );
        cout << " another random genotype:\n    " << g << endl;
        cout << " read back:\n    ";
        int i = 0;
        ss >> g >> i;
        cout << g << endl;
        cout << " next number after genotype: " << i << endl;
}
