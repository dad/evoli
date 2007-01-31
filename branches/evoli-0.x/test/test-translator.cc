#include "translator.hh"
#include "genetic-code.hh"


void printProtein( int* v, int length )
{
        for ( int i=0; i<length; i++ )
        {
                //	GeneticCodeUtil::printResidue( cout, v[i] );
                //	cout << " ";
                cout << v[i] << " ";
        }
        cout << endl;
}


int main()
{
        Genotype g;
        int length = 250;
        int protein[250];
        g.resize( length );


        for ( int i=0; i<length; i++ )
        {
                g[i] = 0;
                protein[i] = -1;
        }

        for ( int i=0; i<100000; i++ )
        {
                Translator t( .1, length );
                t.translate ( g, protein );
                printProtein( protein, length );
                if ( i==30)
                {
                        protein[44056604]=0;
                }

        }
}




