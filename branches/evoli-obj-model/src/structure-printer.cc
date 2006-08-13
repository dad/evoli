#include "protein-folder.hh"
#include "genotype-util.hh"

#include <fstream>

int main( int ac, char **av) {
        if ( ac != 2 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <sequence file>" << endl;
                exit (-1);
        }

        ifstream in( av[1] );

        if ( !in )
        {
                cerr << "Error opening " << av[1] << endl;
                exit(-1);
        }

        Genotype g;
        in >> g;

        // size of the lattice proteins is hardcoded
        const int size = 5;

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();
        // first, get structure
        pair<double, int> fp = GenotypeUtil::translateAndFold( b, g );
        // then, print
        cout << "Sequence encodes ";
        GenotypeUtil::printProtein(cout, g);
        cout << endl << "Structure " << fp.second << ": " << endl;
        b.printStructure( fp.second );

}




