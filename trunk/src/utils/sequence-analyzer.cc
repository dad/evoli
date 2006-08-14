#include "protein-folder.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <algorithm>
#include <fstream>
#include <string>


struct Parameters
{
        double free_energy_cutoff;
        double tr_cost;
        double ca_cost;
        double transl_error_rate;
        string filename;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   Sequence file: " << p.filename << endl;
        return s;
}


Parameters getParams( int ac, char **av )
{
        if ( ac != 2 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <sequence file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.filename = av[i++];

        return p;
}

void analyzeSequences( istream &in, ProteinFolder &b, const Parameters &p )
{
        Genotype g;

        int count = 0;
        const char* tab = "\t";

        cout << p;
 		pair<double, int> fp;
        cout << "# <index> <free energy> <structure ID>" << endl;
        while ( !in.eof() ) {
			in >> g;
                fp = GenotypeUtil::translateAndFold( b, g );
                count += 1;
		        cout << count << tab << fp.first << tab << fp.second << endl;
        }

}



int main( int ac, char **av)
{
        Parameters p = getParams( ac, av );

        ifstream in( p.filename.c_str() );

        if ( !in )
        {
                cerr << "Error opening " << av[1] << endl;
                exit(-1);
        }

        // size of the lattice proteins is hardcoded
        const int size = 5;

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();


        analyzeSequences( in, b, p );
}




