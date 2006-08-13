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
        s << "#   sequence file: " << p.filename << endl;
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

        // currently hardcoded
        p.transl_error_rate = 0.002;

        return p;
}

void analyzeSequences( ifstream &in, ProteinFolder &b, const Parameters &p )
{

        Genotype g;

		string s;
        double Gs1 = 0;
        double Gs2 = 0;
        double tr = -1;
        int count = 0;
        int rep = 0;
        int lastrep = 0;
        double lasttr = 0.0;

        cout << p;

        getline( in, s);
        getline(in, s);

        while ( !in.eof() ) {
				in >> tr >> rep >> g;
                if (rep < lastrep) {
					cout << "# <tr cost> <mean free energy> <std err>" << endl;
					cout << lasttr << " " << Gs1/(double) count << " " << ( (count*Gs2-Gs1*Gs1)/(double)(count*(count-1)) )/sqrt( (double) count ) << endl;
					count = 0;
					Gs1 = 0.0;
					Gs2 = 0.0;
				}
                pair<double, int> fp = GenotypeUtil::translateAndFold( b, g );
                double nu = GenotypeUtil::calcNeutrality( b, g, -5.0 );
                //cout << "# " << fp.first << " " << fp.second << endl;
                cout << tr << " " << rep << " " << fp.first << " " << nu << endl;
                Gs1 += fp.first;
                Gs2 += fp.first*fp.first;
                count += 1;
				lastrep = rep;
				lasttr = tr;
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




