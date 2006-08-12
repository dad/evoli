#include "protein-folder.hh"
#include "genotype-util.hh"

#include <fstream>
#include <sstream>

std::string
trim(std::string& s,const std::string& drop = " ") {
	std::string r=s.erase(s.find_last_not_of(drop)+1);
	return r.erase(0,r.find_first_not_of(drop));
}

int
main( int ac, char **av) {
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

	// size of the lattice proteins is hardcoded
	const int size = 5;

	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();

	Genotype g;
	unsigned int buflen = 500;
	char buf[buflen];
	do {
		in.getline(buf, buflen);
		string s = buf;
		if (s[0] != '#' && trim(s) != "") { // Skip comments
			stringstream ss(buf);
			ss >> g;
			// first, get structure
			pair<double, int> fp = GenotypeUtil::translateAndFold( b, g );
			// then, print
			cout << "Sequence encodes ";
			GenotypeUtil::printProtein(cout, g);
			cout << endl << "Structure " << fp.second << ": " << endl;
			b.printStructure( fp.second );
		}
	} while (!in.eof());
}




