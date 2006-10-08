/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


#include "compact-lattice-folder.hh"
#include "gene-util.hh"

#include <fstream>
#include <sstream>

std::string trim( std::string& s, const std::string& drop = " " )
{
	std::string r=s.erase(s.find_last_not_of(drop)+1);
	return r.erase(0,r.find_first_not_of(drop));
}

int main( int ac, char **av )
{
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
	CompactLatticeFolder b(size);

	unsigned int buflen = 500;
	char buf[buflen];
	do {
		in.getline(buf, buflen);
		string s = buf;
		if (s[0] != '#' && trim(s) != "") { // Skip comments
			Gene g(s);
			// first, get structure
			Protein p = g.translate();
			FoldInfo fp = b.fold(p);
			// then, print
			cout << "Sequence encodes: " << p << endl << "Structure " << fp.getStructure() << ": " << endl;
			b.printStructure( fp.getStructure(), cout, "" );
		}
	} while (!in.eof());
}




