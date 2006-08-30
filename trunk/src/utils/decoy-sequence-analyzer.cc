#include "decoy-contact-folder.hh"
#include "translator.hh"
#include "gene-util.hh"
#include "tools.hh"

#include <fstream>


struct Parameters
{
	string structure_file;
	string structure_dir;
	string sequence_file;
	int protein_length;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   structure file: " << p.structure_file << endl;
	s << "#   structure dir: " << p.structure_dir << endl;
	s << "#   sequence file: " << p.sequence_file << endl;
	s << "#   protein length: " << p.protein_length << endl;
	s << "#" << endl;
	return s;
}

Parameters getParams( int ac, char **av )
{
	if ( ac != 5 )
	{
		cout << "Start program like this:" << endl;
		cout << "  " << av[0] << " <struct list file> <struct dir> <seq list file> <prot length>" << endl;
		exit (-1);
	}

	Parameters p;
	int i = 1;
	p.structure_file = av[i++];
	p.structure_dir = av[i++];
	p.sequence_file = av[i++];
	p.protein_length = atoi( av[i++] );

	return p;
}

int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );

	string path = (p.structure_dir+p.structure_file);
	ifstream fin(path.c_str());
	if (!fin.good()) {// if we can't read the contact maps file, bail out
		cerr << "ERROR: can't read contact maps from " << path << endl;
		return 1;
	}
	double log_nconf = 160.0*log(10.0);
	DecoyContactFolder folder(p.protein_length, log_nconf, fin, p.structure_dir);
	if (!folder.good()) {// if we can't read the contact maps file, bail out
		cerr << "ERROR: couldn't initialize folder." << endl;
		return 1;
	}
	fin.close();

	fin.open(p.sequence_file.c_str());
	cout << p;
	cout << "# <sequence> <free energy> <structure id>" << endl;

	while (!fin.eof()) {
		string seq;
		fin >> seq;
		if (seq[0] != '#') {
			Protein p(seq);
			FoldInfo fi = folder.fold(p);
			cout << seq << " " << fi.getFreeEnergy() << " " << fi.getStructure() << endl;
		}
	}

	return 0;
}




