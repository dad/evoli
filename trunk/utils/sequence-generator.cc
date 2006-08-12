#include "protein-folder.hh"
#include "translator.hh"
#include "genotype-util.hh"
#include "tools.hh"

#include <fstream>


struct Parameters
{
	double free_energy_cutoff;
	int repetitions;
	int random_seed;
	int struct_id;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   free energy cutoff: " << p.free_energy_cutoff << endl;
	s << "#   repetitions: " << p.repetitions << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   target structure id: " << p.struct_id << endl;
	s << "#" << endl;
	return s;
}

Parameters getParams( int ac, char **av )
{
	if ( ac != 5 )
	{
		cout << "Start program like this:" << endl;
		cout << "  " << av[0] << " <free_energy_cutoff> <repetitions> <random seed> [<struct id>|-1]" << endl;
		exit (-1);
	}

	Parameters p;
	int i = 1;
	p.free_energy_cutoff = atof( av[i++] );
	p.repetitions = atoi( av[i++] );
	p.random_seed = atoi( av[i++] );
	p.struct_id = atoi( av[i++] );

	return p;
}


// finds a random sequence with folding energy smaller than cutoff.
void getSequence( ProteinFolder &b, const Parameters &p, ostream &s )
{
	// find a random sequence with folding energy smaller than cutoff
	double G;
	int id;
	Genotype g, g2;
	pair<double, int> fdata;
	bool found = false;

	g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
	fdata = GenotypeUtil::translateAndFold(b,g);

	G = fdata.first;
	id = fdata.second;
	int fail_count = 0;

	do	{
//		 cout << "# " << fail_count << " " << G << " " << id << endl;
		g2 = g;
		bool changed = false;
		do {
			changed = GenotypeUtil::mutateGenotype( g2, 0.02 );
		} while (!changed);
		fdata = GenotypeUtil::translateAndFold( b, g2 );
		if ( fdata.first < G && fdata.second >= 0 ) {
			g = g2;
			G = fdata.first;
			id = fdata.second;
			fail_count = 0;
		}
		else {
			fail_count++;
		}

		if ( fail_count > 50000 )
		{ // start again with random genotype if search is not successful after 30000 iterations
			found = false;
			g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
			fdata = GenotypeUtil::translateAndFold( b, g );
			G = fdata.first;
			id = fdata.second;
			fail_count = 0;
		}
	}
	while( G > p.free_energy_cutoff );

	s << g << " " << G << " " << id << " " << GenotypeUtil::calcNeutrality( b, g, p.free_energy_cutoff ) << endl;
}

// finds a random sequence with folding energy smaller than cutoff and structure given by struct_id
void getSequenceTargeted( ProteinFolder &b, const Parameters &p, const int struct_id, ostream &s )
{
	// find a random sequence with folding energy smaller than cutoff
	double G;
	Genotype g, g2;
	pair<double, int> fdata;
	bool found = false;
	Translator t(0, b.getProteinLength());
	int *seq = new int[b.getProteinLength()];
	double min_free_energy_for_starting = max(0.0, p.free_energy_cutoff);
	double eps = 1e-4;

	// find sequence that encodes our target
	do	{
		g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
		found = t.translateErrorFree(g, seq) && b.isFoldedBelowThreshold(seq, struct_id, min_free_energy_for_starting);
	}
	while ( !found );

	int fail_count = 0;
	int total_fail_count = 0;
	G = fdata.first;

	do {
		// cout << "# " << fail_count << " " << G << " " << endl;
		g2 = g;
		bool changed = false;
		do {
			changed = GenotypeUtil::mutateGenotype( g2, 0.02 );
		} while (!changed);

		if (t.translateErrorFree(g2, seq) && b.isFoldedBelowThreshold(seq, struct_id, G-eps)) {
			g = g2;
			G = b.foldProtein(seq);
			fail_count = 0;
		}
		else {
			fail_count++;
			total_fail_count++;
		}

		if ( fail_count > 50000 || total_fail_count > 1e6 )
		{ // start again with random genotype if search is not successful after some # of iterations
			found = false;
			do  {
				g = GenotypeUtil::createRandomGenotypeNoStops( b.getProteinLength() );
				found = t.translateErrorFree(g, seq) && b.isFoldedBelowThreshold(seq, struct_id, min_free_energy_for_starting);
			}
			while ( !found );
			G = fdata.first;
			fail_count = 0;
			total_fail_count = 0;
		}
	}
	while( G > p.free_energy_cutoff );

	fdata = GenotypeUtil::translateAndFold( b, g);
	s << g << " " << fdata.first << " " << fdata.second << " " << GenotypeUtil::calcNeutrality( b, g, p.free_energy_cutoff )
	  << endl;
	delete [] seq;
}

int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );

	// set random seed
	srand48( p.random_seed );

	// size of the lattice proteins is hardcoded
	const int size = 5;

	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();

	cout << p;
	cout << "# <sequence> <free energy> <structure id> <neutrality>" << endl;

	for ( int i=0; i<p.repetitions; i++ )
	{
		if (p.struct_id < 0) {
			getSequence( b, p, cout );
		}
		else {
			getSequenceTargeted( b, p, p.struct_id, cout );
		}
	}
}




