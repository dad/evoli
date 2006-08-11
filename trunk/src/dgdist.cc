#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>



const char* tab = "\t";
typedef unsigned int uint16;
const int size = 5;

struct Stats
{
	double sum;
	double sumSquared;
	uint16 samples;

	void addValue(double val)
	{
		sum += val;
		sumSquared += val*val;
		samples++;
	}

	void reset(void)
	{
		sum = 0.0;
		sumSquared = 0.0;
		samples = 0;
	}

	double getMean(void) { return sum/samples; }
	double getVariance(void) { return sumSquared/samples - (sum/samples)*(sum/samples); }
	double getSampleVariance(void) { return sumSquared/(samples-1) - (sum/samples)*(sum/samples)*(double)samples/(samples-1); }
	double getStandardDeviation(void) { return sqrt(getVariance()); }
	double getStandardError(void) { return sqrt(getVariance()/samples); }
	double getZscore(double val) { return (val - getMean())/getStandardDeviation(); }
};

struct RunRecord
{
	double cost;
	int runNumber;
	Genotype gene;
};

struct Parameters
{
	double free_energy_cutoff;
	double free_energy_target;
	int structure_ID;
	Genotype g;
};


Genotype getSequence( ProteinFolder &b, const Parameters &p, ostream &s );
Genotype findSequenceForTargetStructure( ProteinFolder &folder, double free_energy, int structure_ID);
Genotype hillclimb( ProteinFolder &b, double free_energy, Genotype g);
void getddGDistribution(ProteinFolder &b, Genotype& g, double* arr);
void getddGDistributionAA(ProteinFolder &b, Genotype& g, double* arr);
bool isPointMutation(int codon1, int codon2);

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   free energy cutoff: " << p.free_energy_cutoff << endl;
	s << "#   free energy target: " << p.free_energy_target << endl;
	s << "#   structure ID: " << p.structure_ID << endl;
	s << "#" << endl;
	return s;
}

Parameters getParams( int ac, char **av )
{
	if ( ac != 4 )
	{
		cout << "Start program like this:" << endl;
		cout << "  " << av[0] << " <free energy cutoff> <free energy target> <structure ID>" << endl;
		exit (-1);
	}

	Parameters p;
	int i = 1;
	p.free_energy_cutoff = atof( av[i++] );
	p.free_energy_target = atof( av[i++] );
	p.structure_ID = atoi( av[i++] );

	return p;
}

int getStructureID( ProteinFolder &b, const Genotype &g )
{
	int l = b.getProteinLength();
	Translator t( 0, l );
	int *seq = new int[l];

	if ( t.translateErrorFree( g, seq ) )
	{
		b.foldProtein( seq );
		return b.getLastFoldedProteinStructureID();
	}
	else
		return -1;
	delete [] seq;
}

Genotype getSequence( ProteinFolder &b, const Parameters &p, ostream &s )
{
	// find a random sequence with folding energy smaller than cutoff
	double G;
	int id;
	Genotype g, g2;
	pair<double, int> fdata;

	// find sequence that translates
	do
	{
		g = GenotypeUtil::createRandomGenotype( b.getProteinLength() );
		fdata = GenotypeUtil::translateAndFold( b, g );
	}
	while ( fdata.second < 0 );

	G = fdata.first;
	id = fdata.second;
	int count = 1;
	do {
//		 cout << "# " << count << " " << G << " " << id << endl;
		g2 = g;
		GenotypeUtil::mutateGenotype( g2, 0.02 );
		fdata = GenotypeUtil::translateAndFold( b, g2 );
		if ( fdata.first < G && fdata.second >= 0 )
		{
			g = g2;
			G = fdata.first;
			id = fdata.second;
		}
		count += 1;
		if ( count > 30000 )
		{ // start again with random genotype if search is not successful after 50000 iterations
			do
			{
				g = GenotypeUtil::createRandomGenotype( b.getProteinLength() );
				fdata = GenotypeUtil::translateAndFold( b, g );
			}
			while ( fdata.second < 0 );
			count = 1;
			G = fdata.first;
			id = fdata.second;
		}
	}
	while( G > p.free_energy_cutoff );

	return g;
}

Genotype findSequenceForTargetStructure( ProteinFolder &b, double free_energy, int structure_ID){
	uint16 length = size*size;
	// find a random sequence with folding energy smaller than cutoff
	double G;
	int id=-1;
	Genotype g, g2;
	pair<double, int> fdata(0.0, -1);
	ProteinStructureFitness fe(&b, structure_ID, free_energy);
	while (fdata.second != structure_ID) {
		g = GenotypeUtil::createRandomGenotype( length );
		fdata = GenotypeUtil::translateAndFold( b, g );
	}

	G = fdata.first;
	id = fdata.second;
	int count = 1;
	do {
		g2 = g;
		GenotypeUtil::mutateGenotype( g2, 0.02 );
		fdata = GenotypeUtil::translateAndFold( b, g2 );
		if ( (fdata.first < G) && (fdata.second == structure_ID) ){
			g = g2;
			G = fdata.first;
			id = fdata.second;
		}
		count += 1;
		if ( count > 30000 ){ // start again with random genotype
			// if search is not successful
			while (fdata.second != structure_ID) {
				g = GenotypeUtil::createRandomGenotype( length );
				fdata = GenotypeUtil::translateAndFold( b, g2 );
			}
			G = fdata.first;
			id = fdata.second;
			count = 1;
		}
		//cout << "# " << count << " " << G << endl;
	}
	while( G > free_energy );
	cout << "# Found initial sequence, structure ID " << id << endl;
	cout << "# with free energy " << G << endl;
	return g;
}

Genotype hillclimb( ProteinFolder &b, double free_energy, Genotype g) {
	// find a random sequence with folding energy smaller than cutoff
	Genotype g2;
	pair<double, int> fdata;

	fdata = GenotypeUtil::translateAndFold( b, g );

	if ( fdata.second < 0 ) {
		cout << "Initial sequence does not translate. Nothing to be done." << endl;
		return g;
	}

	double G = fdata.first;
	int id = fdata.second;
	int count = 0;

	do {
		g2 = g;
		bool changed = GenotypeUtil::mutateGenotype( g2, 0.01 );
		fdata = GenotypeUtil::translateAndFold( b, g2 );
		if ( changed && fdata.second==id && fdata.first < G) {
			g = g2;
			G = fdata.first;
			//cout << -1*fdata.first << endl;
			count += 1;
		}
	}
	while( count < 1000 && G > free_energy);
	return g;
}

bool isPointMutation(int codon1, int codon2) {
	int a1, b1, c1, a2, b2, c2;
	CodonUtil::codonToLetters(a1, b1, c1, codon1);
	CodonUtil::codonToLetters(a2, b2, c2, codon2);

	bool res = false;
	if (a1 == a2 && b1 == b2 && c1 != c2) {
		res = true;
	}
	else if (a1 == a2 && b1 != b2 && c1 == c2) {
		res = true;
	}
	else if (a1 != a2 && b1 == b2 && c1 == c2) {
		res = true;
	}
	return res;
}

void getddGDistribution(ProteinFolder &b, Genotype& g, double* arr) {
	Genotype gcopy = g;
	pair<double, int> fdata;
	fdata = GenotypeUtil::translateAndFold( b, g );

	// Got a valid protein
	double G = fdata.first;
	int structure_ID = fdata.second;

	int i=0;
	for (unsigned int j=0; j<g.size(); j++) {
		int codon = (g[j] + 1) % 64;
		while (codon != g[j]) {
			bool point = isPointMutation(g[j], codon);
			if (point) {
				// This is a single-substitution mutation
				// Fold it.
				gcopy[j] = codon;
				fdata = GenotypeUtil::translateAndFold( b, gcopy );
				double ddG;
				if (fdata.second != structure_ID) {
					ddG = 25.0;
				}
				else {
					ddG = fdata.first - G;
				}
				arr[i++] = ddG;
				gcopy[j] = g[j];
			}
			codon = (codon + 1) % 64;
		}
	}
}

// Get ddG distribution considering all possible aa substitutions, not
// constrained by nucleotide mutations.
void getddGDistributionAA(ProteinFolder &b, Genotype& g, double* arr) {
	pair<double, int> fdata = GenotypeUtil::translateAndFold( b, g );
	int L = b.getProteinLength();
	Translator t(0, L);

	// Got a valid protein
	double G = fdata.first;
	int structure_ID = fdata.second;

	int* prot = new int[L];
	t.translateErrorFree(g, prot);

	int i=0;
	for (int j=0; j<L; j++) {
		int aaj = prot[j];
		int aa = (prot[j] + 1) % 20;
		while (aa != aaj) {
			// Fold it.
			prot[j] = aa;
			double ddG = b.foldProtein(prot) - G;
			int id = b.getLastFoldedProteinStructureID();
			if (id != structure_ID) {
				ddG = 25.0;
			}
			arr[i++] = ddG;
			aa = (aa + 1) % 20;
		}
		prot[j] = aaj;
	}
	delete [] prot;
}


void dGDistExperiment(Parameters& p)
{
	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();
	//FitnessEvaluator *feG = new ProteinFreeEnergyFitness(&b);

	pair <double, int> fdata;
	Genotype g = findSequenceForTargetStructure(b, p.free_energy_cutoff, p.structure_ID);
	fdata = GenotypeUtil::translateAndFold( b, g );

	bool ntDist = false;
	int n = 0;
	if (ntDist) {
		n = g.size()*9;
	}
	else {
		n = g.size()*19;
	}
	double* arr = new double[n];
	for (int j=1; j<=3; j++) {
		Genotype g2 = hillclimb(b, p.free_energy_cutoff - 0.5*j, g);
		fdata = GenotypeUtil::translateAndFold( b, g2 );
		cout << ">>>" << tab << fdata.first << endl;
		if (ntDist) {
			getddGDistribution(b, g2, arr);
		}
		else {
			getddGDistributionAA(b, g2, arr);
		}			
		for (int i=0; i<n; i++) {
			cout << arr[i] << endl;
		}
	}
	delete [] arr;
}

int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );
	dGDistExperiment(p);
}




