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
        double ca_cost;
        double transl_error_rate;
        int num_to_fold;
        char *genotype_file;
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
        s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
        s << "#   transl. error rate: " << p.transl_error_rate << endl;
        s << "#   num. proteins to attempt to fold: " << p.num_to_fold << endl;
        s << "#   genotype file: " << p.genotype_file << endl;
        s << "#" << endl;
        return s;
}

Parameters getParams( int ac, char **av )
{
        if ( ac != 6 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <ca cost> <transl. error rate> <free energy cutoff> <num. to fold> <genotype file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.ca_cost = atof( av[i++] );
        p.transl_error_rate = atof( av[i++] );
        p.free_energy_cutoff = atof( av[i++] );
        p.num_to_fold = atoi( av[i++] );
        p.genotype_file = av[i];

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


void evolStabExperiment(Parameters& p) {
	// Goal: read in each evolved genotype and compute its dG distribution

	// size of the lattice proteins is hardcoded
	const int size = 5;

	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();

	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(p.genotype_file);
	while (!fin.eof())
	{
		RunRecord rec;
		fin >> rec.cost >> rec.runNumber >> rec.gene;
		runResults.push_back(rec);
        //cout << "# " << rec.cost << tab << rec.runNumber << tab << rec.gene << endl;
	}
	fin.close();

	bool ntDist = false;
	int n = 0;
	if (ntDist) {
		n = runResults[0].gene.size()*9;
	}
	else {
		n = runResults[0].gene.size()*19;
	}
	double* arr = new double[n];

	vector<RunRecord>::iterator it = runResults.begin();
	int k = 0;
	for (; it != runResults.end() && k < 1000; it++, k++)
	{
		RunRecord& rec = *it;
		// find initial structure
		int structure_ID = getStructureID( b, rec.gene );
		if ( structure_ID < 0 ) {
			cerr << "Input sequence does not translate!" << endl;
			continue;
		}
		if (ntDist) {
			getddGDistribution(b, rec.gene, arr);
		}
		else {
			getddGDistributionAA(b, rec.gene, arr);
		}
		pair<double, int> fdata = GenotypeUtil::translateAndFold( b, rec.gene );
		cout << ">>>" << tab << fdata.first << tab << rec.cost << tab << rec.runNumber << endl;
		for (int i=0; i<n; i++) {
			cout << arr[i] << endl;
		}
		//cout << endl << endl;
	}
	delete [] arr;
}


int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );
	evolStabExperiment(p);
}




