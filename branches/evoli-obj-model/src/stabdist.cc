#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <set>



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
	int num_to_fold;
	int structure_ID;
	int rand_seed;
	Genotype g;
};


Genotype findSequenceForTargetStructure( ProteinFolder &folder, double free_energy, int structure_ID);
void stabDistExperiment(Parameters& p);
void structureDrawingExperiment(Parameters& params);

ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	//s << "#   free energy cutoff: " << p.free_energy_cutoff << endl;
	s << "#   num. to fold: " << p.num_to_fold << endl;
	s << "#    random seed: " << p.rand_seed << endl;
	//s << "#   structure ID: " << p.structure_ID << endl;
	s << "#" << endl;
	return s;
}

Parameters getParams( int ac, char **av )
{
	if ( ac != 3 )
	{
		cout << "Start program like this:" << endl;
		cout << "  " << av[0] << " <num to fold> <rand seed>" << endl;
		exit (-1);
	}

	Parameters p;
	int i = 1;
	//p.free_energy_cutoff = atof( av[i++] );
	p.num_to_fold = atoi( av[i++] );
	p.rand_seed = atoi( av[i++] );
	//p.structure_ID = atoi( av[i++] );

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

void structureDrawingExperiment(Parameters& params) {
	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();
	int L = b.getProteinLength();
	int *p = new int[L];


	int arr[10] = {415, 414, 820, 873, 350, 19, 200, 55, 300, 1080};
	set<int, less<int> > v(arr,arr+10);
	/*for (int j=0; j<10; j++) {
		v.push_back(arr[j]);
	}*/

	for (int i=0; i<params.num_to_fold; i++) {
		for (int j=0; j<L; j++) {
			p[j] = (int)(myRand()*20);
		}
		b.foldProtein(p);
		int structID = b.getLastFoldedProteinStructureID();
		set<int, less<int> >::iterator ind = v.find(structID);
		if (ind != v.end()) {
			v.erase(ind);
			cout << "ID = " << structID << endl;
			b.printStructure(structID);
		}
	}
	delete [] p;
}

void stabDistExperiment(Parameters& params) {
	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();
	int L = b.getProteinLength();
	int *p = new int[L];

	for (int i=0; i<params.num_to_fold; i++) {
		for (int j=0; j<L; j++) {
			p[j] = (int)(myRand()*20);
		}
		double dG = b.foldProtein(p);
		int structID = b.getLastFoldedProteinStructureID();
		cout << i << tab << dG << tab << structID << endl;
	}
	delete [] p;
}

void stabHistExperiment(Parameters& params) {
	// initialize the protein folder
	ProteinFolder b(size);
	b.enumerateStructures();
	int L = b.getProteinLength();
	int *p = new int[L];

	for (int i=0; i<params.num_to_fold; i++) {
		for (int j=0; j<L; j++) {
			p[j] = (int)(myRand()*20);
		}
		double dG = b.foldProtein(p);
		int structID = b.getLastFoldedProteinStructureID();
		cout << i << tab << dG << tab << structID << endl;
	}
	delete [] p;
}

int main( int ac, char **av)
{
	Parameters p = getParams( ac, av );
	srand48(p.rand_seed);
	stabDistExperiment(p);
	//structureDrawingExperiment(p);
}




