#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "translator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <ctime>
#include <iomanip>

typedef unsigned int uint16;

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
	string cost_id;
	int runNumber;
	double fitness;
	double s;
	int nonsynonymous;
	Genotype gene;
};

ostream& operator<<(ostream& os, const RunRecord& rec) {
	os << rec.cost_id << tab << rec.cost << tab << rec.runNumber << tab << rec.gene << endl;
	return os;
}

class Parameters {
public:
	string genotype_file_name;
	bool valid;


	Parameters( int ac, char **av ) {
		if ( ac != 2 )	{
			valid = false;
			return;
		}

		int i = 1;
		genotype_file_name = av[i++];

		valid = true;
	}
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   genotype file: " << p.genotype_file_name << endl;
	s << "#" << endl;
	return s;
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

void misfoldDistExperiment(Parameters& p)
{
	// size of the lattice proteins is hardcoded
	const int size = 5;

	// Seed random number generator
	long seconds = (long)time(NULL);
	srand48(seconds);

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(p.genotype_file_name.c_str());
	char buf[100];
	// Skip two lines
	fin.getline(buf,100);
	fin.getline(buf,100);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.cost_id >> rec.runNumber >> rec.fitness >> rec.s >> rec.nonsynonymous >> rec.gene;
		if (rec.gene.size() > 0) {
			rec.cost = pow(10.0,atof(rec.cost_id.c_str()));
			runResults.push_back(rec);
			//cout << rec;
		}
	}
	fin.close();

	cout << "# Read " << runResults.size() << " results" << endl;
	cout << p;
	cout << "tr\texpr\trun\tdG\tfitness\ts\tns\tgene\tsnpcfacc\tsnpcfrob\tsnpcftrunc\tsnpcffold" << endl;
	vector<RunRecord>::iterator it = runResults.begin();

	for (; it != runResults.end(); it++) {
		RunRecord& rec = *it;
		//cout << rec << endl;
		if (rec.cost < 0.01 || rec.cost > 100) {
			continue;
		}
		pair<double,int> fdata = GenotypeUtil::translateAndFold(folder, rec.gene);
		if (fdata.first <= -5) {
			cout << rec.cost_id << tab << rec.cost << tab << rec.runNumber << tab << fdata.first << tab
				 << rec.fitness << tab << rec.s << tab << rec.nonsynonymous << tab << rec.gene << endl;
		}
	}
}

int main( int ac, char **av)
{
	Parameters p( ac, av );
	if (p.valid) {
		misfoldDistExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <gene file name>" << endl;
	}
}




