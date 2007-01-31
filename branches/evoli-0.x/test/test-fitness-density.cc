#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "genetic-code.hh"
#include "genotype-util.hh"
#include "codon.hh"
#include <fstream>

const char* tab = "\t";

struct RunRecord
{
	double cost;
	int runNumber;
	Genotype gene;
};

ostream& operator<<(ostream& os, const RunRecord& rec) {
	os << rec.cost << tab << rec.runNumber << tab << rec.gene << endl;
	return os;
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

int main(int argc, char**argv) {
	// size of the lattice proteins is hardcoded
	const int size = 5;

	// Seed random number generator
	long seconds = (long)time(NULL);
	srand48(seconds);

	// Read the results.
	vector<RunRecord> runResults;
	ifstream fin(argv[1]);
	char buf[500];
	// Skip two lines
	fin.getline(buf,500);
	fin.getline(buf,500);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.cost >> rec.runNumber >> rec.gene;
		if (rec.gene.size() > 1 && rec.cost >= 0.01 && rec.cost <= 100) {
			runResults.push_back(rec);
			//cout << rec;
		}
	}
	fin.close();

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	FitnessDensityEvaluator fde;

	cout << "tr\trun\tfdensns\tfdenss" << endl;
	vector<RunRecord>::iterator it = runResults.begin();
	int structure_ID = -1;
	structure_ID = getStructureID( folder, runResults[0].gene );
	if ( structure_ID < 0 )	{
		cerr << "Input sequence does not translate!" << endl;
	}

	for (; it != runResults.end(); it++) {
		RunRecord& rec = *it;
		if (rec.cost < 0.01 || rec.cost > 100 || rec.runNumber > 0) {
			continue;
		}
		// initialize evaluator
		ErrorproneTranslation ept;
		ept.init( &folder, structure_ID, -5, rec.cost, 6, 0.002 );
		double fitnessDensity = 0.0;
		fitnessDensity = fde.getFitnessDensity(rec.gene, ept, 1000);
		cout << rec.cost << tab << rec.runNumber << tab << fitnessDensity << endl; //tab << rec.gene << endl;
		//double fitnessDensityN = fde.getFitnessDensityNonsyn(rec.gene, ept, 1000);
		//double fitnessDensityS = fde.getFitnessDensitySyn(rec.gene, ept, 1000);
		//cout << rec.cost << tab << rec.runNumber << tab << fitnessDensityN << tab << fitnessDensityS << endl; //tab << rec.gene << endl;
	}
	return 0;
}

