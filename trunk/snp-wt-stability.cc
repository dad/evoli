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
	string eval_type;
	int N;
	double ca_cost;
	double error_rate;
	double accuracy_weight;
	double error_weight;
	mutable int structure_ID;
	double free_energy_cutoff;
	double free_energy_minimum;
	int random_seed;
	string genotype_file_name;
	string snp_file_name;
	bool valid;


	Parameters( int ac, char **av ) {
		if ( ac != 13 )	{
			valid = false;
			return;
		}

		int i = 1;
		eval_type = av[i++];
		N = atoi( av[i++] );
		ca_cost = atof( av[i++] );
		error_rate = atof( av[i++] );
		accuracy_weight = atof( av[i++] );
		error_weight = atof( av[i++] );
		structure_ID = atoi( av[i++] );
		free_energy_cutoff = atof( av[i++] );
		free_energy_minimum = atof( av[i++] );
		random_seed = atoi( av[i++] );
		snp_file_name = av[i++];
		genotype_file_name = av[i++];

		valid = true;
	}
};


ostream & operator<<( ostream &s, const Parameters &p )
{
	s << "# Parameters:" << endl;
	s << "#   evaluation type: " << p.eval_type << endl;
	s << "#   structure id: " << p.structure_ID << endl;
	s << "#   free energy maximum: " << p.free_energy_cutoff << endl;
	s << "#   free energy minimum: " << p.free_energy_minimum << endl;
	s << "#   codon adaptation cost factor: " << p.ca_cost << endl;
	s << "#   transl. error rate per codon: " << p.error_rate << endl;
	s << "#   transl. accuracy weight: " << p.accuracy_weight << endl;
	s << "#   transl. error weight: " << p.error_weight << endl;
	s << "#   population size N: " << p.N << endl;
	s << "#   random seed: " << p.random_seed << endl;
	s << "#   SNP file: " << p.snp_file_name << endl;
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

void SNPAssayExperiment(Parameters& p)
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
	ifstream fin(p.snp_file_name.c_str());
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

	// Wild-type results
	vector<RunRecord> wtRunResults;
	ifstream fin(p.genotype_file_name.c_str());
	char buf[100];
	// Skip two lines
	fin.getline(buf,100);
	fin.getline(buf,100);
	while (!fin.eof()) {
		RunRecord rec;
		fin >> rec.cost_id >> rec.runNumber >> rec.gene;
		if (rec.gene.size() > 0) {
			rec.cost = pow(10.0,atof(rec.cost_id.c_str()));
			runResults.push_back(rec);
			cout << rec;
		}
	}
	fin.close();

	cout << "# Read " << runResults.size() << " results" << endl;
	cout << p;
	cout << "tr\texpr\trun\tdG\tfitness\ts\tns\tsnpcfacc\tsnpcfrob\tsnpcftrunc\tsnpcffold\tgene" << endl;
	vector<RunRecord>::iterator it = runResults.begin();

	double last_cost = -1;
	for (; it != runResults.end(); it++) {
		RunRecord rec = *it;

		if (rec.cost < 0.01 || rec.cost > 100) {
			continue;
		}

		ErrorproneTranslation* fe = NULL;
		ErrorproneTranslation* ept = new ErrorproneTranslation();
		ept->init( &folder, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;

		if (rec.cost != last_cost) {
			delete fe;
			if (p.eval_type == "tr") {
				ErrorproneTranslation* ept = new ErrorproneTranslation();
				ept->init( &folder, p.structure_ID, p.free_energy_cutoff, rec.cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
				fe = ept;
			}
			else if (p.eval_type == "acc") {
				AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation();
				afe->init( &folder, p.structure_ID, p.free_energy_cutoff, rec.cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
				fe = afe;
			}
			else if (p.eval_type == "rob") {
				RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation();
				rob->init( &folder, p.structure_ID, p.free_energy_cutoff, rec.cost, p.ca_cost, p.error_rate );
				fe = rob;
			}
			else if (p.eval_type == "con") {
				StabilityConstraint* sc = new StabilityConstraint();
				sc->init( &folder, p.structure_ID, p.free_energy_cutoff, p.free_energy_minimum, rec.cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
				fe = sc;
			}
			if (!fe) {
				cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
				exit(1);
			}
		}

		//cout << rec << endl;
		pair<double,int> fdata = GenotypeUtil::translateAndFold(folder, rec.gene);
		double cffold, cfacc, cfrob, cftrunc;
		fe->calcOutcomes(rec.gene, cfacc, cfrob, cftrunc, cffold);
		if (fdata.first <= -5) {
			cout << rec.cost_id << tab << rec.cost << tab << rec.runNumber << tab << fdata.first << tab
				 << rec.fitness << tab << rec.s << tab << rec.nonsynonymous << tab 
				 << cfacc << tab << cfrob << tab << cftrunc << tab << cffold << tab << rec.gene << endl;
		}
	}
}

int main( int ac, char **av)
{
	Parameters p( ac, av );
	if (p.valid) {
		SNPAssayExperiment(p);
	}
	else {
		cout << "Start program like this:" << endl;
		cout << "\t" << av[0] << " <eval type> <pop size> <ca cost> <error rate> <accuracy weight> <error weight> <structure id> <free energy cutoff> <free energy minimum> <random seed> <SNP file name> <gene file name>" << endl;
	}
}




