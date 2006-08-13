#include "genotype.hh"
#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "translator.hh"
#include "codon.hh"
#include "genotype-util.hh"
#include "tools.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <ctime>

const int size = 5;

ostream& operator<<(ostream& os, const pair<int,int> p) {
	os << p.first << tab << p.second;
	return os;
}

int designability( int ac, char **av) {
	srand48(0);
	ProteinFolder folder( size );
	folder.enumerateStructures();

	vector<int> counts(folder.getNumStructures(), 0);
	vector<int> folded_counts(folder.getNumStructures(), 0);
	long start_time = clock();
	int tot_to_fold = atoi(av[1]);
	int folded_count = 0;
	int total_count = 0;

	int length = size*size;
	Translator t(0, length);
	int* seq = new int[length];
	double cutoff = 0.0;
	//for (int i=0; i<tot_to_fold; i++) {
	while (total_count < tot_to_fold) {
		Genotype g = GenotypeUtil::createRandomGenotypeNoStops(length);
		t.translate(g, seq);
		double dg = folder.foldProtein(seq);
		int struct_id = folder.getLastFoldedProteinStructureID();
		if (dg <= cutoff) {
			folded_counts[struct_id]++;
			folded_count++;
		}
		counts[struct_id]++;
		total_count++;
	}

	vector<pair<int,int> > designabilities;
	for (unsigned int i=0; i<counts.size(); i++) {
		pair<int,int> p(counts[i], i);
		designabilities.push_back(p);
	}
	sort(designabilities.begin(), designabilities.end(), greater<pair<int,int> >());
	int rank = 0;
	cout << "rank\tsid\tcount\tfcount" << endl;
	for (vector<pair<int,int> >::iterator it=designabilities.begin(); it!=designabilities.end(); it++) {
		rank +=1;
		cout << rank << tab << (*it).second << tab << (*it).first << tab << folded_counts[(*it).second] << endl;
	}
	//copy(designabilities.begin(), designabilities.end(), ostream_iterator<pair<int,int> >(cout, "\n"));

	cout << "# Enumeration of " << total_count << " random proteins took " << (double)(clock()-start_time)/CLK_TCK << " seconds." << endl;
	cout << "# Found " << folded_count << " proteins below cutoff of dG = " << cutoff << "." << endl;

	delete [] seq;
	return 0;
}

int thresh_comparison( int ac, char **av) {
	srand48(0);
	ProteinFolder folder( size );
	folder.enumerateStructures();

	vector<int> counts(folder.getNumStructures(), 0);
	vector<int> folded_counts(folder.getNumStructures(), 0);
	long start_time;
	long run_start_time = clock();
	double tot_secs_full=0.0, tot_secs_approx=0.0;
	double tdiff = 0.0;
	bool below_cutoff_full, below_cutoff_approx;

	int tot_to_fold = atoi(av[1]);
	int target_struct_id = atoi(av[2]);
	double cutoff = atof(av[3]);

	int folded_count = 0;
	int total_count = 0;
	int failures = 0;

	int length = size*size;
	Translator t(0, length);
	int* seq = new int[length];
	double dg = 0.0;
	int struct_id = -1;
	//for (int i=0; i<tot_to_fold; i++) {
	while (folded_count < tot_to_fold) {
		Genotype g = GenotypeUtil::createRandomGenotypeNoStops(length);
		t.translate(g, seq);
		start_time = clock();
		below_cutoff_approx = folder.isFoldedBelowThreshold(seq, target_struct_id, cutoff);
		tdiff = (clock()-start_time)/(double)CLK_TCK;
		if (below_cutoff_approx) {
			tot_secs_approx += tdiff;
			start_time = clock();
			dg = folder.foldProtein(seq);
			struct_id = folder.getLastFoldedProteinStructureID();
			below_cutoff_full = (dg <= cutoff) && struct_id == target_struct_id;
			tdiff = (clock()-start_time)/(double)CLK_TCK;
			tot_secs_full += tdiff;
			if (below_cutoff_approx != below_cutoff_full) {
				failures++;
			}
			folded_count++;
		}
		total_count++;
	}
	double total_time = (clock()-run_start_time)/(double)CLK_TCK;

	cout << "# Assaying " << total_count << " random proteins (" << folded_count << " folded) took:" << endl;
	cout << "#     Full: " << tot_secs_full << " seconds." << endl;
	cout << "#   Approx: " << tot_secs_approx << " seconds." << endl;
	cout << "#  Overall: " << total_time << " seconds." << endl;
	cout << "# Savings: " << tot_secs_full/tot_secs_approx << " fold." << endl;
	cout << "# Failed " << failures << " times."  << endl;

	delete [] seq;
	return 0;
}

// finds a random sequence with folding energy smaller than cutoff.
Genotype getSequence( ProteinFolder &b, double free_energy_cutoff)
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
	while( G > free_energy_cutoff );
	return g;
}

void getRandomWeights(int ac, char**av) {
	srand48(clock());
	ProteinFolder folder( size );
	folder.enumerateStructures();

	AccuracyOnlyTranslation ept;
	int struct_id = 574;
	ept.init( &folder, struct_id, -5, 1, 6, 0.01, 10, 10);
	double error_rate;
	double accuracy_weight;
	double error_weight;
	//istringstream gene_str("UGCAAGAAAUGUGAUAUGUUGUGCGUCGACGAUGAUAAGGACUGUAAGGAUAAACGUGUAAAAAUCUGCCCCAAG");
	Genotype seed_genotype = GenotypeUtil::getSequenceForStructure(folder, -5, struct_id);
	//gene_str >> seed_genotype;
	cout << "Got seed." << endl;
	for (int i=0; i<5; i++) {
		ept.getWeightsForTargetAccuracy(seed_genotype, 0.85, error_rate, accuracy_weight, error_weight, 10000, 100000);
		cout << struct_id << tab << error_rate << tab << accuracy_weight << tab << error_weight << endl;
	}
}

int get_transitions(const vector<int>& v) {
	int prev = v[0];
	int num_changes = 0;
	for (unsigned int i=0; i<v.size(); i++) {
		if (v[i] != prev) {
			num_changes++;
		}
		prev = v[i];
	}
	return num_changes;
}

void surface(int ac, char**av) {
	// Read in designabilities file
	// Print surface and core for each structure ID, according to rank
	// Spit out the number of surface/core changes.
	srand48(clock());
	ProteinFolder folder( size );
	folder.enumerateStructures();

	ifstream fin("../analysis/designabilities.txt", ifstream::in);
	char buf[500];
	fin.getline(buf,500);
	int rank, sid, count, fcount;
	while (!fin.eof()) {
		//rank	sid	count	fcount
		fin >> rank >> sid >> count >> fcount;
		vector<int> surf = folder.getSurface(sid);
		int changes = get_transitions(surf);
		cout << rank << tab << sid << tab << changes << tab;
		copy(surf.begin(), surf.end(), ostream_iterator<int>(cout, ""));
		cout << endl;
	}
}

int main(int ac, char**av) {
	//designability(ac,av);
	//thresh_comparison(ac,av);
	//surface(ac, av);
	//getRandomWeights(ac, av);

	ProteinFolder folder( size );
	folder.enumerateStructures();
	srand48(time(NULL));
	for (int i=0; i<10; i++) {
		Genotype g = GenotypeUtil::getSequence(folder, -5);
		pair<double,int> fdata = GenotypeUtil::translateAndFold(folder, g);
		cout << fdata.first << tab << fdata.second << endl;
	}
	return 0;
}


