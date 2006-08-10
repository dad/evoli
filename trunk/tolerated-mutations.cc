#include "genotype.hh"
#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "translator.hh"
#include "codon.hh"
#include "genotype-util.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>

const int size = 5;

class Parameters {
public:
	double free_energy_cutoff;
	double free_energy_minimum;
	int random_seed;
	int structure_ID;
	string file_name;

	bool valid;

	Parameters() {
		valid = false;
	}

	Parameters(int ac, char**av, const string& usage) {
		if ( ac != 5 )	{
			cout << usage << endl;
			valid = false;
		}
		else {
			int i = 1;
			free_energy_cutoff = atof( av[i++] );
			free_energy_minimum = atof( av[i++] );
			random_seed = atoi( av[i++] );
			file_name = av[i++];
			valid = true;
		}
	}

};

int getStructureID( ProteinFolder &b, const Genotype &g ) {
	int l = b.getProteinLength();
	Translator t( 0, l );
	int *seq = new int[l];

	if ( t.translateErrorFree( g, seq ) ) {
		b.foldProtein( seq );
		return b.getLastFoldedProteinStructureID();
	}
	else
		return -1;
	delete [] seq;
}

/**
 * Calculates the neutrality of the given genotype.
 **/
void doMutagenesis( ProteinFolder &folder, const Genotype &gorig, int structure_ID, double max_free_energy, double tr, int rep ) {
	Genotype g = gorig;
	int l = folder.getProteinLength();
	Translator t( 0, l );
	int *seq = new int[l];

	int num_codons = g.size();
	// go through all codons in the gene
	for ( int aapos=0; aapos<num_codons; aapos++ ) {
		// go through all possible point mutations
		// (avoid operator %, which can be very slow)
		int from_codon = g[aapos];
		char from_aa = GeneticCodeUtil::residueLetter(from_codon);
		int nts[3];
		char ntcodes[3];
		CodonUtil::codonToLetters( nts[0], nts[1], nts[2], from_codon );
		for (int codonpos=0; codonpos<3; codonpos++) {
			ntcodes[codonpos] = CodonUtil::intToBase(nts[codonpos]);
		}
		// go through all positions in the codon
		for (int codonpos=0; codonpos<3; codonpos++) {
			// go through all nucleotides
			for (int nt=0; nt<4; nt++) {
				char to_nt = CodonUtil::intToBase(nt);
				char from_nt = ntcodes[codonpos];
				ntcodes[codonpos] = to_nt;
				int to_codon = CodonUtil::lettersToCodon(ntcodes[0], ntcodes[1], ntcodes[2]);
				if (to_codon != from_codon) {
					g[aapos] = to_codon;
					char to_aa = GeneticCodeUtil::residueLetter(to_codon);
					bool folded = false;
					if ( t.translateErrorFree( g, seq ) ) {
						double dG = folder.foldProtein( seq );
						int sid = folder.getLastFoldedProteinStructureID();
						folded = (sid == structure_ID && dG <= max_free_energy);
					}
					cout << tr << tab << rep << tab << (3*aapos+codonpos) << tab;
					CodonUtil::printCodon(cout, from_codon);
					cout << tab;
					CodonUtil::printCodon(cout, to_codon);
					cout << tab << aapos << tab << from_aa << tab << to_aa << tab << folded << endl;
					// Reset the genotype
					g[aapos] = gorig[aapos];
				}
				ntcodes[codonpos] = from_nt;
			}
		}
	}
	delete [] seq;
}

/**
 * Calculates the neutrality of the given genotype, checking all codon substitutions
 **/
void doMutagenesisFull( ProteinFolder &folder, const Genotype &gorig, int structure_ID, double max_free_energy, double tr, int rep ) {
	Genotype g = gorig;
	int l = folder.getProteinLength();
	Translator t( 0, l );
	int *seq = new int[l];

	int num_codons = g.size();
	// go through all codons in the gene
	for ( int aapos=0; aapos<num_codons; aapos++ ) {
		// go through all possible point mutations
		// (avoid operator %, which can be very slow)
		int from_codon = g[aapos];
		char from_aa = GeneticCodeUtil::residueLetter(from_codon);
		// go through all positions in the codon
		for (int to_codon=0; to_codon<64; to_codon++) {
			if (to_codon != from_codon) {
				g[aapos] = to_codon;
				char to_aa = GeneticCodeUtil::residueLetter(to_codon);
				bool folded = false;
				if ( t.translateErrorFree( g, seq ) ) {
					double dG = folder.foldProtein( seq );
					int sid = folder.getLastFoldedProteinStructureID();
					folded = (sid == structure_ID && dG <= max_free_energy);
				}
				cout << tr << tab << rep << tab << aapos << tab;
				CodonUtil::printCodon(cout, from_codon);
				cout << tab;
				CodonUtil::printCodon(cout, to_codon);
				cout << tab << aapos << tab << from_aa << tab << to_aa << tab << folded << endl;
				// Reset the genotype
				g[aapos] = gorig[aapos];
			}
		}
	}
	delete [] seq;
}



int main( int ac, char **av) {
	string prog_name = av[0];
	stringstream usage;
	usage << "Start program like this:\n  " << prog_name  << " <free energy cutoff> <free energy minimum> <random seed> <file>";
	Parameters p( ac, av, usage.str() );

	if (!p.valid) {
		return 0;
	}

	// initialize the protein folder
	ProteinFolder folder(size);
	folder.enumerateStructures();

	Genotype g;
	double tr = -1.0;
	int rep = -1;

	// Read file of genotypes
	ifstream fin;
	fin.open(p.file_name.c_str(), ifstream::in);
	char buf[1000];
	fin.getline(buf, 1000);
	fin.getline(buf, 1000);
	ostream& os = cout;
	os << "tr\trep\tntpos\tfrom_codon\tto_codon\taapos\tfrom_aa\tto_aa\tfolded" << endl;
	while (!fin.eof()) {
		fin >> tr >> rep >> g;
		double expr = pow(10.0,tr);
		int structure_ID = getStructureID(folder, g);
		if (expr < 0.01 || expr > 100) {
			continue;
		}
		os << "# " << tr << tab << rep << tab << g << tab;
		GenotypeUtil::printProtein(os, g);
		// Loop over all nucleotide mutants and output mutagenesis stats.
		doMutagenesisFull(folder, g, structure_ID, p.free_energy_cutoff, expr, rep);
	}
	fin.close();
}




