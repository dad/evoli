/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006, 2007 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/


/** \page tr-driver tr-driver
The program \c tr-driver is used to run actual translational robustness experiments. It takes a number of command-line parameters. A typical call is:
\verbatim
   ./tr-driver tr 25 1000  0.00 3 0.0112 33.4 57.7 599 -5 -100 0.00001 100000 50000 50 111
             0  1  2    3     4 5      6    7    8   9 10   11      12     13    14 15  16
\endverbatim
(The numbers under the call count the parameters.)

The parameters are, in order (the zeroth parameter is the program name):
-#  The type of %FitnessEvaluator (tr=standard translation)
-#  %Protein length in amino acids
-#  %Population size
-#  Base-10 logarithm of the translational cost factor
-#  Codon adaptation cost -- average accuracy ratio between optimal and non-optimal codons
-#  Translational error rate
-#  Accuracy weight
-#  Error weight
-#  Structure ID -- index of structure
-#  Maximum delta free energy of folding
-#  Minimum delta free energy of folding
-#  Mutation rate per nucleotide
-#  Window time -- number of generations for which evolutionary data is being collected
-#  Equilibration time -- evolution time in generations before collection of evolutionary data begins
-#  Number of replicates
-#  %Random number seed -- best to use an odd number

The error rate (parameter 6) and weights (parameters 7 and 8) are
chosen such that, given a particular structure and codon adaptation
cost (which affects translational accuracy), a random gene encoding a
protein that folds to that structure is accurately translated a
particular fraction of the time (say, 85% of the time).  These
parameters can be obtained by running the \ref get-weights program
as follows:
\verbatim
   ./get-weights 3 -5 111 stable-sequence.txt 10000 10 0.85 5
               0 1  2   3                   4     5  6    7 8
\endverbatim
-#  Codon adaptation cost -- average accuracy ratio between optimal and non-optimal codons
-#  Maximum free energy of folding
-#  %Random number seed -- best to use an odd number
-#  Seed gene sequence file containing nucleotide sequence that folds into the target structure with less than maximum free energy
-#  Equilibration time -- number of generations to equilibrate seed sequence by mutational drift
-#  Window time -- number of generations for which to measure random sequence statistics
-#  Target translational accuracy -- fraction of time a random gene is accurately translated
-#  Number of replicates -- averaged data will also be provided

A seed gene sequence may be generated using the \ref sequence-generator
program.  In the file \c stable-sequence.txt given as an argument to \c get-weights, the sequence must appear on its own line, preceded (if at all) only by comment lines beginning with \c # .

A suitable seed gene encoding a 25-aa protein folding into structure 599 with a stability of -5.106 kcal/mol is:
\verbatim
CUUGUCCUAAGGAGACCAUGCAACCGGAUUAACAGUUCAAUGCCGGACAUUUGGUUUCUAGCUCUGGACAAGAAG
\endverbatim



*/


#include "translation-experiment.hh"
#include "compact-lattice-folder.hh"

int main( int ac, char **av)
{
	Parameters p( ac, av );

	if (!p.valid) {
		exit(1);
	}

	// seed the random number generator
	Random::seed(p.random_seed);

	// initialize the protein folder
	int side_length = (int)(sqrt(float(p.protein_length)));
	assert(side_length*side_length == p.protein_length);
	CompactLatticeFolder folder(side_length);

	cout << p;
	// Create Polymerase based on input parameter p.mutation_rate
	double GCtoAT = 69.;
	double ATtoGC = 41.;
	double GCtoTA = 88.;
	double GCtoCG = 52.;
	double ATtoCG = 5.;
	double ATtoTA = 11.;
	//Polymerase poly(p.u, GCtoAT, ATtoGC, GCtoTA, GCtoCG, ATtoCG, ATtoTA );
	Polymerase poly(p.u);
	// Choose the FitnessEvaluator based on input parameters (p.eval_type).
	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation( &folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation( &folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation( &folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = rob;
	}
	else if (p.eval_type == "nu") {
		FoldingOnlyFitness* fof = new FoldingOnlyFitness( &folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = fof;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(1);
	}

	cout << setprecision(4);
	evolutionExperiment( p, *fe, poly);
	//evolutionTest( p, *fe );
	delete fe;

	return 0;
}




