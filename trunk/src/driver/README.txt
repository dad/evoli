
Quick start
===========
How to run the tr-driver script:

             0  1  2    3     4 5      6    7    8   9 10   11      12     13    14 15  16
   ./tr-driver tr 25 1000  0.00 3 0.0112 33.4 57.7 599 -5 -100 0.00001 100000 50000 50 111

The parameters are, in order:
0)  The name of the program
1)  The type of FitnessEvaluator (tr=standard translation)
2)  Protein length in amino acids
3)  Population size
4)  Base-10 logarithm of the translational cost factor
5)  Codon adaptation cost -- average accuracy ratio between optimal 
	and non-optimal codons
6)  Translational error rate
7)  Accuracy weight
8)  Error weight
9)  Structure ID -- index of structure
10) Maximum free energy of folding
11) Minimum free energy of folding
12) Mutation rate per nucleotide
13) Window time -- number of generations for which evolutionary 
	data is being collected
14) Equilibration time -- evolution time in generations before 
	collection of evolutionary data begins
15) Number of replicates
16) Random number seed -- best to use an odd number

The error rate (parameter 6) and weights (parameters 7 and 8) are
chosen such that, given a particular structure and codon adaptation
cost (which affects translational accuracy), a random gene encoding a
protein that folds to that structure is accurately translated a
particular ffraction of the time (say, 85% of the time).  These
parameters can be obtained by running the get-weights program
(src/utils/get-weights) as follows:

               0 1  2   3                   4     5  6    7 8
   ./get-weights 3 -5 111 stable-sequence.txt 10000 10 0.85 5

0)  The name of the program
1)  Codon adaptation cost -- average accuracy ratio between optimal and non-optimal codons
2)  Maximum free energy of folding
3)  Random number seed -- best to use an odd number
4)  Seed gene sequence file containing nucleotide sequence that folds into the target structure with less than maximum free energy
5)  Equilibration time -- number of generations to equilibrate seed sequence by mutational drift
6)  Window time -- number of generations for which to measure random sequence statistics
7)  Target translational accuracy -- fraction of time a random gene is accurately translated
8)  Number of replicates -- averaged data will also be provided

A seed gene sequence may be generated using the sequence-generator
program (src/utils/sequence-generator).  The sequence should appear on
its own line, preceded (if at all) only by comment lines beginning
with "#".

A suitable seed gene encoding a 25-aa protein folding into structure 599 with a stability of -5.106 kcal/mol is:

CUUGUCCUAAGGAGACCAUGCAACCGGAUUAACAGUUCAAUGCCGGACAUUUGGUUUCUAGCUCUGGACAAGAAG


