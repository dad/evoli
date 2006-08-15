#include "translation-experiment.hh"
#include "compact-lattice-folder.hh"

int main( int ac, char **av)
{
	Parameters p( ac, av );

	if (!p.valid) {
		exit(1);
	}

	// seed the random number generator
	srand48(p.random_seed);

	// initialize the protein folder
	int side_length = (int)(sqrt(p.protein_length));
	CompactLatticeFolder* folder = new CompactLatticeFolder(side_length);
	folder->enumerateStructures();

	cout << p;

	// Choose the FitnessEvaluator based on input parameters (p.eval_type).
	ErrorproneTranslation* fe = NULL;
	if (p.eval_type == "tr") {
		ErrorproneTranslation* ept = new ErrorproneTranslation();
		ept->init( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = ept;
	}
	else if (p.eval_type == "acc") {
		AccuracyOnlyTranslation* afe = new AccuracyOnlyTranslation();
		afe->init( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate, p.accuracy_weight, p.error_weight );
		fe = afe;
	}
	else if (p.eval_type == "rob") {
		RobustnessOnlyTranslation* rob = new RobustnessOnlyTranslation();
		rob->init( folder, p.protein_length, p.structure_ID, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate );
		fe = rob;
	}
	if (!fe) {
		cerr << "ERROR: unknown fitness evaluator type '" << p.eval_type << "'.  Exiting..." << endl;
		exit(1);
	}

	cout << setprecision(4);
	evolutionExperiment( p, *fe );
	//evolutionTest( p, *fe );
	delete fe;
	delete folder;

	return 0;
}




