#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "genotype-util.hh"

#include <cmath>

int main()
{
        const int N=50;
        const int size = 4;
        const int length = size*size;

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();

        // initialize the fitness evaluator
        ProteinFreeEnergyFitness fe( &b );

        // initialize the population
        Population p( N );
        p.init( GenotypeUtil::createRandomGenotype( length ), &fe, 0.01 );

        int equil_time = 100;
        int window_size = 30;

        for ( int i=0; ; i++ )
        {
                cout << i*100+1 << endl;
                for ( int j=0; j<100; j++ )
                {
                        p.evolve();
                        //cout << i+1 << " " << log( p.getAveFitness() ) << endl;
                }
                p.prepareCoalescenceCalcs();
                if ( p.calcCoalescenceTime() > equil_time + window_size )
                        break;
        }

        p.printGenebank( cout );
        //  cout << p.calcCoalescenceTime() << endl;

        double ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop;

		p.analyzeDnDs( window_size, ave_dn, ave_ds, ave_N, ave_S, ave_f, ave_fop, ErrorproneTranslation::m_codon_cost );

        cout << ave_dn << " " << ave_ds << " " << ave_N << " " << ave_S << " " << ave_f << endl;
}






