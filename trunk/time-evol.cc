#include "protein-folder.hh"
#include "fitness-evaluator.hh"
#include "population.hh"
#include "genotype-util.hh"

#include <cmath>
#include <cstdio>
#include <fstream>

int N=10;
const int size = 5;
double u = 0.01; // per-site mutation rate
double maxFreeEnergy = -5.0;
int maxTime = 2000;
int repetitions = 250;
int wsize = 300;
//const int repetitions = 1;


void doAnalysis( ProteinFolder &b, ofstream &dataFile )
{
        int length = size*size;
        // find a random sequence with folding energy smaller than cutoff
        double G, G2;
        Genotype g, g2;
        ProteinFreeEnergyFitness fetmp( &b );
        g = GenotypeUtil::createRandomGenotype( length );
        G = -log( fetmp.getFitness( g ) );
        int count = 1;
        do
        {
                g2 = g;
                GenotypeUtil::mutateGenotype( g2, 0.02 );
                G2 = -log( fetmp.getFitness( g2 ) );
                if ( G2 < G )
                {
                        g = g2;
                        G = G2;
                }
                count += 1;
                if ( count > 50000 )
                { // start again with random genotype if search is not successful
                        g = GenotypeUtil::createRandomGenotype( length );
                        G = -log( fetmp.getFitness( g ) );
                        count = 1;
                }
                //    cout << "# " << count << " " << G << endl;
        }
        while( G > maxFreeEnergy );
        int structureID = b.getLastFoldedProteinStructureID();
        //cout << "#Initial sequence found" << endl;

        // initialize the fitness evaluator
        ErrorproneTranslation fe( &b, structureID, maxFreeEnergy, .1, 0, .002 );

        // initialize the population
        Population p( N );
        // we fill the population with the genotype that we found above
        p.init( g, &fe, u );

        for ( int i=0; i<maxTime; i++ )
        {
                p.evolve();
                dataFile << i << " " << p.getAveFitness() << endl;
                cout << i << " " << p.getAveFitness() << endl;
        }
}


int main( int ac, char **av)
{
        if ( ac != 5 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <N> <u> <max time> <max free energy>" << endl;
                exit (-1);
        }

        N = atoi( av[1] );
        u = atof( av[2] );
        maxTime = atoi( av[3] );
        maxFreeEnergy = atof( av[4] );

        char filename[255];
        sprintf( filename, "N%iL%iu%gGmax%g-time-evol.dat", N, size*size, u, maxFreeEnergy );

        ofstream dataFile;
        dataFile.open( filename, ios::out );

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();

        // let the population evolve
        doAnalysis( b, dataFile );
}







