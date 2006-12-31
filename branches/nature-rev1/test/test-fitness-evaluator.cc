#include "fitness-evaluator.hh"
#include "protein-folder.hh"
#include "genetic-code.hh"
#include "genotype-util.hh"
#include "codon.hh"


void printDNA( const vector<int> &v, int length )
{
        for ( int i=0; i<length; i++ )
        {
                CodonUtil::printCodon( cout, v[i] );
                cout << " ";
        }
        cout << endl;
}

void printProtein( const vector<int> &v, int length )
{
        for ( int i=0; i<length; i++ )
        {
                GeneticCodeUtil::printResidue( cout, v[i] );
                cout << " ";
        }
        cout << endl;
}

int main()
{
        const int size = 5;
        int length = size*size;
        double maxFreeEnergy = -1.0;

        ProteinFolder b(size);
        b.enumerateStructures();

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
        cout << "#Initial sequence found" << endl;

        // initialize the fitness evaluator
        ProteinStructureFitness fe( &b, structureID, maxFreeEnergy );

        int sample_size = 1000;
        count = 0;

        for ( int i=0; i<sample_size; i++ )
        {
                g2 = g;
                GenotypeUtil::nonsymPointMutation( g2 );
                if ( fe.getFitness( g2 ) > 0 )
                        count += 1;
        }

        double neutrality = (double) count / (double) sample_size;

        double cost = .1;
        ErrorproneTranslation fe2;
        fe2.init( &b, structureID, maxFreeEnergy, cost, 1, .002 );
        cout << g << endl;
        double f = fe2.getFitness( g );
        double sens = fe2.getLastSensitivity();
        double sens2 = fe2.getLastSensitivityNoStop();
        double neutrl2 = GenotypeUtil::calcNeutrality( b, g, maxFreeEnergy );
        cout << neutrality << " " << neutrl2 << endl << f << " " << sens << " " << sens2 << endl;

}

/*

int main(){
  const int size = 4;
  int length = size*size;
  vector<int> p;

  ProteinFolder b(size);
  b.enumerateStructures();

  ProteinFreeEnergyFitness fe( &b );

  p.reserve(length);

  int sample_size = 1000;
  for ( int i=0; i<sample_size; i++ ){
    for ( int j=0; j<length; j++ )
      p[j] = static_cast<int>( 64*myRand() );
    printDNA( p, length );
    printProtein( p, length );
    cout << fe.getFitness( p ) << endl << endl;
  }
}

*/


