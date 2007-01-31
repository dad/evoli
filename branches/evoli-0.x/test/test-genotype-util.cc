#include "protein-folder.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <fstream>


struct Parameters
{
        double free_energy_cutoff;
        Genotype g;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   Free energy cutoff: " << p.free_energy_cutoff << endl;
        s << "#   initial sequence: " << p.g << endl;
        s << "#" << endl;
        return s;
}


Parameters getParams( int ac, char **av )
{
        if ( ac != 3 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <free energy cutoff> <sequence file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.free_energy_cutoff = atof( av[i++]);

        // read initial sequence
        ifstream in( av[i] );
        in >> p.g;

        return p;
}

void consistencyCheckTranslErrorTable()
{
        bool test_failed = false;
        for ( int c1 = 0; c1 < 64; c1++ )
        {

                if ( GeneticCodeUtil::geneticCode[c1] < 0 )
                        continue;


                for ( int res = 0; res < 20; res++ )
                {

                        bool subst_cand = false;

                        // check if there is at least one codon that corresponds to the given
                        // amino acid and can be converted into the target codon with a single
                        // substitution
                        for ( int c2 = 0; c2 < 64; c2++ )
                        {
                                if ( GeneticCodeUtil::geneticCode[c2] != res )
                                        continue;

                                int l1, l2, l3, l4, l5, l6;
                                CodonUtil::codonToLetters( l1, l2, l3, c1 );
                                CodonUtil::codonToLetters( l4, l5, l6, c2 );

                                int err_count = 0;
                                if ( l1 != l4 )
                                        err_count++;
                                if ( l2 != l5 )
                                        err_count++;
                                if ( l3 != l6 )
                                        err_count++;

                                if ( err_count == 1 )
                                {
                                        subst_cand = true;
                                        break;
                                }
                        }

                        int tentry = GeneticCodeUtil::singleSubstsTranslErrors[c1][res];

                        if ( GeneticCodeUtil::geneticCode[c1] == res )
                        {
                                if ( tentry == 1 )
                                {
                                        CodonUtil::printCodon( cout, c1 );
                                        cout << " " << res << " " << tentry;
                                        cout << " incorrect" << endl;
                                        test_failed = true;
                                }
                        }
                        else if ( ( subst_cand && tentry == 0 ) || ( !subst_cand && tentry == 1 ) )
                        {
                                CodonUtil::printCodon( cout, c1 );
                                cout << " " << res << " " << tentry;
                                cout << " incorrect" << endl;
                                test_failed = true;
                        }
                }
        }

        if ( test_failed )
                cout << "Consistency check failed!" << endl;
}

void analyzeSequence( ProteinFolder &b, const Parameters &p, ostream &s )
{
        int l = b.getProteinLength();

        s << p.g << endl;

        pair<int, int> data;
        for ( int i=0; i<l; i++ )
        {
                CodonUtil::printCodon( s, p.g[i] );

                int count = 0;
                for ( int j=0; j<20; j++ )
                        count += GeneticCodeUtil::singleSubstsTranslErrors[p.g[i]][j];

                data = GenotypeUtil::calcSiteNeutrality( b, p.g, p.free_energy_cutoff, i );

                s << " " << data.first << "/19 " << data.second << "/" << count << " ";

                if ( ErrorproneTranslation::m_codon_cost[p.g[i]] == 0 )
                        cout << "opt" << endl;
                else
                        cout << "not opt" << endl;
        }

}


/*
int main( int ac, char **av)
{
        Parameters p = getParams( ac, av );

        // size of the lattice proteins is hardcoded
        const int size = 5;

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();
        consistencyCheckTranslErrorTable();
        analyzeSequence( b, p, cout );
}
*/

void testGenotypeRandomizer()
{
        int l = 10;
        int *seq = new int[l];
        int *seq2 = new int[l];
        Translator t( 0, l );
        Genotype g = GenotypeUtil::createRandomGenotype( l );

        while ( !t.translateErrorFree( g, seq ) )
                g = GenotypeUtil::createRandomGenotype( l );

        //copy( g.begin(), g.end(), ostream_iterator<int>( cout, " " ) );
        cout << g << " ";
        cout << GenotypeUtil::calcFop( g, ErrorproneTranslation::m_codon_cost ) << endl;
        for ( int i=0; i<10; i++ )
        {
                g = GenotypeUtil::randomizeCodons( g, ErrorproneTranslation::m_codon_cost );
                //copy( g.begin(), g.end(), ostream_iterator<int>( cout, " " ) );
                cout << g << " ";
                cout << GenotypeUtil::calcFop( g, ErrorproneTranslation::m_codon_cost ) << endl;
                t.translateErrorFree( g, seq2 );

                bool equal=true;
                for ( int j=0; j<l; j++ )
                        if ( seq[j] != seq2[j] )
                                equal=false;
                if ( !equal )
                        cout << "randomization has changed amino acid sequence!" << endl;
        }

        delete [] seq;
        delete [] seq2;
}

int main()
{
        for ( int i=0; i<20; i++ )
        {
                testGenotypeRandomizer();
                cout << endl;
        }
}



