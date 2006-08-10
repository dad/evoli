#include "codon.hh"
#include "genetic-code.hh"

#include <iostream>



int main()
{
        double U = .5;

        int c = 0, nc;
        double sdn = 0, sds = 0;

        for ( int i=0; i<10000; i++ )
        {
                nc = CodonUtil::mutateCodon( U, c );
                if ( nc != c )
                {
                        CodonUtil::printCodon(cout, c);
                        cout << " -> ";
                        CodonUtil::printCodon(cout, nc);
                        double dn, ds;
                        GeneticCodeUtil::calcDnDs( dn, ds, c, nc );
                        sdn += dn;
                        sds += ds;
                        cout << " " << dn << " " << ds << endl;
                }
                c = nc;
        }

        cout << "Total nonsynonymous subst.: " << sdn << " Total synonymous subst.: " << sds << endl;

        double syn = 0;
        double non_syn = 0;
        for ( int i=0; i<64; i++ )
        {
                if ( GeneticCodeUtil::geneticCode[i] >= 0 )
                {
                        double s = GeneticCodeUtil::calcSynonymousSites( i );
                        syn += s;
                        non_syn += 3-s;
                }
        }

        cout << "\n\nTotal nonsynonymous sites in all codons: " << non_syn << endl << "Total synonymous sites in all codons: " << syn << endl;
        cout << "Total number of sites: " << syn + non_syn << endl;
}


