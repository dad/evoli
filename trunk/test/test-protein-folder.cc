#include "protein-folder.hh"
#include "tools.hh"



int main()
{
        const int size = 5;
        int *p, *counts;
        int bins=20; // number of bins in free energy histogram

        ProteinFolder b(size);
        b.enumerateStructures();

        b.printStructure( 415 );
        vector<int> v = b.getSurface( 415 );
        copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
        cout << endl;
        b.printStructure( 524 );
        v = b.getSurface( 524 );
        copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
        cout << endl;
        b.printStructure( 335 );
        v = b.getSurface( 335 );
        copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
        cout << endl;
        b.printStructure( 848 );
        v = b.getSurface( 848 );
        copy( v.begin(), v.end(), ostream_iterator<int>( cout, " " ) );
        cout << endl;

        /*
        p = new int[size*size];
        counts = new int[4*bins];
        for ( int i=0; i<4*bins; i++ )
                counts[i]=0;
        int total=0;
        int sample_size = 1000;
        for ( int i=0; i<sample_size; i++ )
        {
                for ( int j=0; j<size*size; j++ )
                        p[j] = static_cast<int>( 20*myRand() );
                double G = b.foldProtein( p );
                int index = (int) ((G+1.)*bins);
                //    cout << G << " " << index << endl;
                if ( G>=-1. && G<3. )
                {
                        counts[index] += 1;
                        total += 1;
                }
        }
        for ( int i=0; i<4*bins; i++ )
                cout << (double) i/(double) bins - 1. << " " << (double) counts[i] / (double) total << " " <<counts[i] << endl;
        //  cout << total << endl;

       delete [] p;
       delete [] counts;
       */
}




