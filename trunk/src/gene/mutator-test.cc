#include <iostream>

#include "mutator.hh"
#include "nucleotide-sequence.hh"

using namespace std;

void testPolymerase( const Polymerase& p, int reps, int L )
{
	vector<int> muts_per_site( L );
	vector<int> AtoCGT( 3 );
	vector<int> CtoGTA( 3 );
	vector<int> GtoTAC( 3 );
	vector<int> TtoACG( 3 );
	NucleotideSequence s( L );
	
	int count = 0;

	for ( int i=0; i<reps; i++ )
	{
		// test base 'A'
		s = NucleotideSequence( L, 'A' );
		p.mutate( s );
		for ( int k=0; k<L; k++ )
		{
			int j=-1;
			switch( s[k] ){
			case 'C': j=0; break;
			case 'G': j=1; break;
			case 'T': j=2; break;
			}
			if ( j >= 0 )
			{
				AtoCGT[j] += 1;
				muts_per_site[k] += 1;
				count += 1;
			}
		}
		// test base 'C'
		s = NucleotideSequence( L, 'C' );
		p.mutate( s );
		for ( int k=0; k<L; k++ )
		{
			int j=-1;
			switch( s[k] ){
			case 'G': j=0; break;
			case 'T': j=1; break;
			case 'A': j=2; break;
			}
			if ( j >= 0 )
			{
				CtoGTA[j] += 1;
				muts_per_site[k] += 1;
				count += 1;
			}
		}
		// test base 'G'
		s = NucleotideSequence( L, 'G' );
		p.mutate( s );
		for ( int k=0; k<L; k++ )
		{
			int j=-1;
			switch( s[k] ){
			case 'T': j=0; break;
			case 'A': j=1; break;
			case 'C': j=2; break;
			}
			if ( j >= 0 )
			{
				GtoTAC[j] += 1;
				muts_per_site[k] += 1;
				count += 1;
			}
		}
		// test base 'T'
		s = NucleotideSequence( L, 'T' );
		p.mutate( s );
		for ( int k=0; k<L; k++ )
		{
			int j=-1;
			switch( s[k] ){
			case 'A': j=0; break;
			case 'C': j=1; break;
			case 'G': j=2; break;
			}
			if ( j >= 0 )
			{
				TtoACG[j] += 1;
				muts_per_site[k] += 1;
				count += 1;
			}
		}
	}

	double m = static_cast<double>( count )/(L*reps*4.);
	cout << "Overall mutation rate: " << m << endl;
	cout << "\nRelative per-site mutation frequencies (should be close to 1.):" << endl;
	cout << "\tsite\tmutation frequency" << endl;
	for ( int i=0; i<L; i++ )
		cout << "\t" << i << "\t" << muts_per_site[i]*L/static_cast<double>( count ) << endl;
	cout << "\nRelative base-to-base mutation frequencies (sum to 1):" << endl;
	cout << "\tfrom\tto\tmutation frequency" << endl;
	cout << "\tA\tC\t" << AtoCGT[0]/static_cast<double>( count ) << endl;
	cout << "\tA\tG\t" << AtoCGT[1]/static_cast<double>( count ) << endl;
	cout << "\tA\tT\t" << AtoCGT[2]/static_cast<double>( count )<< endl;
	cout << "\tC\tG\t" << CtoGTA[0]/static_cast<double>( count ) << endl;
	cout << "\tC\tT\t" << CtoGTA[1]/static_cast<double>( count ) << endl;
	cout << "\tC\tA\t" << CtoGTA[2]/static_cast<double>( count ) << endl;
	cout << "\tG\tT\t" << GtoTAC[0]/static_cast<double>( count ) << endl;
	cout << "\tG\tA\t" << GtoTAC[1]/static_cast<double>( count ) << endl;
	cout << "\tG\tC\t" << GtoTAC[2]/static_cast<double>( count ) << endl;
	cout << "\tT\tA\t" << TtoACG[0]/static_cast<double>( count ) << endl;
	cout << "\tT\tC\t" << TtoACG[1]/static_cast<double>( count ) << endl;
	cout << "\tT\tG\t" << TtoACG[2]/static_cast<double>( count ) << endl;
}

int main()
{
	cout << "=== 1. Polymerase with uniform mutation frequency ===\n" << endl;

	Polymerase p( .1 );
	testPolymerase( p, 1000, 25 );

	cout << "\n\n=== 2. Polymerase with mutation frequencies of Bloom et al. 2006 ===\n" << endl;

	vector<double> AtoCGT;
	vector<double> CtoGTA;
	vector<double> GtoTAC;
	vector<double> TtoACG;
	AtoCGT.push_back( 6. );	// A -> C
	AtoCGT.push_back( 53. );	// A -> G
	AtoCGT.push_back( 25. );	// A -> T
	CtoGTA.push_back( 0. );	// C -> G
	CtoGTA.push_back( 10. );	// C -> T
	CtoGTA.push_back( 3. );	// C -> A
	GtoTAC.push_back( 3. );	// G -> T
	GtoTAC.push_back( 10. );	// G -> A
	GtoTAC.push_back( 0. );	// G -> C
	TtoACG.push_back( 25. );	// T -> A
	TtoACG.push_back( 53. );	// T -> C
	TtoACG.push_back( 6. );	// T -> G

	Polymerase p2( .1, AtoCGT, CtoGTA, GtoTAC, TtoACG );
	testPolymerase( p2, 1000, 25 );
}
