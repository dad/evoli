#include "mutator.hh"


int main()
{
	SimpleMutator m( .5 );
	CodingDNA s( "AAAGGG" );
	cout << s << endl;
	for( int i=0; i<20; i++ ){
		m.mutate( s );
		cout << s << endl;
	}
}
