#include "mutator.hh"
#include "genetic-code.hh"
#include "random.hh"

#include <iostream>
using namespace std;

////////////////////
// SimpleMutator methods
////////////////////

/*SimpleMutator::SimpleMutator(double mutation_rate) {
	m_mutation_rate = mutation_rate;
	}*/

SimpleMutator::SimpleMutator(double mutation_rate) {
	m_mutation_rate = mutation_rate;
}

SimpleMutator::~SimpleMutator() {
}

bool SimpleMutator::mutate(NucleotideSequence& seq) const {
	bool changed = false;
	const char* mutA = "CGT";
	const char* mutC = "GTA";
	const char* mutG = "TAC";
	const char* mutT = "ACG";
	
	if(seq.length() > 300) {
	  double mean = m_mutation_rate*seq.length();
	  int k = Random::rpois(mean);
	  if (k>0){
	    changed = true;
	  }
	  for (int l=0; l<k; l++) {
	    int i = Random::rint(seq.length());
	    int j = Random::rint( 3 );
	    switch( seq[i] ){
	    case 'A':
	      seq[i] = mutA[j]; break;
	    case 'C':
	      seq[i] = mutC[j]; break;
	    case 'G':
	      seq[i] = mutG[j]; break;
	    case 'T':
	      seq[i] = mutT[j]; break;
	    default:
	      assert( false ); // should never get here
	    }
	  }	
	}
	else {
	  for (int i=0; i<seq.length(); i++) {
	    if (Random::runif() < m_mutation_rate) {
	      changed = true;
	      int j = Random::rint( 3 );
	      switch( seq[i] ){
	      case 'A':
		seq[i] = mutA[j]; break;
	      case 'C':
		seq[i] = mutC[j]; break;
	      case 'G':
		seq[i] = mutG[j]; break;
	      case 'T':
		seq[i] = mutT[j]; break;
	      default:
		assert( false ); // should never get here
	      }
	    }
	  }
	}
	
	return changed;
}





////////////////////
// Polymerase methods
////////////////////

Polymerase::Polymerase(double mutation_rate) {
	m_mutation_rate = mutation_rate;
	// default is all mutations are equally likely
	m_AtoCGT = vector<double>( 3, 1. );
	m_CtoGTA = vector<double>( 3, 1. );
	m_GtoTAC = vector<double>( 3, 1. );
	m_TtoACG = vector<double>( 3, 1. );
	prepareMutationMatrix();
}

Polymerase::Polymerase(double mutation_rate, const vector<double>& AtoCGT, const vector<double>& CtoGTA,  const vector<double>& GtoTAC,  const vector<double>& TtoACG ) {
	m_mutation_rate = mutation_rate;
	m_AtoCGT = AtoCGT;
	m_CtoGTA = CtoGTA;
	m_GtoTAC = GtoTAC;
	m_TtoACG = TtoACG;
	prepareMutationMatrix();
}

Polymerase::Polymerase(double mutation_rate, double GCtoAT, double ATtoGC, double GCtoTA, double GCtoCG, double ATtoCG, double ATtoTA ) {
	m_mutation_rate = mutation_rate;
	m_AtoCGT.push_back(ATtoCG);
	m_AtoCGT.push_back(ATtoGC);
	m_AtoCGT.push_back(ATtoTA);
	m_CtoGTA.push_back(GCtoCG);
	m_CtoGTA.push_back(GCtoAT);
	m_CtoGTA.push_back(GCtoTA);
	m_GtoTAC.push_back(GCtoTA);
	m_GtoTAC.push_back(GCtoAT);
	m_GtoTAC.push_back(GCtoCG);
	m_TtoACG.push_back(ATtoTA);
	m_TtoACG.push_back(ATtoCG);
	m_TtoACG.push_back(ATtoGC);
	prepareMutationMatrix();
}

Polymerase::~Polymerase() {
}

void Polymerase::normalize( vector<double> & p )
{
	double p1, p2, p3;
	assert( p.size()==3 );
	p1 = p[0];
	p2 = p[1];
	p3 = p[2];
	double sum=p1+p2+p3;
	p[0]/=sum;
	p[1]=(p1+p2)/sum;
	p[2]=1.;
}

void Polymerase::prepareMutationMatrix()
{
	assert( m_AtoCGT.size() == 3 );
	assert( m_CtoGTA.size() == 3 );
	assert( m_GtoTAC.size() == 3 );
	assert( m_TtoACG.size() == 3 );
	m_mut_rate_ACGT.resize( 4 );

	m_mut_rate_ACGT[0] = m_AtoCGT[0] + m_AtoCGT[1] +m_AtoCGT[2];
	m_mut_rate_ACGT[1] = m_CtoGTA[0] + m_CtoGTA[1] +m_CtoGTA[2];
	m_mut_rate_ACGT[2] = m_GtoTAC[0] + m_GtoTAC[1] +m_GtoTAC[2];
	m_mut_rate_ACGT[3] = m_TtoACG[0] + m_TtoACG[1] +m_TtoACG[2];

	// first normalize "to" rates
	normalize( m_AtoCGT );
	normalize( m_CtoGTA );
	normalize( m_GtoTAC );
	normalize( m_TtoACG );

	// now calculate "from" rates
	double sum = m_mut_rate_ACGT[0] + m_mut_rate_ACGT[1] + m_mut_rate_ACGT[2] + m_mut_rate_ACGT[3];
	m_mut_rate_ACGT[0] *= 4*m_mutation_rate/sum;
	m_mut_rate_ACGT[1] *= 4*m_mutation_rate/sum;
	m_mut_rate_ACGT[2] *= 4*m_mutation_rate/sum;
	m_mut_rate_ACGT[3] *= 4*m_mutation_rate/sum;
}

int Polymerase::chooseMutation( const vector<double> & p ) const
{
	assert( p.size() == 3 );
	double r = Random::runif();
	if ( r<p[0])
		return 0;
	else if ( r<p[1] )
		return 1;
	else return 2;
}


bool Polymerase::mutate(NucleotideSequence& seq) const {
	bool changed = false;
	const char* mutA = "CGT";
	const char* mutC = "GTA";
	const char* mutG = "TAC";
	const char* mutT = "ACG";
	for (int i=0; i<seq.length(); i++) {
		switch( seq[i] ){
		case 'A':
			if (Random::runif() < m_mut_rate_ACGT[0] )
			{
				changed = true;
				seq[i] = mutA[ chooseMutation( m_AtoCGT ) ];
			}
			break;
		case 'C':
			if (Random::runif() < m_mut_rate_ACGT[1] )
			{
				changed = true;
				seq[i] = mutC[ chooseMutation( m_CtoGTA ) ];
			}
			break;
		case 'G':
			if (Random::runif() < m_mut_rate_ACGT[2] )
			{
				changed = true;
				seq[i] = mutG[ chooseMutation( m_GtoTAC ) ];
			}
			break;
		case 'T':
			if (Random::runif() < m_mut_rate_ACGT[3] )
			{
				changed = true;
				seq[i] = mutT[ chooseMutation( m_TtoACG ) ];
			}
			break;
		default:
			assert( false ); // should never get here
		}
	}
	return changed;

}
