#include "compact-lattice-folder.hh"
#include "genetic-code.hh"

#include <cassert>
#include <iostream>
#include <cmath>
#include <iterator>

// this defines the spacing in drawing protein structures
#define CHARS_PER_SITE 2


LatticeStructure::LatticeStructure()
		: m_size(0), m_structure(0)
{}

LatticeStructure::LatticeStructure( const char *structure, int size )
		: m_size(size), m_structure(0)
{
	m_structure = new char[3*size*size];
	strcpy( m_structure, structure );
	calcInteractingPairs();
}

LatticeStructure::LatticeStructure( const LatticeStructure &rhs )
		: m_size(rhs.m_size), m_structure(0)
{
	m_structure = new char[3*m_size*m_size];
	strcpy( m_structure, rhs.m_structure );
	calcInteractingPairs();
}

void LatticeStructure::calcInteractingPairs()
{
	m_interacting_pairs.clear();

	// take first the horizontal bonds
	for ( int i=0; i<m_size; i++ )
		for ( int j=0; j<m_size-1; j++ )
			if ( m_structure[i*(3*m_size-1)+2*j+1] == ' ' )
				m_interacting_pairs.push_back( pair<int,int> ( (int) m_structure[i*(3*m_size-1)+2*j], (int) m_structure[i*(3*m_size-1)+2*j+2] ) );
	// now the vertical bonds
	for ( int i=0; i<m_size-1; i++ )
		for ( int j=0; j<m_size; j++ )
			if ( m_structure[i*(3*m_size-1)+2*m_size-1+j] == ' ' )
				m_interacting_pairs.push_back( pair<int,int> ( (int) m_structure[i*(3*m_size-1)+2*j], (int) m_structure[i*(3*m_size-1)+3*m_size-1+2*j] ) );
}


LatticeStructure::~LatticeStructure()
{
	delete [] m_structure;
}


vector<int> LatticeStructure::getSurface() const
{
	int l = m_size*m_size;
	vector<int> v;
	v.resize( l );
	vector<int>::iterator it = v.begin(), e = v.end();
	for ( ; it != e; it++ )
		*it = 0;

	// first the top row
	for ( int i=0; i<m_size; i++ )
		v[m_structure[2*i]-1] = 1;
	// now the sides
	for ( int i=1; i<m_size; i++ )
	{
		v[m_structure[i*(3*m_size-1)]-1] = 1;
		v[m_structure[i*(3*m_size-1)+2*(m_size-1)]-1] = 1;
	}
	// finally the bottom row
	for ( int i=0; i<m_size; i++ )
		v[m_structure[(m_size-1)*(3*m_size-1)+2*i]-1] = 1;
	return v;
}

void LatticeStructure::draw(ostream& os, const char* prefix) const
{
	StructureUtil::drawStructure( os, m_structure, m_size, prefix );
	vector<pair<int,int> >::const_iterator it,e;
	it = m_interacting_pairs.begin();
	e = m_interacting_pairs.end();
	os << prefix;
	for ( ; it != e; it++ )
		os << "(" << (*it).first << ", " << (*it).second << ") ";
	os << endl << prefix;

	vector<int> v = getSurface();
	copy( v.begin(), v.end(), ostream_iterator<int>( os, " " ) );
	os << endl;
}

void StructureUtil::drawSite( ostream &s, int site )
{
	s.width( CHARS_PER_SITE );
	if ( site > 0 )
	{
		s.fill('0');
		s << site;
	}
	else
	{
		s.fill(' ');
		s << "+ ";
	}
}

void StructureUtil::drawSite( ostream &s, char site )
{
	s.width( CHARS_PER_SITE );
	if ( site > 0 )
	{
		s.fill('0');
		s << (int) site;
	}
	else
	{
		s.fill(' ');
		s << "+ ";
	}
}


void StructureUtil::drawVBond( ostream &s, int bond )
{
	s.width( CHARS_PER_SITE + 2 );
	s.fill(' ');
	if ( bond > 0 )
		s << "|   ";
	else
		s << "    ";
}

void StructureUtil::drawVBond( ostream &s, char bond )
{
	for ( int i=0; i<CHARS_PER_SITE-2; i++ )
		s << " ";
	s << bond << "   ";
}

void StructureUtil::drawHBond( ostream &s, int bond )
{
	if ( bond > 0 )
		s << "--";
	else
		s << "  ";
}

void StructureUtil::drawHBond( ostream &s, char bond )
{
	s << bond << bond;
}


void StructureUtil::drawStructure( ostream& os, const char *s, int size, const char *prefix )
{

	int k=0;
	os << prefix;
	for ( int i=0; i<size-1; i++ )
	{
		for ( int j=0; j<size-1; j++ )
		{
			drawSite( os, s[k++] );
			drawHBond( os, s[k++] );
		}
		drawSite( os, s[k++] );
		os << endl << prefix;
		for ( int j=0; j<size; j++ )
			drawVBond( os, s[k++] );
		os << endl << prefix;
		k-=size;
		for ( int j=0; j<size; j++ )
			drawVBond( os, s[k++] );
		os << endl << prefix;
	}

	for ( int j=0; j<size-1; j++ )
	{
		drawSite( os, s[k++] );
		drawHBond( os, s[k++] );
	}
	drawSite( os, s[k] );
	os << endl;
}


void StructureUtil::flipLeftRight( const char* s, char* d, int size )
{
	int di=0;

	int k=0;
	for ( int j=0; j<size-1; j++ )
	{
		for ( int i=2*size-2; i>=0; i-- )
			d[di++]=s[k+i];
		k+=2*size-1;
		for ( int i=size-1; i>=0; i-- )
			d[di++]=s[k+i];
		k+=size;
	}
	for ( int i=2*size-2; i>=0; i-- )
		d[di++]=s[k+i];

	d[di]=0;
}


void StructureUtil::rotate90( const char * s, char *d, int size )
{
	int di=0;

	for ( int i=0; i<size-1; i++ )
	{
		for ( int j=size-2; j>=0; j-- )
		{
			d[di++]=s[j*(3*size-1)+3*size-1+2*i];
			if ( s[j*(3*size-1)+2*size-1+i] == '|' )
				d[di++]='-';
			else
				d[di++]=' ';
		}
		d[di++]=s[2*i];
		for ( int j=size-1; j>=0; j-- )
		{
			if ( s[j*(3*size-1)+1+2*i] == '-' )
				d[di++]='|';
			else
				d[di++]=' ';
		}
	}
	for ( int j=size-2; j>=0; j-- )
	{
		d[di++]=s[j*(3*size-1)+3*size-1+2*size-2];
		if ( s[j*(3*size-1)+2*size-1+size-1] == '|' )
			d[di++]='-';
		else
			d[di++]=' ';
	}
	d[di++]=s[2*size-2];
	d[di]=0;
}


SelfAvoidingWalk::SelfAvoidingWalk( int n ) :  m_size(n)
{
	m_vbonds.resize(n);
	m_hbonds.resize(n);
	m_sites.resize(n);
	for ( int i=0; i<n; i++ )
	{
		m_vbonds[i].resize(n);
		m_hbonds[i].resize(n);
		m_sites[i].resize(n);
	}
	m_tmp_structure = new char[2*n*n];
	clear();
}

SelfAvoidingWalk::~SelfAvoidingWalk()
{
	delete [] m_tmp_structure;
}

void SelfAvoidingWalk::clear()
{
	for ( int i=0; i<m_size; i++ )
		for ( int j=0; j<m_size; j++ )
		{
			m_vbonds[i][j]=0;
			m_hbonds[i][j]=0;
			m_sites[i][j]=0;
		}
	m_start_x = 0;
	m_start_y = 0;
	m_cur_x = 0;
	m_cur_y = 0;
	m_dx = 1;
	m_dy = 0;
	m_length = 0;
	m_walk.clear();
}


void SelfAvoidingWalk::setStart( int x, int y )
{
	if ( m_length > 0 )
		cout << "Cannot start again. Walk is in progress!" << endl;
	else
	{
		m_start_x = m_cur_x = x;
		m_start_y = m_cur_y = y;
		m_sites[x][y]=++m_length;
	}
}


void SelfAvoidingWalk::eraseLastMove()
{
	m_sites[m_cur_x][m_cur_y] = 0;

	if ( m_dx == 1 )
	{
		m_cur_x--;
		m_hbonds[m_cur_x][m_cur_y] = 0;
	}
	else if ( m_dx == -1 )
	{
		m_hbonds[m_cur_x][m_cur_y] = 0;
		m_cur_x++;
	}
	else if ( m_dy == 1 )
	{
		m_cur_y--;
		m_vbonds[m_cur_x][m_cur_y] = 0;
	}
	else if ( m_dy == -1 )
	{
		m_vbonds[m_cur_x][m_cur_y] = 0;
		m_cur_y++;
	}

	int dtmp;

	m_length--;
	switch ( m_walk[m_length-1] )
	{
	case left:
		dtmp = m_dx;
		m_dx = -m_dy;
		m_dy = dtmp;
		break;
	case right:
		dtmp = m_dx;
		m_dx = m_dy;
		m_dy = -dtmp;
		break;
	case forward:
		break;
	}
	m_walk.pop_back();
}

bool SelfAvoidingWalk::doMove( Direction d )
{
	int dx, dy;

	bool oob_error = false;
	bool intersect_error = false;

	switch ( d )
	{
	case left:
		dx = m_dy;
		dy = -m_dx;
		break;
	case right:
		dx = -m_dy;
		dy = m_dx;
		break;
	case forward:
		dx = m_dx;
		dy = m_dy;
		break;
	default:
		dx = m_dx;
		dy = m_dy;
	}

	if ( dx == 1 )
	{
		if ( m_cur_x + 1 >= m_size ) oob_error = true;
		else if ( m_sites[m_cur_x+1][m_cur_y] > 0 ) intersect_error = true;
		else
		{
			m_hbonds[m_cur_x][m_cur_y] = 1;
			m_cur_x++;
			m_sites[m_cur_x][m_cur_y] = ++m_length;
		}
	}
	else if ( dx == -1 )
	{
		if ( m_cur_x - 1 < 0 ) oob_error = true;
		else if ( m_sites[m_cur_x-1][m_cur_y] > 0 ) intersect_error = true;
		else
		{
			m_cur_x--;
			m_hbonds[m_cur_x][m_cur_y] = 1;
			m_sites[m_cur_x][m_cur_y] = ++m_length;
		}
	}
	else if ( dy == 1 )
	{
		if ( m_cur_y + 1 >= m_size ) oob_error = true;
		else if ( m_sites[m_cur_x][m_cur_y+1] > 0 ) intersect_error = true;
		else
		{
			m_vbonds[m_cur_x][m_cur_y] = 1;
			m_cur_y++;
			m_sites[m_cur_x][m_cur_y] = ++m_length;
		}
	}
	else if ( dy == -1 )
	{
		if ( m_cur_y - 1 < 0 ) oob_error = true;
		else if ( m_sites[m_cur_x][m_cur_y-1] > 0 ) intersect_error = true;
		else
		{
			m_cur_y--;
			m_vbonds[m_cur_x][m_cur_y] = 1;
			m_sites[m_cur_x][m_cur_y] = ++m_length;
		}
	}

	if ( oob_error )
	{
		//    cout << "walk out of bounds!" << endl;
		return false;
	}
	if ( intersect_error )
	{
		//cout << "walk intersects itself!" << endl;
		return false;
	}

	m_dx = dx;
	m_dy = dy;
	m_walk.push_back(d);
	return true;
}


bool SelfAvoidingWalk::setString( int x, int y, const char * string )
{
	clear();
	setStart( x, y );
	Direction d;
	bool end = false;

	for ( int i=0; i<m_size*m_size; i++ )
	{
		switch ( string[i] )
		{
		case 'l':
			d = left;
			break;
		case 'r':
			d = right;
			break;
		case 'f':
			d = forward;
			break;
		default:
			d = forward;
			end = true;
		}
		if ( end ) return false;

		if ( !doMove( d ) )
			return false;
	}
	return true; // returns true only if no error occured.
}



void SelfAvoidingWalk::draw(ostream& os) const
{
	os << "Start: ("<< m_start_x << ", " << m_start_y << "); Length: " << m_length << endl;
	os << "Pos: ("<< m_cur_x << ", " << m_cur_y << "); Direction: (" << m_dx << ", " << m_dy << ")" << endl;
	for ( int i=0; i<m_size-1; i++ )
	{
		for ( int j=0; j<m_size-1; j++ )
		{
			StructureUtil::drawSite( os, m_sites[j][i] );
			StructureUtil::drawHBond( os, m_hbonds[j][i] );
		}
		StructureUtil::drawSite( os, m_sites[m_size-1][i] );
		os << endl;
		// draw the vbonds twice
		for ( int j=0; j<m_size; j++ )
			StructureUtil::drawVBond( os, m_vbonds[j][i] );
		os << endl;
		for ( int j=0; j<m_size; j++ )
			StructureUtil::drawVBond( os, m_vbonds[j][i] );
		os << endl;
	}
	for ( int j=0; j<m_size-1; j++ )
	{
		StructureUtil::drawSite( os, m_sites[j][m_size-1] );
		StructureUtil::drawHBond( os, m_hbonds[j][m_size-1] );
	}
	StructureUtil::drawSite( os, m_sites[m_size-1][m_size-1] );
	os << endl;
}



void SelfAvoidingWalk::getStructure( char * d ) const
{
	int di = 0;

	for ( int i=0; i<m_size-1; i++ )
	{
		for ( int j=0; j<m_size-1; j++ )
		{
			d[di++] = (char) m_sites[j][i];
			if ( m_hbonds[j][i] == 1 )
				d[di++]='-';
			else
				d[di++]=' ';
		}
		d[di++] = (char) m_sites[m_size-1][i];
		for ( int j=0; j<m_size; j++ )
		{
			if ( m_vbonds[j][i] == 1 )
				d[di++]='|';
			else
				d[di++]=' ';
		}
	}
	for ( int j=0; j<m_size-1; j++ )
	{
		d[di++] = (char) m_sites[j][m_size-1];
		if ( m_hbonds[j][m_size-1] == 1 )
			d[di++]='-';
		else
			d[di++]=' ';
	}
	d[di++] = (char) m_sites[m_size-1][m_size-1];

	d[di]=0;
}


CompactLatticeFolder::CompactLatticeFolder( int size )
	: m_size( size ), m_num_structures( 0 ), m_num_folded( 0 )
{
	if ( m_size > 15 )
	{
		cout << "Maximally suported size: 15 (because char is being used to enumerate sites" << endl;
		exit(-1);
	}
	m_last_folded_structure = -1;

	m_ffw_struct = new char[3*m_size*m_size];
	m_ss_struct = new char[3*m_size*m_size];
	m_ss_struct2 = new char[3*m_size*m_size];
	
	enumerateStructures();
}

CompactLatticeFolder::~CompactLatticeFolder()
{
	delete [] m_ffw_struct;
	delete [] m_ss_struct;
	delete [] m_ss_struct2;

	vector<LatticeStructure *>::iterator it = m_structures.begin();

	for ( ; it != m_structures.end(); it++ )
		delete (*it);
}

void CompactLatticeFolder::findFillingWalks( SelfAvoidingWalk &w, int &moves )
{
	if ( w.length() == w.maxLength() )
	{
		w.getStructure( m_ffw_struct );
		storeStructure( m_ffw_struct );
		return;
	}

	moves++;
	if ( w.doMove(SelfAvoidingWalk::forward) )
	{
		findFillingWalks( w, moves );
		w.eraseLastMove();
	}

	moves++;
	if ( w.doMove(SelfAvoidingWalk::left) )
	{
		findFillingWalks( w, moves );
		w.eraseLastMove();
	}

	moves++;
	if ( w.doMove(SelfAvoidingWalk::right) )
	{
		findFillingWalks( w, moves );
		w.eraseLastMove();
	}
}


bool CompactLatticeFolder::findStructure( const char *s )
{
	//  cout << "Searching for structure:" << endl;
	//StructureUtil::drawStructure( s, m_size );
	StructureMapConstIterator it = m_structure_map.find( s );
	bool r = (it != m_structure_map.end());
	if ( r )
	{
		//    cout << "Found structure no "<< (*it).second<< ":" << endl;
		//    StructureUtil::drawStructure( s, m_size );
	}
	return r;
}


void CompactLatticeFolder::storeStructure( const char *s )
{
	//  cout << "Searching for structure:" << endl;
	//  StructureUtil::drawStructure( s, m_size );

	if ( findStructure( s ) )
	{
		return;
	}

	//  cout << "Trying 90 degree rotation" << endl;
	StructureUtil::rotate90( s, m_ss_struct, m_size );
	if ( findStructure( m_ss_struct ) )
	{
		return;
	}

	//  cout << "Trying 180 degree rotation" << endl;
	StructureUtil::rotate90( m_ss_struct, m_ss_struct2, m_size );
	if ( findStructure( m_ss_struct2 ) )
	{
		return;
	}

	//  cout << "Trying 270 degree rotation" << endl;
	StructureUtil::rotate90( m_ss_struct2, m_ss_struct, m_size );
	if ( findStructure( m_ss_struct ) )
	{
		return;
	}

	//  cout << "Trying r-l flip" << endl;
	StructureUtil::flipLeftRight( s, m_ss_struct, m_size );
	if ( findStructure( m_ss_struct ) )
	{
		return;
	}

	//  cout << "Trying r-l flip + 90 degree rotation" << endl;
	StructureUtil::rotate90( m_ss_struct, m_ss_struct2, m_size );
	if ( findStructure( m_ss_struct2 ) )
	{
		return;
	}

	//  cout << "Trying r-l flip + 180 degree rotation" << endl;
	StructureUtil::rotate90( m_ss_struct2, m_ss_struct, m_size );
	if ( findStructure( m_ss_struct ) )
	{
		return;
	}

	//  cout << "Trying r-l flip + 270 degree rotation" << endl;
	StructureUtil::rotate90( m_ss_struct, m_ss_struct2, m_size );
	if ( findStructure( m_ss_struct2 ) )
	{
		return;
	}
	//  cout << "Structure does not exist yet." << endl;
	//  cout << "Adding structure no "<< m_num_structures << endl;

	LatticeStructure *structure = new LatticeStructure( s, m_size );

	m_structures.push_back( structure );
	m_structure_map[structure->getStructure()]=m_num_structures++;
}



void CompactLatticeFolder::enumerateStructures()
{
	if ( m_num_structures > 0 )
	{
	//	cout << "Enumeration of structures has already been performed.\nWe are not going to do this again..." << endl;
		return;
	}

	SelfAvoidingWalk w(m_size);

	//cout << "#Enumerating all possible structures on " << m_size << "x" << m_size << " lattice" << endl;

	int moves = 0;
	for ( int i=0; i<m_size; i++ )
		for ( int j=0; j<m_size-1; j++ )
		{
//			 cout << "#Structures starting at (" << j << ", " << i << ")" << endl;
			w.clear();
			w.setStart( j, i );
			w.doMove( SelfAvoidingWalk::forward );
			moves++;
			findFillingWalks( w, moves );
		}

	vector<vector<int> > pairs;
	pairs.resize(m_size*m_size);
	for ( int i=0; i<m_size*m_size; i++ )
	{
		pairs[i].resize(m_size*m_size);
	}

	vector<LatticeStructure *>::const_iterator it = m_structures.begin();
	for ( ; it!=m_structures.end(); it++ )
	{
		//(*it)->draw();
		const vector<pair<int,int> > &p=(*it)->getInteractingPairs();
		vector<pair<int,int> >::const_iterator it2=p.begin();
		for ( ; it2!=p.end(); it2++ )
		{
			int a = (*it2).first;
			int b = (*it2).second;
			if ( a>b )
			{
				a=b; b=(*it2).first;
			}
			pairs[a-1][b-1] += 1;
		}
	}

	int count = 0;
	for ( int i=0; i<m_size*m_size; i++ )
	{
		for ( int j=0; j<m_size*m_size; j++ )
		{
			//cout.width(3);
			//cout.fill(' ');
			//cout << pairs[i][j] << " ";
			if ( pairs[i][j] > 0 )
				count++;
		}
		//cout << endl;
	}

//	 cout << "#Total number of potentially interacting pairs: " << count << endl;
}

/**
 * Fold the sequence and return information about the result (structure, free energy).
 */
FoldInfo CompactLatticeFolder::fold( const Sequence& s )
{
	assert( m_num_structures > 0 );

	double kT = 0.6;
	double minE = 1e50;
	int minIndex = 0;
	double Z = 0;
	double G;

	for ( int i=0; i<m_num_structures; i++ )
	{
		double E = 0;

		// calculate binding energy of this fold
		const vector<Contact> &pair_list = m_structures[i]->getInteractingPairs();
		vector<pair<int,int> >::const_iterator it=pair_list.begin();
		for ( ; it!=pair_list.end(); it++ )
		{
			double deltaE = contactEnergies[s[(*it).first-1]][s[(*it).second-1]];
			E += deltaE;
		}
		// check if binding energy is lower than any previously calculated one
		if ( E < minE )
		{
			minE = E;
			minIndex = i;
		}
		// add energy to partition sum
		Z +=  exp(-E/kT);

	}


	// calculate free energy of folding
	G = minE + kT * log( Z - exp(-minE/kT) );

	//cout <<  Z << " " << exp(-minE/kT) << " " << G << endl;
	//cout << "Folding energy: " << minE << endl;
	//cout << "Folding free energy: " << G << endl;

	// record the structure into which this protein folds
	m_last_folded_structure = minIndex;
	// increment folded count
	m_num_folded += 1;

	return FoldInfo(G, minIndex);
}

inline double CompactLatticeFolder::getEnergy(const Sequence& p, const int structID) const {
	const vector<Contact> &pair_list = m_structures[structID]->getInteractingPairs();
	vector<Contact>::const_iterator it=pair_list.begin();
	double E = 0.0;
	for ( ; it!=pair_list.end(); it++ ){
		double deltaE = contactEnergies[p[(*it).first-1]][p[(*it).second-1]];
		E += deltaE;
    }
	return E;
}

void CompactLatticeFolder::getMinMaxPartitionContributions(const Sequence& p, const int ci, double& cmin, double& cmax) const {
	double kT = 0.6;
	double min_cont = 1e5;
	double max_cont = -1e5;
	for (int i=0; i<m_size*m_size-3; i++) {
		for (int j=i+3; j<m_size*m_size; j++) {
			double e = contactEnergies[p[i]][p[j]];
			if (e < min_cont) {
				min_cont = e;
			}
			else if (e > max_cont) {
				max_cont = e;
			}
		}
	}
	int num_contacts = 2*(m_size-1)*m_size - m_size*m_size;
	cmin = exp(-max_cont*num_contacts/kT);
	cmax = exp(-min_cont*num_contacts/kT);
	//string tab = "\t";
	//cout << cmin << tab << cmax << tab << num_contacts << endl;
}

bool CompactLatticeFolder::isFoldedBelowThreshold( const Sequence& p, const int structID, double cutoff) const
{
	assert( m_num_structures > 0 );

	double kT = 0.6;
	double Zu = 0;

	double Ef = getEnergy(p, structID);
	double Zf = exp(-Ef/kT);
	double Zcutoff = Zf*exp(cutoff/kT);

	double cmin = 0.0;
	double cmax = 0.0;

	int rep_count = 0;
	int reps_per_minmax = 20;
	getMinMaxPartitionContributions(p, 0, cmin, cmax);
	for ( int ci=0; ci<m_num_structures; ci++ ) {
		if (ci != structID) {
			rep_count++;
			// calculate binding energy of this fold
			double E = getEnergy(p, ci);
			// bail out if structID is not the minimum-energy structure.
			if (E<=Ef) {
				//cout << "E <= Ef" << endl;
				return false;
			}
			Zu += exp(-E/kT);
			// bail out if the unfolded partition function is too large.
			if (Zu >= Zcutoff) {
				//cout << "Zu >= Zcutoff" << endl;
				return false;
			}
			if (true && rep_count > reps_per_minmax) {

				rep_count = 0;
				if (Zu + (m_num_structures-(ci+1))*cmin >= Zcutoff) {
					//cout << "Zu + (m_num_structures-ci)*cmin >= Zcutoff" << endl;
					return false;
				}
				if (Zu + (m_num_structures-(ci+1))*cmax < Zcutoff) {
					//cout << "Zu + (m_num_structures-ci)*cmax < Zcutoff" << endl;
					return true;
				}
			}
		}
	}


	// calculate free energy of folding
	double G = Ef + kT * log( Zu );

	return (G <= cutoff);
}

void CompactLatticeFolder::printContactEnergyTable( ostream &s ) const
{
	s << "Contact energies:" << endl;
	s << "      ";
	for ( int i=0; i<20; i++ )
		s << GeneticCodeUtil::residues[i] << "   ";
	s << endl;
	for ( int i=0; i<20; i++ )
	{
		s << GeneticCodeUtil::residues[i] << " ";
		for ( int j=0; j<20; j++ )
		{
			s.width(5);
			s.fill(' ');
			s << contactEnergies[i][j] << " ";
		}
		s << endl;
	}
}

void CompactLatticeFolder::printStructure( int id, ostream& os, const char* prefix ) const
{
	if ( id < 0 || id >= m_num_structures )
		return;

	m_structures[id]->draw(os, prefix);
}

// Contact energies according to Miyazawa and Jernigan, Residue-Residue Potentials
// with a Favorable Contact Pair Term and an Unfavorable High Packing Density Term,
// for Simulation and Threading. J. Mol. Biol. (1996) 256:623-644
// Table III upper triangle, values e_{ij}.
const double CompactLatticeFolder::contactEnergies[20][20] =
	{
		// CYS
		{-5.44, -4.99, -5.80, -5.50, -5.83, -4.96, -4.95, -4.16, -3.57, -3.16, -3.11, -2.86, -2.59, -2.85, -2.41, -2.27, -3.60, -2.57, -1.95, -3.07},
		// MET
		{-4.99, -5.46, -6.56, -6.02, -6.41, -5.32, -5.55, -4.91, -3.94, -3.39, -3.51, -3.03, -2.95, -3.30, -2.57, -2.89, -3.98, -3.12, -2.48, -3.45},
		// PHE
		{-5.80, -6.56, -7.26, -6.84, -7.28, -6.29, -6.16, -5.66, -4.81, -4.13, -4.28, -4.02, -3.75, -4.10, -3.48, -3.56, -4.77, -3.98, -3.36, -4.25},
		// ILE
		{-5.50, -6.02, -6.84, -6.54, -7.04, -6.05, -5.78, -5.25, -4.58, -3.78, -4.03, -3.52, -3.24, -3.67, -3.17, -3.27, -4.14, -3.63, -3.01, -3.76},
		// LEU
		{-5.83, -6.41, -7.28, -7.04, -7.37, -6.48, -6.14, -5.67, -4.91, -4.16, -4.34, -3.92, -3.74, -4.04, -3.40, -3.59, -4.54, -4.03, -3.37, -4.20},
		// VAL
		{-4.96, -5.32, -6.29, -6.05, -6.48, -5.52, -5.18, -4.62, -4.04, -3.38, -3.46, -3.05, -2.83, -3.07, -2.48, -2.67, -3.58, -3.07, -2.49, -3.32},
		// TRP
		{-4.95, -5.55, -6.16, -5.78, -6.14, -5.18, -5.06, -4.66, -3.82, -3.42, -3.22, -2.99, -3.07, -3.11, -2.84, -2.99, -3.98, -3.41, -2.69, -3.73},
		// TYR
		{-4.16, -4.91, -5.66, -5.25, -5.67, -4.62, -4.66, -4.17, -3.36, -3.01, -3.01, -2.78, -2.76, -2.97, -2.76, -2.79, -3.52, -3.16, -2.60, -3.19},
		// ALA
		{-3.57, -3.94, -4.81, -4.58, -4.91, -4.04, -3.82, -3.36, -2.72, -2.31, -2.32, -2.01, -1.84, -1.89, -1.70, -1.51, -2.41, -1.83, -1.31, -2.03},
		// GLY
		{-3.16, -3.39, -4.13, -3.78, -4.16, -3.38, -3.42, -3.01, -2.31, -2.24, -2.08, -1.82, -1.74, -1.66, -1.59, -1.22, -2.15, -1.72, -1.15, -1.87},
		// THR
		{-3.11, -3.51, -4.28, -4.03, -4.34, -3.46, -3.22, -3.01, -2.32, -2.08, -2.12, -1.96, -1.88, -1.90, -1.80, -1.74, -2.42, -1.90, -1.31, -1.90},
		// SER
		{-2.86, -3.03, -4.02, -3.52, -3.92, -3.05, -2.99, -2.78, -2.01, -1.82, -1.96, -1.67, -1.58, -1.49, -1.63, -1.48, -2.11, -1.62, -1.05, -1.57},
		// GLN
		{-2.59, -2.95, -3.75, -3.24, -3.74, -2.83, -3.07, -2.76, -1.84, -1.74, -1.88, -1.58, -1.68, -1.71, -1.68, -1.51, -2.08, -1.64, -1.21, -1.53},
		// ASN
		{-2.85, -3.30, -4.10, -3.67, -4.04, -3.07, -3.11, -2.97, -1.89, -1.66, -1.90, -1.49, -1.71, -1.54, -1.46, -1.42, -1.98, -1.80, -1.29, -1.73},
		// GLU
		{-2.41, -2.57, -3.48, -3.17, -3.40, -2.48, -2.84, -2.76, -1.70, -1.59, -1.80, -1.63, -1.68, -1.46, -1.21, -1.02, -2.32, -2.29, -1.68, -1.33},
		// ASP
		{-2.27, -2.89, -3.56, -3.27, -3.59, -2.67, -2.99, -2.79, -1.51, -1.22, -1.74, -1.48, -1.51, -1.42, -1.02, -0.91, -2.15, -2.27, -1.80, -1.26},
		// HIS
		{-3.60, -3.98, -4.77, -4.14, -4.54, -3.58, -3.98, -3.52, -2.41, -2.15, -2.42, -2.11, -2.08, -1.98, -2.32, -2.15, -3.05, -2.16, -1.35, -2.25},
		// ARG
		{-2.57, -3.12, -3.98, -3.63, -4.03, -3.07, -3.41, -3.16, -1.83, -1.72, -1.90, -1.62, -1.64, -1.80, -2.29, -2.27, -2.16, -1.55, -0.59, -1.70},
		// LYS
		{-1.95, -2.48, -3.36, -3.01, -3.37, -2.49, -2.69, -2.60, -1.31, -1.15, -1.31, -1.05, -1.21, -1.29, -1.68, -1.80, -1.35, -0.59, -0.12, -0.97},
		// PRO
		{-3.07, -3.45, -4.25, -3.76, -4.20, -3.32, -3.73, -3.19, -2.03, -1.87, -1.90, -1.57, -1.53, -1.73, -1.33, -1.26, -2.25, -1.70, -0.97, -1.75}
	};




/*
// Contact energies according to Miyazawa and Jernigan, Estimation of
// effective interresidue contact energies from protein crystal
// structures: Quasi-chemical approximation. Macromolecules 18:534-552
// (1985).
//Table V, values e_{ij}.
const double CompactLatticeFolder::contactEnergies[20][20] =
	{
		// CYS
		{ -5.44, -5.05, -5.63, -5.03, -5.03, -4.46, -4.76, -3.89, -3.38, -3.16, -2.88, -2.86, -2.73, -2.59, -2.08, -2.66, -3.63, -2.70, -1.54, -2.92 },
		// MET
		{ -5.05, -6.06, -6.68, -6.33, -6.01, -5.52, -6.37, -4.92, -3.99, -3.75, -3.73, -3.55, -3.17, -3.50, -3.19, -2.90, -3.31, -3.49, -3.11, -4.11 },
		// PHE
		{ -5.63, -6.68, -6.85, -6.39, -6.26, -5.75, -6.02, -4.95, -4.36, -3.72, -3.76, -3.56, -3.30, -3.55, -3.51, -3.31, -4.61, -3.54, -2.83, -3.73 },
		// ILE
		{ -5.03, -6.33, -6.39, -6.22, -6.17, -5.58, -5.64, -4.63, -4.41, -3.65, -3.74, -3.43, -3.22, -2.99, -3.23, -2.91, -3.76, -3.33, -2.70, -3.47 },
		// LEU
		{ -5.03, -6.01, -6.26, -6.17, -5.79, -5.38, -5.50, -4.26, -3.96, -3.43, -3.43, -3.16, -3.09, -2.99, -2.91, -2.59, -3.84, -3.15, -2.63, -3.06 },
		// VAL
		{ -4.46, -5.52, -5.75, -5.58, -5.38, -4.94, -5.05, -4.05, -3.62, -3.06, -2.95, -2.79, -2.67, -2.36, -2.56, -2.25, -3.38, -2.78, -1.95, -2.96 },
		// TRP
		{ -4.76, -6.37, -6.02, -5.64, -5.50, -5.05, -5.42, -4.44, -3.93, -3.37, -3.31, -2.95, -3.16, -3.11, -2.94, -2.91, -4.02, -3.56, -2.49, -3.66 },
		// TYR
		{ -3.89, -4.92, -4.95, -4.63, -4.26, -4.05, -4.44, -3.55, -2.85, -2.50, -2.48, -2.30, -2.53, -2.47, -2.42, -2.25, -3.33, -2.75, -2.01, -2.80 },
		// ALA
		{ -3.38, -3.99, -4.36, -4.41, -3.96, -3.62, -3.93, -2.85, -2.51, -2.15, -2.15, -1.89, -1.70, -1.44, -1.51, -1.57, -2.09, -1.50, -1.10, -1.81 },
		// GLY
		{ -3.16, -3.75, -3.72, -3.65, -3.43, -3.06, -3.37, -2.50, -2.15, -2.17, -2.03, -1.70, -1.54, -1.56, -1.22, -1.62, -1.94, -1.68, -0.84, -1.72 },
		// THR
		{ -2.88, -3.73, -3.76, -3.74, -3.43, -2.95, -3.31, -2.48, -2.15, -2.03, -1.72, -1.59, -1.59, -1.51, -1.45, -1.66, -2.31, -1.97, -1.02, -1.66 },
		// SER
		{ -2.86, -3.55, -3.56, -3.43, -3.16, -2.79, -2.95, -2.30, -1.89, -1.70, -1.59, -1.48, -1.37, -1.31, -1.48, -1.46, -1.94, -1.22, -0.83, -1.35 },
		// GLN
		{ -2.73, -3.17, -3.30, -3.22, -3.09, -2.67, -3.16, -2.53, -1.70, -1.54, -1.59, -1.37, -0.89, -1.36, -1.33, -1.26, -1.85, -1.85, -1.02, -1.73 },
		// ASN
		{ -2.59, -3.50, -3.55, -2.99, -2.99, -2.36, -3.11, -2.47, -1.44, -1.56, -1.51, -1.31, -1.36, -1.59, -1.43, -1.33, -2.01, -1.41, -0.91, -1.43 },
		// GLU
		{ -2.08, -3.19, -3.51, -3.23, -2.91, -2.56, -2.94, -2.42, -1.51, -1.22, -1.45, -1.48, -1.33, -1.43, -1.18, -1.23, -2.27, -2.07, -1.60, -1.40 },
		// ASP
		{ -2.66, -2.90, -3.31, -2.91, -2.59, -2.25, -2.91, -2.25, -1.57, -1.62, -1.66, -1.46, -1.26, -1.33, -1.23, -0.96, -2.14, -1.98, -1.32, -1.19 },
		// HIS
		{ -3.63, -3.31, -4.61, -3.76, -3.84, -3.38, -4.02, -3.33, -2.09, -1.94, -2.31, -1.94, -1.85, -2.01, -2.27, -2.14, -2.78, -2.12, -1.09, -2.17 },
		// ARG
		{ -2.70, -3.49, -3.54, -3.33, -3.15, -2.78, -3.56, -2.75, -1.50, -1.68, -1.97, -1.22, -1.85, -1.41, -2.07, -1.98, -2.12, -1.39, -0.06, -1.85 },
		// LYS
		{ -1.54, -3.11, -2.83, -2.70, -2.63, -1.95, -2.49, -2.01, -1.10, -0.84, -1.02, -0.83, -1.02, -0.91, -1.60, -1.32, -1.09, -0.06,  0.13, -0.67 },
		// PRO
		{ -2.92, -4.11, -3.73, -3.47, -3.06, -2.96, -3.66, -2.80, -1.81, -1.72, -1.66, -1.35, -1.73, -1.43, -1.40, -1.19, -2.17, -1.85, -0.67, -1.18 }
	};

// Table VI, values e_{ij}+e_{rr}-e_{ir}-e_{jr}
const double CompactLatticeFolder::contactEnergies[20][20] =
      {
	      // CYS
	      { -1.06,  0.19, -0.23,  0.16, -0.08,  0.06,  0.08,  0.04,  0.00, -0.08,  0.19, -0.02,  0.05,  0.13,  0.69,  0.03, -0.19,  0.24,  0.71,  0.00 },
	      // MET
	      {  0.19,  0.04, -0.42, -0.28, -0.20, -0.14, -0.67, -0.13,  0.25,  0.19,  0.19,  0.14,  0.46,  0.08,  0.44,  0.65,  0.99,  0.31,  0.00, -0.34 },
	      // PHE
	      { -0.23, -0.42, -0.44, -0.19, -0.30, -0.22, -0.16,  0.00,  0.03,  0.38,  0.31,  0.29,  0.49,  0.18,  0.27,  0.39, -0.16,  0.41,  0.44,  0.20 },
	      // ILE
	      {  0.16, -0.28, -0.19, -0.22, -0.41, -0.25,  0.02,  0.11, -0.22,  0.25,  0.14,  0.21,  0.36,  0.53,  0.35,  0.59,  0.49,  0.42,  0.36,  0.25 },
	      // LEU
	      { -0.08, -0.20, -0.30, -0.41, -0.27, -0.29, -0.09,  0.24, -0.01,  0.23,  0.20,  0.25,  0.26,  0.30,  0.43,  0.67,  0.16,  0.35,  0.19,  0.42 },
	      // VAL
	      {  0.06, -0.14, -0.22, -0.25, -0.29, -0.29, -0.07,  0.02, -0.10,  0.16,  0.25,  0.18,  0.24,  0.50,  0.34,  0.58,  0.19,  0.30,  0.44,  0.09 },
	      // TRP
	      {  0.08, -0.67, -0.16,  0.02, -0.09, -0.07, -0.12, -0.04, -0.09,  0.18,  0.22,  0.34,  0.08,  0.06,  0.29,  0.24, -0.12, -0.16,  0.22, -0.28 },
	      // TYR
	      {  0.04, -0.13,  0.00,  0.11,  0.24,  0.02, -0.04, -0.06,  0.09,  0.14,  0.13,  0.09, -0.20, -0.20, -0.10,  0.00, -0.34, -0.25, -0.21, -0.33 },
	      // ALA
	      {  0.00,  0.25,  0.03, -0.22, -0.01, -0.10, -0.09,  0.09, -0.13, -0.07, -0.09, -0.06,  0.08,  0.28,  0.26,  0.12,  0.34,  0.43,  0.14,  0.10 },
	      // GLY
	      { -0.08,  0.19,  0.38,  0.25,  0.23,  0.16,  0.18,  0.14, -0.07, -0.38, -0.26, -0.16, -0.06, -0.14,  0.25, -0.22,  0.20, -0.04,  0.11, -0.11 },
	      // THR
	      {  0.19,  0.19,  0.31,  0.14,  0.20,  0.25,  0.22,  0.13, -0.09, -0.26,  0.03, -0.08, -0.14, -0.11,  0.00, -0.29, -0.19, -0.35, -0.09, -0.07 },
	      // SER
	      { -0.02,  0.14,  0.29,  0.21,  0.25,  0.18,  0.34,  0.09, -0.06, -0.16, -0.08, -0.20, -0.14, -0.14, -0.26, -0.31, -0.05,  0.17, -0.13,  0.01 },
	      // GLN
	      {  0.05,  0.46,  0.49,  0.36,  0.26,  0.24,  0.08, -0.20,  0.08, -0.06, -0.14, -0.14,  0.29, -0.25, -0.17, -0.17, -0.02, -0.52, -0.38, -0.42 },
	      // ASN
	      {  0.13,  0.08,  0.18,  0.53,  0.30,  0.50,  0.06, -0.20,  0.28, -0.14, -0.11, -0.14, -0.25, -0.53, -0.32, -0.30, -0.24, -0.14, -0.33, -0.18 },
	      // GLU
	      {  0.69,  0.44,  0.27,  0.35,  0.43,  0.34,  0.29, -0.10,  0.26,  0.25,  0.00, -0.26, -0.17, -0.32, -0.03, -0.15, -0.45, -0.74, -0.97, -0.10 },
	      // ASP
	      {  0.03,  0.65,  0.39,  0.59,  0.67,  0.58,  0.24,  0.00,  0.12, -0.22, -0.29, -0.31, -0.17, -0.30, -0.15,  0.04, -0.39, -0.72, -0.76,  0.04 },
	      // HIS
	      { -0.19,  0.99, -0.16,  0.49,  0.16,  0.19, -0.12, -0.34,  0.34,  0.20, -0.19, -0.05, -0.02, -0.24, -0.45, -0.39, -0.29, -0.12,  0.22, -0.21 },
	      // ARG
	      {  0.24,  0.31,  0.41,  0.42,  0.35,  0.30, -0.16, -0.25,  0.43, -0.04, -0.35,  0.17, -0.52, -0.14, -0.74, -0.72, -0.12,  0.11,  0.75, -0.38 },
	      // LYS
	      {  0.71,  0.00,  0.44,  0.36,  0.19,  0.44,  0.22, -0.21,  0.14,  0.11, -0.09, -0.13, -0.38, -0.33, -0.97, -0.76,  0.22,  0.75,  0.25,  0.11 },
	      // PRO
	      {  0.00, -0.34,  0.20,  0.25,  0.42,  0.09, -0.28, -0.33,  0.10, -0.11, -0.07,  0.01, -0.42, -0.18, -0.10,  0.04, -0.21, -0.38,  0.11,  0.26 }
      };
*/



