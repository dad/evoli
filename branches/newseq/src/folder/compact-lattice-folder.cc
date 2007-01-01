/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1
*/

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
		cout << "Maximally suported size: 15 (because char is being used to enumerate sites)" << endl;
		exit(-1);
	}

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
FoldInfo* CompactLatticeFolder::fold( const Sequence& s ) const
{
	assert( m_num_structures > 0 );

	double kT = 0.6;
	double minE = 1e50;
	int minIndex = 0;
	double Z = 0;
	double G;

	for ( int i=0; i<m_num_structures; i++ ) {
		double E = 0;

		// calculate binding energy of this fold
		const vector<Contact> &pair_list = m_structures[i]->getInteractingPairs();
		vector<pair<int,int> >::const_iterator it=pair_list.begin();
		for ( ; it!=pair_list.end(); it++ )
		{
			double deltaE = contactEnergy(s[(*it).first-1], s[(*it).second-1]);
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

	// increment folded count
	m_num_folded += 1;

	return new FoldInfo(G, minIndex);
}

double CompactLatticeFolder::getEnergy(const Sequence& p, StructureID sid) const {
	const vector<Contact> &pair_list = m_structures[sid]->getInteractingPairs();
	vector<Contact>::const_iterator it=pair_list.begin();
	double E = 0.0;
	for ( ; it!=pair_list.end(); it++ ){
		double deltaE = contactEnergy(p[(*it).first-1], p[(*it).second-1]);
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
			double e = contactEnergy(p[i], p[j]);
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
			s << contactEnergy(i, j) << " ";
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
