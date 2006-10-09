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


#ifndef COMPACT_LATTICE_FOLDER_HH
#define COMPACT_LATTICE_FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>

#include "folder.hh"
#include "protein-contact-energies.hh"

// uncomment next line for older versions of gcc
//#include <hash_map>
// comment out next two lines for older versions of gcc
#include <ext/hash_map>
using namespace __gnu_cxx;

using namespace std;

class LatticeStructure : public ContactStructure
{
private:
	int m_size;
	char *m_structure;
	vector<Contact> m_interacting_pairs;

	void calcInteractingPairs();
	const LatticeStructure & operator=( const LatticeStructure & );
public:
	LatticeStructure();
	LatticeStructure( const LatticeStructure & );
	LatticeStructure( const char *structure, int size );
	~LatticeStructure();

	char * getStructure() const	{
		return m_structure;
	}
	const vector<Contact>& getInteractingPairs() const {
		return m_interacting_pairs;
	}
	const vector<Contact>& getContacts() const {
		return m_interacting_pairs;
	}

	vector<int> getSurface() const;
	void draw(ostream& os, const char* prefix) const;
};

class StructureUtil
{
public:
	/**
	* Draws the structure specified in the string s to the specified ostream.
	**/
	static void drawStructure( ostream& os, const char *s, int size, const char *prefix );
	static void flipLeftRight( const char* s, char* d, int size );
	static void rotate90( const char* s, char* d, int size );

	static void drawSite( ostream &s, int site );
	static void drawSite( ostream &s, char site );
	static void drawHBond( ostream &s, int bond );
	static void drawHBond( ostream &s, char bond );
	static void drawVBond( ostream &s, int bond );
	static void drawVBond( ostream &s, char bond );
};


class SelfAvoidingWalk
{
public:
	enum Direction {right, left, forward};
private:
	vector<vector<int> > m_vbonds;
	vector<vector<int> > m_hbonds;
	vector<vector<int> > m_sites;
	vector<Direction> m_walk;
	const int m_size;

	int m_start_x;
	int m_start_y;
	int m_cur_x;
	int m_cur_y;
	int m_dx;
	int m_dy;
	int m_length;

	char *m_tmp_structure;

	SelfAvoidingWalk();
	SelfAvoidingWalk( const SelfAvoidingWalk & );
	const SelfAvoidingWalk & operator=( const SelfAvoidingWalk & );
public:
	SelfAvoidingWalk( int n );
	~SelfAvoidingWalk();

	// modifiers
	void clear();
	bool doMove( Direction d );
	void eraseLastMove();
	bool setString( int x, int y, const char* string );
	void setStart( int x, int y );

	// accessors
	void draw(ostream& os) const;

	/**
	 * Saves the current structure in the provided string d. d must
	 * be a sufficiently long array.
	 **/
	void getStructure( char *d ) const;
	int startX() const {
		return m_start_x;
	}
	int startY() const {
		return m_start_y;
	}
	int length() const {
		return m_length;
	}
	int maxLength() const {
		return m_size*m_size;
	}
};


/**
 * Needed for the structure map in StructureBank.
 **/
struct eqstr {
	bool operator()(const char* s1, const char* s2) const {
		return (strcmp(s1, s2) == 0);
	}
};

/**
 * Needed for the structure map in StructureBank.
 **/
struct ltstr {
	bool operator()(const char* s1, const char* s2) const {
		return (strcmp(s1, s2) < 0);
	}
};


class CompactLatticeFolder : public Folder {
private:
	// some useful typedefs
	typedef hash_map<const char*, int, hash<const char*>, eqstr> StructureMap;
	typedef hash_map<const char*, int, hash<const char*>, eqstr>::iterator StructureMapIterator;
	typedef hash_map<const char*, int, hash<const char*>, eqstr>::const_iterator StructureMapConstIterator;

	// the contact energies between residues
//	static const double contactEnergies[20][20];
	// the size of the square lattice (protein length is size^2)
	const int m_size;
	// the number of proteins folded
	mutable int m_num_folded;

	int m_num_structures; // total number of structures
	vector<LatticeStructure *> m_structures; // list of all potential structures
	StructureMap m_structure_map; // lookup table for structures

	char * m_ffw_struct;  // variable used by findFillingWalk();
	char * m_ss_struct; // variable used by storeStructure();
	char * m_ss_struct2; // variable used by storeStructure();

	CompactLatticeFolder();
	CompactLatticeFolder( const CompactLatticeFolder & );
	const CompactLatticeFolder & operator=( const CompactLatticeFolder & );
protected:
	void findFillingWalks( SelfAvoidingWalk &w, int &moves );
	bool findStructure( const char*s );
	void storeStructure( const char* s );
	void enumerateStructures();
	/**
	* Wrapper function to encapsulate the lookup of the
	* contact energy from a table.
	*
	* @param residue1 The first of the two contacting residues.
	* Mapping from int to residue type is as in class @ref ProteinContactEnergies.
	* @param residue2 The second of the two contacting residues.
	* @return The contact energy between the two residues.
	**/
	double contactEnergy( int residue1, int residue2 ) const {
		return ProteinContactEnergies::MJ96TableIII[residue1][residue2]; }
	//return ProteinContactEnergies::MJ85TableVI[residue1][residue2]; }

public:
	CompactLatticeFolder( int size );
	virtual ~CompactLatticeFolder();

	virtual bool good() const { return m_structures.size() > 0; }
	virtual FoldInfo* fold( const Sequence& s ) const;
	bool isFoldedBelowThreshold( const Sequence&s, const int structID, double cutoff) const;
	void getMinMaxPartitionContributions(const Sequence& s, const int ci, double& cmin, double& cmax) const;
	/**
	 * @param s The sequence whose energy is sought.
	 * @param sid The structure ID of the target conformation.
	 * @return The contact energy of a sequence in the target conformation.
	 **/ 
	virtual double getEnergy(const Sequence& s, StructureID sid) const;

	void printContactEnergyTable( ostream &s ) const;
	void printStructure( int id, ostream& os, const char* prefix ) const;
	vector<int> getSurface( int id ) const
	{
		return m_structures[id]->getSurface();
	}

	LatticeStructure* getStructure(const int id) const {
		return m_structures[id];
	}
	unsigned int getNumFolded() const {
		return m_num_folded;
	}
};


#endif //COMPACT_LATTICE_FOLDER_HH






