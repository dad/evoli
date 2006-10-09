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


#ifndef FOLDER_HH
#define FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>
#include "sequence.hh"

using namespace std;
typedef int StructureID;

typedef pair<int,int> Contact;

/**
 * Stores a contact structure: a list of contacts.
 **/
class ContactStructure {
private:
protected:
public:
	virtual ~ContactStructure() {}
	virtual const vector<Contact>& getContacts() const = 0;
};

class FoldInfo {
protected:
	double m_free_energy;
	StructureID m_structure_id;
public:
	FoldInfo() {
		m_free_energy = 0.0;
		m_structure_id = (StructureID)-1;
	}
		
	FoldInfo(double fe, StructureID sid) {
		m_free_energy = fe;
		m_structure_id = sid;
	}
	virtual ~FoldInfo() {}
	
	StructureID getStructure() const { return m_structure_id; }
	double getFreeEnergy() const { return m_free_energy; }
};

class Folder {
private:
protected:
public:
	virtual ~Folder() {}
	/**
	 * Fold sequence and return folding information.
	 * @param s The sequence whose energy is sought.
	 * @return Folding information (e.g. free energy, structure identifier).
	 **/
	virtual FoldInfo* fold(const Sequence& s) const = 0;
	/**
	 * @param s The sequence whose energy is sought.
	 * @param sid The structure ID of the target conformation.
	 * @return The contact energy of a sequence in the target conformation.
	 **/ 
	virtual double getEnergy(const Sequence& s, StructureID sid) const = 0;
	/**
	 * Has the folder been properly initialized?
	 **/
	virtual bool good() const = 0;
};

#endif // FOLDER_HH






