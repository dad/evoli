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
	Fold sequence and return folding information in a \ref FoldInfo
	object. Note that the function fold() returns a pointer to a \ref FoldInfo
	object, not an instance of this object. The ownership of this object is
	transferred to the function issuing the call. It is strongly recommended
	to store this pointer in an auto-pointer, as in this example:
	\code
// assume the variable f holds a Folder object:
auto_ptr<FoldInfo> fi( f.fold( s ) ); // fold sequence s
StructureID sid = fi->getStructure(); // assign structure ID to variable sid
	\endcode
	\warning Do not write code of the form f.fold()->getStructure(),
	as this statement would result in a memory leak (\ref FoldInfo object
	is not deleted).
	@param s The protein sequence to fold.
	@return A pointer to a \ref FoldInfo object containing all
	the folding information.
	*/
	virtual FoldInfo* fold(const Sequence& s) const = 0;
	/**
	 * @param s The sequence whose energy is sought.
	 * @param sid The structure ID of the target conformation.
	 * @return The contact energy of a sequence in the target conformation.
	 **/ 
	virtual double getEnergy(const Sequence& s, StructureID sid) const = 0;
	/**
	This function assesses whether the folder has been properly initialized.
	@return True if the folder is in good working order, False otherwise.
	 **/
	virtual bool good() const = 0;
};

#endif // FOLDER_HH






