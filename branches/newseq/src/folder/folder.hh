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
#include "genetic-code.hh" // this is possibly a bad dependence

using namespace std;
typedef int StructureID;

typedef pair<int,int> Contact;

/**
\brief This class stores a contact structure: a list of contacts.
*/
class ContactStructure {
private:
protected:
public:
	virtual ~ContactStructure() {}
	virtual const vector<Contact>& getContacts() const = 0;
};


/**
\brief A \ref FoldInfo object contains data generated during the folding of a protein.

At a minimum, a \ref FoldInfo object contains the \ref StructureID of the minimum free energy structure of the protein sequence, as well as a bool indicating whether the sequence folds stably into the minimum free energy structure or not.
*/
class FoldInfo {
protected:
	double m_free_energy;
	StructureID m_structure_id; ///< \ref StructureID of the minimum free energy structure of the folded protein sequence.
public:
	FoldInfo() : m_free_energy( 0.0 ),
		m_structure_id( static_cast<StructureID>( -1 ) ) { ; }
		
	FoldInfo(double fe, StructureID sid) : m_free_energy( fe ),
		m_structure_id( sid ) { ; }

	virtual ~FoldInfo() {}

	/**
	\return The \ref StructureID of the minimum free energy structure of the folded protein sequence.
	*/
	StructureID getStructure() const { return m_structure_id; }
	double getDeltaG() const { return m_free_energy; }
};

/**
\brief Abstract basis class for a class that can fold a protein sequence into a (three-dimensional) structure.
*/
class Folder {
private:
protected:
	/**
	 * Turn amino acid characters into indices to enable fast lookup during folding.
	 *
	 * @return Whether all indices were valid (>= 0).
	 **/
	bool getAminoAcidIndices(const Sequence& s, vector<int>& aa_indices) const
	{
		int index = 0;
		for (int i=0; i<s.size() && index >= 0; i++) {
			index = GeneticCodeUtil::aminoAcidLetterToIndex(s[i]);
			aa_indices[i] = index;
		}
		return index >= 0;
	}


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
	\brief This function assesses whether the folder has been properly initialized.
	@return True if the folder is in good working order, False otherwise.
	 **/
	virtual bool good() const = 0;

	/**
	@return The number of proteins that have been folded so far with this Folder instance.
	*/
	virtual uint getNumFolded() const = 0;

};

#endif // FOLDER_HH






