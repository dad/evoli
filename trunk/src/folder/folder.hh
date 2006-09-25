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
	virtual FoldInfo fold(const Sequence& s) = 0;
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






