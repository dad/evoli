#ifndef PROTEIN_FOLDER_HH
#define PROTEIN_FOLDER_HH

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
	virtual const vector<Contact>& getContacts() const = 0;
};

class FoldInfo {
private:
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
	StructureID getStructure() const { return m_structure_id; }
	double getFreeEnergy() const { return m_free_energy; }
};

class ProteinFolder {
private:
protected:
public:
	/**
	 * foldProtein(): core interface
	 * Pure virtual to force overriding.
	 **/
	virtual FoldInfo fold(const Sequence& s) = 0;
};

#endif // PROTEIN_FOLDER_HH






