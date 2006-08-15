#ifndef DECOY_CONTACT_FOLDER_HH
#define DECOY_CONTACT_FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>

#include "protein-folder.hh"


using namespace std;

/**
 * Stores a contact structure: a list of contacts.
 **/
class DecoyContactStructure : public ContactStructure {
protected:
	vector<Contact> m_contacts;
public:
	DecoyContactStructure() {}

	/**
	 * Read from a file.
	 **/
	void read(ifstream& fin);
	virtual const vector<Contact>& getContacts() const { return m_contacts; }
};

class DecoyContactFolder : public ProteinFolder {
private:
	DecoyContactFolder();
protected:
	int m_length;
	double m_log_num_conformations;
	vector<DecoyContactStructure*> m_structures;
	static const double DecoyContactFolder::contactEnergies [20][20];
	int m_num_folded;

public:
	DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs);
	/**
	 * foldProtein(): core interface
	 **/
	virtual FoldInfo foldProtein(Protein& p);

	
};

#endif // DECOY_CONTACT_FOLDER_HH






