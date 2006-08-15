#ifndef DECOY_CONTACT_FOLDER_HH
#define DECOY_CONTACT_FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>

#include "protein-folder.hh"


using namespace std;

typedef pair<int,int> Contact;

/**
 * Stores a contact structure: a list of contacts.
 **/
class DecoyContactStructure {
protected:
	vector<Contact> m_contacts;
public:
	DecoyContactStructure() {}

	/**
	 * Read from a file.
	 **/
	void read(ifstream& fin);
};

class DecoyContactFolder : public ProteinFolder {
private:
	DecoyContactFolder();
protected:
	int m_length;
	vector<DecoyContactStructure> m_structures;
public:
	DecoyContactFolder(int length, vector<DecoyContactStructure>& structs);
	/**
	 * foldProtein(): core interface
	 **/
	virtual FoldInfo foldProtein(Protein& p);

	
};

#endif // DECOY_CONTACT_FOLDER_HH






