#ifndef DECOY_CONTACT_FOLDER_HH
#define DECOY_CONTACT_FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

#include "folder.hh"

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
	int getMaxResidueNumber();
};

class DecoyContactFolder : public Folder {
private:
	DecoyContactFolder();
protected:
	int m_length;
	double m_log_num_conformations;
	vector<DecoyContactStructure*> m_structures;
	static const double DecoyContactFolder::contactEnergies [20][20];
	int m_num_folded;

public:
	/**
	 * Create DCF using structures provided in structs.
	 * log_num_confs is a numerical fudge-factor; use log(10^160).
	 **/
	DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs);
	/**
	 * Create DCF; obtain structures from files listed in fin and stored in directory dir.
	 * log_num_confs is a numerical fudge-factor; use log(10^160).
	 **/
	DecoyContactFolder(int length, double log_num_confs, ifstream& fin, const string& dir);

	~DecoyContactFolder();

	/**
	 * Fold the protein and return folding information (structure, free energy).
	 **/
	virtual FoldInfo fold(const Sequence& s);
	/**
	 * Returns the number of proteins that have been folded since initialization.
	 **/
	int getNumFolded() { return m_num_folded; }

	/**
	 * Has folder been properly initialized?
	 **/
	virtual bool good();
	
};

struct ContactMapUtil {
	static void readContactMapsFromFile(ifstream& fin, const string& dir, vector<DecoyContactStructure*>& structs);
};

#endif // DECOY_CONTACT_FOLDER_HH






