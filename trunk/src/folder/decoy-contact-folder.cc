#include "decoy-contact-folder.hh"
#include <fstream>
#include <cmath>
#include "genetic-code.hh"

void DecoyContactStructure::read(ifstream& fin) {
	int r1, r2;
	string r1aa, r2aa;
	while (!fin.eof()) {
		fin >> r1 >> r1aa >> r2 >> r2aa;
		if (r1 >= 0 && r2 >= 0) {
// COW debugging: what happens if we leave out trivial contacts
//			if ( r1-r2 > 0 || r1-r2 < -0 )
				m_contacts.push_back(Contact(r1,r2));
		}
	}
}

int DecoyContactStructure::getMaxResidueNumber() {
	int max_res = -1;
	vector<Contact>::const_iterator it=m_contacts.begin();
	for ( ; it!=m_contacts.end(); it++ )	{
		if ((*it).first > max_res)
			max_res = (*it).first;
		if ((*it).second > max_res)
			max_res = (*it).second;
	}
	return max_res;
}

DecoyContactFolder::DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs) {
	m_length = length;
	m_structures = structs;
	m_log_num_conformations = log_num_confs;
	m_num_folded = 0;
}

DecoyContactFolder::DecoyContactFolder(int length, double log_num_confs, ifstream& fin, const string& dir): m_structures(0) {
	m_length = length;
	ContactMapUtil::readContactMapsFromFile(fin, dir, m_structures);
	m_log_num_conformations = log_num_confs;
	m_num_folded = 0;
}

DecoyContactFolder::~DecoyContactFolder() {
	vector<DecoyContactStructure*>::iterator it = m_structures.begin();
	for ( ; it != m_structures.end(); it++) {
		delete *it;
	}
}

bool DecoyContactFolder::good() {
	return m_structures.size() > 0;
}

/**
 * Fold the protein and return folding information (structure, free energy).
 **/
FoldInfo DecoyContactFolder::fold(const Sequence& s) {
	double kT = 0.6;
	double minG = 1e50;
	int minIndex = 0;

	double sumG = 0.0;
	double sumsqG = 0.0;

	for ( unsigned int sid = 0; sid < m_structures.size(); sid++) {
		double G = 0;

		// calculate binding energy of this fold
		const vector<Contact> &pair_list = m_structures[sid]->getContacts();
		vector<Contact>::const_iterator it=pair_list.begin();
		for ( ; it!=pair_list.end(); it++ )	{
			int s1 = (*it).first;
			int s2 = (*it).second;
			if (s1 < m_length && s2 < m_length) {
				double contact_G = contactEnergy( s[s1], s[s2] );
				G += contact_G;
				//cout << "(" << s1 << ", " << s2 << ") -> " << GeneticCodeUtil::residues[s[s1]] 
				//	 << ":" << GeneticCodeUtil::residues[s[s2]] << " " << contact_G << " " << G << endl << flush;
			}
		}
		// check if binding energy is lower than any previously calculated one
		if ( G < minG )
		{
			minG = G;
			minIndex = sid;
		}
		// add energy to partition sum
		sumG += G;
		sumsqG += G*G;
	}

	// remove min. energy
	sumG -= minG;
	sumsqG -= minG*minG;

	unsigned int num_confs = m_structures.size() -1;
	double mean_G = sumG/num_confs;
	double var_G = (sumsqG - (sumG*sumG)/num_confs)/(num_confs-1.0);
	// calculate free energy of folding
	double dG = minG + (var_G - 2*kT*mean_G)/(2*kT) + kT * m_log_num_conformations;

//	cout << "minG:" << minG << endl;
//	cout << "mean_G:" << mean_G << endl;
//	cout << "var_G:" << var_G << endl;
//	cout << "dG:" << dG << endl;
//	cout << "(var_G - 2*kT*mean_G)/(2.0*kT): " << ((var_G - 2*kT*mean_G)/(2.0*kT)) << endl;
//	cout << "kT ln N: " << kT * m_log_num_conformations << endl;

	// increment folded count
	m_num_folded += 1;

	return FoldInfo(dG, minIndex);
}

void ContactMapUtil::readContactMapsFromFile(ifstream& fin, const string& dir, vector<DecoyContactStructure*>& structs) {
	string filename;
	char buf[100];
	while (!fin.eof()) {
		fin.getline(buf, 100);
		string path = dir + string(buf);
		ifstream cfile(path.c_str());
		if (!fin.good())
			continue;
		DecoyContactStructure* cstruct = new DecoyContactStructure();
		cstruct->read(cfile);
		cfile.close();
		structs.push_back(cstruct);
	}
}
