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


#include "decoy-contact-folder.hh"
#include <fstream>
#include <cmath>
#include "genetic-code.hh"

double DecoyContactFolder::BAD_ENERGY = 999999.0;

void DecoyContactStructure::read(istream& fin) {
	int r1, r2;
	string r1aa, r2aa;
	do {
		fin >> r1 >> r1aa >> r2 >> r2aa;
		//cout << r1 << r1aa << r2 << r2aa << endl;
		if (!fin.eof() && r1 >= 0 && r2 >= 0) {
			// COW debugging: what happens if we leave out trivial contacts
			if ( abs(r1-r2) != 1) {
				m_contacts.push_back(Contact(r1,r2));
			}
		}
	}while (!fin.eof());
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

DecoyContactFolder::DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs, double deltaGCutoff, StructureID targetSID)
	: DGCutoffFolder( deltaGCutoff, targetSID )
{
	m_length = length;
	m_structures = structs;
	m_log_num_conformations = log_num_confs;
	m_num_folded = 0;
}

DecoyContactFolder::DecoyContactFolder(int length, double log_num_confs, ifstream& fin, const string& dir, double deltaGCutoff, StructureID targetSID )
	: DGCutoffFolder( deltaGCutoff, targetSID ), m_structures(0)
{
	m_length = length;
	m_log_num_conformations = log_num_confs;
	ContactMapUtil::readContactMapsFromFile(fin, dir, m_structures);
	//cout << "# structures = " << m_structures.size() << endl;
	//cout << "# lognumconfs = " << log_num_confs << endl;
	m_num_folded = 0;
}

DecoyContactFolder::~DecoyContactFolder() {
	vector<DecoyContactStructure*>::iterator it = m_structures.begin();
	for ( ; it != m_structures.end(); it++) {
		delete *it;
	}
}

bool DecoyContactFolder::good() const {
	return m_structures.size() > 0;
}



double DecoyContactFolder::getEnergy(const Protein& s, StructureID sid) const {
	double G = 0;
	if (sid >= m_structures.size()) {
		return DecoyContactFolder::BAD_ENERGY;
	}
	vector<unsigned int> aa_indices(s.size());
	getAminoAcidIndices(s, aa_indices);
	const vector<Contact> &pair_list = m_structures[sid]->getContacts();
	vector<Contact>::const_iterator it=pair_list.begin();
	for ( ; it!=pair_list.end(); it++ )	{
		int s1 = (*it).first;
		int s2 = (*it).second;
		if (s1 < m_length && s2 < m_length) {
			double contact_G = contactEnergy( aa_indices[s1], aa_indices[s2] );
			G += contact_G;
			//cout << "(" << s1 << ", " << s2 << ") -> " << GeneticCodeUtil::residues[s[s1]] 
			//	 << ":" << GeneticCodeUtil::residues[s[s2]] << " " << contact_G << " " << G << endl << flush;
		}
	}
	return G;
}

/**
 * Fold the protein and return folding information (structure, free energy).
 **/
DecoyFoldInfo* DecoyContactFolder::fold(const Protein& s) const {
	double kT = 0.6;
	double minG = 1e50;
	int minIndex = -1;

	double sumG = 0.0;
	double sumsqG = 0.0;

	vector<unsigned int> aa_indices(s.size());
	bool valid = getAminoAcidIndices(s, aa_indices);
	if (!valid) {
		return new DecoyFoldInfo(false, false, 9999, -1, 9999, 9999, 9999);
	}

	for ( unsigned int sid = 0; sid < m_structures.size(); sid++) {
		double G = 0;
		// calculate binding energy of this fold
		const vector<Contact> &pair_list = m_structures[sid]->getContacts();
		vector<Contact>::const_iterator it=pair_list.begin();
		int num_contacts = 0;
		for ( ; it!=pair_list.end(); it++ )	{
			int s1 = (*it).first;
			int s2 = (*it).second;
			if (s1 < m_length && s2 < m_length) {
				double contact_G = contactEnergy( aa_indices[s1], aa_indices[s2] );
				// DAD: debugging
				/*if (sid == 24) {
					cout << num_contacts << "\t" << s1 << "\t" << s2 << "\t" << GeneticCodeUtil::residueLetters[s[s1]+1] << "\t" << GeneticCodeUtil::residueLetters[s[s2]+1] << "\t" << contact_G << endl;
					}*/
				G += contact_G;
				//cout << "(" << s1 << ", " << s2 << ") -> " << GeneticCodeUtil::residues[s[s1]] 
				//	 << ":" << GeneticCodeUtil::residues[s[s2]] << " " << contact_G << " " << G << endl << flush;
				num_contacts++;
			}
		}
		// check if binding energy is lower than any previously calculated one
		if ( G < minG )
		{
			minG = G;
			minIndex = sid;
		}
		// add energy to partition sum
        // DAD: debugging
		//cout << G << " energy for str. " << sid << " (" << num_contacts << " contacts of " << pair_list.size() << ")" << endl;
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

	/*cout << "minG:" << minG << endl;
	cout << "mean_G:" << mean_G << endl;
	cout << "var_G:" << var_G << endl;
	cout << "dG:" << dG << endl;
	cout << "(var_G - 2*kT*mean_G)/(2.0*kT): " << ((var_G - 2*kT*mean_G)/(2.0*kT)) << endl;
	cout << "kT ln N: " << kT * m_log_num_conformations << endl;
	*/
	
	// increment folded count
	m_num_folded += 1;
	return new DecoyFoldInfo(dG<m_deltaG_cutoff, minIndex==m_target_sid, dG, minIndex, mean_G, var_G, minG);
}


void ContactMapUtil::readContactMapsFromFile(ifstream& fin, const string& dir, vector<DecoyContactStructure*>& structs) {
	string filename;
	if ( !fin.good() ){
		cout << "# Warning: cannot read from stream in ContactMapUtil::readContactMapsFromFile" << endl;
		return;
	}
	do {
		fin >> filename;
		if (!fin.eof()) {
			string path = dir + filename;
			ifstream cfile(path.c_str());
			//cout << path << endl;
			if (cfile.good()) {
				DecoyContactStructure* cstruct = new DecoyContactStructure();
				cstruct->read(cfile);
				cfile.close();
				structs.push_back(cstruct);
			}
			else {
				cout << "# Warning: bad file in ContactMapUtil::readContactMapsFromFile: " << path << endl;
			}
		}
	} while (!fin.eof());
}
