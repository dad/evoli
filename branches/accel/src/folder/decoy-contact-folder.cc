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
#include "random.hh"

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
	} while (!fin.eof());
	init();
}

void DecoyContactStructure::init() {
	// create 
	for (const_iterator it=m_contacts.begin(); it != m_contacts.end(); it++) {
		const Contact& c = *it;
		m_contacts_for_site[c.first].push_back(c);
		m_contacts_for_site[c.second].push_back(c);
	}

	/* DAD debugging
	   uint max_res = getMaxResidueNumber();
	   for (uint i=0; i<max_res; i++) {
	   const ContactSet& cs = m_contacts_for_site[i];
	   if (cs.size() > 0) {
	   cout << "init site " << i << ": ";
	   for (const_iterator it=cs.begin(); it != cs.end(); it++) {
	   const Contact& c = *it;
	   cout << c << ",";
	   }
	   cout << endl;
	   }
	   }
	*/
}

const DecoyContactStructure::ContactSet& DecoyContactStructure::contactsForSite(uint site) const {
	ContactSetMap::const_iterator it = m_contacts_for_site.find(site);
	if (it != m_contacts_for_site.end()) {
		return it->second;
	}
	return nullSet();
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
	initializeStructuresForResidues();
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
	initializeStructuresForResidues();
}


bool inList(const StructureID& sid, const vector<StructureID>& list) {

	//void DecoyContactStructure::inList(StructureID sid, const vector<Contact> m_structures_for_residue[uint residue_number]) {
	//(s1 < m_length) {
	//cout << "inList: " << sid << " " << list.size() << endl;
	for (uint i=0; i<list.size(); i++) {
		if (list[i] == sid) {
			return true;
		}
	}
	return false;
}

/*******************************Beginning of Structures For Residue Implementation*************************************************/

//vector<vector<StructureID> >& lookup_table = folder.getStructureLookupTable();
// compare structures against contact lookups

// Return true if the given StructureID is in the list, false otherwise.
bool structureInList(const vector<StructureID>& list, const StructureID& sid) {
	bool in_list = false;
	for (uint i=0; i<list.size() && !in_list; i++) {
		in_list = (list[i] == sid);
	}
	return in_list;   
}

void DecoyContactFolder::initializeStructuresForResidues() {
	// Set up cache of structures.
	// Allocate list of structure lists
	m_structures_for_residue.resize(m_length);
	// For each structure ID,
	// get the structure
	// run through the contacts
	// for both the first and second residue in the contact,
	// add the structure ID to the corresponding location.
  
	for (unsigned int sid = 0; sid < m_structures.size(); sid++) {
		const vector<Contact> &pair_list = m_structures[sid]->getContacts();
		vector<Contact>::const_iterator it = pair_list.begin();
		int num_contacts = 0;
		for ( ; it!=pair_list.end(); it++ )	{
			int s1 = (*it).first;
			int s2 = (*it).second;
			if (s1 < m_length && !inList(sid, m_structures_for_residue[s1])){
				m_structures_for_residue[s1].push_back(sid);
			}	   
	  
			if (s2 < m_length && !inList(sid, m_structures_for_residue[s2])){
				m_structures_for_residue[s2].push_back(sid);
			}	   
		}
	}
}

/*******************************End of Structures For Residue Implementation*************************************************/



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
	vector<unsigned int> aa_indices(s.size());
	if (sid >= m_structures.size()) {
		return 999999;
	}
	bool valid = getAminoAcidIndices(s, aa_indices);
	const vector<Contact> &pair_list = m_structures[sid]->getContacts();
	vector<Contact>::const_iterator it=pair_list.begin();
	for ( ; it!=pair_list.end(); it++ )	{
		int s1 = (*it).first;
		int s2 = (*it).second;
		if (s1 < m_length && s2 < m_length) {
			double contact_G = contactEnergy( aa_indices[s1], aa_indices[s2] );
			G += contact_G;
		}
	}
	return G;
}
//getStructuresWithResidueContact(uint residue_number)
double DecoyContactFolder::getEnergy(const vector<uint>& aa_indices, StructureID sid) const {
	double G = 0.0;
	if (sid >= m_structures.size()) {
		return 999999;
	}
	const vector<Contact> &pair_list = m_structures[sid]->getContacts();
	vector<Contact>::const_iterator it = pair_list.begin();
	for ( ; it!=pair_list.end(); it++ )	{
		int s1 = (*it).first;
		int s2 = (*it).second;
		if (s1 < m_length && s2 < m_length) {
			double contact_G = contactEnergy( aa_indices[s1], aa_indices[s2] );
			G += contact_G;
		}
	}
	return G;
}

// Problem with this function is that, due to rounding, values drift from reality rather quickly.
double DecoyContactFolder::updateEnergy(const DecoyContactStructure* structure, double old_energy, const Protein& old_protein, const Protein& new_protein, const vector<uint>& sites_changed) const {
	//triplet of (index, from_aa, to_aa)) {
	// take old energy
	double new_energy = old_energy;
	for (vector<uint>::const_iterator sit= sites_changed.begin(); sit != sites_changed.end(); sit++ ) {
		const DecoyContactStructure::ContactSet& set = structure->contactsForSite( *sit );
		for (DecoyContactStructure::ContactSet::const_iterator it=set.begin(); it != set.end(); it++) {
			const Contact& contact = *it;
			// subtract the energy of old contact
			uint from_index1 = GeneticCodeUtil::aminoAcidLetterToIndex(old_protein[contact.first]);
			uint from_index2 = GeneticCodeUtil::aminoAcidLetterToIndex(old_protein[contact.second]);
			new_energy -= contactEnergy( from_index1, from_index2 );
			// get energy of new contact, add to old_energy
			uint to_index1 = GeneticCodeUtil::aminoAcidLetterToIndex(new_protein[contact.first]);
			uint to_index2 = GeneticCodeUtil::aminoAcidLetterToIndex(new_protein[contact.second]);
			new_energy += contactEnergy( to_index1, to_index2 );
		}
	}
	return new_energy;
}

/**
 * Fold the protein and return folding information (structure, free energy).
 **/
DecoyFoldInfo* DecoyContactFolder::fold(const Protein& s) const {
	double kT = 0.6;
	double minG = 1e50;
	int minIndex = -1;
	vector<unsigned int> aa_indices(s.size());
	double sumG = 0.0;
	double sumsqG = 0.0;
  
	bool valid = getAminoAcidIndices(s, aa_indices);
	if (!valid) {
		return new DecoyFoldInfo(false, false, 9999, -1, 9999, 9999, 9999);
	}
	for ( unsigned int sid = 0; sid < m_structures.size(); sid++) {
		double G = 0;
		// calculate binding energy of this fold
    
		/* change into a function*/	
    
		const vector<Contact> &pair_list = m_structures[sid]->getContacts();
		vector<Contact>::const_iterator it = pair_list.begin();
		const vector<Contact>::const_iterator& end = pair_list.end();
    
		int num_contacts = 0;
    
		for ( ; it!=end; it++ ) {
			int s1 = (*it).first;
			int s2 = (*it).second;
			if (s1 < m_length && s2 < m_length) {
				double contact_G = contactEnergy( aa_indices[s1], aa_indices[s2] );
				//	cout << "Fold G for contacts: " << s1 << "and " << s2 
				//   << "for SID: " << sid << " is: " << G << endl;
				G += contact_G;
				num_contacts++;
			}
		}
    

		/*end of function*/

		//cout << "fold: " << G << " sid = " << sid << endl;

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
	
	// increment folded count
	m_num_folded += 1;
	return new DecoyFoldInfo(dG<m_deltaG_cutoff, minIndex==m_target_sid, dG, minIndex, mean_G, var_G, minG);
}

DecoyHistoryFoldInfo* DecoyContactFolder::foldWithHistory(const Protein & p, const DecoyHistoryFoldInfo* history) const {
	
	double kT = 0.6;
	double minG = 1e50;
	int minIndex = -1;
	vector<unsigned int> aa_indices(p.size());
	double sumG = 0.0;
	double sumsqG = 0.0;

	// history = NULL;
	// What if there's no history?
	// Make a new DHFI
	if (history == NULL || Random::runif()<0.5) {
		bool valid = getAminoAcidIndices(p, aa_indices);
		vector<double> energies(m_structures.size(), 0.0);
		for (unsigned int sid = 0; sid < m_structures.size(); sid++) {
			double G = getEnergy(p, sid);
			//cout << "null: sid=" << sid << " energy=" << G << endl; 
			energies[sid] = G;

			if ( G < minG ) {
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
	
		// compute statistics
		unsigned int num_confs = m_structures.size() - 1;
		double mean_G = sumG/num_confs;
		double var_G = (sumsqG - (sumG*sumG)/num_confs)/(num_confs-1.0);
		// calculate free energy of folding
		double dG = minG + (var_G - 2*kT*mean_G)/(2*kT) + kT * m_log_num_conformations;

		return new DecoyHistoryFoldInfo(dG<m_deltaG_cutoff, minIndex==m_target_sid, dG, minIndex, mean_G, var_G, minG, p, energies);
	}
	// Now actually implement interesting logic.
	else {
		// First find the residue at which there's a difference
		// Note: what do you do if there's more than one difference?
		const Protein& q = history->getProtein();
		//	  double G = getEnergy(q, sid);
		// vector<uint> getDifferences(const Sequence& p, const Sequence& q) {...}
		vector<double> energies = history->getEnergies();
		vector<uint> diff_indices = p.getDifferences(q);
		uint diffs = diff_indices.size();
		//cout << "Diffs are: " << diffs << endl;
		if (diffs == 0) {
			DecoyHistoryFoldInfo* new_hist = new DecoyHistoryFoldInfo(*history);
			//cout << new_hist->getDeltaG() << " " << history->getDeltaG() << endl;
			return new_hist;
		}
		else {
			// do some interesting logic
			vector<uint>::iterator site_changed = diff_indices.begin();
			for (unsigned int sid = 0; sid < m_structures.size(); sid++) {
				// compute energy for structure sid
				// replace energies[sid] with that new energy.
				energies[sid] = updateEnergy(m_structures[sid], energies[sid], q, p, diff_indices);
				//cout << "new: ind=" << changed_index << " sid=" << sid << " energy=" << energies[sid] << endl;
			}
			// Now we're ready to finish up.
			for (unsigned int sid = 0; sid < m_structures.size(); sid++) {
				double G = energies[sid];
				//cout << "fwh: " << energies[sid] << " sid = " << sid << endl;
				if ( G < minG ) {
					minG = G;
					minIndex = sid;
				}
				// add energy to partition sum
				sumG += G;
				sumsqG += G*G;
			}
			// remove min. energy
			sumG -= minG;
			// cout << "MinG = " << minG << endl;
			sumsqG -= minG*minG;
			// cout << "SumsqG = " << sumsqG << endl;
	  
			// compute statistics
			unsigned int num_confs = m_structures.size() - 1;
			double mean_G = sumG/num_confs;
			double var_G = (sumsqG - (sumG*sumG)/num_confs)/(num_confs-1.0);
			// calculate free energy of folding
			double dG = minG + (var_G - 2*kT*mean_G)/(2*kT) + kT * m_log_num_conformations;
			//cout << dG << " " << var_G << endl;
			// update the history
			return new DecoyHistoryFoldInfo(dG<m_deltaG_cutoff, minIndex==m_target_sid, dG, minIndex, mean_G, var_G, minG, p, energies);
		}
	}
	return NULL;
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

const vector<StructureID> DecoyContactFolder::getStructuresWithResidueContact(uint residue_number) const {
	//cout << residue_number << " " << m_structures_for_residue.size() << endl;
    return m_structures_for_residue[residue_number];
}
