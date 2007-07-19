/*
This file is part of the E.voli project.
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


#ifndef DECOY_CONTACT_FOLDER_HH
#define DECOY_CONTACT_FOLDER_HH

#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>

#include "folder.hh"
#include "protein-contact-energies.hh"

using namespace std;

/**
 * Stores folding information.
 **/
 
class DecoyFoldInfo : public FoldInfo {
protected:
	double m_var_G; 
	double m_mean_G;
	double m_min_G;
	
public:
  DecoyFoldInfo( const FoldInfo& fi, double mean_G, double var_G, double min_G)
    : FoldInfo( fi.foldIsStable(), fi.foldMatchesTarget(), fi.getDeltaG(), fi.getStructure()) {
    m_var_G = var_G;
    m_mean_G = mean_G;
    m_min_G = min_G;
  }
  DecoyFoldInfo( bool fold_is_stable, bool fold_is_target, double fe, StructureID sid, double mean_G, double var_G, double min_G)
    : FoldInfo( fold_is_stable, fold_is_target, fe, sid)
  {
    m_var_G = var_G;
    m_mean_G = mean_G;
    m_min_G = min_G;
  }

  virtual ~DecoyFoldInfo() {}
  
  double getUnfoldedDeltaGMean() const { return m_mean_G; }
  double getUnfoldedDeltaGVariance() const { return m_var_G; }
  double getMinEnergy() const { return m_min_G; }
};

/*************************** Beginning of DecoyHistoryFoldInfo Implementation********************/

class DecoyHistoryFoldInfo : public DecoyFoldInfo {
protected:
  Protein m_p;
  vector<double> m_energy_list;

public:
  DecoyHistoryFoldInfo( const DecoyFoldInfo& fi, const Protein& p, vector<double>& energy_list) : DecoyFoldInfo( fi.foldIsStable(), fi.foldMatchesTarget(), fi.getDeltaG(), fi.getStructure(), fi.getUnfoldedDeltaGMean(), fi.getUnfoldedDeltaGVariance(), fi.getMinEnergy()), m_p(p) {
    m_energy_list = energy_list;
  }
  DecoyHistoryFoldInfo( bool fold_is_stable, bool fold_is_target, double fe, StructureID sid, double mean_G, double var_G, double min_G, const Protein& p, vector<double>& energy_list)
    : DecoyFoldInfo( fold_is_stable, fold_is_target, fe, sid, mean_G, var_G, min_G), m_p(p)
	{
	  m_energy_list = energy_list;
	}

  virtual ~DecoyHistoryFoldInfo() {}

  const Protein getProtein() const { return m_p; }
  vector<double> getEnergies() const { return m_energy_list; }
};


/******************* End of DecoyHistoryFoldInfo Implementation*************************/


/**
 * Stores a contact structure: a list of contacts.
 **/

class DecoyContactStructure : public ContactStructure {
protected:
	vector<Contact> m_contacts;
public:
	DecoyContactStructure() {}
	virtual ~DecoyContactStructure() {}

	/**
	 * Read from a file.
	 **/
	void read(istream& fin);
	virtual const vector<Contact>& getContacts() const { return m_contacts; }
	int getMaxResidueNumber();
};

/** \brief  The DecoyContactFolder implements the protein folding model of P. D. Williams
 * et al., PLoS Comp. Biol. 2:e69.
 *
 * This folder implementation folds proteins using a set of decoy structures given as
 * contact maps. The folder can either read the contact maps from disk or accept
 * ready-made objects of type @ref ContactStructure.
 *
 * Here is an example of the typical usage for this folder:
 * \code
// the sequence to fold
Protein p( "STLRFVAVGDWGGVPNAPFHTAREMANAKEIARTVQIMGADFIMSLGDNFYFTGVHDANDKRFQETFEDVFSDRALRNIPWYVLAGNHDHLGNVSAQIAYSKISKRWNFPSPYYRLRFKVPRSNITVAIFMLDTVMLCGNSDDFVSQQPEMPRDLGVARTQLSWLKKQLAAAKEDYVLVAGHYPIWSIAEHGPTRCLVKNLRPLLAAYGVTAYLCGHDHNLQYLQDENGVGYVLSGAGNFMDPSVRHQRKVPNGYLRFHYGSEDSLGGFTYVEIGSKEMSITYVEASGKSLFKTSLPRRP" );

ifstream fin("test/data/williams_contact_maps/maps.txt"); // file which contains the names of the contact maps

// initialize the folder with protein length 300, fudge factor log(10^160),
// in directory "test/data/williams_contact_maps/":
DecoyContactFolder folder( 300, 160.0*log(10.0), fin, "test/data/williams_contact_maps/");
if ( folder.good() )
{	// we store the result from the fold() function in an auto_ptr,
	// so that we don't have to worry about memory management
	auto_ptr<FoldInfo> fi( folder.fold( p ) );
	cout << fi->getStructure() << " " << fi->getDeltaG() << endl;
}
\endcode

 **/

class DecoyContactFolder : public DGCutoffFolder {
private:
	DecoyContactFolder(); ///< Prohibited constructor.
	DecoyContactFolder( const DecoyContactFolder & ); ///< Prohibited constructor.
	DecoyContactFolder &operator=( const DecoyContactFolder & ); ///< Prohibited assignment operator.

protected:
	int m_length; ///< Length of the proteins to fold.
	double m_log_num_conformations; ///< Fudge factor for the folding process.
	vector<DecoyContactStructure*> m_structures; ///< Vector of the contact maps used as decoys.
//	static const double DecoyContactFolder::contactEnergies [20][20]; ///< Table of contact energies.
	mutable int m_num_folded; ///< Number of proteins folded since creation of the folder object.

  // List of structure lists.  Entry i is a list of structures with contacts involving residue i.
  vector<vector<StructureID> > m_structures_for_residue;

  /**
   * Fetch the structures for a given residue;
   **/
  /*  const vector<StructureID> getStructuresWithResidueContact(uint residue_number) const {
    return m_structures_for_residue[residue_number];
    }*/

  void initializeStructuresForResidues();

	/**
	* Wrapper function to encapsulate the lookup of the
	* contact energy from a table.
	*
	* @param residue1 The first of the two contacting residues.
	* Mapping from int to residue type is as in class @ref ProteinContactEnergies.
	* @param residue2 The second of the two contacting residues.
	* @return The contact energy between the two residues.
	**/
	double contactEnergy( int residue1, int residue2 ) const {
		return ProteinContactEnergies::ProteinContactEnergies::WilliamsPLoSCB2006[residue1][residue2]; }
	//=MJ85TableVI[residue1][residue2]; }

public:
	/**
	 * Create DecoyContactFolder using already existing structure objects.
	 *
	 * @param length Length of the proteins to fold.
	 * @param log_num_confs A numerical fudge-factor; use log(10^160).
	 * @param structs A vector of pointers to contact structures. The folder object
	 * assumes ownership over these structures and will delete them upon exit.
	 * @param deltaG_cutoff The DeltaG cutoff, as in \ref DGCutoffFolder.
	 * @param target_sid The target structure ID, as in \ref DGCutoffFolder.
	 **/
	DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs, double deltaG_cutoff = 0., StructureID target_sid = -1 );
	/**
	 * Create DecoyContactFolder; structures are read from disk.
	 *
	 * @param length Length of the proteins to fold.
	 * @param log_num_confs A numerical fudge-factor; use log(10^160).
	 * @param fin An ifstream that contains filenames of the contact map files to be read.
	 * @param dir Directory in which contact maps are stored (must end with "/" or platform-appropriate directory separator)
	 * @param deltaG_cutoff The DeltaG cutoff, as in \ref DGCutoffFolder.
	 * @param target_sid The target structure ID, as in \ref DGCutoffFolder.
	 **/
	DecoyContactFolder(int length, double log_num_confs, ifstream& fin, const string& dir, double deltaG_cutoff = 0., StructureID target_sid = -1 );

	~DecoyContactFolder();

	/**
	 * Folds a protein. See Folder::fold() for details.
	 *
	 * @param s The sequence to be folded.
	 * @return The folding information (of type DecoyFoldInfo).
	 **/
	virtual DecoyFoldInfo* fold(const Protein& p) const;

	/**
	 **/
	virtual DecoyHistoryFoldInfo* foldWithHistory(const Protein& p, const DecoyHistoryFoldInfo* history) const;

	/**
	 * @param s The sequence whose energy is sought.
	 * @param sid The structure ID of the target conformation.
	 * @return The contact energy of a sequence in the target conformation.
	 **/ 
	double getEnergy(const Protein& s, StructureID sid) const;

	/**
	@return The number of proteins that have been folded so far with this Folder instance.
	*/
	uint getNumFolded() const { return m_num_folded; }

	/**
	This function assesses whether the folder has been properly initialized.
	@return True if the folder is in good working order, False otherwise.
	 **/
  virtual bool good() const;

  const DecoyContactStructure* getStructure(StructureID sid) const { return m_structures[sid]; }

  const vector<StructureID> getStructuresWithResidueContact(uint residue_number) const {
    return m_structures_for_residue[residue_number];
  }


  /*vector <Contact>& getContacts() const { return m_contacts; }
  virtual const vector<Contact>& getContacts() const { return m_contacts; }
  */
};

  /*
  const DecoyContactStructure* getStructure(StructureID sid) const {
    return m_structures[sid];
    //  int k = Random::rint(100);
    }*/

struct ContactMapUtil {
	static void readContactMapsFromFile(ifstream& fin, const string& dir, vector<DecoyContactStructure*>& structs);
};

#endif // DECOY_CONTACT_FOLDER_HH





