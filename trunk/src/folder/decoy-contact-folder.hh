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
	DecoyFoldInfo(double fe, StructureID sid, double mean_G, double var_G, double min_G) : FoldInfo(fe, sid) {
		m_var_G = var_G;
		m_mean_G = mean_G;
		m_min_G = min_G;
	}
	virtual ~DecoyFoldInfo() {}

	double getUnfoldedFreeEnergyMean() const { return m_mean_G; }
	double getUnfoldedFreeEnergyVariance() const { return m_var_G; }
	double getMinEnergy() const { return m_min_G; }
};


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
if ( folder->good() )
{
	FoldInfo fi = folder->fold( p );
	cout << fi.getStructure() << " " << fi.getFreeEnergy() << endl;
}
\endcode

 **/


class DecoyContactFolder : public Folder {
private:
	DecoyContactFolder(); ///< Prohibited constructor.
	DecoyContactFolder( const DecoyContactFolder & ); ///< Prohibited constructor.
	DecoyContactFolder &operator=( const DecoyContactFolder & ); ///< Prohibited assignment operator.

protected:
	int m_length; ///< Length of the proteins to fold.
	double m_log_num_conformations; ///< Fudge factor for the folding process.
	vector<DecoyContactStructure*> m_structures; ///< Vector of the contact maps used as decoys.
//	static const double DecoyContactFolder::contactEnergies [20][20]; ///< Table of contact energies.
	int m_num_folded; ///< Number of proteins folded since creation of the folder object.

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
	 **/
	DecoyContactFolder(int length, double log_num_confs, vector<DecoyContactStructure*>& structs);
	/**
	 * Create DecoyContactFolder; structures are read from disk.
	 *
	 * @param length Length of the proteins to fold.
	 * @param log_num_confs A numerical fudge-factor; use log(10^160).
	 * @param fin An ifstream that contains filenames of the contact map files to be read.
	 * @param dir Directory in which contact maps are stored (must end with "/" or platform-appropriate directory separator)
	 **/
	DecoyContactFolder(int length, double log_num_confs, ifstream& fin, const string& dir);

	~DecoyContactFolder();

	/**
	 * Folds a protein.
	 *
	 * @param s The sequence to be folded.
	 * @return The folding information (of type DecoyFoldInfo).
	 **/
	virtual FoldInfo fold(const Sequence& s);
	/**
	 * Folds a protein.
	 *
	 * @param s The sequence to be folded.
	 * @return The folding information (of type DecoyFoldInfo).
	 **/
	DecoyFoldInfo foldStats(const Sequence& s);
	/**
	 * @param s The sequence whose energy is sought.
	 * @param sid The structure ID of the target conformation.
	 * @return The contact energy of a sequence in the target conformation.
	 **/ 
	double getEnergy(const Sequence& s, StructureID sid) const;
	/**
	 * @return The number of proteins that have been folded since initialization.
	 **/
	int getNumFolded() { return m_num_folded; }

	/**
	 * @return A bool indicating whether the folder is in a usable state.
	 **/
	virtual bool good() const;
	
};

struct ContactMapUtil {
	static void readContactMapsFromFile(ifstream& fin, const string& dir, vector<DecoyContactStructure*>& structs);
};

#endif // DECOY_CONTACT_FOLDER_HH






