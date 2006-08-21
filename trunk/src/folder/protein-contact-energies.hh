#ifndef PROTEIN_CONTACT_ENERGIES_HH
#define PROTEIN_CONTACT_ENERGIES_HH



/** \brief A static class containing protein contact energies.
*
* This class provides several energy tables published in the literature. The mapping from 
* int to amino acid in all tables is, starting counting at 0 and incrementing by one for each
* residue, CYS, MET, PHE, ILE, LEU, VAL, TRP, TYR, ALA, GLY, THR, SER, GLN, ASN, GLU, 
* ASP, HIS, ARG, LYS, PRO.
**/

class ProteinContactEnergies {
private:
	ProteinContactEnergies();
	ProteinContactEnergies( const ProteinContactEnergies & );
	const ProteinContactEnergies & operator=( const ProteinContactEnergies & );

public:
	/**
	* Contact energies according to Miyazawa and Jernigan, Residue-residue potentials
	* with a favorable contact pair term and an unfavorable high packing density term,
	* for simulation and threading. J. Mol. Biol. (1996) 256:623-644, Table III upper triangle,
	* values e_{ij}.
	*/
	static const double MJ96TableIII[20][20];

	/**
	* Contact energies according to Miyazawa and Jernigan, Estimation of
	* effective interresidue contact energies from protein crystal
	* structures: Quasi-chemical approximation. Macromolecules 18:534-552,
	* (1985), Table V, values e_{ij}.
	*/
	static const double MJ85TableV[20][20];


	/**
	* Contact energies according to Miyazawa and Jernigan, Estimation of
	* effective interresidue contact energies from protein crystal
	* structures: Quasi-chemical approximation. Macromolecules 18:534-552,
	* (1985), Table VI, values e_{ij}+e_{rr}-e_{ir}-e_{jr}.
	*/
	static const double MJ85TableVI[20][20];
};


#endif //PROTEIN_CONTACT_ENERGIES_HH






