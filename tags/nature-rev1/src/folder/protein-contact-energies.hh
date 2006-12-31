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

	/**
	 * Contact energies according to Williams et al. PLoS Computational
	 * Biology 2006.  Energies are MJ85TableVI with a CYS-CYS energy equal
	 * to the SER-SER energy.
	 */
	static const double WilliamsPLoSCB2006[20][20];
};


#endif //PROTEIN_CONTACT_ENERGIES_HH






