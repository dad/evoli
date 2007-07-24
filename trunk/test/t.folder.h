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


#ifndef _T_FOLDER_H__
#define _T_FOLDER_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "coding-sequence.hh"
#include "protein.hh"
#include "folder-util.hh"
#include "random.hh"

#include <fstream>
#include <cmath>
#include <memory>

struct TEST_CLASS( folder_basic )
{
	const static int side_length = 5;
	const static int gene_length = side_length*side_length*3;

	void TEST_FUNCTION( init_lattice )
	{
		CompactLatticeFolder folder(side_length);
		Protein p("CSVMQGGKTVFQMPIIERVMQAYNI"); //Gene::createRandomNoStops(gene_length);
		auto_ptr<FoldInfo> fi( folder.fold(p) );
		TEST_ASSERT( fabs(fi->getDeltaG()-0.564) < 1e-2 );
		TEST_ASSERT( fi->getStructure() == (StructureID)225 );
		return;
	}
	void TEST_FUNCTION( init_decoy )
	{
		ifstream fin("test/data/williams_contact_maps/maps.txt");
		int protein_length = 300;
		double log_nconf = 160.0*log(10.0);
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
		TEST_ASSERT( folder.good() );
		if (!folder.good() )
			return;
		int num_to_fold = 1;
		for (int j=0; j<num_to_fold; j++) {
			CodingDNA g = CodingDNA::createRandomNoStops(protein_length*3);
			Protein p = g.translate();
			auto_ptr<FoldInfo> fi( folder.fold(p) );
			TEST_ASSERT(fi->getStructure()>-1);
			//cout << "folded:" << tab << fi.getStructure() << tab << fi.getDeltaG() << endl;
		}
		// Clean up
		return;
	}

	void TEST_FUNCTION( decoy_known_stability ) {
		ifstream fin("test/data/rand_contact_maps/maps.txt");
		int protein_length = 300;
		double log_nconf = 160.0*log(10.0);
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/rand_contact_maps/");
		TEST_ASSERT( folder.good() );
		if (!folder.good() )
			return;
		Protein p("PRPEEEKKKREREEKRRKEDKLERIRDLPRKILKMIVEPKRRKKGETEDDDEKESKRREEMEKFKREFFTICIKLLECEEEMARRREKRREEEDIDSLRELMKDCRRFIDDPRRVEQQSQRLDFRSRRKLEDEKDDEDKRKPDFLFEFEMCEEDMRRRPLDRVKDICRVCCEMDEEEEIREEEEFFRPEEEDMKLKSFRESFKDVRRCILRKFEKSRREKSAEFLRHEIPMFSSEDEEDRKKKDRRRQRPMMRHFMKRIKEKEEERKKREFKEQEEPKPKSFKWKTEEEMEELGEQEKRV");
		auto_ptr<FoldInfo> fi( folder.fold(p) );
		//cout << endl << fi->getDeltaG() << " " << fi->getStructure() << endl;
		//cout << p << endl;
		TEST_ASSERT( fabs(fi->getDeltaG()-0.00732496)<1e-4 );
		TEST_ASSERT( fi->getStructure() == (StructureID)0 );
	}


	void TEST_FUNCTION( contact_reader ) {
		string seq = "STLRFVAVGDWGGVPNAPFHTAREMANAKEIARTVQIMGADFIMSLGDNFYFTGVHDANDKRFQETFEDVFSDRALRNIPWYVLAGNHDHLGNVSAQIAYSKISKRWNFPSPYYRLRFKVPRSNITVAIFMLDTVMLCGNSDDFVSQQPEMPRDLGVARTQLSWLKKQLAAAKEDYVLVAGHYPIWSIAEHGPTRCLVKNLRPLLAAYGVTAYLCGHDHNLQYLQDENGVGYVLSGAGNFMDPSVRHQRKVPNGYLRFHYGSEDSLGGFTYVEIGSKEMSITYVEASGKSLFKTSLPRRP";

		const char* fname = "1qhwA_6_CB.cmap";
		ifstream fin;
		DecoyContactStructure structure;
		string filename = string("test/data/williams_contact_maps/")+fname;
		fin.open(filename.c_str());
		TEST_ASSERT( fin.good() );
		if (!fin.good()) // if we can't read the contact maps, bail out
			return;
		structure.read(fin);
		fin.close();
		fin.open(filename.c_str());
		vector<Contact> contacts = structure.getContacts();
		int r1, r2;
		char r1aa, r2aa;
		for (int i=0; i<contacts.size() && !fin.eof();) {
			fin >> r1 >> r1aa >> r2 >> r2aa;
			if ( abs(r1-r2)>1 ) {
				TEST_ASSERT(contacts[i].first == r1);
				TEST_ASSERT(contacts[i].second == r2);
				TEST_ASSERT(seq[r1] == r1aa);
				TEST_ASSERT(seq[r2] == r2aa);
				i++;
			}
			//cout << contacts[i].first << " " << r1 << " " << contacts[i].second << " " << r2 << " " << seq[r1] << " " << r1aa << " " << seq[r2] << " " << r2aa << endl;
		}
		fin.close();
	}

	void TEST_FUNCTION( williams_test ) {
		// Compare to 
		// Load the real sequence
		// Compute free energy
		string native_1qhw_seq = "STLRFVAVGDWGGVPNAPFHTAREMANAKEIARTVQIMGADFIMSLGDNFYFTGVHDANDKRFQETFEDVFSDRALRNIPWYVLAGNHDHLGNVSAQIAYSKISKRWNFPSPYYRLRFKVPRSNITVAIFMLDTVMLCGNSDDFVSQQPEMPRDLGVARTQLSWLKKQLAAAKEDYVLVAGHYPIWSIAEHGPTRCLVKNLRPLLAAYGVTAYLCGHDHNLQYLQDENGVGYVLSGAGNFMDPSVRHQRKVPNGYLRFHYGSEDSLGGFTYVEIGSKEMSITYVEASGKSLFKTSLPRRP";
		string stable_seq = "IREDEWEVRRKKKDVRWDMKKQEEDKKKWEMMRCFMCCIHKKRKTERWEDEWMPEEEEKRMRELWEHCIEMIMCWWCCDEEMREDPWMRFWWKEEKRMKEMCRECKKWWRVTEEICMDRHMLCECWKICIKKNKMCEEEFDMCLCIRIKIKKKRCDCERERDKCHNACMWKINMFPLCLEEEEEMEWEFCWCRKIEPWIKRPVQFPGWIFCCKRKRKMRFEKGKGWCWCMCECEEEHEEEECMCKHREMEKSCIEKGGIKFKKGDKKEMDMREQDCCDCKTWKWKEEKEMEGMAECRMMA";

		int protein_length = stable_seq.size();
		double log_nconf = 160.0*log(10.0);
		ifstream fin("test/data/williams_contact_maps/maps.txt");
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
		TEST_ASSERT( folder.good() );
		if (!folder.good() )
			return;
		Protein p(stable_seq);
		auto_ptr<FoldInfo> fi( folder.fold( p ) );
		TEST_ASSERT(fi->getStructure()==34);
		// Clean up
	}

	void TEST_FUNCTION( sequence_for_structure )
	{
		CompactLatticeFolder folder(side_length);
		double max_dg = -1;
		int sid = 574;
		Random::seed(11);
		int nfolded = folder.getNumFolded();
		CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
		//cout << "num folded: " << (folder.getNumFolded() - nfolded) << endl;
		Protein p = g.translate();
		auto_ptr<FoldInfo> fi( folder.fold(p) );
		TEST_ASSERT( fi->getDeltaG() <= max_dg );
		TEST_ASSERT( fi->getStructure() == (StructureID)sid );
		return;
	}

	bool getAminoAcidIndices(const Protein& p, vector<unsigned int>& aa_indices)
	{
		int index = 0;
		for (unsigned int i=0; i<p.size() && index >= 0; i++) {
			index = GeneticCodeUtil::aminoAcidLetterToIndex(p[i]);
			aa_indices[i] = index;
			//cout << index << " " << i << " " << p[i] << endl;
		}
		return index >= 0;
	}

	void TEST_FUNCTION( get_energy ) {
	  int protein_length = 500;
		double log_nconf = 160.0*log(10.0);
		ifstream fin("test/data/williams_contact_maps/maps.txt");
		TEST_ASSERT( fin.good() );
		if (!fin.good()) // if we can't read the contact maps, bail out
			return;
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
		TEST_ASSERT(folder.good());
		if (!folder.good())
			return;
		// Test whether energies yield same value as fold()

		for (int i = 0; i<100; i++) {

		  CodingDNA g = CodingDNA::createRandomNoStops(protein_length*3);
		  Protein p = g.translate();

		  auto_ptr<FoldInfo> fi( folder.fold( p ) );
		  double realdG = fi->getDeltaG();

		  double kT = 0.6;
		  double minG = 1e50;
		  int minIndex = -1;
		  vector<unsigned int> aa_indices(p.size());
		  double sumG = 0.0;
		  double sumsqG = 0.0;

		  bool valid = getAminoAcidIndices(p, aa_indices);
		  for ( unsigned int sid = 0; sid < folder.getNumStructures(); sid++) {
			double G = folder.getEnergy(p, sid);
			// check if binding energy is lower than any previously calculated one
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
		
		  unsigned int num_confs = folder.getNumStructures() -1;
		  double mean_G = sumG/num_confs;
		  double var_G = (sumsqG - (sumG*sumG)/num_confs)/(num_confs-1.0);
		  // calculate free energy of folding
		  double dG = minG + (var_G - 2*kT*mean_G)/(2*kT) + kT * log_nconf;
		  TEST_ASSERT(abs(dG - realdG) < 1e-6);
		}
	}

};


#endif // _T_FOLDER_H__
