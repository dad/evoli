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
		int protein_length = 500;
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
		Protein p("PRPEEEKKKREREEKRRKEDKLERIRDLPRKILKMIVEPKRRKKGETEDDDEKESKRREEMEKFKREFFTICIKLLECEEEMARRREKRREEEDIDSLRELMKDCRRFIDDPRRVEQQSQRLDFRSRRKLEDEKDDEDKRKPDFLFEFEMCEEDMRRRPLDRVKDICRVCCEMDEEEEIREEEEFFRPEEEDMKLKSFRESFKDVRRCILRKFEKSRREKSAEFLRHEIPMFSSEDEEDRKKKDRRRQRPMMRHFMKRIKEKEEERKKREFKEQEEPKPKSFKWKTEEEMEELGEQEKRV");
		int protein_length = p.length();
		double log_nconf = 160.0*log(10.0);
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/rand_contact_maps/");
		TEST_ASSERT( folder.good() );
		if (!folder.good() )
			return;
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
		//cout << "xx" << p << endl;
		auto_ptr<FoldInfo> fi( folder.fold(p) );
		TEST_ASSERT( fi->getDeltaG() <= max_dg );
		TEST_ASSERT( fi->getStructure() == (StructureID)sid );
		return;
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
		//cout << "Williams:" << endl;
		//cout << "folded:" << tab << fi.getStructure() << tab << fi.getDeltaG() << endl;
		// Clean up
	}

 	void TEST_FUNCTION( fold_with_history )
	{
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
		return;

		//FoldInfo* real_fi = folder.fold(p);
		//auto_ptr<FoldInfo> fi( real_fi );

		DecoyHistoryFoldInfo *dhfi = NULL;
		dhfi = folder.foldWithHistory(p, dhfi);
		auto_ptr<DecoyHistoryFoldInfo> auto_dhfi(dhfi);
		
		TEST_ASSERT(fi->getDeltaG() == auto_dhfi->getDeltaG());

		TEST_ASSERT(auto_dhfi->getProtein() == p);
	}

 	void TEST_FUNCTION( with_history_twice )
	{
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

		//FoldInfo* real_fi = folder.fold(p);
		//auto_ptr<FoldInfo> fi( real_fi );

		DecoyHistoryFoldInfo *dhfi = NULL;
		dhfi = folder.foldWithHistory(p, dhfi);

		dhfi = folder.foldWithHistory(p, dhfi);
	
		auto_ptr<DecoyHistoryFoldInfo> auto_dhfi(dhfi);
		
		TEST_ASSERT(fi->getDeltaG() == auto_dhfi->getDeltaG());

		TEST_ASSERT(auto_dhfi->getProtein() == p);
	}


	
	void TEST_FUNCTION( compare_contacts )
	{
	  int protein_length = 500;
	  CodingDNA g = CodingDNA::createRandomNoStops(protein_length*3);
	  Protein p = g.translate();
		double log_nconf = 160.0*log(10.0);
		ifstream fin("test/data/williams_contact_maps/maps.txt");
		TEST_ASSERT( fin.good() );
		if (!fin.good()) // if we can't read the contact maps, bail out
			return;
		DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
		TEST_ASSERT(folder.good());
		if (!folder.good())
			return;

		auto_ptr<FoldInfo> fi( folder.fold( p ) );
		DecoyHistoryFoldInfo *dhfi = folder.foldWithHistory(p, NULL);
		auto_ptr<DecoyHistoryFoldInfo> auto_dhfi(dhfi);
		
		//vector<vector<StructureID> >& lookup_table = folder.getStructureLookupTable();
		// compare structures against contact lookups
		TEST_ASSERT(abs(fi->getDeltaG()- auto_dhfi->getDeltaG()) < 1e-6);
		//vector<DecoyContactStructure*>& structures = folder.getStructures();
		
		// Assume existence of a lookup table
		// Pick 1000 structure IDs at random
		for ( int i = 0; i < 1000; i++) {
		  StructureID sid = (StructureID)Random::rint(folder.getNumStructures());
		  // cout << sid << endl;
		  const DecoyContactStructure* structure = folder.getStructure(sid);
		  if (structure == NULL) {
			i--;
			continue;
		  }
		  // Within each structure, pick 10 contacts at random
		  for (int j = 0; j < 10; j++) {
		    const vector<Contact> contacts = structure->getContacts();
		    int cid = Random::rint(contacts.size());
			//cout << cid << " " << contacts.size() << endl;
		    int first_aa_index = contacts[cid].first;
		    // Success means that the table's entry at first_aa_index contains structure ID sid.
		    const vector<StructureID>& structures = folder.getStructuresWithResidueContact(first_aa_index);
			//cout << "got here" << endl;
		    bool found = false;
		    for (int k = 0; k < structures.size() && !found; k++) {
		      found = (structures[k] == sid);
		    }
		    TEST_ASSERT(found);
		  }
		}
	
		// Take first contact residue, look up in table
		// Confirm that the structure ID exists in the table.
		
		
		
		
		//FoldInfo* real_fi = folder.fold(p);
		//auto_ptr<FoldInfo> fi( real_fi );
	}

	void TEST_FUNCTION( no_history )
	{
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

		for (int i=0; i<10; i++) {
		  CodingDNA g = CodingDNA::createRandomNoStops(protein_length*3);
		  Protein p = g.translate();
		  auto_ptr<FoldInfo> fi( folder.fold( p ) );
		  DecoyHistoryFoldInfo *dhfi = folder.foldWithHistory(p, NULL);
		  auto_ptr<DecoyHistoryFoldInfo> auto_dhfi(dhfi);
		  TEST_ASSERT(dhfi->getProtein() == p);
		  const vector<double>& energies = dhfi->getEnergies();
		  TEST_ASSERT(energies.size() == folder.getNumStructures());
		  // cout << fi->getDeltaG() << " " << auto_dhfi->getDeltaG() << endl;
		  TEST_ASSERT(abs(fi->getDeltaG()-auto_dhfi->getDeltaG()) < 1e-6);
		}
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
	void TEST_FUNCTION( with_some_history )
	{
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

		for (int i=0; i<10; i++) {
		  CodingDNA g = CodingDNA::createRandomNoStops(protein_length*3);
		  Protein p = g.translate();
		  DecoyHistoryFoldInfo *dhfi1 = folder.foldWithHistory(p, NULL);
		  //cout << dhfi1 <<"\n\n"<< endl;//Error checking
		  auto_ptr<DecoyHistoryFoldInfo> auto_dhfi1(dhfi1);	
		  //cout << dhfi1 <<"\n\n"<< endl;//Error checking
		  DecoyHistoryFoldInfo *dhfi2 = folder.foldWithHistory(p, dhfi1);
		  //cout << dhfi2 <<"\n\n"<< endl;//Error checking
		  TEST_ASSERT(dhfi2 != NULL);		 
		  if (dhfi2 != NULL){
			auto_ptr<DecoyHistoryFoldInfo> auto_dhfi2(dhfi2);
			const vector<double>& energies1 = dhfi1->getEnergies();
			const vector<double>& energies2 = dhfi2->getEnergies();
			TEST_ASSERT(energies1.size() == energies2.size());
			for (int j=0; j<energies1.size(); j++) {
			  TEST_ASSERT(abs(energies1[j] - energies2[j]) < 1e-6);
			}
			TEST_ASSERT(abs(dhfi1->getDeltaG() - dhfi2->getDeltaG()) < 1e-6);
		  }
		
		}
	}
 
// New test: make sure invalidated structures are added only once.

// New test: inList function.  

//New test to ensure that dhfi and dfi are the same. 
};


#endif // _T_FOLDER_H__
