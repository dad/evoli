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
#include "protein.hh"
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
			Gene g = Gene::createRandomNoStops(protein_length*3);
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

		int protein_length = 300;
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

};


#endif // _T_FOLDER_H__
