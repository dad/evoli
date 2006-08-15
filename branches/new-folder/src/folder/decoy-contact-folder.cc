#include "decoy-contact-folder.hh"
#include <cmath>


void DecoyContactStructure::read(ifstream& fin) {
	
}


DecoyContactFolder::DecoyContactFolder(int length, vector<DecoyContactStructure>& structs) {
	m_length = length;
	m_structures = structs;
}

FoldInfo DecoyContactFolder::foldProtein(Protein& p) {
	return FoldInfo(0, 1e6);
}
