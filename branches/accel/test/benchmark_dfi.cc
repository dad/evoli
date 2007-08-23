#include "protein.hh"
#include "random.hh"
#include "folder.hh"
#include "decoy-contact-folder.hh"
#include "fitness-evaluator.hh"
#include "coding-sequence.hh"
#include "compact-lattice-folder.hh"
#include "genetic-code.hh"
#include "folder-util.hh"

#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <time.h>

/*char calc_frac(Protein x, int y){
  string aa = "ACDEFGHIKLMNPQRSTVWY";
 
  for (int i = 0; i < 20; i++){
  int c = 0;  
  for (unsigned int j = 0; j < y; j++){   
  if(aa[i] == x[j]) c++;
  }
  cout << aa[i] << ": has " << c << "amino acids and a fraction of: " 
  << static_cast<float> (c)/y << endl;
  }
  cout << "\n"; 
  }
*/
int main(){
	long t;
	time(&t);
	clock_t start, start2;
	int MAX = 1000;

	int protein_length = 500;
	int gene_length = protein_length*3;
	double log_nconf = 100.0*log(10.0);
	double max_dg2 = 0;
	double sum_dg2;
	double mean_dg2;
	double sum_sq_dg2;  
	double variance_dg2;
  
	vector<Contact> m_contacts;
	double array[MAX];
	int sid = 0; 
	Random::seed(11);
  
	ifstream fin("test/data/williams_contact_maps/maps.txt");
       
	DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
 
	if (!folder.good() ){
		cout<< "Error: folder is not good!" << endl;
		return 1;
	}
	CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg2, sid);
	CodingDNA orig_g = g;
	CodingDNA g2 = g;
	Protein p = g.translate();
	SimpleMutator mut(0.0003);
	Random::seed(111);

	cout << endl << endl << " Starting DFI performance test..." << endl;
	start2 = clock();
    
	for (int i=0; i< MAX;) {    
		g2 = g;
		bool changed = mut.mutate(g2);
		if(changed && g2.encodesFullLength()){
			Protein p = g2.translate();
			auto_ptr<DecoyFoldInfo> dfi(folder.fold(p));
		
			if ( dfi->getDeltaG() <= max_dg2 && dfi->getStructure() == (StructureID)sid ) {	  //array[i]= dfi->getDeltaG();
				sum_dg2 += dfi->getDeltaG();
				sum_sq_dg2 += dfi->getDeltaG()*dfi->getDeltaG();
				array[i] = dfi->getDeltaG();
				i++;
				g = g2; 
			}
		}
	}
  

	clock_t duration2 = clock() - start2;

	cout << "/*********************************DECOY FOLD INFO STATS************************************/" << endl << endl;


	double seconds2 = static_cast<double> (duration2)/CLOCKS_PER_SEC;
	cout << " That took: " << seconds2 << " seconds." << endl;
	double seconds_per_protein2 = seconds2/static_cast<double>(MAX);
	cout << " That was: " << seconds_per_protein2 << " seconds per protein." << endl;
	double proteins_per_second2 = static_cast<double>(MAX)/seconds2;
	cout << " That was: " << proteins_per_second2 << " proteins per second." << endl;
	cout << " The sum of the delta Gs is: " << sum_dg2 << endl;
	mean_dg2 = sum_dg2/MAX;
	cout << " The mean of the delta Gs is: " << mean_dg2 << endl;
	double sum_sqdev_dg2 = 0.0;
	for (int i = 0; i <= MAX; i++){
		sum_sqdev_dg2 += pow(array[i]- mean_dg2,2);
	}
	//    cout << " The sum of the squared deviations from the mean is: " << sum_sqdev_dg2 << endl;
	variance_dg2 = sum_sqdev_dg2/(MAX-1);
	cout << " The variance of the delta Gs is: " << variance_dg2 << endl;
	cout << " The standard deviation is: " << sqrt(variance_dg2) << endl << endl << endl;

	return 0;
}

	
