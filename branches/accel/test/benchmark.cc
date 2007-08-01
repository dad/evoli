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
  clock_t start;
  //Random::seed(t);
  // int size = 500;
  int MAX = 1000;

  int protein_length = 500;
  int gene_length = protein_length*3;
  double log_nconf = 100.0*log(10.0);
  double max_dg = 0;
  double sum_dg;
  double mean_dg;
  double sum_sq_dg;  
  double variance_dg;
  
  vector<Contact> m_contacts;
  double array[MAX];
  int sid = 0; 
  Random::seed(11);
  
  ifstream fin("test/data/williams_contact_maps/maps.txt");
       
  DecoyContactFolder folder(protein_length, log_nconf, fin, "test/data/williams_contact_maps/");
 
  if (!folder.good() ){
      cout<< "Error: folder is not good!" << endl;
  } else {
    //CodingDNA g = Gene::createRandomNoStops(gene_length);
	CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
	CodingDNA g2 = g;
    Protein p = g.translate();
    auto_ptr<DecoyFoldInfo> dfi(folder.fold(p));
    DecoyHistoryFoldInfo *dhfi = NULL;
    //auto_ptr<DecoyHistoryFoldInfo> dhfi(folder.foldWithHistory(p, dhfi));
    cout << "Starting dG = " << dfi->getDeltaG() << endl;
    SimpleMutator mut(0.0001);
    cout << " Starting performance test..." << endl;
    start = clock();
    
    for (int i=0; i< MAX;) {    
      g2 = g;
      bool changed = mut.mutate(g2);
      if(changed && g2.encodesFullLength()){
	Protein p = g2.translate();
		
	auto_ptr<DecoyFoldInfo> dfi(folder.fold(p));
	cout << "New dG = " << dfi->getDeltaG() << endl;

	/***************************Error Line Begins***********************/
		
	DecoyHistoryFoldInfo* new_dhfi = folder.foldWithHistory(p, dhfi);
       	if (new_dhfi != NULL) {
	  delete dhfi;
	  dhfi = new_dhfi;
	}
	//	DecoyHistoryFoldInfo* dhfi = folder.foldWithHistory(p, dhfi);
	cout << "FWH dG = " << dhfi->getDeltaG() << endl;
	//	cout << "NFWH dG = " << new_dhfi->getDeltaG() << endl;
	      
		
		/***************************Error Line Ends***********************/
		
		
		//auto_ptr<DecoyHistoryFoldInfo> auto_dhfi(dhfi);
		//cout << dhfi->getDeltaG() << endl;
		if ( dhfi->getDeltaG() <= max_dg && dhfi->getStructure() == (StructureID)sid ) {		 //array[i]= dhfi->getDeltaG();
		  sum_dg += dhfi->getDeltaG();
		  sum_sq_dg += dhfi->getDeltaG()*dhfi->getDeltaG();
		  array[i] = dhfi->getDeltaG();
		  i++;
		  g = g2; 
		}
	  }
    }
  

    clock_t duration = clock() - start;
    cout << " That took: " <<static_cast<double>(duration)/CLOCKS_PER_SEC << " seconds." << endl;
    cout << " That was: " <<static_cast<double>((duration)/CLOCKS_PER_SEC)/MAX << " seconds per protein." << endl;
    cout << " That was: " <<static_cast<double> (MAX/(duration/CLOCKS_PER_SEC))<< " proteins per second." << endl;
    cout << " The sum of the delta Gs is: " << sum_dg << endl;
    mean_dg = sum_dg/MAX;
    cout << " The mean of the delta Gs is: " << mean_dg << endl;
	double sum_sqdev_dg = 0.0;
    for (int i = 0; i <= MAX; i++){
      sum_sqdev_dg += pow(array[i]- mean_dg,2);
	}
    //    cout << " The sum of the squared deviations from the mean is: " << sum_sqdev_dg << endl;
    variance_dg = sum_sqdev_dg/(MAX-1);
    cout << " The variance of the delta Gs is: " << variance_dg << endl;
    cout << " The standard deviation is: " << sqrt(variance_dg) << endl;
  }

  return 0;
}

	
	
