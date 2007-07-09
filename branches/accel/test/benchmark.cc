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

char calc_frac(Protein x, int y){
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


int main(){
  long t;
  time(&t);
  clock_t start;
  Random::seed(t);
  int size = 300;
  int MAX = 10000;

  int protein_length = 300;
  int gene_length = protein_length*3;
  double log_nconf = 100.0*log(10.0);
  double max_dg = 0;
  double sum_dg;
  double mean_dg;
  double sum_sqdev_dg;  
  double variance_dg;
  
  double array[MAX];
  int sid = 0; //previously 574;
  Random::seed(11);
  
  ifstream fin("data/williams_contact_maps/maps.txt");
  
  cout << " Starting performance test" << endl;
  DecoyContactFolder folder(protein_length, log_nconf, fin, "data/williams_contact_maps/");
  //
  if (!folder.good() )
      cout<< "G, DAD we have a problem!!" << endl;
  else {
    // CodingDNA g = FolderUtil::getSequenceForStructure( folder, gene_length, max_dg, sid);
    CodingDNA g("TTCCGTTTTTTCCGTGTATGTGACCCGATCATATCCTTACTTATGCATATCTGGCCAGGCCCAATTCGAACACGCCGGCGTTCAAAAATATGGATAAAGGAATTTTGTATCAAGTGTGAACTATGTCCGGAAACATTCCGTCACAAAATGATAGTGTGCGAAAAAAATCCCCCGAGGCTGATATACGCAAAACTCCTACATAGACGCATGAGACGCTTAATTGATATACTAGGTATGATACACCATCGTAAAATCCGTGAGGAAGAGCAGCCGTGCCCTGCTCAGGTGATACTCCCGTTTTTGAGACTAGAGGTGAAGCACAATGTGGCCGCGATTCGCCGATGGACCGCGCCGACTGTGCTAAAAAAGCTCAAATTGGACATAATCCGCTTGGTTCAGATCGTGCCCAAGCTCTGTGAATGCTGTTCATTCGATTCCATCGAAGACTGCCGACGGAAAAGATTTTCGATTAAGGATAGATCGATAACAATAAAGACGGAGGCATGCAAAATCTTTTTTAGGCGGCGCAGATTTCGTCTGCGTACAATCGTAGATGAGCAAGTTCTCTCCCGTTCATCAGCACGGCGACGGAGACGCAAGATCGCCGCCTTCCTTTTTAGCGACCGTTTGCCAATCAGACTCGAACGACGTGATCGGAAGAAGAAAAGGCAAATATGCGAGGAGGAAGAGCCGTTTGACTGCGAGCCGGAGGCTCGCGACATACGATTGAAATTTAGGAGACTGCTAAGACCTAAGTCCTTTATTCGGGACGCGGAACGGGAAACCGATAAGCAAATTCGTTCGAATAGATGCATCGTTTCGGAAGTTTTAAAAGCCATCTCGCACGAGTCTATAGTCCTCAGGATCAAGCCCACACCGTTATTACGTATGAATGCCTTT");
    //cout << g << endl;
   // int nfolded = folder.getNumFolded();  
  
    //   cout << "num folded: " << (folder.getNumFolded()) << endl;
   start = clock();
    for (int i=0; i< MAX;) {    
      SimpleMutator mut(0.001);
      CodingDNA g2 = g;
      bool changed = mut.mutate(g2);
      if(changed){
	//	cout << "Ahhh, it has been mutated!" << endl;
	Protein p = g2.translate();
        auto_ptr<FoldInfo> fi( folder.fold(p) );
	if ( fi->getDeltaG() <= max_dg && fi->getStructure() == (StructureID)sid ) {
	  array[i]= fi->getDeltaG();
	  sum_dg +=  fi->getDeltaG();
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

	
	
