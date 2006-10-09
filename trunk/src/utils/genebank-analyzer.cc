#include "protein.hh"
#include "compact-lattice-folder.hh"
#include "folder.hh"
#include "translator.hh"
#include "population.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "gene-util.hh"

#include <algorithm>
#include <fstream>
#include <string>

struct GenebankData
{
        int id;
        int parent_id;
        int count;
        int birth_time;
        int struct_id;
        char coalescent;
        double free_energy;
        double w_saved;
        double w_new;
        double sensitivity;
        double sensitivity_ns;
        double nu;
        double fop;
        Gene g;
};


/* need to add more fields here to increase the number of parameters 
	in the new genebank output files
*/
struct Params
{
		
		string eval_type;
		int structure_ID;
		int free_energy_cutoff;
		int free_energy_min;
		int protein_length;
		double u;
        double tr_cost;
        double ca_cost;
        double transl_error_rate;
        int transl_acc_wt;
        double transl_error_wt;
		int N;
        int window_size;
        int equilibration_time;
        int coalescent_time;
};


void readLine( istream &in )
{
        char c;

        do
        {
			in.get( c );
        }
        while( in && c != '\n' );
}



void readGenebank( istream &in, Params &p, vector<GenebankData> &v )
{
        char c;
        string s;

        do
        {
                in.get( c );
                if ( c == '#' )
                {
                        in >> s;
						/* get the evaluation type */
						if ( s == "evaluation" ) {
							in >> s;
							in >> p.eval_type;
						}
						/* get the structure ID*/
						else if ( s == "structure" ) {
							in >> s;
							in >> p.structure_ID;
						}
                        // get free energy cutoff
                        else if ( s == "free"  )
                        {
                                in >> s;
                                in >> s;
								if ( s == "maximum:" ) {
									in >> p.free_energy_cutoff;
								}
								else
									in >> p.free_energy_min;
                        }
                        // get mutation rate
                        else if ( s == "per-site" )
                        {
                                in >> s;
                                in >> s;
                                in >> s;
                                in >> p.u;
                        }
                        // get tr cost factor
                        else if ( s == "transl." )
                        {
                                in >> s;
                                if ( s == "robust." )
                                {
                                        in >> s;
                                        in >> s;
                                        in >> p.tr_cost;
                                }
                                else if ( s == "error" )
                                {
                                        in >> s;
										if ( s == "rate" )
										{
											in >> s;
											in >> s;
											in >> p.transl_error_rate;
										}
										else
										{
											in >> p.transl_error_wt;
										}
                                }
								else if ( s == "accuracy" )
								{
										in >> s;
										in >> p.transl_acc_wt;
								}
								
                        }
                        // get ca cost factor
                        else if ( s == "codon" )
                        {
                                in >> s;
                                in >> s;
                                in >> s;
                                in >> p.ca_cost;
                        }
                        // get population size
                        else if ( s == "population" )
                        {
                                in >> s;
                                in >> s;
                                in >> p.N;
                        }
                        // get window size
                        else if ( s == "window" )
                        {
                                in >> s;
                                in >> s;
                                in >> p.window_size;
                        }
                        // get equilibration time
                        else if ( s == "equilibration" )
                        {
                                in >> s;
                                in >> p.equilibration_time;
                        }
                        readLine( in );
                }
                else
                {
                        in.putback( c );
                        break;
                }
        }
        while ( 1 );

        GenebankData d;
        v.clear();

        while ( in >> d.id )
        {
                in >> d.parent_id >> d.count >> d.birth_time >> d.coalescent >> d.w_saved >> d.g;

                if ( d.id == 0 )  {
					p.protein_length = (int)d.g.codonLength();
					continue;
				}
				
                v.push_back( d );

                if ( d.count > 1 )
                        break;
        }

        // read one more if possible
        if ( in >> d.id )
        {
                in >> d.parent_id >> d.count >> d.birth_time >> d.coalescent >> d.w_saved >> d.g;
                // this one gives us the coalescent time:
                p.coalescent_time = d.birth_time - 1;
        }
        else
                p.coalescent_time = d.birth_time;

        if ( p.coalescent_time < p.window_size + p.equilibration_time )
                p.coalescent_time = p.window_size + p.equilibration_time;
}

void printGenebank( CompactLatticeFolder &b, const Params p, const vector<GenebankData> &v )
{
        cout << "# free energy cutoff: " << p.free_energy_cutoff << endl;
        cout << "# per-site mutation rate: " << p.u << endl;
        cout << "# transl. robustness cost factor: " << p.tr_cost << endl;
        cout << "# codon adaptation cost factor: " << p.ca_cost << endl;
        cout << "# transl. error rate: " << p.transl_error_rate << endl;
        cout << "#\n#<time> <fitness (saved)> <fitness (recalc.)> <free energy> <structure id> <neutrality> <sensitivity> <sensitivity (no stop)> <fop>" << endl;

        vector<GenebankData>::const_iterator cit, ce;
        cit = v.begin();
        ce = v.end();

        while( cit != ce )
        {
                GenebankData d = *cit;
                FoldInfo fi = b.fold( d.g.translate() );
                d.free_energy = fi.getFreeEnergy();
                d.struct_id = fi.getStructure();
				
                ErrorproneTranslation *ept = new ErrorproneTranslation();
				ept->init( &b, p.protein_length, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost,
									p.transl_error_rate, p.transl_acc_wt, p.transl_error_wt );
									
                d.w_new = ept->getFitness( d.g );
                
				d.nu = GeneUtil::calcNeutrality( b, d.g.translate(), p.free_energy_cutoff );
                d.fop = GeneUtil::calcFop( d.g, ept->getOptimalCodons(false) );

                cout << d.birth_time << " " << d.w_saved << " " << d.w_new << " ";
                cout << d.free_energy << " " << d.struct_id << " " << d.nu << " ";
                cout << d.fop << endl;
				
                cit++;
                if ( cit == ce )
                        break;

                // output intermediate steps (good for plotting)
                /*
                                for ( int i=d.birth_time+1; i<(*cit).birth_time; i++ )
                                {
                                        cout << i << " " << d.w_saved << " " << d.w_new << " ";
                                        cout << d.free_energy << " " << d.struct_id << " " << d.nu << " ";
                                        cout << d.sensitivity << " " << d.fop << endl;
                                }
                */
        }

}

void analyzeSurfaceCore( CompactLatticeFolder &b, Params p, const vector<GenebankData> &v )
{
        // first, get structure
        //	pair<double, int> fp = GeneUtil::translateAndFold( b, v[0].g);
		FoldInfo fi = b.fold( v[0].g.translate() );

        //b.printStructure( fp.second );
        vector<int> surface = b.getSurface( fi.getStructure() );

        // now, do analysis
        int start = p.coalescent_time - p.window_size + 1;

        vector<GenebankData>::const_iterator cit, tmp, ce;
        cit = tmp = v.begin();
        ce = v.end();
        if ( cit == ce )
                return;

        while( 1 )
        {
                tmp++;
                if ( tmp == ce )
                        break;
                if ( (*tmp).birth_time < start )
                        cit = tmp;
                else
                        break;
        }

        cout << "# free energy cutoff: " << p.free_energy_cutoff << endl;
        cout << "# per-site mutation rate: " << p.u << endl;
        cout << "# transl. robustness cost factor: " << p.tr_cost << endl;
        cout << "# codon adaptation cost factor: " << p.ca_cost << endl;
        cout << "# window size: " << p.window_size << endl;
        cout << "# analysis start time: " << start << endl;
        cout << "# analysis end: " << p.coalescent_time << endl;

        GenebankData d = (*cit);
        cit++;
        double nn = 0, ns = 0, nnSurf = 0, nnCore = 0, nsSurf = 0, nsCore = 0;
        double dn = 0, ds = 0, dnSurf = 0, dnCore = 0, dsSurf = 0, dsCore = 0;
        double Fop = 0, FopSurf = 0, FopCore = 0;
        //cout << "#\n#<time> <nn> <ns> <nnSurf> <nnCore> <nsSurf> <nsCore> <dn> <ds> <dnSurf> <dnCore> <dsSurf> <dsCore> <sequence>" << endl;
        //cout << start << " 0 0 0 0 0 0 0 0 0 0 0 0 " << d.g << endl;
        while( cit != ce )
        {
                double n, s, nSurf, nCore, sSurf, sCore;
                double N, S, NSurf, NCore, SSurf, SCore;
                GeneUtil::calcDnDs( n, s, d.g, (*cit).g );
                GeneUtil::calcDnDsSurfaceCore( nSurf, nCore, sSurf, sCore, d.g, (*cit).g, surface );
                S = GeneUtil::calcSynonymousSites( d.g );
                N = GeneUtil::calcTotalSites( d.g ) - S;
                GeneUtil::calcSNSitesSurfaceCore( NSurf, NCore, SSurf, SCore, d.g, surface );
                nn += n;
                ns += s;
                nnSurf += nSurf;
                nnCore += nCore;
                nsSurf += sSurf;
                nsCore += sCore;
                dn += n/N;
                ds += s/S;
                dnSurf += nSurf/NSurf;
                dnCore += nCore/NCore;
                dsSurf += sSurf/SSurf;
                dsCore += sCore/SCore;
                //cout << d.birth_time << " " << n << " " << s << " " << nSurf << " " << nCore << " " << sSurf << " " << sCore <<  " " << n/N << " " << s/S << " " << nSurf/NSurf << " " << nCore/NCore << " " << sSurf/SSurf << " " << sCore/SCore << " " << (*cit).g;
                //cout << endl;
                d = (*cit);
                cit++;
                // calculate the fraction of optimal codons from the last sequence
                if ( cit == ce )
                {
						ErrorproneTranslation* ept = new ErrorproneTranslation();
						ept->init( &b, p.protein_length, p.structure_ID, p.free_energy_cutoff, 1, p.ca_cost,
									p.transl_error_rate, p.transl_acc_wt, p.transl_error_wt );
                        Fop = GeneUtil::calcFop( d.g, ept->getOptimalCodons(false) );
                        GeneUtil::calcFopSurfaceCore( FopSurf, FopCore, d.g, ErrorproneTranslation::m_codon_cost, surface );
                        break; // we're done
                }
        }

        cout << "# <nn> <ns> <nnSurf> <nnCore> <nsSurf> <nsCore> <dn> <ds> <dnSurf> <dnCore> <dsSurf> <dsCore> <Fop> <FopSurf> <FopCore>" << endl;
        cout << nn << " " << ns << " " << nnSurf << " " << nnCore << " " << nsSurf << " " << nsCore <<  " " << dn << " " << ds << " " << dnSurf << " " << dnCore << " " << dsSurf << " " << dsCore << " " << Fop << " " << FopSurf << " " << FopCore << endl;
}

/*
void analyzeMutations(Params p, const vector<GenebankData> &v )
{
        int start = p.coalescent_time - p.window_size + 1;

        vector<GenebankData>::const_iterator cit, tmp, ce;
        cit = tmp = v.begin();
        ce = v.end();
        if ( cit == ce )
                return;

        while( 1 )
        {
                tmp++;
                if ( tmp == ce )
                        break;
                if ( (*tmp).birth_time < start )
                        cit = tmp;
                else
                        break;
        }

        cout << "# free energy cutoff: " << p.cutoff << endl;
        cout << "# per-site mutation rate: " << p.u << endl;
        cout << "# transl. robustness cost factor: " << p.tr_cost << endl;
        cout << "# codon adaptation cost factor: " << p.ca_cost << endl;
        cout << "# window size: " << p.window_size << endl;
        cout << "#\n# analysis start time: " << start << endl;
        cout << "# analysis end: " << p.coalescent_time << endl;
        cout << "#\n#<time> <dn> <ds> <sequence>" << endl;



        GenebankData d = (*cit);
        cout << start << " 0 0 " << d.g << endl;
        cit++;
        while( cit != ce )
        {
                double dn, ds;
                GenotypeUtil::calcDnDs( dn, ds, d.g, (*cit).g );
                cout << (*cit).birth_time << " " << dn << " " << ds << " " << (*cit).g;
                if ( dn + ds > 1.01 )
                        cout << " multiple mutations";
                cout << endl;
                d = (*cit);
                cit++;
        }

}
*/

/* 
 * debugging function to print the parameters read from the Genebank file.
 * this should be useful if the format of the genebank files every change
*/
void printParams( Params p )
{
	
  cout << "evaluation type: " << p.eval_type << endl;
  cout << "structure id: " << p.structure_ID << endl;
  cout << "free energy maximum:" << p.free_energy_cutoff << endl;
  cout << "free energy minimum: " << 	p.free_energy_min << endl;
  cout << "protein length " << p.protein_length << endl;
  cout << "per-site mutation rate u: " << p.u << endl;
  cout << "transl. robust. cost factor: " << p.tr_cost << endl;
  cout << "codon adaptation cost factor: " << p.ca_cost << endl;
  cout << "transl. error rate per codon: " << p.transl_error_rate << endl;
  cout << "transl. accuracy weight: " << p.transl_acc_wt << endl;
  cout << "transl. error weight: " << p.transl_error_wt << endl;
  cout << "population size N: " << p.N << endl;
  cout << "window size tau: " << p.window_size << endl;
  cout << "equilibration time: " << p.equilibration_time << endl;
  cout << "coalescent time: " << p.coalescent_time << endl;
  
}

int main( int ac, char **av)
{

        if ( ac != 2 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <genebank file>" << endl;
                exit (-1);
        }

        ifstream in( av[1] );

        if ( !in )
        {
                cerr << "Error opening " << av[1] << endl;
                exit(-1);
        }

        // size of the lattice proteins is hardcoded
        const int size = 5;

        // initialize the protein folder
        CompactLatticeFolder b(size);

        Params p;
        vector<GenebankData> v;

        readGenebank( in, p, v );
	// printParams( p );		
	printGenebank( b, p, v );
	// analyzeMutations( p, v );
	analyzeSurfaceCore( b, p, v );

}




