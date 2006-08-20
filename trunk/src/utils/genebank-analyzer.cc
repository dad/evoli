#include "protein-folder.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "genotype-util.hh"

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
        Genotype g;
};

struct Params
{
        double cutoff;
        double tr_cost;
        double ca_cost;
        double transl_error_rate;
        double u;
        int window_size;
        int equilibration_time;
        int coalescent_time;
        int N;
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
                        // get free energy cutoff
                        if ( s == "free"  )
                        {
                                in >> s;
                                in >> s;
                                in >> p.cutoff;
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
                                else
                                {
                                        in >> s;
                                        in >> p.transl_error_rate;
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

                if ( d.id == 0 )
                        continue;

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

void printGenebank( ProteinFolder &b, const Params p, const vector<GenebankData> &v )
{
        cout << "# free energy cutoff: " << p.cutoff << endl;
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
                pair<double, int> fp = GenotypeUtil::translateAndFold( b, d.g );
                d.free_energy = fp.first;
                d.struct_id = fp.second;
                ErrorproneTranslation fe( &b, d.struct_id, p.cutoff, p.tr_cost, p.ca_cost, p.transl_error_rate );
                d.w_new = fe.getFitness( d.g );
                d.sensitivity = fe.getLastSensitivity();
                d.sensitivity_ns = fe.getLastSensitivityNoStop();
                d.nu = GenotypeUtil::calcNeutrality( b, d.g, p.cutoff );
                d.fop = GenotypeUtil::calcFop( d.g, ErrorproneTranslation::m_codon_cost );

                cout << d.birth_time << " " << d.w_saved << " " << d.w_new << " ";
                cout << d.free_energy << " " << d.struct_id << " " << d.nu << " ";
                cout << d.sensitivity << " " << d.sensitivity_ns << " " << d.fop << endl;
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

void analyzeSurfaceCore( ProteinFolder &b, Params p, const vector<GenebankData> &v )
{
        // first, get structure
        pair<double, int> fp = GenotypeUtil::translateAndFold( b, v[0].g );

        //b.printStructure( fp.second );
        vector<int> surface = b.getSurface( fp.second );

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

        cout << "# free energy cutoff: " << p.cutoff << endl;
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
                GenotypeUtil::calcDnDs( n, s, d.g, (*cit).g );
                GenotypeUtil::calcDnDsSurfaceCore( nSurf, nCore, sSurf, sCore, d.g, (*cit).g, surface );
                S = GenotypeUtil::calcSynonymousSites( d.g );
                N = GenotypeUtil::calcTotalSites( d.g ) - S;
                GenotypeUtil::calcSNSitesSurfaceCore( NSurf, NCore, SSurf, SCore, d.g, surface );
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
                        Fop = GenotypeUtil::calcFop( d.g, ErrorproneTranslation::m_codon_cost );
                        GenotypeUtil::calcFopSurfaceCore( FopSurf, FopCore, d.g, ErrorproneTranslation::m_codon_cost, surface );
                        break; // we're done
                }
        }

        cout << "# <nn> <ns> <nnSurf> <nnCore> <nsSurf> <nsCore> <dn> <ds> <dnSurf> <dnCore> <dsSurf> <dsCore> <Fop> <FopSurf> <FopCore>" << endl;
        cout << nn << " " << ns << " " << nnSurf << " " << nnCore << " " << nsSurf << " " << nsCore <<  " " << dn << " " << ds << " " << dnSurf << " " << dnCore << " " << dsSurf << " " << dsCore << " " << Fop << " " << FopSurf << " " << FopCore << endl;
}

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
        ProteinFolder b(size);

        Params p;
        vector<GenebankData> v;

        readGenebank( in, p, v );
  //      printGenebank( b, p, v );
   //   analyzeMutations( p, v );
        analyzeSurfaceCore( b, p, v );
}




