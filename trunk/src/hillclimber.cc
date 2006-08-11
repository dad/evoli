#include "protein-folder.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <fstream>


struct Parameters
{
        double free_energy_cutoff;
        double tr_cost;
        double ca_cost;
        double error_rate;
        int max_steps;
        int random_seed;
        Genotype g;
};


ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   Free energy cutoff: " << p.free_energy_cutoff << endl;
        s << "#   transl. robustness cost: " << p.tr_cost << endl;
        s << "#   codon adaptation cost: " << p.ca_cost << endl;
        s << "#   error rate: " << p.error_rate << endl;
        s << "#   max. steps: " << p.max_steps << endl;
        s << "#   random seed: " << p.random_seed << endl;
        s << "#   initial sequence: " << p.g << endl;
        s << "#" << endl;
        return s;
}


Parameters getParams( int ac, char **av )
{
        if ( ac != 8 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <free energy cutoff> <tr cost> <ca cost> <error rate> <max. steps> <random seed> <sequence file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.free_energy_cutoff = atof( av[i++]);
        p.tr_cost = atof( av[i++] );
        p.ca_cost = atof( av[i++] );
        p.error_rate = atof( av[i++] );
        p.max_steps = atoi( av[i++] );
        p.random_seed = atoi( av[i++] );


        // read initial sequence
        ifstream in( av[i] );
        in >> p.g;

        return p;
}


void hillclimb( ProteinFolder &b, const Parameters &p, ostream &s )
{
        // find a random sequence with folding energy smaller than cutoff
        Genotype g, g2;
        pair<double, int> fdata;

        g = p.g;
        fdata = GenotypeUtil::translateAndFold( b, p.g );

        if ( fdata.second < 0 )
        {
                s << "Initial sequence does not translate. Nothing to be done." << endl;
                return;
        }


        double G = fdata.first;
        int id = fdata.second;
        ErrorproneTranslation fe( &b, id, p.free_energy_cutoff, p.tr_cost, p.ca_cost, p.error_rate );
        double w = fe.getFitness( g );
        double sens = fe.getLastSensitivity();
        double sens2 = fe.getLastSensitivityNoStop();
        double nu = GenotypeUtil::calcNeutrality( b, g, p.free_energy_cutoff );
        double fop = GenotypeUtil::calcFop( g, ErrorproneTranslation::m_codon_cost );

        int count = 1;
        s << "# <count> <free energy> <neutrality> <fitness> <sensitivity> <sensitivity (no stop)> <fop> <sequence>" << endl;
        s << count << " " << G << " " << nu << " " << w << " " << sens << " " << sens2 << " " << fop << " " << g << endl;

        double fop_ave = 0;
        do
        {
                fop_ave += fop;
                g2 = g;
                bool changed = GenotypeUtil::mutateGenotype( g2, 0.01 );
                double w2 = fe.getFitness( g2 );
                if ( changed && w2 >= w - 0.000001 )
                        //cout << w2 << " " << w << " " << w2 - w << endl;
                        //    if ( w2 - w >= -.1 && w2 > 0 )
                {

                        g = g2;
                        w = w2;
                        fdata = GenotypeUtil::translateAndFold( b, g2 );
                        G = fdata.first;
                        id = fdata.second;
                        sens = fe.getLastSensitivity();
                        sens2 = fe.getLastSensitivityNoStop();
                        nu = GenotypeUtil::calcNeutrality( b, g, p.free_energy_cutoff );
                        fop = GenotypeUtil::calcFop( g, ErrorproneTranslation::m_codon_cost );
                        s << count << " " << G << " " << nu << " " << w << " " << sens << " " << sens2 << " " << fop << " " << g << endl;
                }

                count += 1;
        }
        while( count < p.max_steps );

        cout << "ave fop: " << fop_ave / (double) p.max_steps << endl;
}

void hillclimb2( ProteinFolder &b, const Parameters &p, ostream &s )
{
        // find a random sequence with folding energy smaller than cutoff
        Genotype g, g2;
        pair<double, int> fdata;

        g = p.g;
        fdata = GenotypeUtil::translateAndFold( b, p.g );

        if ( fdata.second < 0 )
        {
                s << "Initial sequence does not translate. Nothing to be done." << endl;
                return;
        }


        int id = fdata.second;

        int count = 0;

        do
        {
                g2 = g;
                bool changed = GenotypeUtil::mutateGenotype( g2, 0.01 );
                fdata = GenotypeUtil::translateAndFold( b, g2 );
                if ( changed && fdata.second==id && fdata.first < p.free_energy_cutoff )
                {
                        g = g2;
                        cout << -1*fdata.first << endl;
                        count += 1;
                }
        }
        while( count < p.max_steps );
}


int main( int ac, char **av)
{
        Parameters p = getParams( ac, av );

        // size of the lattice proteins is hardcoded
        const int size = 5;

        // initialize the protein folder
        ProteinFolder b(size);
        b.enumerateStructures();

        // set random seed
        srand48( p.random_seed );

        hillclimb2( b, p, cout );
}




