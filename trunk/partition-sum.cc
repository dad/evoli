#include "protein-folder.hh"
#include "translator.hh"
#include "fitness-evaluator.hh"
#include "tools.hh"
#include "genotype-util.hh"

#include <fstream>
#include <string>
#include <cstdio>


struct Parameters
{
        double free_energy_cutoff;
        double ca_cost;
        double error_rate;
        double u;
        double min_free_energy;
        double sampling_bias;
        int samples;
        int random_seed;
        int window_size;
        string call;
        Genotype g;
};

struct PSumEntry
{
        double sensitivity;
        double free_energy;
        double fop;
        double fitness;
        vector<double> syn_mut_sensitivities;
        vector<double> nonsyn_mut_sensitivities;
};

ostream & operator<<( ostream &s, const Parameters &p )
{
        s << "# Parameters:" << endl;
        s << "#   Free energy cutoff: " << p.free_energy_cutoff << endl;
        s << "#   codon adaptation cost: " << p.ca_cost << endl;
        s << "#   transl. error rate: " << p.error_rate << endl;
        s << "#   mutation rate: " << p.u << endl;
        s << "#   window size: " << p.window_size << endl;
        s << "#   number of samples: " << p.samples << endl;
        s << "#   random seed: " << p.random_seed << endl;
        s << "#   initial sequence: " << p.g << endl;
        s << "#   program call: " << p.call << endl;
        s << "#" << endl;
        return s;
}


Parameters getParams( int ac, char **av )
{
        if ( ac != 9 )
        {
                cout << "Start program like this:" << endl;
                cout << "  " << av[0] << " <ca cost> <transl. error rate> <free energy cutoff> <per-site mutation rate> <window size> <samples> <random seed> <sequence file>" << endl;
                exit (-1);
        }

        Parameters p;
        int i = 1;
        p.ca_cost = atof( av[i++] );
        p.error_rate = atof( av[i++] );
        p.free_energy_cutoff = atof( av[i++]);
        p.u = atof( av[i++] );
        p.window_size = atoi( av[i++] );
        p.samples = atoi( av[i++] );
        p.random_seed = atoi( av[i++] );

        p.min_free_energy = -7.5;
        p.sampling_bias = 4.2;

        // read initial sequence
        ifstream in( av[i] );
        in >> p.g;

        // record programm call
        for ( int i=0; i<ac; i++ )
                p.call = p.call + av[i] + " ";
        return p;
}

void calcSinglePointMutants( ProteinFolder &b, Genotype g, ErrorproneTranslation *fe, PSumEntry *e )
{
        e->syn_mut_sensitivities.clear();
        e->nonsyn_mut_sensitivities.clear();

        GenotypeIterator gi = g.begin(), ge = g.end();
        for ( ; gi != ge; gi++ )
        {
                int c_orig = ( *gi );
                CodonPointMutationIterator ci( c_orig );
                do
                {
                        ( *gi ) = ci.codon();
                        // exclude stop codons
                        if ( GeneticCodeUtil::geneticCode[ *gi ] >= 0 )
                        {
                                double dn, ds;
                                GeneticCodeUtil::calcDnDs( dn, ds, ( *gi ), c_orig );
                                double w = fe->getFitness( g );
                                if ( w > 0 )
                                {
                                        if ( dn > 0 )
                                                e->nonsyn_mut_sensitivities.push_back( fe->getLastSensitivity() );
                                        else
                                                e->syn_mut_sensitivities.push_back( fe->getLastSensitivity() );
                                }

                        }
                        //cout << g << " " << nu_syn << " " << nu_nonsyn << " " << count << endl;

                }
                while( ci.next() );
                ( *gi ) = c_orig; // set back to original codon
        }
        //cout << "Nu: " << GenotypeUtil::calcNeutrality( b, g, free_energy_cutoff ) << " Nu nonsyn: " << nu_nonsyn << " Nu syn: " << nu_syn << endl;
}

double calcSamplingProb( const Parameters& p, double free_energy )
{
        //return 1.;
        double prob = exp( p.sampling_bias * (p.min_free_energy - free_energy) );
        if ( prob > 1 )
                prob = 1;
        else if ( prob < 1e-4 )
                prob = 1e-4;
        return prob;
}

void recordEntry( ProteinFolder &b, const Parameters &p, vector<PSumEntry> *entries, ErrorproneTranslation *fe, const Genotype &g )
{
        PSumEntry e;
        e.fitness = fe->getFitness( g );
        e.free_energy = fe->getLastFreeEnergy();
        e.sensitivity = fe->getLastSensitivity();
        e.fop = GenotypeUtil::calcFop( g, ErrorproneTranslation::m_codon_cost );
        calcSinglePointMutants( b, g, fe, &e );
        entries->push_back( e );
}

vector<PSumEntry> getPSumEntries( ProteinFolder &b, const Parameters &p )
{
        Genotype g = p.g, g2;
        vector<PSumEntry> entries;


        // check if initial sequence is sane
        pair<double, int> fdata = GenotypeUtil::translateAndFold( b, p.g );
        int struct_id = fdata.second;
        if ( struct_id < 0 )
        {
                cerr << "Initial sequence does not translate. Nothing to be done." << endl;
                return entries;
        }
        // set up fitness evaluator.
        // tr_cost is irrelevant at this point, set it to .01
        ErrorproneTranslation fe( &b, struct_id, p.free_energy_cutoff, .01, p.ca_cost, p.error_rate );

        // record the first entry
        recordEntry( b, p, &entries, &fe, g );

        // now do all the remaining ones
        int count = 1;
        do
        {
                g2 = g;
                GenotypeUtil::pointMutation( g2 ); // do a single point mutation
                fdata = GenotypeUtil::translateAndFold( b, g2 ); // try to fold the sequence
                if ( fdata.second == struct_id && fdata.first < p.free_energy_cutoff
                                && myRand() < calcSamplingProb( p, fdata.first ) )
                {
                        g = g2;
                        recordEntry( b, p, &entries, &fe, g );
                        count += 1;
                        if ( count % 100 == 0 )
                                cout << "[" << count << "/" << p.samples <<"] " << flush;
                }
        }
        while( count < p.samples );
        return entries;
}

double calcFitness( double sensitivity, double tr_cost, double transl_error_rate )
{
        double f = sensitivity * transl_error_rate;
        return exp( - tr_cost * f / ( 1 - f ) );
}

void updateFitnesses( vector<PSumEntry> &entries, double tr_cost, double transl_error_rate )
{
        vector<PSumEntry>::iterator it = entries.begin(), e = entries.end();
        for ( ; it != e; it++ )
        {
                (*it).fitness = calcFitness( (*it).sensitivity, tr_cost, transl_error_rate );
        }
}

double calcFixProb( double fitness_or, double fitness_tar, int N )
{
        double dx = log( fitness_tar/fitness_or );
        double fixprob;
        if ( dx > 1e-10 )
                fixprob = (1-exp(-2*dx))/(1-exp(-2*N*dx));
        else if ( dx < -1e-10 )
                fixprob = exp((2*N-2)*dx)*(1-exp(2*dx))/(1-exp(2*N*dx));
        else
                fixprob = 1./N;
        //cout << dx << " " << fitness_or << " " << fitness_tar << " " << fixprob << endl;
        return fixprob;
}

double calcSynTransitionProb( const PSumEntry &e, double tr_cost, double transl_error_rate, int N )
{
        vector<double>::const_iterator it = e.syn_mut_sensitivities.begin(), end = e.syn_mut_sensitivities.end();
        double result = 0;

        for ( ; it != end; ++it )
        {
                result += calcFixProb( e.fitness, calcFitness( (*it), tr_cost, transl_error_rate ), N );
        }
        return result;
}

double calcNonsynTransitionProb( const PSumEntry &e, double tr_cost, double transl_error_rate, int N )
{
        vector<double>::const_iterator it = e.nonsyn_mut_sensitivities.begin(),
                                            end = e.nonsyn_mut_sensitivities.end();
        double result = 0;

        for ( ; it != end; ++it )
        {
                result += calcFixProb( e.fitness, calcFitness( (*it), tr_cost, transl_error_rate ), N );
        }
        return result;
}

double calcPContrib( const PSumEntry &e, double beta, double scale, const Parameters &p )
{
        double x = log( e.fitness );
        return exp( beta * x - scale ) / calcSamplingProb( p, e.free_energy );
}

double findMinSens( const vector<PSumEntry> &entries )
{
        vector<PSumEntry>::const_iterator it = entries.begin(), e = entries.end();
        double min_sens = (*it).sensitivity;
        it++;
        for ( ; it != e; it++ )
        {
                if ( (*it).sensitivity < min_sens )
                        min_sens = (*it).sensitivity;
        }
        return min_sens;
}

void doAnalysis( ProteinFolder &b, const Parameters &p )
{
        vector<PSumEntry> entries = getPSumEntries( b, p );
        double sens_min = findMinSens( entries );
        double f_min = sens_min * p.error_rate;

        int Narray[5] = {100, 300, 1000, 3000, 10000};
        for ( int i=0; i<5; i++ )
        {
                int N = Narray[i];
                double beta = 2*N - 2;
                double evol_rate_scale = p.u * p.window_size * N / 3.; // have to divide by 3 because there are 3 possible mutations at each site

                // setup output file
                char filename[255];
                sprintf( filename, "N%ica%gGmax%g-pred.dat", N, p.ca_cost, p.free_energy_cutoff );
                ofstream data_file( filename, ios::out );
                data_file << p;
                data_file << "# <expr. level> <dn> <ds> <fitness> <fop>" << endl;

                for ( double tr_cost = .001; tr_cost < 10000; tr_cost *=1.2 )
                {
                        updateFitnesses( entries, tr_cost, p.error_rate );
                        double scale = beta * ( - tr_cost * f_min / ( 1 - f_min ) );
                        double psum = 0;
                        double fave = 0;
                        double fopave = 0;
                        double dnave = 0;
                        double dsave = 0;
                        vector<PSumEntry>::const_iterator it = entries.begin(), e = entries.end();
                        for ( ; it != e; it++ )
                        {
                                double contrib = calcPContrib( (*it), beta, scale, p );
                                fave += (*it).fitness * contrib;
                                fopave += (*it).fop * contrib;
                                dnave += calcNonsynTransitionProb( (*it), tr_cost, p.error_rate, N ) * evol_rate_scale * contrib;
                                dsave += calcSynTransitionProb( (*it), tr_cost, p.error_rate, N ) * evol_rate_scale * contrib;
                                psum += contrib;
                        }
                        data_file << tr_cost << " " << dnave/psum << " " << dsave/psum << " " << fave/psum << " " << fopave/psum << endl;
                }
        }


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

        doAnalysis( b, p );
}





