#include "mutator.hh"

Mutator::Mutator( ProteinFolder &folder, double cutoff ) : m_folder( folder ), m_cutoff( cutoff )
{
        m_L = m_folder.getProteinLength();
        m_sequence = new int[m_L];

}

Mutator::~Mutator()
{
        delete m_sequence;
}

