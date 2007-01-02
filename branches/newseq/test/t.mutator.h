/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <drummond@alumni.princeton.edu>

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


#ifndef _T_MUTATOR_H__
#define _T_MUTATOR_H__
#include "cutee.h"
#include "decoy-contact-folder.hh"
#include "compact-lattice-folder.hh"
#include "fitness-evaluator.hh"
#include "protein.hh"
#include "gene-util.hh"
#include "tools.hh"
#include <fstream>
#include <cmath>

struct TEST_CLASS( mutator_basic )
{
	void TEST_FUNCTION( is_mutated ) {
		CodingDNA g = CodingDNA::createRandomNoStops(75);
		SimpleMutator mut(0.01);
		for (int i=0; i<100; i++) {
			CodingDNA g2 = g;
			bool changed = mut.mutate(g2);
			TEST_ASSERT( changed == (g2 != g) );
		}
	}
};

#endif // _T_MUTATOR_H__
