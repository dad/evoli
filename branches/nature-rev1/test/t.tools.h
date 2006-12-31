/*
This file is part of the evoli project.
Copyright (C) 2004, 2005, 2006 Claus Wilke <cwilke@mail.utexas.edu>,
Allan Drummond <dadrummond@gmail.com>

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


#ifndef _T_TOOLS_H__
#define _T_TOOLS_H__
#include "cutee.h"
#include "tools.hh"
#include "random.hh"
#include <fstream>
#include <cstdlib>

#define MAXINT 0x7FFFFFFF

struct TEST_CLASS( tools_basic )
{
	void TEST_FUNCTION( itoa_decimal ) {
		srand48(0);
		for (int i=0; i<1000; i++) {
			int target = Random::rint( MAXINT/2 );
			if ( Random::runif() < 0.5)
				target = -target;
			string s = itoa(target, 10);
			int res = atoi(s.c_str());
			//cout << target << " " << s << " " << res << endl;
			TEST_ASSERT(res==target);
		}
		return;
	}

};


#endif // _T_TOOLS_H__
