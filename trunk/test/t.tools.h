#ifndef _T_TOOLS_H__
#define _T_TOOLS_H__
#include "cutee.h"
#include "tools.hh"
#include <fstream>
#include <cstdlib>

struct TEST_CLASS( tools_basic )
{
	void TEST_FUNCTION( itoa_decimal ) {
		srand48(0);
		for (int i=0; i<1000; i++) {
			TEST_ASSERT(atoi(itoa(i,10).c_str())==i);
		}
		return;
	}

};


#endif // _T_TOOLS_H__
