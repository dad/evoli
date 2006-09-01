#ifndef _T_TOOLS_H__
#define _T_TOOLS_H__
#include "cutee.h"
#include "tools.hh"
#include <fstream>
#include <cstdlib>

#define MAXINT 0x7FFFFFFF

struct TEST_CLASS( tools_basic )
{
	void TEST_FUNCTION( itoa_decimal ) {
		srand48(0);
		for (int i=0; i<1000; i++) {
			int target = (int)((MAXINT/2.0)*drand48());
			if (drand48() < 0.5)
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
