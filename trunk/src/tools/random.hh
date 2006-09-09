#ifndef RANDOM_HH
#define RANDOM_HH

typedef unsigned int uint;

/** \brief Static class encapsulating the random number generator.

The random number implementation is the "Mersenne Twister" random number generator MT19937, which generates pseudorandom integers uniformly distributed in 0..(2^32 - 1) starting from any odd seed in 
0..(2^32 - 1).  This version is a recode by Shawn Cokus
(Cokus@math.washington.edu) on March 8, 1998 of a version by Takuji
Nishimura (who had suggestions from Topher Cooper and Marc Rieffel in
July-August 1997).

Effectiveness of the recoding (on Goedel2.math.washington.edu, a DEC Alpha
running OSF/1) using GCC -O3 as a compiler: before recoding: 51.6 sec. to
generate 300 million random numbers; after recoding: 24.0 sec. for the same
(i.e., 46.5% of original time), so speed is now about 12.5 million random
number generations per second on this machine.

According to the URL http://www.math.keio.ac.jp/~matumoto/emt.html
(and paraphrasing a bit in places), the Mersenne Twister is "designed
with consideration of the flaws of various existing generators," has
a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
equidistributed, and "has passed many stringent tests, including the
die-hard test of G. Marsaglia and the load test of P. Hellekalek and
S. Wegenkittl."  It is efficient in memory usage (typically using 2506
to 5012 bytes of static data, depending on data type sizes, and the code
is quite short as well).  It generates random numbers in batches of 624
at a time, so the caching and pipelining of modern systems is exploited.
It is also divide- and mod-free.

The code is at least 20-30% faster than the standard drand48()
random number generator supplied with the GNU libc. Calculating the sum
of 10^7 random doubles on a Pentium 4 takes on the order of 0.6 seconds
for this random number generator, and on the order of 0.9 seconds for
drand48().
*/

class Random
{
public:
	/**
	* @return A random 32 bit unsigned integer.
	*/
	static uint rint();

	/**
	* @return A random double chosen from uniform distribution on the interval [0, 1).
	*/
	static double runif();

	/**
	* This function returns random normal variates. It is based on the method by
	* J. H. Ahrens and U. Dieter, Extensions of Forsythe's method for random sampling
	* from the normal distribution. Math. Comput. 27,124 (OCT. 1973), 927 - 937.
	*
	* @param mean Mean of the normal distribution, default is 0.
	* @param sd Standard deviation of the normal distribution, default is 1.
	* @return A random double chosen from the normal distribution with given mean and standard deviation.
	*/
	static double rnorm( double mean = 0, double sd = 1 );

	/**
	* Seeds the random number generator. Any 32-bit unsigned integer
	* will work, but the best seeds are odd numbers in the interval 
	* 0..(2^32 - 1).
	*
	* @param i Unsigned integer to use as seed.
	*/
	static void seed( uint i );
};



#endif
