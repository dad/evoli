#include "expressible-gene.hh"

ostream& operator<<(ostream& os, const ExpressibleGene& g) {
	g.print(os);
	return os;
}
