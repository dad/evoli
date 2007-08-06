#include "sequence.hh"

ostream& operator<<(ostream& os, const Contact& c) {
	os << "(" << c.first << "," << c.second << ")";
	return os;
}
