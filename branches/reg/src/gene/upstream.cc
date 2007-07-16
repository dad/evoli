#include "upstream.hh"

ostream& operator<<(ostream& os, const BindingInteraction& bi) {
	bi.print(os);
	return os;
}
