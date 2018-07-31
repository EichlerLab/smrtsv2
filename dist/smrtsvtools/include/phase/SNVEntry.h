#ifndef SNVENTRY_H_
#define SNVENTRY_H_

#include <string>

using namespace std;

namespace phase {

class SNVEntry {
public:

	// CTOR/DTOR
	SNVEntry(int pos, int ps, char allele0, char allele1) :
		pos(pos), ps(ps), allele0(allele0), allele1(allele1) {};

	// Data members
	const int pos;
	const int ps;
	const char allele0;
	const char allele1;
};

} /* namespace phase */

#endif /* SNVENTRY_H_ */
