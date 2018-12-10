#ifndef SRC_PHASE_READPHASER_H_
#define SRC_PHASE_READPHASER_H_

#include <iostream>

#include "phase/PhaseTable.h"

using namespace std;

namespace phase {

/**
 * Tags long reads in an alignment with phase information.
 */
class ReadPhaser {
public:

	ReadPhaser(const PhaseTable *phaseTable) :
		phaseTable(phaseTable), verbose(false) {
			return;
	}

	void phase(const vector<string> *alignFileList, const string outFileName);

	void setVerbose(bool verbose) {
		this->verbose = verbose;

		return;
	}

private:
	bool verbose;
	const PhaseTable *phaseTable;  // Table of phase information (SNVs)
};

} /* namespace phase */

#endif /* SRC_PHASE_READPHASER_H_ */
