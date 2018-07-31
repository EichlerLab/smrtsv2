/*
 * PhaseTable.h
 *
 *  Created on: Aug 14, 2017
 *      Author: paudano
 */

#ifndef PHASETABLE_H_
#define PHASETABLE_H_

#include <unordered_map>
#include <string>
#include <vector>

#include "phase/SNVEntry.h"

using namespace std;

namespace phase {

class PhaseTable {
public:

	// CTOR/DTOR
	PhaseTable();
	virtual ~PhaseTable();

	// Add a single entry
	void add(const string &chr, const int pos, const int ps, const char allele0, const char allele1);

	// Load entries from a file
	void load(const string &snvTableFileName);

	// Get a list of SNV entries
	const vector<SNVEntry> *getEntry(const char *templateName) const;

	// Get number of entries
	unsigned int size() const;

	// Clear all entries
	void clear();

private:
	unordered_map<string, int> chromHash;  /** Maps chromosome names to an integer. */

	vector<vector<SNVEntry>*> snvVector;  /** Vector SNVs (inner vector) for each chromosome (outer vector). */

	vector<int> lastPos;  /** Last position in snvVector where an entry was added. */

	string lastChr;    /** The last added SNV was an this chromosome (init as empty string). */
	int lastChrIndex;  /** Index of the last chromosome. */
};

} /* namespace phase */

#endif /* PHASETABLE_H_ */
