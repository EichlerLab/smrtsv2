#include <string>

#include "phase/StringUtil.h"

namespace phase {

/**
 * Determine if a string is a single-letter nucleotide in all caps.
 *
 * @param s String to check.
 *
 * @return TRUE if s is a nucleotide character, and FALSE otherwise.
 */
bool isNucleotide(const string s) {
	return find(NUCLEOTIDE_LIST.begin(), NUCLEOTIDE_LIST.end(), s) != NUCLEOTIDE_LIST.end();
}

} /* namespace phase */
