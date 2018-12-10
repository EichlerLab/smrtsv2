#ifndef INCLUDE_PHASE_STRINGUTIL_H_
#define INCLUDE_PHASE_STRINGUTIL_H_

#include <vector>
#include <string>
#include <algorithm>
#include <cstdarg>
#include <stdexcept>

#include "boost/static_assert.hpp"

using namespace std;

namespace phase {

const vector<string> NUCLEOTIDE_LIST{"A", "C", "G", "T"};

/**
 * Determine if a string is a single-letter nucleotide in all caps.
 *
 * @param s String to check.
 *
 * @return TRUE if s is a nucleotide character, and FALSE otherwise.
 */
bool isNucleotide(const string s);

/**
 * Throw a runtime error.
 *
 * @param message Message with formatting characters.
 * @parma ... Format arguments.
 */
template<typename T>
void throw_error(const char *message, ...) {

	BOOST_STATIC_ASSERT_MSG(is_base_of<runtime_error, T>::value, "Must throw exceptions of type runtime_error");

	char *msgBuf = new char[1024];

	// Get message
	va_list args;
	va_start(args, message);

	vsprintf(msgBuf, message, args);

	va_end(args);

	// Generate exception
	T rtErr(msgBuf);

	delete(msgBuf);

	throw rtErr;

	return;
}

} /* namespace phase */

#endif /* INCLUDE_PHASE_STRINGUTIL_H_ */
