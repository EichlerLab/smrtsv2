/**
 * Error reporting utilities.
 */

#include <cstdarg>
#include <iostream>

// Constants
const int ERR_BUF_SIZE = 2048;

extern char *progName;

using namespace std;

void err(const char *message, ...) {

	char *msgBuf = new char[ERR_BUF_SIZE];

	// Get message
	va_list args;
	va_start(args, message);

	vsprintf(msgBuf, message, args);

	va_end(args);

	// Print error
	cerr << progName << ": Error: " << string(msgBuf) << endl;

	// Clean up
	delete(msgBuf);

	return;
}

void err(const std::exception &ex) {

	// Print error
	cerr << progName << ": Error: " << ex.what() << endl;

	return;
}

void warn(const char *message, ...) {

	char *msgBuf = new char[ERR_BUF_SIZE];

	// Get message
	va_list args;
	va_start(args, message);

	vsprintf(msgBuf, message, args);

	va_end(args);

	// Print error
	cerr << progName << ": Warning: " << string(msgBuf) << endl;

	// Clean up
	delete(msgBuf);

	return;
}

void warn(const std::exception &ex) {

	// Print error
	cerr << progName << ": Warning: " << ex.what() << endl;

	return;
}
