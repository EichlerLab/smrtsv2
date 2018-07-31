#include <iostream>
#include <fstream>

#include "phase/PhaseTable.h"
#include "phase/StringUtil.h"

#include "boost/algorithm/string.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string/trim.hpp"
#include "boost/lexical_cast.hpp"

using namespace std;

namespace phase {

PhaseTable::PhaseTable() {
	
	lastChr = "";  // No chromosome added
	lastChrIndex = -1;  // Dummy value
	
	return;
}

void PhaseTable::add(const string &chr, const int pos, const int ps, const char allele0, const char allele1) {

	int chrIndex;
	unordered_map<string, int>::const_iterator mapIter;

	bool newLastChr;  // Set to true if chr != lastChr

	if (chr.empty())
		throw invalid_argument("Empty chromosome name");

	// Get index
	if (chr.compare(lastChr) == 0) {
		// Last chromosome added
		chrIndex = lastChrIndex;
		newLastChr = false;

	} else if ((mapIter = chromHash.find(chr)) != chromHash.end()) {
		// Not last chromosome, but chromosome exists

		chrIndex = chromHash[chr];
		newLastChr = true;

	} else {
		// New chromosome
		chrIndex = snvVector.size();
		snvVector.push_back(new vector<SNVEntry>);

		chromHash[chr] = chrIndex;

		newLastChr = true;
	}

	// Add entry
	snvVector[chrIndex]->push_back(SNVEntry(pos, ps, allele0, allele1));

	// Update last chr index state
	if (newLastChr) {
		lastChrIndex = chrIndex;
		lastChr = chr;
	}

	return;
}

void PhaseTable::load(const string &snvTableFileName) {

	ifstream inFile(snvTableFileName);  // Input file
	string line;  // Line read from SNV table
	unsigned int lineCount;  // Line count

	vector<string> lineTok;  // line tokenized by tab
	vector<string>::iterator tokIter;  // Iterator over lineTok

	string chr;
	int pos;
	int ps;
	char allele0;
	char allele1;

	// Check file
	if (! inFile.is_open()) {
		string errString("Error opening SNV table: ");
		errString.append(snvTableFileName);

		throw runtime_error(errString);
	}

	// Read each line
	lineCount = 0;

	while(getline(inFile, line)) {
		lineCount += 1;

		boost::trim(line);  // Strip leading/trailing whitespace

		// Skip empty lines and comments
		if (line.size() == 0 || line[0] == '#')
			continue;

		// Tokenize BED record on tabs
		boost::split(lineTok, line, boost::is_any_of("\t"));

		if (lineTok.size() != 10)
			throw_error<runtime_error>("Expected 4 fields on line %d: Received %d", lineCount, lineTok.size());

		// Skip non-SNV entries
		if (lineTok[7] == "." || lineTok[8] == "." || lineTok[9] == ".")
			continue;

		// Get fields
		chr = lineTok[0];

		try {
			pos = boost::lexical_cast<int>(lineTok[1]);

		} catch (boost::bad_lexical_cast &ex) {
			throw_error<runtime_error>("POS field (2) on line %d is not an integer: %s", lineCount, lineTok[1].c_str());
		}

		try {
			ps = boost::lexical_cast<int>(lineTok[7]);

		} catch (boost::bad_lexical_cast &ex) {
			throw_error<runtime_error>("PS (phase-set) field (8) on line %d is not an integer: %s", lineCount, lineTok[7].c_str());
		}

		if (! isNucleotide(lineTok[8]))
			throw_error<runtime_error>("ALLELE_0 field (9) on line %1% is not an valid nucleotide: %2%", lineCount, lineTok[8].c_str());

		if (! isNucleotide(lineTok[9]))
			throw_error<runtime_error>("ALLELE_1 field (10) on line %1% is not an valid nucleotide: %2%", lineCount, lineTok[9].c_str());

		allele0 = lineTok[8][0];
		allele1 = lineTok[9][0];

		// Add allele
		add(chr, pos, ps, allele0, allele1);
	}

	return;
}

const vector<SNVEntry> *PhaseTable::getEntry(const char *templateName) const {
	return snvVector[chromHash.at(templateName)];  // throws std::out_of_range
}

unsigned int PhaseTable::size() const {

	unsigned int size = 0;

	for (unsigned int index = 0; index < snvVector.size(); ++index)
		size += snvVector[index]->size();

	return size;
}

void PhaseTable::clear() {

	for (unsigned int index = 0; index < snvVector.size(); ++index)
		delete(snvVector[index]);

	chromHash.clear();
	snvVector.clear();
	lastPos.clear();

	lastChr = "";  // No chromosome added
	lastChrIndex = -1;  // Dummy value

	return;
}

PhaseTable::~PhaseTable() {

	// Free SNV Vector entries.
	for (unsigned int index = 0; index < snvVector.size(); ++index)
		delete(snvVector[index]);

	return;
}

} /* namespace phase */
