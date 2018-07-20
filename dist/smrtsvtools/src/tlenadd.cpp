/**
 * Calculate template length (tlen) for each read. This field is not set by newer
 * versions of PacBio BLASR, but some tools rely on it.
 */

#include <iostream>

#include "util/err.h"
#include "util/constants.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "htslib/sam.h"

using namespace std;

namespace po = boost::program_options;

const uint32_t CIGAR_MATCH = 0x0F;

char *progName;


// Calculate tLen (not done by newer versions of BLASR)
int getReferenceLength(const bam1_t *alignRecordIn);


//
// Main
//
int main(int argc, char *argv[]) {

	// Set program name
	progName = argv[0];

	//
	// Declarations
	//

	// Define types and modes
	enum ALIGN_FILE_TYPE {SAM, BAM, CRAM};
	const char *ALIGN_WRITE_MODE[] = {"w", "wb", "wc"};

	// Options
	string inFileName;
	string outFileName;
	string refFileName;

	// File types
	ALIGN_FILE_TYPE outFileType;

	// BAM Objects
	samFile *inFile;      // Alignment input file
	bam_hdr_t *inHeader;  // Header for input file

	samFile *outFile;      // Output file

	bam1_t *alignRecordIn = bam_init1();  // Alignment record

	int recordCount;  // Number of records processed

	// Record processing objects
	string qname;  // Query name (of the read)
	string chr;    // Chromosome read aligns to
	int chrId;     // Chromosome ID
	int pos;       // Chromosome position


	//
	// Parse Options
	//

	// Setup options
	po::options_description prog_opts("Fix local alignments and place them back in reference context");

	prog_opts.add_options()
			("help,h",
					"Print help."
			)

			("in,i", po::value<string>(&inFileName),
					"Input SAM."
			)

			("out,o", po::value<string>(&outFileName),
					"Output CRAM."
			)

			("ref,r", po::value<string>(&refFileName),
					"Reference FASTA the output and header files are aligned to. Option required if output is a CRAM file."
			)
			;

	// Process options
	po::variables_map vm;

	try {
		po::store(po::command_line_parser(argc, argv).options(prog_opts).run(), vm);
		po::notify(vm);

	} catch (const std::exception &ex) {
		err(ex);
		return ERR_USAGE;
	}

	//
	// Get file types
	//

	// Get output file type
	if (boost::algorithm::iends_with(outFileName, ".sam"))
		outFileType = SAM;
	else if (boost::algorithm::iends_with(outFileName, ".bam"))
		outFileType = BAM;
	else if (boost::algorithm::iends_with(outFileName, ".cram"))
		outFileType = CRAM;
	else {
		err("Unrecognized type for output file: %s: Must be sam, bam, or cram", outFileName.c_str());
		return ERR_USAGE;
	}

	// Check for reference (if output is CRAM)
	if (outFileType == CRAM && refFileName == "") {
		err("Reference file is required when output file is in cram format");
		return ERR_USAGE;
	}


	//
	// Open alignment files
	//

	inFile = nullptr;
	outFile = nullptr;

	// Open input file
	inFile = hts_open(inFileName.c_str(), "r");

	if (! inFile) {
		err("Error opening input file: %s", inFileName.c_str());
		return ERR_IO;
	}

	inHeader = sam_hdr_read(inFile);

	// Open output file
	outFile = hts_open(outFileName.c_str(), ALIGN_WRITE_MODE[outFileType]);

	if (! outFile) {
		err("Error opening input file: %s", outFileName.c_str());
		return ERR_IO;
	}

	// Write header
	if (refFileName != "") {
		if (hts_set_opt(outFile, CRAM_OPT_REFERENCE, refFileName.c_str()) < 0) {
			err("Error setting reference file name for output: %s", refFileName.c_str());
			return ERR_IO;
		}
	}

	if (sam_hdr_write(outFile, inHeader) < 0) {
		err("Error writing headers to output file %s", outFileName.c_str());
	}


	//
	// Process reads
	//
	recordCount = 0;

	while (sam_read1(inFile, inHeader, alignRecordIn) >= 0) {
		++recordCount;

		//cout << "Record: " << recordCount << endl;  // DBGTMP

		alignRecordIn->core.isize = getReferenceLength(alignRecordIn);

		// cout << "\tSize = " << alignRecordIn->core.isize << endl;

		// Write
		// cout << "\tWriting" << endl;

		if (sam_write1(outFile, inHeader, alignRecordIn) < 0) {
			err("Error writing record %d to output file %s", recordCount, outFileName.c_str());
		}

		// cout << "\tDone writing" << endl;
	}


	//
	// Close
	//

	bam_hdr_destroy(inHeader);

	sam_close(inFile);
	sam_close(outFile);


	//
	// Exit
	//

	return ERR_NONE;
}



/**
 * Get the number of reference bases covered by an alignment record.
 */
int getReferenceLength(const bam1_t *alignRecordIn) {

	// Init
	const uint32_t *cigar_array = bam_get_cigar(alignRecordIn);  // Get CIGAR array (uint32 array)
	const uint32_t * const cigar_array_end = cigar_array + alignRecordIn->core.n_cigar;

	int tlen = 0;  // Calculated template length

	uint32_t cigarOp;
	uint32_t cigarLen;

	// Process CIGAR records
	while (cigar_array < cigar_array_end) {
		cigarOp = *cigar_array & CIGAR_MATCH;
		cigarLen = *cigar_array >> 4;


		switch (cigarOp) {
		case BAM_CEQUAL:
		case BAM_CDIFF:
		case BAM_CDEL:
		case BAM_CMATCH:
		case BAM_CREF_SKIP:
			tlen += cigarLen;
			break;
		}

		// Next record
		cigar_array += 1;
	}

	return tlen;
}
