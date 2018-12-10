/*
 * alignfixup.cpp
 *
 * Reads are extracted from a region, assembled, and aligned back to that region (which is
 * extracted from the reference). This tool takes the aligned assembly and finalizes them
 * for merging with other assemblies.
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <cstring>

#include "util/err.h"
#include "util/constants.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "htslib/sam.h"

using namespace std;

namespace po = boost::program_options;


// Declarations
char *progName;

const int ALIGN_BUF_LEN = 512;  // Length of buffers holding alignment data

const int FLAG_NON_PRIMARY = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;  // Flags if read is not primary

// Set chrN and pos from a string in the format "chr:pos-end" where coordinates are 1-based.
bool setAlignPos(const char *chrNameHdr, int &chrId, int &pos, bam_hdr_t *outHeader, int recordCount);

// Set contig name on the output record
bool setContigName(bam1_t *alignRecordIn, bam1_t *alignRecordOut, const uint8_t *regionNameC, int regionNameLen);

// Calculate number of reference bases covered. Older versions of BLASR stored this in the TLEN field, but newer
// versions do not. Calculate it here for tools that use the field (e.g. inversion calling).
const uint32_t CIGAR_MATCH = 0x0F;

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
	string headFileName;
	string refFileName;
	string regionName;

	uint8_t *regionNameC;
	size_t regionNameLen;

	// File types
	ALIGN_FILE_TYPE outFileType;

	// BAM Objects
	samFile *inFile;      // Alignment input file
	bam_hdr_t *inHeader;  // Header for input file

	samFile *outFile;      // Output file
	bam_hdr_t *outHeader;  // Output header

	bam1_t *alignRecordIn = bam_init1();  // Alignment record
	bam1_t *alignRecordOut = bam_init1();  // Alignment record

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

			("head,d", po::value<string>(&headFileName),
					"Extract headers from this alignment file."
			)

			("ref,r", po::value<string>(&refFileName),
					"Reference FASTA the output and header files are aligned to. Option required if output is a CRAM file."
			)

			("region,g", po::value<string>(&regionName),
					"Name of the region contigs belong to. This is pre-pended to contigs."
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

	// Get region name as a c-string (append separator)
	regionName = regionName + "|";
	regionNameC = (uint8_t *) regionName.c_str();
	regionNameLen = regionName.length();


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

	// Get headers (re-use inFile to read headers)
	inFile = hts_open(headFileName.c_str(), "r");

	if (! inFile) {
		err("Error opening header file: %s", headFileName.c_str());
		return ERR_IO;
	}

	outHeader = sam_hdr_read(inFile);
	hts_close(inFile);
	inFile = nullptr;

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

	if (sam_hdr_write(outFile, outHeader) < 0) {
		err("Error writing headers to output file %s", outFileName.c_str());
	}


	//
	// Process reads
	//
	recordCount = 0;

	while (sam_read1(inFile, inHeader, alignRecordIn) >= 0) {
		++recordCount;

		// Skip supplementary
		if (alignRecordIn->core.flag & FLAG_NON_PRIMARY)
			continue;

		// Get chromosome and position
		if (! setAlignPos(inHeader->target_name[alignRecordIn->core.tid], chrId, pos, outHeader, recordCount))
			return ERR_FORMAT;  // Err message already output

		alignRecordOut->core.tid = chrId;
		alignRecordOut->core.pos = pos + alignRecordIn->core.pos;

		// Make globally-unique chr name
		if (! setContigName(alignRecordIn, alignRecordOut, regionNameC, regionNameLen))
			return ERR_FORMAT;  // Err message already output

		// Set other core data
		alignRecordOut->core.qual = alignRecordIn->core.qual;
		alignRecordOut->core.flag = alignRecordIn->core.flag;
		alignRecordOut->core.n_cigar = alignRecordIn->core.n_cigar;
		alignRecordOut->core.l_qseq = alignRecordIn->core.l_qseq;

		alignRecordOut->core.mtid = -1;
		alignRecordOut->core.mpos = -1;

		alignRecordOut->core.bin = alignRecordIn->core.bin;
		alignRecordOut->core.unused1 = alignRecordIn->core.unused1;
		alignRecordOut->core.isize = getReferenceLength(alignRecordIn);

		if (alignRecordIn->core.mtid >= 0)
			warn("Found split alignment in record %d: %s: RNEXT and PNEXT have been cleared", recordCount, (char *) alignRecordIn->data);

		// Write
		if (sam_write1(outFile, outHeader, alignRecordOut) < 0) {
			err("Error writing record %d to output file %s", recordCount, outFileName.c_str());
		}
	}


	//
	// Close
	//

	bam_hdr_destroy(outHeader);
	bam_hdr_destroy(inHeader);

	sam_close(inFile);
	sam_close(outFile);


	//
	// Exit
	//

	return ERR_NONE;
}

/**
 * Sets chromosome ID and position. Not thread-safe.
 *
 * @param chrNameHdr: Chromosome name from the input header.
 * @param chrId: Value is set to the chromosome ID.
 * @param pos: Value is set to the alignment position.
 * @param outHeader: Header for the output file.
 * @param recordCount: Record number for error reporting.
 */
bool setAlignPos(const char *chrNameHdr, int &chrId, int &pos, bam_hdr_t *outHeader, int recordCount) {

	char *chrEnd;
	char *posEnd;

	static char chrName[ALIGN_BUF_LEN];

	// Copy chromosome name
	strncpy(chrName, chrNameHdr, ALIGN_BUF_LEN);

	// Find end of chr name
	chrEnd = chrName;

	while (*chrEnd != '\0' && *chrEnd != ':')
		++chrEnd;

	if (*chrEnd == '\0' || chrEnd == chrName) {
		err("Unrecognized reference format for input contig in record %d: \"%s\": Missing chr: Expected \"chr:start-end\"",
				recordCount, chrName
		);

		return false;
	}

	// Find end of the position string
	posEnd = chrEnd + 1;

	while (*posEnd != '\0' && *posEnd != '-')
		++posEnd;

	if (*posEnd == '\0' || posEnd == chrEnd + 1) {
		err("Unrecognized reference format for input contig in record %d: \"%s\": Missing pos: Expected \"chr:start-end\"",
				recordCount, chrName
		);

		return false;
	}

	// Get chr id from header
	*chrEnd = '\0';

	if ((chrId = bam_name2id(outHeader, chrName)) < 0) {
		err("Contig name \"%d\" in record %d does not exist in headers",
				chrName, recordCount
		);

		return false;
	}

	++chrEnd;  // Point to beginning of position

	// Get position
	*posEnd = '\0';

	if ((pos = atoi((char *) chrEnd)) == 0) {
		err("Position \"%s\" in record %d is not an integer",
				(char *) chrEnd, recordCount);

		return false;
	}

	--pos;  // To 0-based coordinates

	return true;
}

/**
 * Set contig name in output record.
 */
bool setContigName(bam1_t *alignRecordIn, bam1_t *alignRecordOut, const uint8_t *regionNameC, int regionNameLen) {

	int qnameLen;
	int newQnameLen;
	int extraNul;

	uint8_t *dataNoQname;
	uint8_t *dataPtr;

	// Set new qname length
	qnameLen = alignRecordIn->core.l_qname - alignRecordIn->core.l_extranul;  // Name length without null-padding, but with null terminator
	newQnameLen = qnameLen + regionNameLen;
	extraNul = 4 - (newQnameLen - 1) % 4 - 1;  // Number of nulls to pad to reach 32-bit boundary

	// Setup transfer
	alignRecordOut->l_data = newQnameLen + 1 + extraNul + (alignRecordIn->l_data - alignRecordIn->core.l_qname);

	// Expand data buffer
	if (alignRecordOut->l_data > alignRecordOut->m_data) {
		alignRecordOut->m_data = alignRecordOut->l_data;
		kroundup32(alignRecordOut->m_data);
		alignRecordOut->data = (uint8_t *) realloc(alignRecordOut->data, alignRecordOut->m_data);
	}

	// Copy prefix
	dataPtr = alignRecordOut->data;
	memcpy(dataPtr, regionNameC, regionNameLen);
	dataPtr += regionNameLen;

	// Copy qname and pad
	memcpy(dataPtr, alignRecordIn->data, qnameLen);
	dataPtr += qnameLen;

	for (int count = 0; count < extraNul; ++count) {
		*dataPtr = '\0';
		++dataPtr;
	}

	// Copy remaining line
	dataNoQname = alignRecordIn->data + alignRecordIn->core.l_qname;

	memcpy(dataPtr, dataNoQname, alignRecordIn->l_data - alignRecordIn->core.l_qname);

	// Set record core
	alignRecordOut->core.l_qname = newQnameLen + extraNul;
	alignRecordOut->core.l_extranul = extraNul;

	return true;
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
