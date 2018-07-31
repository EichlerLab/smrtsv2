#include <iostream>
#include <fstream>

#include "phase/ReadPhaser.h"
#include "phase/PhaseWriter.h"

#include "htslib/sam.h"

using namespace std;

namespace phase {

const char SEQI_TO_CHAR[] {
	'*', 'A', 'C', '*',
	'G', '*', '*', '*',
	'T', '*', '*', '*',
	'*', '*', '*', 'N',
};

void ReadPhaser::phase(const vector<string> *alignFileList, const string outFileName) {

	// Declarations
	PhaseWriter phaseWriter(outFileName);  // Phase writer

	string *alignFileName;  // Name of an alignment input file

	samFile *bamIn;                     // Alignment input file
	bam_hdr_t *bamHdr;                  // Alignment file header
	bam1_t *alignRecord = bam_init1();  // Alignment record

	uint32_t lastTid;             // Last template ID - ID snvVector contains records for
	const vector<SNVEntry> *snvVector;  // Vector of SNVEntry elements describing SNV sites

	vector<SNVEntry>::const_iterator snvStart;  // Position in SNV array where reads begin
	vector<SNVEntry>::const_iterator snvNext;   // Next SNV to process in the read
	vector<SNVEntry>::const_iterator snvEnd;    // End position of the SNV

	int32_t posTemplate;  // Template (reference) position
	int32_t posQuery;     // Query (read) position

	const uint32_t *cigar;  // Array of CIGAR operations
	uint32_t cigarIndex;    // Current index in cigar
	uint32_t nCigar;        // Number of CIGAR operations
	uint32_t cigarOp;       // CIGAR operation code
	uint32_t cigarLen;      // CIGAR operation length

	int count0;      // Number of SNVs in haplotype 0
	int count1;      // Number of SNVs in haplotype 1
	int countObase;  // Number of SNVs that were not consistent with haplotype 0 or 1
	int countGap;    // Number of SNVs aligning to a gap in the read

	int ps;  // Phase set - detects phase set changes

	uint8_t *seq;  // Sequence of bases - encoded 2 bases per byte
	char base;     // Base at a SNV position in a read

	// Init state
	lastTid = -1;  // Out of range, will load template in the first sequence read
	snvVector = nullptr;

	// Read each file
	for (vector<string>::const_iterator inIter = alignFileList->begin(); inIter != alignFileList->end(); ++inIter) {

		// Report input file
		if (verbose)
			cout << "Opening: " << *inIter << endl;

		// Open alignment file
		bamIn = hts_open(inIter->c_str(), "r");
		bamHdr = sam_hdr_read(bamIn);

		// Read each record
		while (sam_read1(bamIn, bamHdr, alignRecord) >= 0) {

			// Skip unmapped reads
			if (alignRecord->core.flag & BAM_FUNMAP)
				continue;

			// Get SNV table for the template
			if (alignRecord->core.tid != lastTid) {
				snvVector = phaseTable->getEntry(bamHdr->target_name[alignRecord->core.tid]);

				snvStart = snvVector->begin();
				snvEnd = snvVector->end();
			}

			// Seek to start to the first position within this read
			while (snvStart != snvEnd && snvStart->pos < alignRecord->core.pos)
				++snvStart;

			// Setup to traverse SNVs in this read
			snvNext = snvStart;

			nCigar = alignRecord->core.n_cigar;
			cigar = bam_get_cigar(alignRecord);
			seq = bam_get_seq(alignRecord);

			posTemplate = alignRecord->core.pos;
			posQuery  = 0;

			count0 = 0;
			count1 = 0;
			countObase = 0;
			countGap = 0;

			ps = snvNext->ps;

			// Parse CIGAR - find SNV sites
			cigarIndex = 0;

			while (cigarIndex < nCigar && snvNext != snvEnd) {
				cigarOp = cigar[cigarIndex] & BAM_CIGAR_MASK;
				cigarLen = cigar[cigarIndex] >> BAM_CIGAR_SHIFT;

				// Match CIGAR op
				switch(cigarOp) {

				// Matching bases
				case BAM_CMATCH:
				case BAM_CEQUAL:
				case BAM_CDIFF:

					if (posTemplate + cigarLen > snvNext->pos) {
						base = SEQI_TO_CHAR[bam_seqi(seq, posQuery + (snvNext->pos - posTemplate))];

						// Increment base count
						if (base == snvNext->allele0)
							count0 += 1;
						else if (base == snvNext->allele1)
							count1 += 1;
						else
							countObase += 1;

						// Next SNV
						snvNext += 1;

						// Check phase-set
						if (snvNext != snvEnd && ps != snvNext->ps) {
							phaseWriter.addPhase(ps, count0, count1, countObase, countGap);

							ps = snvNext->ps;

							count0 = 0;
							count1 = 0;
							countObase = 0;
							countGap = 0;
						}

						// Continue - Next SNV may be in the same CIGAR op
						continue;
					}

					posTemplate += cigarLen;
					posQuery += cigarLen;

					break;

				// Inserted bases (skip)
				case BAM_CINS:

					posQuery += cigarLen;

					break;

				// Deleted and skipped bases
				case BAM_CDEL:
				case BAM_CREF_SKIP:

					if (posTemplate + cigarLen > snvNext->pos) {

						// Increment count gap
						countGap += 1;

						// Next SNV
						snvNext += 1;

						// Check phase-set
						if (snvNext != snvEnd && ps != snvNext->ps) {
							phaseWriter.addPhase(ps, count0, count1, countObase, countGap);

							ps = snvNext->ps;

							count0 = 0;
							count1 = 0;
							countObase = 0;
							countGap = 0;
						}

						// Continue - Next SNV may be in the same CIGAR op
						continue;
					}

					posTemplate += cigarLen;

					break;

				// Clipped bases
				case BAM_CSOFT_CLIP:

					posQuery += cigarLen;

					break;

				// Ignore other operations (hard clipping, padding)
				}

				// Next CIGAR op
				cigarIndex += 1;

			}  // loop: CIGAR ops

			// Add phase data
			phaseWriter.addPhase(ps, count0, count1, countObase, countGap);

			// Write record
			phaseWriter.writePhase(alignRecord, bamHdr);

		}  // loop: SAM records

	}  // loop: Alignment input files

	// Free resources
	bam_hdr_destroy(bamHdr);


	return;
}

} /* namespace phase */

