#ifndef SRC_PHASE_PHASEWRITER_H_
#define SRC_PHASE_PHASEWRITER_H_

#include <vector>
#include <string>
#include <fstream>

#include "phase/StringUtil.h"

#include "htslib/sam.h"

using namespace std;

namespace phase {

/**
 * Saves information for one read and one phase set.
 */
class PhaseWriterNode {
public:
	PhaseWriterNode(int phaseSet, int count0, int count1, int countObase, int countGap) :
		phaseSet(phaseSet),
		count0(count0), count1(count1),
		countObase(countObase), countGap(countGap) {

		return;
	}

	/** Phase set. */
	const int phaseSet;

	/** Number of bases supporting haplotype 0. */
	const int count0;

	/** Number of bases supporting haplotype 1. */
	const int count1;

	/** Number of bases not matching haplotype 0 or 1. */
	const int countObase;

	/** Number of HET sites aligning to a gap in the read alignment. */
	const int countGap;
};


/**
 * Tracks information for all phase sets for a read and writes to output files.
 */
class PhaseWriter {

public:

	/**
	 * Create a new writer.
	 *
	 * @throws runtime_error If <code>outFileName</code> is not "" and cannot be opened.
	 */
	PhaseWriter(string outFileName);

	/** Frees writer nodes. */
	virtual ~PhaseWriter();

	/**
	 * Add phase information for one phase set.
	 *
	 * @param phaseSet Phase set.
	 * @param count0 Number of HETs called haplotype 0.
	 * @param count1 Number of HETs called haplotype 1.
	 * @param countObase Number of HETs that were neither haplotype 0 or 1, but were a base.
	 * @param countGap Number of HETs with a gap at the locus.
	 */
	void addPhase(int phaseSet, int count0, int count1, int countObase, int countGap);

	/**
	 * Write information for one read. addPhase must have already have been called or this call is ignored.
	 * Phase information is written and cleared for the next sequence read.
	 *
	 * @param alignRecord BAM alignment record ptr for one read.
	 * @param bamHdr BAM header ptr.
	 */
	void writePhase(const bam1_t *alignRecord, const bam_hdr_t *bamHdr);

private:

	/** Linked-list of phase information for one read. */
	vector<PhaseWriterNode> phaseNodeList;

	/** File of tabular output per read/phase-set. Shows counts (hap0, hap1, obase, and gap). */
	ofstream phaseTableFile;
};

} /* namespace phase */

#endif /* SRC_PHASE_PHASEWRITER_H_ */
