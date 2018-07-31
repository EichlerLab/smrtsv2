#include <stdexcept>

#include "phase/PhaseWriter.h"
#include "phase/StringUtil.h"

#include "htslib/sam.h"

using namespace std;

namespace phase {

PhaseWriter::PhaseWriter(string outFileName) {

	// Open table file
	if (outFileName == "")
		throw_error<runtime_error>("Cannot create phase output writer: Output file name is empty");

	phaseTableFile.open(outFileName);

	if (! phaseTableFile.is_open())
		throw_error<runtime_error>("Failed opening output file: %s", outFileName.c_str());

	phaseTableFile << "#CHROM\tPOS\tEND\tID\tPS\tHAP0\tHAP1\tOTHER\tGAP" << endl;

	return;
}

PhaseWriter::~PhaseWriter() {
	phaseNodeList.clear();

	// Close table file
	phaseTableFile.close();

	return;
}

void PhaseWriter::addPhase(int phaseSet, int count0, int count1, int countObase, int countGap) {

	if (count0 > 0 || count1 > 0 || countObase > 0 || countGap > 0)
		phaseNodeList.push_back(PhaseWriterNode(phaseSet, count0, count1, countObase, countGap));

	return;
}

void PhaseWriter::writePhase(const bam1_t *alignRecord, const bam_hdr_t *bamHdr) {

	// Write phase table entries
	for (vector<PhaseWriterNode>::const_iterator phaseNode = phaseNodeList.begin(); phaseNode != phaseNodeList.end(); ++phaseNode) {
		phaseTableFile <<
				bamHdr->target_name[alignRecord->core.tid] << '\t' <<
				alignRecord->core.pos << '\t' <<
				alignRecord->core.pos +  bam_cigar2rlen(alignRecord->core.n_cigar, bam_get_cigar(alignRecord)) << '\t' <<
				(char *) bam_get_qname(alignRecord) << '\t' <<
				phaseNode->phaseSet << '\t' <<
				phaseNode->count0 << '\t' <<
				phaseNode->count1 << '\t' <<
				phaseNode->countObase << '\t' <<
				phaseNode->countGap << '\n'
				;
	}

	// Clear phase state
	phaseNodeList.clear();

	return;
}

} /* namespace phase */
