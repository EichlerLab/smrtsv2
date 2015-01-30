#include "sam.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <numeric>

using namespace std;
typedef struct {  
    int beg, end;  
    samfile_t *in;  
} tmpstruct_t;  
  

class Output {
public:
	ofstream *outFilePtr;
	int minClipping;
	int minq;
	string chrom;
};

int GetTLen(const bam1_t *b) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int tlen = 0;
	int i;
	for (i = 0; i < len; i++) {
		int op = cigar[i] & 0xf;
		int oplen = cigar[i] >> 4;
		if (op == BAM_CMATCH or op == BAM_CDEL) {
			tlen += oplen;
		}
	}
	return tlen;
}

int GetStrand(const bam1_t *b) {
	if (b->core.flag & 16) {
		return 1;
	}
	else {
		return 0;
	}
}

int GetClipping(const bam1_t *b, int &leftClip, int &rightClip) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int qlen = 0;
	int i;
	leftClip = rightClip = 0;
	
	if (len > 0) {
		int first = 0;
		int op = cigar[first] & 0xf;
		int oplen = cigar[first] >> 4;
		if (op == BAM_CHARD_CLIP) {
			first++;
			op = cigar[first] & 0xf;
			oplen = cigar[first] >> 4;
		}
		if (op == BAM_CSOFT_CLIP) {
			leftClip = oplen;
		}
		int last = len - 1;
		if (last > first) {
			op = cigar[last] & 0xf;
			oplen = cigar[last] >> 4;
			if ( op == BAM_CHARD_CLIP) {
				last--;
			}
			if (last > first) {
				op = cigar[last] & 0xf;
				oplen = cigar[last] >> 4;
				if (op == BAM_CSOFT_CLIP) {
					rightClip = oplen;
				}
			}
		}
	}
	return 0;
}

int numProcessed = 0;

const int MIN_ALIGNED_LENGTH = 500;
int numFound = 0;
int WriteHardStop(const bam1_t *b, void *data) {
  Output *output = (Output*) data;
	// 
	// determine overlap with exon
	// 
	if (b->core.qual < output->minq) {
		return 0;
	}
	int leftClipping, rightClipping;
	GetClipping(b, leftClipping, rightClipping);
	int tStart = b->core.pos;
	int tEnd   = b->core.pos + GetTLen(b);
	if (b->core.l_qseq - leftClipping - rightClipping  < MIN_ALIGNED_LENGTH ) {
		return 0;
	}
	if (leftClipping > output->minClipping or rightClipping > output->minClipping) {
		(*output->outFilePtr) << output->chrom << "\t" << tStart << "\t" << tEnd << "\t" << bam1_qname(b) << "\t" << leftClipping << "\t" << rightClipping << "\t" << GetStrand(b) << endl;
	}		
	++numProcessed;
}


int GetQLen(const bam1_t *b) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int qlen = 0;
	int i;
	for (i = 0; i < len; i++) {
		int op = cigar[i] & 0xf;
		int oplen = cigar[i] >> 4;
		if (op == BAM_CMATCH or op == BAM_CINS) {
			qlen += oplen;
		}
	}
	return qlen;
}


int main(int argc, char* argv[]) {
	

	bam_index_t *idx;
	bam_plbuf_t *buf;

	bamFile bamFile;

	if (argc < 4) {
		cout << "usage: splitreads input.bam minClipping out.bed [region]" << endl;
		exit(1);
	}
	string outFileName;
	string region;
	samfile_t *in;  
	in = samopen(argv[1], "rb", 0);

	int minClippingLength = atoi(argv[2]);
	outFileName = argv[3];

  idx = bam_index_load(argv[1]);

	region = "";
	if (argc == 5) {
		region = argv[4];
	}
	
	Output output;

	ofstream outFile(outFileName.c_str());


	output.outFilePtr = &outFile;
	output.minq = 0;
	output.minClipping = minClippingLength;

	bam1_t *entry = new bam1_t;

	int i;


	tmpstruct_t tmp;  
	int ref;
	
	if (region != "") {
		ifstream regionFile(region.c_str());
		vector<string> regions;
		if (regionFile) {
			while (regionFile) {
				regionFile >> region;
				regions.push_back(region);
			}
		}
		else {
			regions.push_back(region);
		}			
		int r;
		for (r = 0; r < regions.size(); r++) {
			region = regions[r];
			cout << "Parsing region " << region << endl;
			bam_parse_region(in->header, region.c_str(), &ref,  
											 &tmp.beg, &tmp.end); // parse the region  
			output.chrom = in->header->target_name[ref];
			bam_fetch(in->x.bam, idx, ref, tmp.beg, tmp.end, &output, WriteHardStop);
		}
	}
	else {

		bam1_t *b = bam_init1();	

		while (bam_read1(in->x.bam, b) > 0) {
			output.chrom = in->header->target_name[b->core.tid];
			WriteHardStop(b, &output);
		}
	}

	outFile.close();
	
}
