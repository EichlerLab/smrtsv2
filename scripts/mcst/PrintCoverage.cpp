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
  

class Coverage {
public:
	map<string, vector<int> > counts;
	int minq;
	int start;
	int end;
	map<int, string> tidKey;
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

int IncrementCoverage(const bam1_t *b, void *data) {
  Coverage *coverage = (Coverage*) data;
	// 
	// determine overlap with exon
	// 
	if (b->core.qual < coverage->minq) {
		return 0;
	}
	int tStart = b->core.pos;
	int tEnd   = b->core.pos + GetTLen(b);
	int ovpStart = max(tStart, coverage->start);
	int ovpEnd   = min(tEnd, coverage->end);
	
	int i;

	for (i = ovpStart; i < ovpEnd; i++) {
		coverage->counts[i-coverage->start] +=1;
	}

}

int main(int argc, char* argv[]) {
	

	bam_index_t *idx;
	bam_plbuf_t *buf;

	bamFile bamFile;

	if (argc != 4) {
		cout << "usage: bamc input.bam region out.txt [options ]" << endl;
		cout << " Prints coverage by base." << endl;
		cout << " Options: " << endl;
		cout << "   -q q   min mapping quality (30)" << endl;
		cout << "   -b b   bin size . " << endl;
		exit(1);
	}
	string outFileName;
	string region;
	samfile_t *in;  
	in = samopen(argv[1], "rb", 0);
		//	bamFile     = bam_open(argv[1], "rb");
	region = argv[2];
	outFileName = argv[3];

  idx = bam_index_load(argv[1]);
	int argi = 4;
	int minQuality = 30;
	while ( argi < argc ) {
		if (strcmp(argv[argi], "-q") == 0) {
			++argi;
			minQuality = atoi(argv[argi]);
		}
		++argi;
	}
		
	ofstream outFile(outFileName.c_str());

	bam1_t *entry = new bam1_t;

	int i;
	//	bam_header_t *header = bam_header_read(bamFile);

	tmpstruct_t tmp;  
	int ref;
	bam_parse_region(in->header, region.c_str(), &ref,  
									 &tmp.beg, &tmp.end); // parse the region  
	Coverage coverage;
	coverage.counts.resize(tmp.end-tmp.beg,0);
	coverage.start = tmp.beg;
	coverage.end = tmp.end;
	coverage.minq = 30;
		
	bam_fetch(in->x.bam, idx, ref, tmp.beg, tmp.end, &coverage, IncrementCoverage);

	for (i = 0; i < coverage.counts.size(); i++) {
		outFile << tmp.beg + i << "\t" << coverage.counts[i] << endl;
	}
	outFile.close();
	
}
