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

int GetNBins(int length, int binSize) {
	return length / binSize + ( (length % binSize) ? 0 : 1);
}
int GetBin(int pos, int binSize) {
	return pos / binSize;
}

void SetMax(int pos, int q, vector<unsigned char> &quals) {
	if (q + 1 > quals[pos]) {
		quals[pos] = q + 1;
	}
}
void StoreMapQV(int pos, int readLen, int binSize, unsigned char mapqv, vector<unsigned char> &quals) {
	int readi = 0;
	for (readi = 0; readi < readLen; readi += binSize) {
		SetMax(GetBin(pos+readi, binSize), mapqv, quals);
	}
}

class Exon {
public:
  vector<int> coverage;
  int start;
  int end;
	int minmapq;
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



int SetMapq(const bam1_t *b, void *data) {
  Exon *exon = (Exon*) data;
	// 
	// determine overlap with exon
	// 
	if (b->core.qual < exon->minmapq) {
		return 0;
	}
	int tStart = b->core.pos;
	int tEnd   = b->core.pos + GetTLen(b);
	int ovpStart = max(tStart, exon->start);
	int ovpEnd   = min(tEnd, exon->end);

	int i;

	for (i = ovpStart; i < ovpEnd; i++) {
		exon->coverage[i-exon->start] = 1;
	}
	return ovpEnd - ovpStart;
}

int main(int argc, char* argv[]) {
	

	bam_index_t *idx;
	bam_plbuf_t *buf;

	bamFile bamFile;

	if (argc != 4) {
		cout << "usage: maxq input.bam binSize out.bed" << endl;
		exit(1);
	}
	string bedFileName, outFileName;
	bamFile     = bam_open(argv[1], "rb");
	bedFileName = argv[2];
	outFileName = argv[3];

  idx = bam_index_load(argv[1]);
  ifstream bedFile(bedFileName.c_str());
	ofstream outFile(outFileName.c_str());

	bam1_t *entry = new bam1_t;

	int i;
	bam_header_t *header = bam_header_read(bamFile);

	map<string, int> chromToTID;
	for (i = 0; i < header->n_targets; i++) {
		chromToTID[header->target_name[i]] = i;
	}

	while (bedFile) {
		string line;
		getline(bedFile, line);
		if (line == "") break;
		stringstream strm(line);
		string chrom;
		int start, end;
		strm >> chrom >> start >> end;
		float fracCov = 0;
		if (end - start > 0) {
			int tid = chromToTID[chrom];
			
			Exon exon;
			exon.coverage.resize(end - start);
			exon.start = start; exon.end = end; exon.minmapq = 10;
			bam_fetch(bamFile, idx, tid, start, end, &exon, SetMapq);
			float nCov = accumulate(exon.coverage.begin(), exon.coverage.end(),0);
			fracCov = nCov / (end - start);
		}
		outFile << line << "\t" << fracCov << endl;
	}
	
	return 0;
	
}
