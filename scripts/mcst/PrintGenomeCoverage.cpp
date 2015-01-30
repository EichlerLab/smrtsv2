#include "bam.h"
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
#include <iomanip>

using namespace std;
typedef struct {
    int beg, end;
    samfile_t *in;
} tmpstruct_t;



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

int main(int argc, char* argv[]) {


	bam_index_t *idx;
	bam_plbuf_t *buf;


	if (argc < 4) {
		cout << "usage: bamc out.bed [-in <filename> -in  <filename> ... ]  ... [options ]" << endl;
		cout << " Prints coverage by base." << endl;
		cout << " Options: " << endl;
		cout << "   -q q   min mapping quality (30)" << endl;
		cout << "   -b b   bin size . " << endl;
		exit(1);
	}
	string outFileName;
	outFileName = argv[1];

	int argi = 2;
	int minQuality = 30;
  int bin = 50;
	vector<string> bamFileNames;
	while ( argi < argc ) {
		if (strcmp(argv[argi], "-q") == 0) {
			++argi;
			minQuality = atoi(argv[argi]);
			++argi;
		}
		if (strcmp(argv[argi], "-b") == 0) {
			++argi;
			bin = atoi(argv[argi]);
			++argi;
		}
		if (strcmp(argv[argi], "-in") == 0) {
			++argi;
			bamFileNames.push_back(argv[argi]);
			++argi;
		}
	}
	ofstream outFile(outFileName.c_str());
	int i;
	int bamI = 0;
	vector<vector<int> > coverage;
	bam_header_t *header;
	int readIndex =0 ;
	long totalNumBases = 0;
	for (bamI = 0; bamI < bamFileNames.size(); bamI++) {

		BGZF *in;
		in = bam_open(bamFileNames[bamI].c_str(), "rb");


		header = bam_header_read(in);
		if (bamI == 0) {


			coverage.resize(header->n_targets);
			for (i = 0; i < header->n_targets; i++) {
				coverage[i].resize(header->target_len[i]/bin + (header->target_len[i]%bin== 0? 0 : 1), 0);
			}
		}

		bam1_t *b =  bam_init1();

		while (bam_read1(in, b) >= 0) {

			int tStart = b->core.pos;
			int tEnd   = b->core.pos + GetTLen(b);
			if (b->core.qual >= minQuality) {
				vector<int>* v = &coverage[b->core.tid];
				for (i = tStart; i < tEnd; i++) {
					(*v)[i/bin]+=1;
				}
			}
			++readIndex;
			totalNumBases += tEnd - tStart;
			if (readIndex % 10000 == 0) {
				cerr << "processed " << readIndex << " reads " << totalNumBases / 1000000000 << "Gb" << endl;
			}

			bam_destroy1(b);
			b = bam_init1();
		}
		bam_header_destroy(header);
		bam_close(in);
	}

	for (i = 0; i < header->n_targets; i++) {
		int p;
		int lastBinIndex = header->target_len[i] / bin + (header->target_len[i] / bin == 0? 0 : 1);
		for (p = 0; p < lastBinIndex - 1; p++) {
			outFile << header->target_name[i] << "\t" << p*bin << "\t" << (p+1)*bin << "\t" << std::setw(4) << float(coverage[i][p]) / bin << endl;
		}
		int lastBinLength = header->target_len[i] - (lastBinIndex -1) * bin;
		if (lastBinLength > 0) {
			outFile << header->target_name[i] << "\t" << p *bin << "\t" << header->target_len[i] << "\t" << float(coverage[i][p]) / lastBinLength << endl;

		}
	}
  outFile.close();
}

