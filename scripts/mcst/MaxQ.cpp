#include "sam.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <assert.h>

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

int main(int argc, char* argv[]) {
	

	bam_index_t *idx;
	bam_plbuf_t *buf;

	samfile_t *samFile;

	if (argc != 4) {
		cout << "usage: maxq input.bam binSize out.bed" << endl;
		exit(1);
	}

	samFile = samopen(argv[1], "rb", NULL);
//	idx = bam_index_load(argv[1]); // load BAM index
	

	
	int binSize = atoi(argv[2]);
	string outFileName = argv[3];
/*
	if (idx == 0) {
		cout << "BAM indexing file is not available." << endl;
		return 1;
	}
*/
	vector<vector< unsigned char> > maxMapQ;
	bam1_t *entry = new bam1_t;
	int nTargets = samFile->header->n_targets;
	maxMapQ.resize(nTargets);
	int i;

	for (i = 0;  i < nTargets; i++) {
		maxMapQ[i].resize(GetNBins(samFile->header->target_len[i], binSize), 0);
	}
	int nread = 0;
	while (samread(samFile, entry) > 0) {
		++nread;
		if (nread % 1000000 == 0) {
			cerr << nread << endl;
		}
		if (entry->core.tid >= 0) {
			assert(entry->core.tid < maxMapQ.size());
			StoreMapQV( entry->core.pos, entry->core.l_qseq, binSize, entry->core.qual, maxMapQ[entry->core.tid]);
		}
	}
	ofstream outFile;
	outFile.open(outFileName.c_str());
	for (i = 0; i < nTargets; i++) {
		int j;
		for (j = 0; j < maxMapQ[i].size(); j++) {
			outFile << samFile->header->target_name[i] << "\t" << j*binSize << "\t" << (j+0)*binSize - 1 << "\t" << (int) maxMapQ[i][j] << endl;
		}
	}
}
