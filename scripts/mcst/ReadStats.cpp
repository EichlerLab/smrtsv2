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


float GetAccuracy(const bam1_t *b) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int tlen = 0;
	int i;
	int nMatch = 0;
	int nDel = 0;
	int nIns = 0;
	for (i = 0; i < len; i++) {
		int op = cigar[i] & 0xf;
		int oplen = cigar[i] >> 4;
		if (op == BAM_CMATCH) {
			nMatch += oplen;
		}
		if (op == BAM_CINS) {
			nIns += oplen;
		}
		if (op == BAM_CDEL) {
			nDel += oplen;
		}
	}
	return (nMatch*1.0)/ (nMatch + nIns + nDel);
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

	if (argc < 3) {
		cout << "usage: bamacc file.bam output_table.txt [-n]" << endl;
		cout << "  -n  prints first column as read names." << endl;
		exit(1);
	}

	samFile = samopen(argv[1], "rb", NULL);
	
	string outFileName = argv[2];

	bool printName = false;
	if (argc == 4 and strcmp(argv[3], "-n") == 0) {
		printName = true;
	}

	bam1_t *entry = new bam1_t;


	int i;

	int nread = 0;
	ofstream outFile;
	outFile.open(outFileName.c_str());
	if (printName) {
		outFile << "name\t";
	}
	outFile << "qLength\ttLength\taccuracy"<< endl;
	while (samread(samFile, entry) > 0) {
		++nread;
		if (nread % 100000 == 0) {
			cerr << nread << endl;
		}
		if (entry->core.qual > 30) {
			int tLen = GetTLen(entry);
			int qLen = GetQLen(entry);
			float accuracy = GetAccuracy(entry);
			if (printName){ 
				outFile << bam1_qname(entry) << "\t";
			}

			outFile << samFile->header->target_name[entry->core.tid] << "\t" << entry->core.pos << "\t" << bam_calend(&entry->core, bam1_cigar(entry)) << "\t";
			outFile << qLen << "\t" << tLen << "\t" << accuracy << endl;
		}
	}
	outFile.close();
}
