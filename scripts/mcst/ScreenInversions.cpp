#include <iostream>
#include <string>
#include <sstream>
#include <map>

using namespace std;
#include "InversionAlign.h"

#include "FASTAReader.h"
#include "FASTASequence.h"
#include <semaphore.h>

class MappingSemaphores {
  public:
	sem_t reader;
	sem_t writer;
	sem_t unaligned;
	sem_t hitCluster;
	MappingSemaphores& operator=(MappingSemaphores &rhs) {
		return *this;
	}

	void InitializeAll() {
		sem_init(&reader, 0, 1);
		sem_init(&writer, 0, 1);
		sem_init(&unaligned, 0, 1);
		sem_init(&hitCluster, 0, 1);
	}
};

ifstream samIn;
vector<FASTASequence> genome;
map<string, FASTASequence*> chromMap;
MappingSemaphores semaphores;
ofstream tableOut;
int seqIndex;
bool isFofn, makeRC, makeDotplot, doSoftClip;
int window = 0;
int wordSize = 11;
void ScreenInversions() {

	while (true) {
		sem_wait(&semaphores.reader);
		if (samIn == false) {
			sem_post(&semaphores.reader);
			break;
		}
		string line;
		getline(samIn, line);

		if (line.size() == 0 or line[0] == '@') {
			sem_post(&semaphores.reader);
			if (line.size() == 0) {
				break;
			}
			else {
				continue;
			}
		}
		//
		// Finished reading line
		//
		sem_post(&semaphores.reader);

		if (isFofn) { 
			ifstream localSam(line.c_str());
			
			getline(localSam, line);
			//
			// Advance to the first alignment
			//

			while (line.size() > 0 and line[0] == '@') {
				getline(localSam, line);
			}
		}
		stringstream lineStrm(line);
		string readName, chrom, cigar, tmpString, seq;
		int flag, mapq, tmp, pos, tLen;
		int pnext;
		lineStrm >> readName >> flag >> chrom >> pos >> mapq >> cigar >> tmpString >> pnext >> tLen >> seq;
		FASTASequence ref;
		int additionalLength = 0;

		if (chromMap.find(chrom) == chromMap.end()) {
			cerr << "skipping:" << readName << " " << chrom << endl;
			continue;
		}
		if (pos + window + tLen > chromMap[chrom]->length) {
			int endWindow = chromMap[chrom]->length - pos - tLen;
			additionalLength += endWindow;
		}
		else {
			additionalLength += window;
		}
		if (pos - window < 0) {
			additionalLength += pos - window;
			pos = 0;
		}
		else {
			additionalLength += window;
			pos = pos - window;
		}
		ref.ReferenceSubstring(*chromMap[chrom], pos, tLen+additionalLength);
		ref.CopyTitle("genome");
		vector<int> lens;
		vector<char> ops;
		stringstream cigarStrm(cigar);
		while (cigarStrm) {
			int l;
			char o;

			if ( (cigarStrm >> l >> o) == 0) {
				break;
			}
			lens.push_back(l);
			ops.push_back(o);
		}
		
		int i= 0;
		int startOp = 0;
		int clipFront = 0;
		while (i < lens.size()) { 
			if (ops[i] == 'M' or ops[i] == 'I' or ops[i] == 'D') {
				startOp = i;
				break;
			}
			if (ops[i] == 'S') {
				clipFront = lens[i];
			}
			i++;
		}
		int clipEnd = 0;
		i = lens.size() - 1;
		while (i > startOp) {
			if (ops[i] == 'M' or ops[i] == 'I' or ops[i] == 'D') {
				startOp = i;
				break;
			}
			if (ops[i] == 'S') {
				clipEnd = lens[i];
			}
			i--;
		}			
		if (doSoftClip == false) {
			clipFront = clipEnd = 0;
		}
		seq = seq.substr(clipFront, (seq.size() - clipEnd) - clipFront);

		if (seq.size() == 0) {
			continue;
		}
		FASTASequence query;
		query.CopyTitle(readName);
		query.seq = (Nucleotide*) seq.c_str();
		query.length = seq.size();
		query.deleteOnExit = false;

		int res;
		int invCoords[4], invCoordsRC[4];

		if (makeRC) {
			FASTASequence rcQuery;
			query.MakeRC(rcQuery);


			int forwardScore, reverseScore;
			forwardScore = SDPAlignLite(query, ref, wordSize);
			reverseScore = SDPAlignLite(rcQuery, ref, wordSize);
			res = -1;
			int resRC = -1;
			float frac = 0.4;
			if (forwardScore > frac * query.length) {
				res = InversionAlign(query, ref, wordSize, invCoords, makeDotplot);
			}
			if (reverseScore > frac * query.length ) {
				resRC = InversionAlign(rcQuery, ref, wordSize, invCoordsRC, makeDotplot);
			}
			if (reverseScore > frac * query.length) {
				if (resRC > res and resRC > forwardScore) {
					memcpy(invCoords, invCoordsRC, sizeof(int)*4);
					res = resRC;
				}
			}
			rcQuery.Free();
		}
		else {
			res = InversionAlign(query, ref, wordSize, invCoords, makeDotplot);
		}
		if (res > 20) {

			float frac = 0.01;

			if (invCoords[0] > frac * query.length and invCoords[2] < (1-frac)*query.length) {
				cout << readName << " " << res << endl;
				sem_wait(&semaphores.writer);
				tableOut << chrom << "\t" << pos + invCoords[2] << "\t" << pos + invCoords[3] + wordSize << "\t" << readName << "\t" << res << endl;
				sem_post(&semaphores.writer);
			}
		}
		
		++seqIndex;
		if (seqIndex % 1000 ==0 ) {
			cerr << "processed " << seqIndex << endl;
		}
		delete ref.title;
		delete query.title;
	}
}
	

int main(int argc, char* argv[]) {
	stringstream helpStrm;
		helpStrm << " Usage: screenInversions sam genome output " << endl;
		helpStrm << "   -j (int) Use n cores." << endl;
		helpStrm << "   -f   Sam is fofn of sam files." << endl;
		helpStrm << "   -w (int) Expand reference by window." << endl;
		helpStrm << "   -r   Allow reverse complement." << endl;
		helpStrm << "   -d   Make dotplot" << endl;
		helpStrm << "   -k (int) Use k for matching ( k <= 15). " << endl;
		helpStrm << "   --noClip  Do not clip" << endl;
		string help = helpStrm.str();

	if (argc < 4) {
		cout << help;
		exit(1);
	}
	string inFileName = argv[1];
	string genomeFileName = argv[2];
	string outputFasta = argv[3];
	int argi = 4;
	int nProc = 1;
	isFofn = false;
	makeRC = false;
	makeDotplot = false;
	doSoftClip = true;
	window = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-j") == 0) {
			nProc = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-f") == 0) {
			isFofn = true;
		}
		else if (strcmp(argv[argi], "-r") == 0) {
			makeRC = true;
		}
		else if (strcmp(argv[argi], "-d") == 0) {
			makeDotplot = true;
		}
		else if (strcmp(argv[argi], "-w") == 0) {
			window = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "--noClip") == 0) {
			doSoftClip = false;
		}
		else if (strcmp(argv[argi], "-k") == 0) {
			wordSize = atoi(argv[++argi]);
		}
		else {
			cout << help << endl;
			cout << "Error with option " << argv[argi] << endl;
			exit(1);
		}
		
		++argi;
	}
	tableOut.open(outputFasta.c_str());
	//
	// Since many alignments will be performed against the genome, 
	// read it all at once, rather than seeking.
	//

	FASTAReader genomeReader;
	genomeReader.Initialize(genomeFileName);
	genomeReader.ReadAllSequences(genome);
	int i;
	for (i = 0; i < genome.size(); i++) {
		chromMap[genome[i].title] = &genome[i];
	}
	samIn.open(inFileName.c_str());
	int seqIndex = 0;

	pthread_t *threads = new pthread_t[nProc];
	pthread_attr_t *threadAttr = new pthread_attr_t[nProc];
	//  MappingSemaphores semaphores;
	//
	// When there are multiple processes running along, sometimes there
	// are semaphores to worry about.
	//

	semaphores.InitializeAll();
	int procIndex;
	for (procIndex = 0; procIndex < nProc; procIndex++ ){
		pthread_attr_init(&threadAttr[procIndex]);
	}

	for (procIndex = 0; procIndex < nProc; procIndex++) {
		pthread_create(&threads[procIndex],&threadAttr[procIndex], (void* (*)(void*))ScreenInversions, NULL);
	}

	for (procIndex = 0; procIndex < nProc; procIndex++) {
		pthread_join(threads[procIndex], NULL);
	}

}

