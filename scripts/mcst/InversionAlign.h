#ifndef INVERSION_ALIGN_H_
#define INVERSION_ALIGN_H_

#include <math.h>
#include "algorithms/alignment/AlignmentUtils.h"
#include "algorithms/alignment/sdp/SparseDynamicProgramming.h"
#include "algorithms/alignment/sdp/SDPFragment.h"
#include "tuples/TupleList.h"
#include "tuples/TupleMatching.h"
#include "tuples/DNATupleList.h"
#include "algorithms/alignment/GraphPaper.h"
#include "FASTASequence.h"
#include <sstream>

template<typename T_Fragment>
class StrandedLexicographicFragmentSort {
 public:
	int operator()(const T_Fragment &a, const T_Fragment &b) const {
		//
		// Sort giving priority to reverse strand.
		//
		if (a.x == b.x && a.y == b.y) {
			return a.strand > b.strand;
		}
		else {
			return a.LessThanYX(b);
		}
	}

};


class StrandedFragment: public Fragment {
 public:
	int strand;
  StrandedFragment() : strand(1), Fragment() {}
  StrandedFragment(const DNALength &xp, const DNALength &yp) :Fragment(xp,yp) {
	}
	StrandedFragment &operator=(const StrandedFragment &rhs) {
  *((Fragment*)this) = (Fragment&)rhs;
    strand = rhs.strand;
  }
};


void StoreMatches(FASTASequence &query, FASTASequence &target, int wordSize, vector<StrandedFragment> &fragmentSet) {
	vector<StrandedFragment> rcFragmentSet;
	TupleList<PositionDNATuple> targetTupleList;

	/*
		Collect a set of matching fragments between query and target.
		Since this function is an inner-loop for alignment, anything to
		speed it up will help.  One way to speed it up is to re-use the
		vectors that contain the sdp matches. 
	*/

	TupleMetrics tm, tmSmall;
	tm.Initialize(wordSize);

	SequenceToTupleList(target, tm, targetTupleList);

  targetTupleList.Sort();
			 
  //
  // Store in fragmentSet the tuples that match between the target
  // and query.
  //

	StoreMatchingPositions(query, tm, targetTupleList, fragmentSet, 100);
	VectorIndex f;  
	for (f = 0; f < fragmentSet.size(); f++) {
		fragmentSet[f].weight = tm.tupleSize;
		fragmentSet[f].length = tm.tupleSize;
		fragmentSet[f].strand = 0;
  }
	
	
	FASTASequence rc;
	query.MakeRC(rc);
	
	StoreMatchingPositions(rc, tm, targetTupleList, rcFragmentSet, 100);

	for (f = 0; f < rcFragmentSet.size(); f++) {
		rcFragmentSet[f].weight = tm.tupleSize;
		rcFragmentSet[f].length = tm.tupleSize;
		rcFragmentSet[f].strand = 1;
		rcFragmentSet[f].x = query.length - rcFragmentSet[f].x - tm.tupleSize;
  }
	fragmentSet.insert(fragmentSet.begin(), rcFragmentSet.begin(), rcFragmentSet.end());
	std::sort(fragmentSet.begin(), fragmentSet.end(), StrandedLexicographicFragmentSort<StrandedFragment>());	

}




int SDPAlignLite(FASTASequence &query, FASTASequence &target, int wordSize) {

	vector<StrandedFragment> fragmentSet;
	TupleList<PositionDNATuple> targetTupleList;

	/*
		Collect a set of matching fragments between query and target.
		Since this function is an inner-loop for alignment, anything to
		speed it up will help.  One way to speed it up is to re-use the
		vectors that contain the sdp matches. 
	*/

	TupleMetrics tm, tmSmall;
	tm.Initialize(wordSize);

	SequenceToTupleList(target, tm, targetTupleList);

  targetTupleList.Sort();
			 
  //
  // Store in fragmentSet the tuples that match between the target
  // and query.
  //

	StoreMatchingPositions(query, tm, targetTupleList, fragmentSet, 100);
	VectorIndex f;  
	for (f = 0; f < fragmentSet.size(); f++) {
		fragmentSet[f].weight = tm.tupleSize;
		fragmentSet[f].length = tm.tupleSize;
		fragmentSet[f].strand = 0;
  }
	
	std::sort(fragmentSet.begin(), fragmentSet.end(), StrandedLexicographicFragmentSort<StrandedFragment>());

	//
	// Assume inversion will be in rc max frament chain set.
	//
	std::vector<int> maxFragmentChain;
	SDPLongestCommonSubsequence(query.length, fragmentSet, tm.tupleSize, 2, 2, -5, maxFragmentChain, Local);
  //GlobalChainFragmnetList(fragmentSet, maxFragmentChain);
	return maxFragmentChain.size();

 }

int RemoveSparseMatches(vector<StrandedFragment> &fragmentSet, int maxGap =100) {
	vector<bool> toRemove(fragmentSet.size(), 0);
	int f;
	if (fragmentSet.size() < 3) {
		fragmentSet.clear();
		return 0;
	}
	for (f = 1; f < fragmentSet.size()-1; f++) {
		if (fragmentSet[f].x - fragmentSet[f-1].x > maxGap or 
				fragmentSet[f].y - fragmentSet[f-1].y > maxGap or 
				fragmentSet[f+1].x - fragmentSet[f].x > maxGap or 
				fragmentSet[f+1].y - fragmentSet[f].y > maxGap) {
			toRemove[f] = true;
		}
	}
	int c = 0;
	for (f = 0; f < fragmentSet.size(); f++) {
		if (toRemove[f] == false) {
			fragmentSet[c] = fragmentSet[f];
			c++;
		}
	}
	int nRemoved = fragmentSet.size() - c;
	fragmentSet.resize(c);
	return nRemoved;
}
		

void RemoveDuplicates(vector<StrandedFragment> &fragmentSet) {
  int f = 0;
  int fCur = 0;
  while (f < fragmentSet.size()) {
		assert(fCur <= f);
    fragmentSet[fCur] = fragmentSet[f];
		int fi = f;
    while (f < fragmentSet.size() and fragmentSet[fi].x == fragmentSet[f].x and fragmentSet[fi].y == fragmentSet[f].y) {
      f++;
    }
    fCur++;

  }
  fragmentSet.resize(fCur);
}


int RemoveColRowMatches(	vector<StrandedFragment> &fragmentSet) {
	//
	// Try removing *all* duplicates for now.
	//
	int nFragments = fragmentSet.size();
	int f = 0;
	int fCur = 0;
	int fIter = 0;
	while (fIter + 1 < fragmentSet.size()) {
		int fNext = fIter + 1;
    while (fNext < fragmentSet.size() and 
      		 (fragmentSet[fIter].x == fragmentSet[fNext].x or fragmentSet[fIter].y == fragmentSet[fNext].y)) {
			fNext++;
		}
		if (fIter + 1 < fNext) {
			//
			// If there is a match, do not copy any of the values.
			//
			fIter = fNext;
			continue;
		}
		else {
			fragmentSet[fCur] = fragmentSet[fIter];
		}
		fIter = fNext;
    fCur++;
  }
  fragmentSet.resize(fCur);
	return nFragments - fragmentSet.size();
}
void ClearFragmentSet(vector<StrandedFragment> &fragmentSet) {
	int f;
	for (f = 0; f < fragmentSet.size(); f++) {
		fragmentSet[f].index = f;
		fragmentSet[f].chainPrev = 0;
		fragmentSet[f].above = -1;
	}
}

void CreateCandidateInversionChains(vector<StrandedFragment> &rcFragmentSet, 
																		int tupleSize, 
																		int queryLength,
																		int minChainLength,
																		vector<vector<StrandedFragment> > &chains) {

	chains.resize(0);
	bool chainWasAdded = true;
	int curChain = 0;
	while (chainWasAdded) {
		chainWasAdded = false;

		std::vector<int> rcMaxFragmentChain;
		ClearFragmentSet(rcFragmentSet);
		SDPLongestCommonSubsequence(queryLength, rcFragmentSet, tupleSize, 2, 2, -5, rcMaxFragmentChain, Global);


		if (rcMaxFragmentChain.size() == 0) {
			break;
		}
	

		//			chains[curChain].push_back(rcMaxFragmentChain[0]);
		int chainStart = 0;
		vector<int> toRemove;
		int f;
		for (f = 1; f < rcMaxFragmentChain.size(); f++) {
			StrandedFragment cur = rcFragmentSet[rcMaxFragmentChain[f]];
			StrandedFragment prev = rcFragmentSet[rcMaxFragmentChain[f-1]];
				
			int drift =  (int) (cur.x-prev.x)-(cur.y-prev.y);   //(int) sqrt(pow(cur.x - prev.x,2) + pow(cur.y - prev.y,2));
			drift = (drift < 0) ? -drift : drift;

			if (drift > tupleSize * 2) {
				//				cout << "candidate: " << f - chainStart << endl;
				if (f - chainStart >= minChainLength) {
					chainWasAdded = true;
					int c;
					chains.push_back(vector<StrandedFragment>() );
					for (c = chainStart; c < f; c++) {
						chains[curChain].push_back(rcFragmentSet[rcMaxFragmentChain[c]]);
						toRemove.push_back(rcMaxFragmentChain[c]);
					}
					curChain++;
				}
				chainStart = f;		
			}
		}
			
		if (toRemove.size() > 0) {
			// Remove this chain
			sort(toRemove.begin(), toRemove.end());
			int i, r;
			f = 0; i = 0; r = 0;
			for (f = 0; f < rcFragmentSet.size(); f++) {
				if (f == toRemove[r]){ 
					++r;
					continue;
				}
				else {
					rcFragmentSet[i] = rcFragmentSet[f];
					i++;
				}
			}
			rcFragmentSet.resize(i);
		}
	}
}	

void GetChainBox(vector<StrandedFragment> &chain, int tupleSize, int b[]) {
	if (chain.size() == 0) {
		b[0] = b[1] = b[2] = b[3] = 0;
		return;
	}
		 
	int last = chain.size() - 1;
	b[0] = chain[0].x + 1;
	b[1] = chain[last].x - tupleSize + 1;
	b[2] = chain[0].y + tupleSize;
	b[3] = chain[last].y;
}

void SwapStrand(vector<StrandedFragment> &chain, int queryLength) {
	int i;
	for (i = 0; i< chain.size(); i++) {
		chain[i].x = queryLength - chain[i].x - 1;
	}
}

//
// Now, for each chain, rotate the coordinates according to the reversal of the 
// aligned sequence.
//
void RotateChain(vector<StrandedFragment> &chain, int tupleSize) {
	int last = chain.size() - 1;
	int chainEnd = chain[0].x + 1;
	int chainBegin = chain[last].x - tupleSize + 1;
	int chainLength = chainEnd - chainBegin;
	int ci;
	for (ci = 0; ci < chain.size(); ci++) {
		chain[ci].x = chainBegin + (chainEnd - chain[ci].x - 1);
	}
}

int RemoveBox(vector<StrandedFragment> &fragmentSet, int b[]) {
	int i, f;
	i = 0; f = 0;
	int nRemoved = 0;
	for (f = 0; f < fragmentSet.size(); f++) {
		if (!(fragmentSet[f].x >= b[0] and fragmentSet[f].x < b[1] and 
					fragmentSet[f].y >= b[2] and fragmentSet[f].y < b[3])) {
			fragmentSet[i]= fragmentSet[f];
			i++;
		}
		else {
			++nRemoved;
		}
	}
	fragmentSet.resize(i);
	return nRemoved;
}


int TestFragmentSet(vector<StrandedFragment> &fragmentSet, 
										int tupleSize, 
										int queryLength,
										vector<StrandedFragment> & rcChain) {
	int chainBegin, chainEnd, chainYBegin, chainYEnd;
	int box[4];

	SwapStrand(rcChain, queryLength);
	RotateChain(rcChain, tupleSize);
	GetChainBox(rcChain, tupleSize, box);
	int removed;
	removed = RemoveBox(fragmentSet, box);
	vector<int> maxFragmentChain;
	fragmentSet.insert(fragmentSet.begin(), rcChain.begin(), rcChain.end());
	ClearFragmentSet(fragmentSet);
	std::sort(fragmentSet.begin(), fragmentSet.end(), StrandedLexicographicFragmentSort<StrandedFragment>());	
	RemoveDuplicates(fragmentSet);
	SDPLongestCommonSubsequence(queryLength, fragmentSet, tupleSize, 5, 5, -5, maxFragmentChain, Global);	
	return maxFragmentChain.size();
}

int InversionAlign(FASTASequence &query, FASTASequence &target, int wordSize, int invCoords[], bool debugDotPlot=false) {

	//
	// Buffers that are useful to keep around.
	//
	vector<StrandedFragment> fragmentSet;
	vector<StrandedFragment> rcFragmentSet;

	TupleList<PositionDNATuple> targetTupleList;
	std::vector<int> maxFragmentChain, rcMaxFragmentChain, maxFragmentChainWithRC;

	/*
		Collect a set of matching fragments between query and target.
		Since this function is an inner-loop for alignment, anything to
		speed it up will help.  One way to speed it up is to re-use the
		vectors that contain the sdp matches. 
	*/

	TupleMetrics tm, tmSmall;
	tm.Initialize(wordSize);

	SequenceToTupleList(target, tm, targetTupleList);

  targetTupleList.Sort();

			 
  //
  // Store in fragmentSet the tuples that match between the target
  // and query.
  //

	StoreMatchingPositions(query, tm, targetTupleList, fragmentSet, 1); 
	
	FASTASequence queryRC;
	query.MakeRC(queryRC);
	StoreMatchingPositions(queryRC, tm, targetTupleList, rcFragmentSet, 1); 	

	sort(fragmentSet.begin(), fragmentSet.end());
	sort(rcFragmentSet.begin(), rcFragmentSet.end());
	int nRem, nRCRem;
	nRem = RemoveSparseMatches(fragmentSet);
	nRCRem = RemoveSparseMatches(rcFragmentSet);
	if (rcFragmentSet.size() < 10) {
		queryRC.Free();
		return 0;
	}

  // 
  // The method to store matching positions is not weight aware.
  // Store the weight here.
  //
	VectorIndex f;  
	for (f = 0; f < fragmentSet.size(); f++) {
		fragmentSet[f].weight = tm.tupleSize;
		fragmentSet[f].length = tm.tupleSize;
		fragmentSet[f].strand = 0;
  }

	for (f = 0; f < rcFragmentSet.size(); f++) {
		
		rcFragmentSet[f].weight = tm.tupleSize;
		rcFragmentSet[f].length = tm.tupleSize;
		rcFragmentSet[f].strand = 1;
		// 
		// The reverse complement fragments are generated with respect to the 5'->3'
		// of the reverse strand.  In order to find the optimal alignment that the
		// forward strand incorporates alignments from the reverse, the coordinates
		// must be swapped.
		//
	}

	std::sort(fragmentSet.begin(), fragmentSet.end(), StrandedLexicographicFragmentSort<StrandedFragment>());


	
  FlatMatrix2D<int> graphScoreMat;
  FlatMatrix2D<Arrow> graphPathMat;
  FlatMatrix2D<int> graphBins;
	
	//
	// If there are way too many fragments, first overlay with graph paper.
	//

	int nOnOpt = fragmentSet.size();

  if (fragmentSet.size() > 100000) {
		int nCol = 50;
		vector<bool> onOptPath(fragmentSet.size(), false);
    nOnOpt = GraphPaper<StrandedFragment>(fragmentSet, nCol, nCol,
																						graphBins, graphScoreMat, graphPathMat,
																						onOptPath);
		int prev = fragmentSet.size();
		RemoveOffOpt(fragmentSet, onOptPath);
  }    
	
	graphScoreMat.Clear();
	graphPathMat.Clear();
	graphBins.Clear();
	
	if (fragmentSet.size() == 0) {
		//
		// This requires at least one seeded tuple to begin an alignment.
		//
		return 0;
	}

	if (debugDotPlot) {
		ofstream dotplot("dotplot.0.tab");
		for (f = 0; f < fragmentSet.size(); f++) {
			dotplot << fragmentSet[f].x << "\t" << fragmentSet[f].y << "\t" << 0 << endl;
		}
		for (f = 0; f < rcFragmentSet.size(); f++) {
			dotplot << query.length - rcFragmentSet[f].x << "\t" << rcFragmentSet[f].y << "\t" << 1 << endl;
		}

	}

	vector<vector<StrandedFragment> > rcChains;
  CreateCandidateInversionChains(rcFragmentSet, tm.tupleSize, query.length, 20, rcChains);
	
	if (debugDotPlot) {

		ofstream dotplot("candidate_boxex.tab");
		int c;
		for (c = 0; c < rcChains.size(); c++) {
			vector<StrandedFragment> tmp;
			tmp = rcChains[c];
			SwapStrand(tmp, query.length);
			int l = tmp.size() - 1;
			dotplot << tmp[0].x << " " << tmp[0].y << " " << tmp[l].x << " " << tmp[l].y << " " << c << endl;
		}
		dotplot.close();
	}
		
	/*
	std::sort(rcFragmentSet.begin(), rcFragmentSet.end(), StrandedLexicographicFragmentSort<StrandedFragment>());	
	RemoveDuplicates(rcFragmentSet);
	*/
	int forwardScore;
	SDPLongestCommonSubsequence(query.length, fragmentSet, tm.tupleSize, 2, 2, -5, maxFragmentChain, Global);	
	forwardScore = maxFragmentChain.size();
	int c;
	int optFragmentSetScore = 0;
	int maxChain = -1;
	for (c = 0; c < rcChains.size(); c++) {
		int invCandidateScore = 0;
		int b[4];
		GetChainBox(rcChains[c], tm.tupleSize, b);
		vector<StrandedFragment> fragmentSetCopy;
		fragmentSetCopy = fragmentSet;
		invCandidateScore = TestFragmentSet(fragmentSetCopy, tm.tupleSize, query.length, rcChains[c]);
		if (debugDotPlot) {
			stringstream namestrm;
			namestrm << "dotplot.candidate." << c << ".tab";
			ofstream dotplot(namestrm.str().c_str());
			for (f = 0; f < fragmentSetCopy.size(); f++) {
				dotplot << fragmentSetCopy[f].x << "\t" << fragmentSetCopy[f].y << "\t" << 0 << endl;
			}
			dotplot.close();

			stringstream boxstrm;
			namestrm << "box.candidate." << c << ".tab";
			ofstream boxplot(boxstrm.str().c_str());
			boxplot << b[0] << " " << b[1]  << " " << b[2] << " " << b[3] << endl;
			boxplot.close();

		}
		if (invCandidateScore > optFragmentSetScore and invCandidateScore - 20 > forwardScore) {
			maxChain = c;
			optFragmentSetScore = invCandidateScore;
		}
	}

	//
	// Now turn the max fragment chain into real a real alignment.
	//
	if (rcChains.size() > 0 and maxChain >= 0) {
		GetChainBox(rcChains[maxChain], tm.tupleSize, invCoords);
		return rcChains[maxChain].size();
	}
	else {
		return 0;
	}
}



#endif

