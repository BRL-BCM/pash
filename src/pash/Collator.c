/*
Copyright (c) 2004-2008 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/

/***************************************************************************
* Collator.c
*
* DESCRIPTION:
* Collates bindump files (the result of hash table inversion) to identify collinear groups of hits that exceed given scoring thresholds
*
* Records within
* a given section are sorted by horizontal HBE.  This allows the collator
* to scan through the files in a single pass and generate the collated
* results in linear time (linear in the size of the input and number of hits).
*
* AUTHOR:
* Cristian Coarfa (coarfa@bcm.edu)
* Ken Kalafus (kkalafus@bcm.tmc.edu)
*
* April 2008 - added whole-chunk/whole-read collation
*            - does not require bindumps anymore
*            - works with hive hash
* March 2007 - added banded collation and scoring
* 		       - account for gaps
* 		       - change the internal representation of matches
* 		       - make sure most of the functions called by the collator are inlined
***************************************************************************/

/**************************************************
* INCLUDES
**************************************************/
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "PashLib.h"
#include "Collator.h"
#include "PashDebug.h"
#include "FixedHashKey.h"
#include "IgnoreList.h"
#include "SAMInfo.h"

#define DEB_PROGRESS               0
#define DEB_SCAN_HORIZONTAL_SEQ		 0
#define DEB_PERFORM_COLLATION 		 0
#define DEB_SET_MATCH_STREAM		   0
#define DEB_RESET_COLLATOR		     0
#define DEB_INSERT_MATCH_STREAM		 0
#define DEB_HEAP_DELETE 					 0
#define DEB_DIAGONAL_COVERAGE		   0
#define DEB_COLL_HEURISTIC_1		   0
#define DEB_KMER_SW			           0
#define DEB_ADVANCE_STREAM         0
#define DEB_SW_CANDIDATES          0
#define DEB_BEST_ANCHORING         0
#define DEB_WEAK_CANDIDATES        0
#define DEB_HWIN    							 0
#define DEB_BSW 									 0
#define DEB_FILTER_OUTPUT					 0
#define DEB_REP_READS   					 0
#define DEB_TRACEBACK 						 0
#define DEB_EVALUATE_ANCHORING     0
#define DEB_BAD_ANCHORING          0
#define DEB_FAIL_SW  							 0
#define DEB_VARIANTS               0
#define DEB_MiX_2S_AND_REG_PASH    0
#define DEB_TOPPERCENT             0
#define DEB_CHH_VARIANTS           0

typedef struct {
  int horizontalStart;
  int horizontalStop;
  int verticalStart;
  int verticalStop;
  int numMatches;
  int numMismatches;
  int numVerticalGaps;
  int numHorizontalGaps;
  int numGapBases;
  int numBlocks;
  unsigned int blockSizes[MAX_NUMBER_BLOCKS];
  unsigned int horizontalBlockStarts[MAX_NUMBER_BLOCKS];
  unsigned int verticalBlockStarts[MAX_NUMBER_BLOCKS];
} AlignmentSummary;

typedef enum {InBlock, InGap} TraceStatus;

unsigned long swCalls;
unsigned long kswCalls;
unsigned long failedSWCalls;
unsigned long reallyPoorAnchorings;
unsigned long predSkelScore, tSkelScore;

inline int bandedSW(int *scoringMatrix, char* verticalSequence, char *horizontalSequence, int sizeVerticalSequence, int band);
int bandedSWAlignmentInfo(int *scoringMatrix,
													char* verticalSequence, char *horizontalSequence,
													int sizeVerticalSequence, int band,
													float targetScore, AlignmentSummary* alignmentSummary);

int bandedSWAlignmentInfoBisulfiteSeq(int *scoringMatrix,
													char* verticalSequence, char *horizontalSequence,
													int sizeVerticalSequence, int band,
													float targetScore, AlignmentSummary* alignmentSummary);


void resetSamInfo(SAMInfo* samInfo);
int callVariants(char* readSequence, char* templateSequence,
								 SAMInfo *samInfo, AlignmentSummary* alignmentSummary, int bisulfiteSequencing);
int outputRegularPashLine(guint32 sequenceId, int swScore, SequenceInfo* sequenceInfo,
                          guint32 chromLength, char* currentSequence,
                          guint32 alignmentHorizontalStart,
                          AlignmentSummary* alignmentSummary,
                          char strand, char *outputLine, SAMInfo* samInfo, const char* qualityScores);
int outputBisulfiteMappingLine(guint32 sequenceId, int swScore, SequenceInfo* sequenceInfo,
                          guint32 chromLength, char* currentSequence,
                          guint32 alignmentHorizontalStart,
                          AlignmentSummary* alignmentSummary,
                          char *outputLine, SAMInfo* samInfo, const char* qualityScores,
                          PashParameters* pp);


static inline char revComplementQuick(char base)
{
	if (base>='a') {
		base = base-32; // ('a'-'A');
	}
	if (base>'G') {
		if (base=='T') {
			return 'A';
		} else { // N
			return 'N';
		}
	} else if (base<'G') {
		if (base=='A') {
			return 'T';
		} else { // C
			return 'G';
		}
	} else { // G
		return 'C';
	}
}


int checkHeap(MatchStream** matchStreams, int numStreams) {

/*
 	int i;


	for(i=1; i<=numStreams/2; i++) {

		int succ;
		succ = 2*i;

		if (succ<=numStreams && !STREAM_LESS_THAN2(matchStreams[i], matchStreams[succ])) {
			fprintf(stderr, "!!!! [%d]	%d %d %d >  %d < %x %x %x\t",
						i,
						matchStreams[i]->verticalSeqID,
						matchStreams[i]->verticalOffset,
						matchStreams[i]->horizontalOffset,
						matchStreams[i]->okey,
						matchStreams[i]->verticalOffset,
						matchStreams[i]->horizontalOffset,
						matchStreams[i]->okey);
			fprintf(stderr, "	%d %d %d > %d < %x %x %x\n",
						matchStreams[succ]->verticalSeqID,
						matchStreams[succ]->verticalOffset,
						matchStreams[succ]->horizontalOffset,
						matchStreams[succ]->okey,
						matchStreams[succ]->verticalOffset,
						matchStreams[succ]->horizontalOffset,
						matchStreams[succ]->okey);
		fflush(stderr);
		return 0;
		}
		succ = 2*i+1;
		if (succ<=numStreams && !STREAM_LESS_THAN2(matchStreams[i], matchStreams[succ])) {
			fprintf(stderr, "!!!! [%d]	%d %d %d >  %d < %x %x %x\t",
						i,
						matchStreams[i]->verticalSeqID,
						matchStreams[i]->verticalOffset,
						matchStreams[i]->horizontalOffset,
						matchStreams[i]->okey,
						matchStreams[i]->verticalOffset,
						matchStreams[i]->horizontalOffset,
						matchStreams[i]->okey);
			fprintf(stderr, "	%d %d %d > %d < %x %x %x\n",
						matchStreams[succ]->verticalSeqID,
						matchStreams[succ]->verticalOffset,
						matchStreams[succ]->horizontalOffset,
						matchStreams[succ]->okey,
						matchStreams[succ]->verticalOffset,
						matchStreams[succ]->horizontalOffset,
						matchStreams[succ]->okey);
			fflush(stderr);
			return 0;
		}
	}
*/
	return 1;
}

void dumpHeap(MatchStream** matchStreams, int numStreams) {
	int i;
	fprintf(stderr, "HEAP ");
	for (i=0; i<=numStreams; i++) {
		fprintf(stderr, "[%d] %p ",
						i, matchStreams[i]);
		fflush(stderr);
		fprintf(stderr, "	%d %d %d >  %d <\t",
						matchStreams[i]->verticalSeqID,
						matchStreams[i]->okey >> 16,
						matchStreams[i]->okey & 0x0000ffff,
						matchStreams[i]->okey);
	}
	fprintf(stderr, "\n");
}

int insertMatchStream(MatchStream** matchStreams, int numStreams, MatchStream *oneMatchStream) {
	int i;
	i = numStreams+1;
	int hole = i;
        int pred = hole >> 1;
        
	while (STREAM_LESS_THAN(oneMatchStream, matchStreams[pred]) ) {
		matchStreams[hole] = matchStreams[pred];
		hole = pred;
		pred >>= 1;
    //int test;
    //test = hole/2;
	  //xDEBUG(DEB_INSERT_MATCH_STREAM, fprintf(stderr, "[%d] %d %d %d\n", i, pred, hole, test));
	}

	matchStreams[hole] = oneMatchStream;
	xDEBUG(DEB_INSERT_MATCH_STREAM, { fprintf(stderr, "inserted match stream at %d: kmer %d: (vseq %d, voff %d, hoff %d, diag %d, okey %d %x)\n",
								 i, 0, matchStreams[hole]->verticalSeqID,
								 matchStreams[hole]->okey >> 16,
								 matchStreams[hole]-> okey & 0x0000ffff,
								 (matchStreams[hole]->okey & 0x0000ffff) -(matchStreams[hole]->okey >> 16),
								 matchStreams[hole]->okey, matchStreams[hole]->okey);
	});
/*
	if (!checkHeap(matchStreams, numStreams)) {
		dumpHeap(matchStreams, numStreams);
		assert(0);
	}
*/
	return 0;
}


int deleteHeapMin(MatchStream** matchStreams, int numStreams) {
  int hole = 1;
  int succ = 2;
	xDEBUG(DEB_HEAP_DELETE, fprintf(stderr, "start heap delete\n"); dumpHeap(matchStreams, numStreams));
  while (succ < numStreams) {
    if (STREAM_LESS_THAN(matchStreams[succ+1],matchStreams[succ])) {
      succ++;
    }
    matchStreams[hole] = matchStreams[succ];
    hole = succ;
    succ <<= 1;
		xDEBUG(DEB_HEAP_DELETE, fprintf(stderr, "heap delete step\n"); dumpHeap(matchStreams, numStreams));
  }

  // bubble up rightmost element
  int pred = hole >> 1;
  while (STREAM_LESS_THAN(matchStreams[numStreams], matchStreams[pred])) { // must terminate since min at root
    matchStreams[hole] = matchStreams[pred];
    xDEBUG(DEB_HEAP_DELETE, fprintf(stderr, "set hole %d to pred %d\n", hole, pred); dumpHeap(matchStreams, numStreams));
		hole = pred;
		pred >>= 1;
  }

  // finally move data to hole
	matchStreams[hole] = matchStreams[numStreams];
  matchStreams[numStreams] = matchStreams[numStreams+1];
	//verticalSeqID = MAX_VERT_ID; // mark as deleted
	xDEBUG(DEB_HEAP_DELETE, fprintf(stderr, "end delete mine\n"); dumpHeap(matchStreams, numStreams));
/*
	if (!checkHeap(matchStreams, numStreams-1)) {
		dumpHeap(matchStreams, numStreams-1);
		assert(0);
	}
*/
	return 0;
}

/** Setup a match stream and insert it into the priority queue.
 * @param matchStream - current match stream
 * @param matchStreamPriorityQueue priority queue of match streams
 * @param numStreams number of streams in the priority queue
 * @return 0 success, 1 otherwise
 */
inline int setMatchStreamAndInsertInQueue(MatchStream* matchStream, MatchStream** matchStreamPriorityQueue,
																					int numStreams,
																					guint32 kmer, guint32 horizontalOffset) {
	matchStream->verticalSeqID = (*matchStream).intListRunner.list[0];
	guint32 verticalOffset = (*matchStream).intListRunner.list[1];
	/*gint64 vkey = (gint64)matchStream->verticalSeqID;
	gint64 vkey16 = vkey << 16;
	gint64 vkey24 = vkey << 24;
	gint64 vkey32 = vkey << 32;

	xDEBUG(DEB_SET_MATCH_STREAM,
				 fprintf(stderr, "vseq %u, vkey=%u %x vkey16=%u %x vkey24=%u %x vkey32 %u %x\n",
								 matchStream->verticalSeqID, vkey, vkey, vkey16, vkey16, vkey32, vkey32));


	guint64 key = (((guint64)matchStream->verticalSeqID)<< 32)|
								(((guint64)matchStream->verticalOffset) <<16)|
							     ((guint64) matchStream->horizontalOffset);
	*/
	matchStream->okey = (((guint32)verticalOffset)<<16)|
											 ((guint32)horizontalOffset);
	/*
	guint16 hOff = (guint16)(matchStream->okey& 0x0000ffff);
	guint16 vOff = (guint16)(matchStream->okey >> 16);
	if (vOff != matchStream->verticalOffset || hOff !=  matchStream->horizontalOffset) {

	 fprintf(stderr, "WRONG match stream for kmer %d: (vseq %d, voff %d, hoff %d, okey=%d, "
								 "\n vseq %x, voff %x, hoff %x, okey=%d okey=%x , diag %d) %d %x %d %x\n",
								 0, matchStream->verticalSeqID, matchStream->verticalOffset,
								 matchStream->horizontalOffset,
								 matchStream->okey,
								 matchStream->verticalSeqID, matchStream->verticalOffset,
								 matchStream->horizontalOffset,
								 matchStream->okey, matchStream->okey,
								 matchStream->horizontalOffset-matchStream->verticalOffset, vOff, vOff, hOff, hOff);
		fflush(stderr);
		assert(0);
	}
	*/
	xDEBUG(DEB_SET_MATCH_STREAM,
				 fprintf(stderr, "set up match stream for kmer %d: (vseq %d, voff %d, hoff %d, okey=%d, "
								 "\n vseq %x, voff %x, hoff %x, okey=%d okey=%x , diag %d)\n",
								 0, matchStream->verticalSeqID,
								 matchStream->okey >> 16,
								 matchStream->okey & 0x0000ffff,
								 matchStream->okey,
								 matchStream->verticalSeqID,
								 matchStream->okey >> 16,
								 matchStream->okey & 0x0000ffff,
								 matchStream->okey, matchStream->okey,
								 (matchStream->okey >> 16) -
								 (matchStream->okey & 0x0000ffff))
				 );
	insertMatchStream(matchStreamPriorityQueue, numStreams, matchStream);
	xDEBUG(DEB_SET_MATCH_STREAM, dumpHeap(matchStreamPriorityQueue, numStreams+1));
	return 0;
}

/** Advance the top match stream and restore the priority queue.
@param matchStreamPriorityQueue priority queue of match streams
@param numStreams number of match streams
@return 0 if the match stream could be advanced, 1 if the match stream reached its limit
*/
int advanceTopMatchStream(MatchStream** matchStreamPriorityQueue, int numStreams, int numDiagonals) {
	xDEBUG(DEB_ADVANCE_STREAM, fprintf(stderr, "START advanceTopMatchStream:  "));
 	MatchStream* topStream = matchStreamPriorityQueue[1];
	guint32 left = topStream->intListRunner.left - 2;
	xDEBUG(DEB_ADVANCE_STREAM, fprintf(stderr, "trying to advance a runner with %d left\n", left));
	guint32 horizontalOffset = topStream->okey & 0x0000ffff;
	while (left>0) {
		topStream->intListRunner.list +=2;

		guint32 verticalOffset = topStream->intListRunner.list[1];
    if ((horizontalOffset>verticalOffset+numDiagonals) || (verticalOffset>numDiagonals) ) {
		//if ((verticalOffset>numDiagonals)) {
			left -= 2;
			continue;
		} else {
			topStream->verticalSeqID  = topStream->intListRunner.list[0];
			topStream->intListRunner.left = left;
			topStream->okey = (((guint32) verticalOffset)<<16)|
												 ((guint32) horizontalOffset);
			deleteHeapMin(matchStreamPriorityQueue, numStreams);
			xDEBUG(DEB_ADVANCE_STREAM, dumpHeap(matchStreamPriorityQueue, numStreams);
					 fprintf(stderr, "deleted  min element"));
			insertMatchStream(matchStreamPriorityQueue, numStreams-1, topStream);
			xDEBUG(DEB_ADVANCE_STREAM, dumpHeap(matchStreamPriorityQueue, numStreams); fprintf(stderr, "RETURN 0 from advanceTopMatchStream\n"));
			return 0;
		}
	}

	xDEBUG(DEB_ADVANCE_STREAM, dumpHeap(matchStreamPriorityQueue, numStreams); fprintf(stderr, "will RETURN 1 from advanceTopMatchStream\n"));
	deleteHeapMin(matchStreamPriorityQueue, numStreams);
	xDEBUG(DEB_ADVANCE_STREAM, dumpHeap(matchStreamPriorityQueue, numStreams); fprintf(stderr, "RETURN 1 from advanceTopMatchStream\n"));
	xDEBUG(DEB_ADVANCE_STREAM,  fprintf(stderr, "STOP advanceTopMatchStream"));
	return 1;
}

inline CollatorControl* initCollatorControl(int numberOfDiagonals) {
	CollatorControl* c = (CollatorControl*) malloc(sizeof(CollatorControl));
	xDieIfNULL(c, fprintf(stderr, "could not allocate memory for the collator control "
																"in %s:%d\n", __FILE__, __LINE__ ), 1);
	c->matchStreams = (MatchStream*) malloc(sizeof(MatchStream)*MAX_MATCH_STREAMS);
	c->matchStreamPtrs = (MatchStream**) malloc(sizeof(MatchStream*)*MAX_MATCH_STREAMS);
	xDieIfNULL(c->matchStreams, fprintf(stderr, "could not allocate memory for the collator control match streams"
																							"in %s:%d\n", __FILE__, __LINE__ ), 1);
	c->maxMatchStreams	= MAX_MATCH_STREAMS;
	c->validMatchStreams = 0;
	c->numberOfDiagonals = numberOfDiagonals;
	c->matchPairsCapacity = 3*MAX_MATCH_STREAMS;
	c->matchPairs = (MatchPair*) malloc(sizeof(MatchPair)*c->matchPairsCapacity);
	xDieIfNULL(c->matchPairs, fprintf(stderr, "could not allocate memory for match pairs at %s:%d\n",
																					__FILE__, __LINE__), 1);
	c->bswMemory = (int*) malloc(MAX_READ_SIZE*(30+3*DEFAULT_BAND)*sizeof(int));
	return c;
}

/** Setup the match stream for a horizontal kmer, querying the hive hash.
@param c collator control data structure
@param kmer horizontal kmer
@parma hh HiveHash
@param horizontalOffset horizontal offset of the kmer in the current horizontal window
*/
inline void addMatchStreamCollatorControl(CollatorControl *c, guint32 kmer, HiveHash* hh, guint32 horizontalOffset) {

	int i = c->validMatchStreams;
/*
	if (i>MAX_MATCH_STREAMS) {
		fprintf(stderr, "number of streams exceeds hard limit of %d at %s:%d\n",
						MAX_MATCH_STREAMS, __FILE__, __LINE__);
		exit(1);
	}
*/

	MatchStream* matchStream = &c->matchStreams[i+1];
	hh->getIntListRunner(kmer, &matchStream->intListRunner);
	if (matchStream->intListRunner.left > 0) {
		setMatchStreamAndInsertInQueue(matchStream, c->matchStreamPtrs, c->validMatchStreams, kmer, horizontalOffset);
		c->validMatchStreams++;
	}
}

void freeCollatorControl(CollatorControl* c) {
	free(c->matchStreams);
	c->matchStreams = NULL;
	free(c);
}

inline void resetCollatorControl(CollatorControl *c) {
	xDEBUG(DEB_RESET_COLLATOR, fprintf(stderr, "resetting collator\n"));
	c->validMatchStreams = 0;
	int i;
        c->matchStreams[0].verticalSeqID  = 0;
	for (i=1; i<=c->numberOfDiagonals+1; i++) {
	  c->matchStreams[i].okey  = G_MAXUINT16 <<16;
		c->matchStreams[i].verticalSeqID  = MAX_VERT_ID;
		c->matchStreamPtrs[i] = & c->matchStreams[i];
	}
	c->matchStreamPtrs[0] = & c->matchStreams[0];
}

static inline void advanceToNextHorizontalSequence(SequenceHash* sequenceHash) {
  sequenceHash->numberOfChunksInCurrentSequence = 0;
  sequenceHash->offsetOfSequenceBufferInRealSequence=0;
  sequenceHash->currentSequenceChunk = 0;
}

/// Scan the horizontal sequence (typically chromosome/genome) agains the hivehash.
int scanHorizontalSequence(PashParameters* pp, SequenceHash* sequenceHash) {
	FastaUtil* fastaUtilHorizontal=pp->fastaUtilHorizontal;
	int i;
	guint32 currentHashSize;
  swCalls=0;
  kswCalls=0;
  failedSWCalls=0;
	reallyPoorAnchorings=0;
	predSkelScore = 0;
	tSkelScore = 0;
	double withinTopPercent = 1 -pp->topPercent;
  guint32 currentForwardChunkStart, currentForwardChunkStop;
  guint32 sequenceLength;
  char currentKmer [MAX_MASK_LEN+1];
  int maskLen, maskWeight, kmerPos, maskPos;
  int startOffset, maxOffset, minOffset, offsetGap;
  Mask mask = pp->mask;
  int currentSequencePos;
  int numberOfDiagonals;
  int minPosition;
  int sequenceToKeep;
  int diagonalChunkSize;
  int radiusChunkStart, radiusChunkStop; // we want the sequence around a certain radius for potential alignment step
	int targetChunkStart, targetChunkStop;
	CollatorControl *cc;
	HiveHash*hiveHash=(HiveHash*)sequenceHash->hiveHash;
	guint32** hashSkeleton=hiveHash->hashSkeleton;
  guint32 forwardKey, reverseKey;
  guint32 memoryTarget;
  maskLen =pp->mask.maskLen;
  maskWeight = pp->mask.keyLen;
  currentKmer[maskWeight] = '\0';
  offsetGap = 1;
  numberOfDiagonals = pp->numberOfDiagonals;
  diagonalChunkSize = 2*numberOfDiagonals;
	printNow();

	char tmpOutputFileName[MAX_FILE_NAME_SIZE];
	sprintf(tmpOutputFileName, "%s.tmp.%d", pp->outputFile, getpid());
  FILE* tmpOutputFilePtr = fopen(tmpOutputFileName, "wt");
  if (tmpOutputFilePtr==NULL) {
    fprintf(stderr, "could not open temporary output file %s\n", tmpOutputFileName);
    fflush(stderr);
    exit(2);
  }
	xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "starting horizontal scanning\n"));
	rewindFastaUtil(fastaUtilHorizontal);
	int useIgnoreList = pp->useIgnoreList;
	IgnoreList ignoreList = pp->ignoreList;
	char currentSequence[MAX_FILE_NAME_SIZE+1];
	cc = initCollatorControl(numberOfDiagonals);
  // if at limit of memory, then stop, because we have enough info to resume the vert hash filling
  while(!fastaUtilHorizontal->parsingDone) {
    if (sequenceHash->numberOfChunksInCurrentSequence==0 ||
        (sequenceHash->numberOfChunksInCurrentSequence>0 &&
         sequenceHash->currentSequenceChunk>=sequenceHash->numberOfChunksInCurrentSequence)) {
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "getting a new sequence: bufferSeqPos=%d\n",
                                            fastaUtilHorizontal->currentSequenceBufferPos));
      // need to read in next sequence
      fastaUtilKeepPartialBuffer(fastaUtilHorizontal, 0);
      nextChunkFastaUtil(fastaUtilHorizontal);
      if (fastaUtilHorizontal->currentSequenceBufferPos<=0) {
        advanceToNextHorizontalSequence(sequenceHash);
        continue;
      }
      sequenceLength = fastaUtilHorizontal->sequencesInformation[fastaUtilHorizontal->currentSequenceIndex].sequenceLength;
      // print current sequence length
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ,
             fprintf(stderr, "got current sequence [%d] %s of length %d: %d characters\n",
                    fastaUtilHorizontal->currentSequenceIndex,
                    fastaUtilHorizontal->sequencesInformation[fastaUtilHorizontal->currentSequenceIndex].sequenceName,
                    sequenceLength,
                    fastaUtilHorizontal->currentSequenceBufferPos));
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, {
        fastaUtilHorizontal->deflineBuffer[fastaUtilHorizontal->currentDeflineBufferPos]='\0';
        fastaUtilHorizontal->sequenceBuffer[fastaUtilHorizontal->currentSequenceBufferPos]='\0';
        fprintf(stderr, "current fasta util state\n");
             fprintf(stderr, "current fasta file [%d]->%s\n",
                     fastaUtilHorizontal->currentFileIndex,
                     fastaUtilHorizontal->fileArray[fastaUtilHorizontal->currentFileIndex]);
             fprintf(stderr, "current line %d\n", fastaUtilHorizontal->currentLine);
             fprintf(stderr, "def: %s \n seq: %s\n",
                     fastaUtilHorizontal->deflineBuffer, fastaUtilHorizontal->sequenceBuffer);
      }	);
      
      if (pp->bisulfiteSequencingMapping) {
        if (strstr(fastaUtilHorizontal->deflineBuffer, "#RC.pash.")!=NULL) {
          pp->reverseStrandDnaMethMapping=1;
          strcpy(pp->actualChromName, &fastaUtilHorizontal->deflineBuffer[strlen("#RC.pash.")]);
          fprintf(stderr, "def: %s ; reverse complement strand of %s\n",
              fastaUtilHorizontal->deflineBuffer, pp->actualChromName);
          pp->reverseComplementSequenceLength = sequenceLength; 
                
        } else {
                fprintf(stderr, "def: %s ; forward strand\n",
                       fastaUtilHorizontal->deflineBuffer);
          pp->reverseStrandDnaMethMapping = 0;
        }
      } else {
        fprintf(stderr, "def: %s\n",
                     fastaUtilHorizontal->deflineBuffer);
      }
			printNow();
			fflush(stderr);
			strncpy(currentSequence, fastaUtilHorizontal->deflineBuffer, MAX_FILE_NAME_SIZE);
      // set # of sequences
      sequenceHash->numberOfChunksInCurrentSequence = (sequenceLength-1)/numberOfDiagonals+1;
      sequenceHash->offsetOfSequenceBufferInRealSequence = 0;
      // set current chunk to 0
      sequenceHash->currentSequenceChunk = 0;
    }

    xDEBUG(DEB_SCAN_HORIZONTAL_SEQ,
             fprintf(stderr, "sequence status (%d, %d, %d)\n",
                     sequenceHash->numberOfChunksInCurrentSequence,
                     sequenceHash->offsetOfSequenceBufferInRealSequence,
                     sequenceHash->currentSequenceChunk));

    currentForwardChunkStart = numberOfDiagonals*sequenceHash->currentSequenceChunk;
    currentForwardChunkStop = currentForwardChunkStart+2*numberOfDiagonals-1;
    if (currentForwardChunkStop>=sequenceLength) {
      currentForwardChunkStop = sequenceLength-1;
    }
		radiusChunkStart = currentForwardChunkStart - numberOfDiagonals-DEFAULT_BAND;
		radiusChunkStop = currentForwardChunkStop + numberOfDiagonals+DEFAULT_BAND;
		targetChunkStart = radiusChunkStart;
		targetChunkStop  = radiusChunkStop;
		if (radiusChunkStart<0) {
			radiusChunkStart = 0;
		}
		if (radiusChunkStop >= sequenceLength) {
			radiusChunkStop = sequenceLength-1;
		}
    xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "f start %d, f stop %d radiusStart %d radiusStop %d\n",
                                          currentForwardChunkStart, currentForwardChunkStop,
																					radiusChunkStart, radiusChunkStop));
    if (radiusChunkStop < sequenceHash->offsetOfSequenceBufferInRealSequence +
				fastaUtilHorizontal->currentSequenceBufferPos  ) {
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "f chunk+radius is available\n"));
      startOffset =0;
      sequenceHash->lastSequenceId ++;
      maxOffset = currentForwardChunkStop - currentForwardChunkStart +1 - maskLen;
		/*	if (maxOffset >= numberOfDiagonals) {
				maxOffset = numberOfDiagonals-1;
			}
		*/
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr,"f st %d f stop %d start offset=%d maxOffset=%d\n",
                                            currentForwardChunkStart, currentForwardChunkStop, startOffset, maxOffset));
      resetCollatorControl(cc);
			for (startOffset = 0; startOffset<=maxOffset; startOffset+= offsetGap) {
        for (currentSequencePos=startOffset+currentForwardChunkStart-sequenceHash->offsetOfSequenceBufferInRealSequence,
             maskPos = 0,kmerPos = 0;
             kmerPos < maskWeight;
             currentSequencePos++, maskPos++) {
          if (mask.mask[maskPos]) {
            currentKmer [kmerPos] = fastaUtilHorizontal->sequenceBuffer[currentSequencePos];
            kmerPos++;
          }
        }
        getKeyForSeq(currentKmer, &forwardKey);
        xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "found forward kmer %s %d %x, h seq id %d, h offset %d\n",
                                              currentKmer, forwardKey, forwardKey, sequenceHash->lastSequenceId,
                                              startOffset));
        // do something
				if (forwardKey != BAD_KEY) {
					if(!useIgnoreList || !isIgnored(forwardKey, ignoreList)) {
            if(hashSkeleton[forwardKey]!=NULL) {
							addMatchStreamCollatorControl(cc, forwardKey,  hiveHash, startOffset);
						}
					}
				}
      }
			xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "do something\n"));
			if (cc->validMatchStreams>0) {
				// setup target alignment template
				xDEBUG(DEB_HWIN, fprintf(stderr,"tStart %d rStart=%d tStop = %d rStop = %d\n",
																 targetChunkStart, radiusChunkStart, radiusChunkStop, targetChunkStop));
				if (targetChunkStart==radiusChunkStart && targetChunkStop==radiusChunkStop) {
					memcpy(&cc->targetTemplate[0],
								 &fastaUtilHorizontal->sequenceBuffer[radiusChunkStart-sequenceHash->offsetOfSequenceBufferInRealSequence],
								radiusChunkStop-radiusChunkStart+1);
				} else {
					int jjj;
					for (jjj=0; jjj<=targetChunkStop-targetChunkStart; jjj++) {
						cc->targetTemplate[jjj]='@';
					}
					memcpy(&cc->targetTemplate[radiusChunkStart-targetChunkStart],
								 &fastaUtilHorizontal->sequenceBuffer[radiusChunkStart-sequenceHash->offsetOfSequenceBufferInRealSequence],
								radiusChunkStop-radiusChunkStart+1);
				}
				cc->targetTemplateStart = targetChunkStart;
                                cc->targetTemplate[targetChunkStop-targetChunkStart+1]='\0';
				performCollation(tmpOutputFilePtr, cc,  sequenceHash, pp,
												 fastaUtilHorizontal->deflineBuffer
												 , currentForwardChunkStart,
												 currentForwardChunkStop,
												 fastaUtilHorizontal->currentSequenceIndex);
			}
      // all input was consumed, move on to the next sequence
      sequenceHash->currentSequenceChunk++;
    } else {
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ,
             fprintf(stderr, "f and reverse chunk are not simultaneously available\n"));
      minPosition = radiusChunkStart;
      if (minPosition< sequenceHash->offsetOfSequenceBufferInRealSequence) {
        fprintf(stderr, "assertion failed at %s:%d\n", __FILE__, __LINE__);
        exit(1);
      }
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "minPosition: %d, f chunkStart %d base offset %d, end of buffer %d\n",
            minPosition,   currentForwardChunkStart, sequenceHash->offsetOfSequenceBufferInRealSequence,
            sequenceHash->offsetOfSequenceBufferInRealSequence+fastaUtilHorizontal->currentSequenceBufferPos));
      if (minPosition>=sequenceHash->offsetOfSequenceBufferInRealSequence+fastaUtilHorizontal->currentSequenceBufferPos) {
        // don't use any of the current sequence
        sequenceHash->offsetOfSequenceBufferInRealSequence += fastaUtilHorizontal->currentSequenceBufferPos;
        sequenceToKeep = 0;
				fprintf(stderr, "should not get here !!!!\n");
				exit(1);
      } else  {
        sequenceToKeep = fastaUtilHorizontal->currentSequenceBufferPos + sequenceHash->offsetOfSequenceBufferInRealSequence - minPosition;
        xDEBUG(DEB_SCAN_HORIZONTAL_SEQ, fprintf(stderr, "keeping %d bases\n", sequenceToKeep));
				sequenceHash->offsetOfSequenceBufferInRealSequence += fastaUtilHorizontal->currentSequenceBufferPos - sequenceToKeep;
			}

      fastaUtilKeepPartialBuffer(fastaUtilHorizontal, sequenceToKeep);
      nextChunkFastaUtil(fastaUtilHorizontal);
      fastaUtilHorizontal->sequenceBuffer[fastaUtilHorizontal->currentSequenceBufferPos]='\0';
      xDEBUG(DEB_SCAN_HORIZONTAL_SEQ,
          fprintf(stderr, "got current sequence [%d] %s of length %d: %d characters %s\n",
                  fastaUtilHorizontal->currentSequenceIndex,
                  fastaUtilHorizontal->sequencesInformation[fastaUtilHorizontal->currentSequenceIndex].sequenceName,
                  sequenceLength,
                  fastaUtilHorizontal->currentSequenceBufferPos,
                  fastaUtilHorizontal->sequenceBuffer));
    }
    xDEBUG(DEB_SCAN_HORIZONTAL_SEQ,
      fprintf(stderr, "++ sequence status (%d, %d, %d)\n",
              sequenceHash->numberOfChunksInCurrentSequence,
              sequenceHash->offsetOfSequenceBufferInRealSequence,
              sequenceHash->currentSequenceChunk));
	}
	fclose(tmpOutputFilePtr);
	
  
	filterOutput(tmpOutputFileName, pp->outputFilePtr,
               pp->verticalSequencesInfos, pp->maxMappings,
               withinTopPercent);
	
  
  
          unlink(tmpOutputFileName);
	fprintf(stderr, "anchorings %ld total sw calls %ld failed calls %ld really poor anchorings %ld predSkelScore %ld tSkelScore %ld\n",
					kswCalls,  swCalls, failedSWCalls, reallyPoorAnchorings,
					predSkelScore, tSkelScore);
	return 0;
}

void performCollation(FILE* outputFilePtr, CollatorControl* c, SequenceHash* sequenceHash,
											PashParameters *pp, char* currentSequence,
											guint32 start, guint32 stop, guint32 currentChrom) {
	guint32 currentVerticalSequenceId;
	guint32 numberOfMatchPairs;
	guint32 chromLength = pp->fastaUtilHorizontal->sequencesInformation[currentChrom].sequenceLength;
	int iii, diagIdx;
	char currentKmer[MAX_MASK_LEN];
	int minDiagonal, maxDiagonal, currentDiagonalCoverage, maxDiagonalCoverage;
	int currentDiagonal, currentMaxDiagonal;
	int numberOfDiagonals = c->numberOfDiagonals;
	int range = 3; // this will become a more adustable constant later
	if (range<numberOfDiagonals) {
		range = numberOfDiagonals;
	}
  double withinTopPercent = 1-pp->topPercent;
	int maxMinDiagonal = numberOfDiagonals - range;
	Mask mask = pp->mask;
	int kmerWeight = mask.keyLen;
	int kmerSpan = mask.maskLen;
	int numberOfMatchRunStarts;
	int extendPreviousRunFlag;
	int matchPairIndex, matchRunStartIndex;
	MatchPair* matchPairs = c->matchPairs;
	int bestExtendScore, bestMatchPair, bestMatchRun;
	int currentExtendScore;
	int matchGain = 2;
	int gapOpenPenalty = -3;
	int gapExtendPenalty = -3;
	int matchPairsCapacity = c->matchPairsCapacity;
	int firstPreviouMatch, previousMatchIndex;
	int overlap;
	int gapBases;
	SequenceInfo *verticalSequenceInfos = pp->verticalSequencesInfos;
	static guint32 iterNo = 0;
  int bisulfiteSequencingMapping = pp->bisulfiteSequencingMapping;
	iterNo +=1;

	xDEBUG(DEB_PROGRESS, fprintf(stderr, "start collation of %s with %d valid streams %d---%d \n",
																				currentSequence, c->validMatchStreams, start, stop));

	guint32 maxReadMappings = pp->maxMappings;
	int numberMatchesCutoff = numberOfDiagonals/pp->wordOffset*5;
/*
while !done
		get all kmers with the same sequence index and sorted by vertical offset
		  collect them in MatchPairs
		determine the maximum diagonal window, depending on read size
		collate the kmers within those diagonals: rough collation, no runs
    dump blocks
    check for collation inconsistencies (e.
    g. do we have parallel runs ?)
	*/
	/* build heap of match streams*/
	while (c->validMatchStreams>0) {
		/*for (diagIdx=0; diagIdx<c->numberOfDiagonals; diagIdx++) {
			c->diagonalCoverage[diagIdx] = 0;
		}
		*/

		// pick top of priority queue
		currentVerticalSequenceId = c->matchStreamPtrs[1]->verticalSeqID;
		matchPairs[0].verticalSeqID = currentVerticalSequenceId;
		matchPairs[0].verticalOffset = c->matchStreamPtrs[1]->okey >> 16;
		matchPairs[0].horizontalOffset = c->matchStreamPtrs[1]->okey & 0x0000ffff;
		//currentDiagonal = c->matchStreamPtrs[1]->diagonal;
		//matchPairs[0].diagonal = currentDiagonal;
		matchPairs[0].diagonal = -matchPairs[0].verticalOffset+matchPairs[0].horizontalOffset;

		if (advanceTopMatchStream(c->matchStreamPtrs, c->validMatchStreams, c->numberOfDiagonals+10) == 1) {
			c->validMatchStreams --;
		}
		numberOfMatchPairs = 1;
		xDEBUG(DEB_PERFORM_COLLATION, fprintf(stderr, "initialize new match set (%d, %d, %d) \n",
																					matchPairs[0].verticalSeqID, matchPairs[0].verticalOffset,
																					matchPairs[0].horizontalOffset));
		if (c->validMatchStreams>0){
		xDEBUG(DEB_PERFORM_COLLATION, fprintf(stderr, "validMatchStreams=%d, top match (vsid %d, kmer %d, hoff %d) \n",
																					c->validMatchStreams, c->matchStreamPtrs[1]->verticalSeqID,
																					0, c->matchStreamPtrs[1]->okey >> 16));
		}
		while (c->validMatchStreams>0 && c->matchStreamPtrs[1]->verticalSeqID == currentVerticalSequenceId) {
			// todo: regrow number of match pairs dynamically
			//assert(numberOfMatchPairs<matchPairsCapacity);
			matchPairs[numberOfMatchPairs].verticalOffset = c->matchStreamPtrs[1]->okey >> 16;
			matchPairs[numberOfMatchPairs].horizontalOffset = c->matchStreamPtrs[1]-> okey & 0x0000ffff; // horizontalOffset;
			//currentDiagonal = c->matchStreamPtrs[1]->diagonal;
			currentDiagonal = matchPairs[numberOfMatchPairs].horizontalOffset - matchPairs[numberOfMatchPairs].verticalOffset;
			if (currentDiagonal>=0) {
				matchPairs[numberOfMatchPairs].verticalSeqID = currentVerticalSequenceId;
				matchPairs[numberOfMatchPairs].diagonal = currentDiagonal;
						numberOfMatchPairs ++;
				xDEBUG(DEB_COLL_HEURISTIC_1, fprintf(stderr, "[000] have %d matches to collate\n", numberOfMatchPairs));
			}
			if (advanceTopMatchStream(c->matchStreamPtrs, c->validMatchStreams, c->numberOfDiagonals+10) == 1) {
				c->validMatchStreams --;
			}
			xDEBUG(DEB_COLL_HEURISTIC_1, fprintf(stderr, "validMatchStreams=%d, top match (vsid %d, kmer %d, hoff %d) \n",
																					c->validMatchStreams, c->matchStreamPtrs[1]->verticalSeqID,
																					0, c->matchStreamPtrs[1]->okey & 0x0000ffff));
		}

		// TODO: resurrect this; mapping based on 1 seed is pointless
		xDEBUG(DEB_COLL_HEURISTIC_1, fprintf(stderr, "have %d matches to collate\n", numberOfMatchPairs));

		// determine cost of kmer sw

		if (numberOfMatchPairs == 1) {
			xDEBUG(DEB_COLL_HEURISTIC_1, fprintf(stderr, "1 match gets ignored\n"));
			continue;
		}
		// don't match sequences with a lot of self repeats
		if (numberOfMatchPairs >=numberMatchesCutoff) {
			continue;
		}

		// I think this is inefficient; one should just start w/ kmer collation; all this diagonal traversal is pointless
		// especially for the slightly longer reads
		xDEBUG(DEB_PERFORM_COLLATION, { fprintf(stderr, "current matches: [%d] ", matchPairs[0].verticalSeqID);
			currentKmer[kmerWeight] = '\0';
			for (iii = 0; iii<numberOfMatchPairs; iii++) {
				getSeqForKey1(matchPairs[iii].kmer, kmerWeight, currentKmer);
				fprintf(stderr, "\t(%d,%d,%d,%d:%s)",
				matchPairs[iii].verticalOffset, matchPairs[iii].horizontalOffset,
				matchPairs[iii].diagonal, matchPairs[iii].kmer, currentKmer); }
				fprintf(stderr, "\n");});

		// now do kmer-level alignment
		numberOfMatchRunStarts = 0;
		kswCalls +=1;
		for (matchPairIndex=0; matchPairIndex<numberOfMatchPairs; matchPairIndex++) {
			currentDiagonal = matchPairs[matchPairIndex].diagonal;
			// attempt to collate current match pair
			xDEBUG(DEB_KMER_SW,
						 fprintf(stderr, "trying to place match pair %d, (%d, %d)[%d]\n",
							matchPairIndex,
							matchPairs[matchPairIndex].verticalOffset,
							matchPairs[matchPairIndex].horizontalOffset,
							matchPairs[matchPairIndex].diagonal));
			bestMatchPair = -1;
			bestMatchRun = -1;
			bestExtendScore = 0;
			xDEBUG(DEB_KMER_SW, fprintf(stderr, "number of match runs %d\n", numberOfMatchRunStarts));
			for (matchRunStartIndex=0;
					 matchRunStartIndex<numberOfMatchRunStarts;
					 matchRunStartIndex++) {
				for(firstPreviouMatch = previousMatchIndex = c->matchRunStarts[matchRunStartIndex];previousMatchIndex>=0;) {
					xDEBUG(DEB_KMER_SW,
						 fprintf(stderr, "trying to extend with match pair %d, (%d, %d)[%d]\n",
							previousMatchIndex,
							matchPairs[previousMatchIndex].verticalOffset,
							matchPairs[previousMatchIndex].horizontalOffset,
							matchPairs[previousMatchIndex].diagonal));
					if (matchPairs[previousMatchIndex].diagonal ==
							matchPairs[matchPairIndex].diagonal) {
						// extend current run on the same diagonal
						xDEBUG(DEB_KMER_SW, fprintf(stderr, "linking with a run on the same diagonal %d\n",
																				matchPairs[matchPairIndex].diagonal));
						overlap = matchPairs[matchPairIndex].verticalOffset -
							matchPairs[previousMatchIndex].verticalOffset;
						if (overlap<0) {
							guint32 sequenceId = currentVerticalSequenceId/2;
							SequenceInfo *sequenceInfo = verticalSequenceInfos+sequenceId;

							fprintf(stderr, "current matches: [%d] %s ",
											matchPairs[0].verticalSeqID,
											sequenceInfo->sequenceName);
							currentKmer[kmerWeight] = '\0';
							for (iii = 0; iii<numberOfMatchPairs; iii++) {
								guint32 okey = (guint32)matchPairs[iii].verticalOffset << 16| ((guint32) matchPairs[iii].horizontalOffset);
								fprintf(stderr, "\t(%d,%d,%d,%x)",
								matchPairs[iii].verticalOffset, matchPairs[iii].horizontalOffset, okey, okey);
							}
							fprintf(stderr, "\n");


							dumpHeap(c->matchStreamPtrs, c->validMatchStreams+1);
							assert(overlap>0);
						}
						if (overlap > kmerSpan) {
							currentExtendScore = matchPairs[previousMatchIndex].bestRunScore + kmerWeight*matchGain;
							xDEBUG(DEB_KMER_SW, fprintf(stderr, "adding a full kmer weight %d, current score = %d\n",
																				kmerWeight, currentExtendScore));
						} else {
							currentExtendScore = matchPairs[previousMatchIndex].bestRunScore + mask.maskOverlapContribution[overlap]*matchGain;
							xDEBUG(DEB_KMER_SW,
										 fprintf(stderr, "adding a overlap of %d kmer contribution %d  matchGain=%d  current score = %d\n",
												overlap, mask.maskOverlapContribution[overlap], matchGain, currentExtendScore));

						}
						if (currentExtendScore>bestExtendScore) {
							xDEBUG(DEB_KMER_SW,
										 fprintf(stderr, "current score = %d exceeds best score of %d, position %d\n",
												currentExtendScore, bestExtendScore, previousMatchIndex));
							bestExtendScore = currentExtendScore;
							bestMatchPair = previousMatchIndex;
							bestMatchRun = matchRunStartIndex;
						}
						break;
					} else {
						// different diagonals
						if (matchPairs[previousMatchIndex].verticalOffset + kmerSpan <= matchPairs[matchPairIndex].verticalOffset) {
						  xDEBUG(DEB_KMER_SW, fprintf(stderr, "linking over a gap w/ diagonal  %d\n",
																				matchPairs[previousMatchIndex].diagonal));
							// gap start
							gapBases = matchPairs[previousMatchIndex].diagonal - matchPairs[matchPairIndex].diagonal;
							if (gapBases<0) {
								gapBases = - gapBases;
							}
							currentExtendScore = matchPairs[previousMatchIndex].bestRunScore +
								(gapBases-1)*gapExtendPenalty+gapOpenPenalty + kmerWeight*matchGain;
							if (currentExtendScore>bestExtendScore) {
								xDEBUG(DEB_KMER_SW,
											 fprintf(stderr, "current score = %d exceeds best score of %d, position %d\n",
													currentExtendScore, bestExtendScore, previousMatchIndex));
								bestExtendScore = currentExtendScore;
								bestMatchPair = previousMatchIndex;
								bestMatchRun = matchRunStartIndex;
							}
							break;
						} else {
							// tough luck, go to previous match pair in current run
							previousMatchIndex = matchPairs[previousMatchIndex].previousMatchPairInRun;
							xDEBUG(DEB_KMER_SW,fprintf(stderr, "tough luck, going to previous matchpair in run %d\n",
													previousMatchIndex));
						}
					}
				}
			}

			if (bestExtendScore > kmerWeight*matchGain) {
				matchPairs[matchPairIndex].bestRunScore = bestExtendScore;
				matchPairs[matchPairIndex].previousMatchPairInRun = bestMatchPair;

				if (c->matchRunStarts[bestMatchRun] == bestMatchPair) {
					// add a new match run start
					c->matchRunStarts[bestMatchRun] = matchPairIndex;
					c->matchRunScores[bestMatchRun] = bestExtendScore;
					xDEBUG(DEB_KMER_SW,
												 fprintf(stderr, "we extended and replaced the top of a previous run for current kmer:"
																 "best score = %d position %d previous matchpair in run %d\n",
														bestExtendScore, matchPairs[matchPairIndex].bestRunScore,
														bestMatchPair, matchPairs[matchPairIndex].previousMatchPairInRun));
				} else {
					xDEBUG(DEB_KMER_SW,
												 fprintf(stderr, "we extended a previous partial run for current kmer: "
																 "best score = %d bestRunScore=%d, best match pair %d "
																 "previous matchpair in run %d new run start %d\n",
														bestExtendScore, matchPairs[matchPairIndex].bestRunScore,
														bestMatchPair, matchPairs[matchPairIndex].previousMatchPairInRun,
														numberOfMatchRunStarts));
					c->matchRunStarts[numberOfMatchRunStarts] = matchPairIndex;
					c->matchRunScores[numberOfMatchRunStarts] = bestExtendScore;
					numberOfMatchRunStarts ++;
				}
			} else {
				matchPairs[matchPairIndex].bestRunScore = kmerWeight*matchGain;
				matchPairs[matchPairIndex].previousMatchPairInRun = -1;
				xDEBUG(DEB_KMER_SW,
											 fprintf(stderr, "could not extended a previous run for current kmer, starting a new one: best score = %d previous matchpair in run %d\n",
													matchPairs[matchPairIndex].bestRunScore,
													matchPairs[matchPairIndex].previousMatchPairInRun));
				/// START NEW RUN !!
				c->matchRunStarts[numberOfMatchRunStarts] = matchPairIndex;
				c->matchRunScores[numberOfMatchRunStarts] = kmerWeight*matchGain;
				numberOfMatchRunStarts ++;
				xDEBUG(DEB_KMER_SW,
											 fprintf(stderr, "added a match pair start %d\n",
													numberOfMatchRunStarts));
			}
		}

		// find the run with best score, backtrace it, and dump the blocks, matching bases, gap, gap bases
		guint32 bestMatchScore, cutoffScore;
		cutoffScore=12;
		bestMatchScore = 0;
		bestMatchRun = 0;
		for (matchRunStartIndex=0; matchRunStartIndex<numberOfMatchRunStarts; matchRunStartIndex++) {
			if (c->matchRunScores[matchRunStartIndex]>bestMatchScore) {
				bestMatchScore = c->matchRunScores[matchRunStartIndex];
			}
		}

		if (bestMatchScore<10) {
			continue;
		}
		guint32 sequenceId;
    if(bisulfiteSequencingMapping) {
      sequenceId=currentVerticalSequenceId;
    } else {
      sequenceId = currentVerticalSequenceId/2;
    }
		SequenceInfo *sequenceInfo = verticalSequenceInfos+sequenceId;
		guint32 bestAnchoringScore = sequenceInfo->bestAnchoringScore ;
		if (bestMatchScore >= bestAnchoringScore*3/4) {
			xDEBUG(DEB_SW_CANDIDATES, fprintf(stderr, "[%d][%s] xAC %d vsq %d\n",
																				sequenceId, sequenceInfo->sequenceName, bestMatchScore, sequenceInfo->bestAnchoringScore*3/4));
			// check if best match score exceeds 20% of best match score
			if (bestMatchScore>sequenceInfo->bestAnchoringScore) {
				sequenceInfo->bestAnchoringScore = bestMatchScore;
				xDEBUG(DEB_BEST_ANCHORING, fprintf(stderr, "[%d][%s] upgraded best anchoring score to %d\n",
																				sequenceId, sequenceInfo->sequenceName, sequenceInfo->bestAnchoringScore));
			}

			// sw
			int sequenceLength = sequenceInfo->sequenceLength;

			if (sequenceInfo->bestSWScore>=sequenceLength && sequenceInfo->bestScoreMappings>maxReadMappings) {
				xDEBUG(DEB_REP_READS, fprintf(stderr, "skip aligning read %s after %d mappings\n",
																			sequenceInfo->sequenceName, sequenceInfo->bestScoreMappings));
				continue;
			}
      if (bisulfiteSequencingMapping) {
        // copy forward read
        const char* fwdSequence = pp->verticalFastqUtil->retrieveSequence(sequenceId); 
        memcpy(&c->readTemplate[0], fwdSequence, sequenceLength+1);
      } else {
        if (currentVerticalSequenceId%2==0) {
          // copy forward read
          const char* fwdSequence = pp->verticalFastqUtil->retrieveSequence(sequenceId); 
          memcpy(&c->readTemplate[0], fwdSequence, sequenceLength+1);
        } else {
          const char* revcSequence = pp->verticalFastqUtil->retrieveRevComplementSequence(sequenceId); 
          memcpy(&c->readTemplate[0], revcSequence, sequenceLength+1);
        }
      }
			c->readTemplate[sequenceLength]='\0';
			xDEBUG(DEB_HWIN, c->readTemplate[sequenceLength]='\0'; fprintf(stderr, "read template %s\n", c->readTemplate));
			// determine banding
			long vStop, hStop, vStart, hStart;

			int iii;
			for (matchRunStartIndex=0; matchRunStartIndex<numberOfMatchRunStarts; matchRunStartIndex++) {
				if (c->matchRunScores[matchRunStartIndex]==bestMatchScore) {
					int bestMatchIndex= c->matchRunStarts[matchRunStartIndex];
					vStop = c->matchPairs[bestMatchIndex].verticalOffset;
					hStop = c->matchPairs[bestMatchIndex].horizontalOffset;
					for (matchPairIndex=bestMatchIndex; c->matchPairs[matchPairIndex].previousMatchPairInRun>=0;
							 matchPairIndex=c->matchPairs[matchPairIndex].previousMatchPairInRun);
					if (matchPairIndex==bestMatchIndex) {
						continue;
					}
					vStart = c->matchPairs[matchPairIndex].verticalOffset;
					hStart = c->matchPairs[matchPairIndex].horizontalOffset;
					xDEBUG(DEB_HWIN, fprintf(stderr, "[%ld,%ld] start %uld %uld %ld bmr %d %d %d \n",
																	start, stop,  matchPairIndex, vStart, hStart,
																	 bestMatchIndex, vStop, hStop
																	));
					xDEBUG(DEB_HWIN,fprintf(stderr, "%s onto %s %d\t%d\t%c\t\n",
											sequenceInfo->sequenceName, currentSequence,
											start+hStart+1, start+hStop+kmerSpan,
											currentVerticalSequenceId%2==0?'+':'-'));
					long stopDiagonal = hStop - vStop;
					long startDiagonal = hStart -vStart;
					long alignmentHorizontalStart,alignmentHorizontalStop;
					int band;
				  long tmp;
					//if (startDiagonal<0 || stopDiagonal<0) {
					//	continue;
					//}
					if (stopDiagonal!=startDiagonal) {
						alignmentHorizontalStart = (long)start + (long)hStart - (long)vStart;
						alignmentHorizontalStop = (long)start +(long)hStop -(long)vStop;
						if (alignmentHorizontalStop<alignmentHorizontalStart) {
							tmp = alignmentHorizontalStart;
							alignmentHorizontalStart = alignmentHorizontalStop;
							alignmentHorizontalStop = tmp;
						}
						alignmentHorizontalStart -= DEFAULT_BAND;
						alignmentHorizontalStop  += DEFAULT_BAND;
						band = alignmentHorizontalStop-alignmentHorizontalStart+2*DEFAULT_BAND+1;
					} else {
						alignmentHorizontalStart = (long)start +(long)hStart -(long)vStart-DEFAULT_BAND;
						band = DEFAULT_BAND+DEFAULT_BAND+1;
					}
					if (band<20) {
						if (alignmentHorizontalStart<c->targetTemplateStart) {
							xDEBUG(DEB_HWIN, printf("%d vs %d\n. ", alignmentHorizontalStart, c->targetTemplateStart));
							assert(0);
						}
						xDEBUG(DEB_HWIN, fprintf(stderr, "using band %d, h start=%d, targetStart=%d, index=%d \n",
																		 band, alignmentHorizontalStart, c->targetTemplateStart,
																		 alignmentHorizontalStart-c->targetTemplateStart));
						/*int swScore = bandedSW(c->bswMemory,
										 c->readTemplate,
										 &c->targetTemplate[alignmentHorizontalStart-c->targetTemplateStart],
										 sequenceLength,
										 band);
						*/
						
						
						xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, "hwin anchoring \n%s\n%s\n",
																									 c->readTemplate, &c->targetTemplate[(long)start +(long)hStart-(long)vStart-c->targetTemplateStart]));
						xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, "start evaluation of anchoring \n%s\n%s\n",
																									 c->readTemplate, &c->targetTemplate[(long)start -c->targetTemplateStart]));
						int anchoringMatches=0;
						int anchoringMismatches=0;
						bestMatchIndex= c->matchRunStarts[matchRunStartIndex];
						for (matchPairIndex=bestMatchIndex; matchPairIndex>=0;) {
							long blockVstart, blockVstop, blockHstart, blockHstop, tmpBlockVstart, tmpBlockHstart;
							blockVstop= c->matchPairs[matchPairIndex].verticalOffset+kmerSpan-1;
							blockHstop= c->matchPairs[matchPairIndex].horizontalOffset+kmerSpan-1;
							blockVstart= c->matchPairs[matchPairIndex].verticalOffset;
							blockHstart= c->matchPairs[matchPairIndex].horizontalOffset;
							xDEBUG(DEB_EVALUATE_ANCHORING,
										 fprintf(stderr, " matchPairIndex=%d, initialize block w/ start (%d,%d) stop (%d,%d) on diagonal %d \n",
														matchPairIndex, blockVstart, blockHstart, blockVstop, blockHstop, blockHstart-blockVstart));
							int blockIndex;
							for (blockIndex=c->matchPairs[matchPairIndex].previousMatchPairInRun;
									 blockIndex>=0;
									 blockIndex=c->matchPairs[blockIndex].previousMatchPairInRun) {
								tmpBlockVstart= c->matchPairs[blockIndex].verticalOffset;
								tmpBlockHstart= c->matchPairs[blockIndex].horizontalOffset;
								if ( (blockHstart-blockVstart)==(tmpBlockHstart-tmpBlockVstart)) {
									blockVstart = tmpBlockVstart;
									blockHstart = tmpBlockHstart;
									xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, " extend block to start (%d,%d) stop (%d,%d) on diagonal %d\n",
																										 blockVstart, blockHstart, blockVstop, blockHstop, blockHstart-blockVstart));
								} else {
									break;
								}
							}
							matchPairIndex=blockIndex;
							
							xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, "now evaluating matches/mismatches\n"));
							char* horizontalBlockSequence =
								&c->targetTemplate[(long)start + (long)blockHstart-c->targetTemplateStart];
							xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, "horiz block seq %s\n", horizontalBlockSequence));
							int hBlockIndex, vBlockIndex;
							for (hBlockIndex=0, vBlockIndex=blockVstart; vBlockIndex<=blockVstop; vBlockIndex++, hBlockIndex++) {
								xDEBUG(DEB_EVALUATE_ANCHORING, fprintf(stderr, "cmp v[%d]=%c to h[%d]=%c\n",
																											 vBlockIndex, c->readTemplate[vBlockIndex],
																											 hBlockIndex, horizontalBlockSequence[hBlockIndex]));
								if ( (c->readTemplate[vBlockIndex]==horizontalBlockSequence[hBlockIndex])
                    ||
                    (bisulfiteSequencingMapping && c->readTemplate[vBlockIndex]=='T'
                     && horizontalBlockSequence[hBlockIndex]=='C')) {
									anchoringMatches++;
								} else {
									anchoringMismatches++;
								}
							}
						}
						guint32 skeletonScore;
						xDEBUG(DEB_FAIL_SW, fprintf(stderr, "%s matches %d mismatches %d\n", sequenceInfo->sequenceName,
																				anchoringMatches, anchoringMismatches));
						if(anchoringMismatches*4>anchoringMatches){
							xDEBUG(DEB_BAD_ANCHORING, fprintf(stderr, "really poor anchoring M %d m %d \n", anchoringMatches, anchoringMismatches));
							reallyPoorAnchorings++;
							continue;
						} else {
							skeletonScore = anchoringMatches-3 * anchoringMismatches;
                                                        xDEBUG(DEB_BAD_ANCHORING, fprintf(stderr,"cmp skel %d vs %d\n", skeletonScore, sequenceInfo->bestSkeletonScore));  
							if (skeletonScore<sequenceInfo->bestSkeletonScore*7/10) {
								tSkelScore ++;
								continue;
							}	
						}
						
						AlignmentSummary alignmentSummary;
						swCalls += 1;
						failedSWCalls += 1;
						int swScore;
            if (bisulfiteSequencingMapping) {
              swScore=bandedSWAlignmentInfoBisulfiteSeq(c->bswMemory,
										 c->readTemplate,
										 &c->targetTemplate[alignmentHorizontalStart-c->targetTemplateStart],
										 sequenceLength,
										 band, sequenceInfo->bestSWScore*withinTopPercent, &alignmentSummary);
            } else {
              swScore=bandedSWAlignmentInfo(c->bswMemory,
										 c->readTemplate,
										 &c->targetTemplate[alignmentHorizontalStart-c->targetTemplateStart],
										 sequenceLength,
										 band, sequenceInfo->bestSWScore*withinTopPercent, &alignmentSummary);
            }

						
						xDEBUG(DEB_HWIN, fprintf(stderr, "got alignment score =%d \n", swScore));
						if (swScore>=sequenceInfo->bestSWScore*withinTopPercent) {
						    xDEBUG(DEB_TOPPERCENT,fprintf(stderr, ">>> about to print %s onto %s %d\t%d\t%c\t%d    exceeds threhold of %d %g %g\n",
												sequenceInfo->sequenceName, currentSequence,
												start+hStart+1, start+hStop+kmerSpan,
												currentVerticalSequenceId%2==0?'+':'-',
												swScore, sequenceInfo->bestSWScore, withinTopPercent, sequenceInfo->bestSWScore*withinTopPercent));
                failedSWCalls -= 1;
								if (skeletonScore>sequenceInfo->bestSkeletonScore) {
									sequenceInfo->bestSkeletonScore=skeletonScore;
								}
							if (swScore>sequenceInfo->bestSWScore) {
								sequenceInfo->bestChrom = currentChrom;
								sequenceInfo->bestStart = start+hStop+kmerSpan;
								sequenceInfo->bestSWScore = swScore;
								sequenceInfo->bestScoreMappings=1;
							} else {
								if (swScore==sequenceInfo->bestSWScore) {
                  if( (sequenceInfo->bestChrom != currentChrom ||
										 (sequenceInfo->bestChrom == currentChrom )&&
										 ( sequenceInfo->bestStart+sequenceLength<start+hStart+1))
                    ) {
                    sequenceInfo->bestChrom = currentChrom;
                    sequenceInfo->bestStart = start+hStop+kmerSpan;
                    sequenceInfo->bestSWScore = swScore;
                    sequenceInfo->bestScoreMappings++;
                  } else {
                  	continue;
                  }
                }
							}
							xDEBUG(DEB_HWIN,fprintf(stderr, "??? about to print\n"));

							if (swScore>=kmerSpan) {
								xDEBUG(DEB_FAIL_SW, fprintf(stderr, "best score M %d m %d\n", anchoringMatches, anchoringMismatches));
								xDEBUG(DEB_HWIN,fprintf(stderr, "### about to print %s onto %s %d\t%d\t%c\t%d\n",
												sequenceInfo->sequenceName, currentSequence,
												start+hStart+1, start+hStop+kmerSpan,
												currentVerticalSequenceId%2==0?'+':'-',
												swScore));
								char outputLine[20*MAX_LINE_LENGTH];
								/*	fprintf(outputFilePtr, "%s\t%d\t%d\t%s\t%c\t%d\t%d\n",
												currentSequence,
												start+hStart+1, start+hStop+kmerSpan,
												sequenceInfo->sequenceName,
												currentVerticalSequenceId%2==0?'+':'-', swScore, sequenceId);
								*/
								SAMInfo samInfo;
								resetSamInfo(&samInfo);
								callVariants(c->readTemplate,
														&c->targetTemplate[alignmentHorizontalStart-c->targetTemplateStart],
														&samInfo, &alignmentSummary, bisulfiteSequencingMapping);
								if (alignmentSummary.numMismatches != samInfo.numberOfBasePairVariants) {
									fprintf(stderr, "incorrect number of mismatches for read %s: %d vs %d\n",
													sequenceInfo->sequenceName,
													alignmentSummary.numMismatches, samInfo.numberOfBasePairVariants);
								}
                const char* qualityScores = pp->verticalFastqUtil->retrieveQualityScores(sequenceId);
								if (bisulfiteSequencingMapping) {
                  outputBisulfiteMappingLine(sequenceId, swScore, sequenceInfo, chromLength, currentSequence,
                                        alignmentHorizontalStart, &alignmentSummary,
                                        outputLine, &samInfo, qualityScores, pp);
                } else {
                  outputRegularPashLine(sequenceId, swScore, sequenceInfo, chromLength, currentSequence,
                                        alignmentHorizontalStart, &alignmentSummary,
                                        currentVerticalSequenceId%2==0?'+':'-',
                                        outputLine, &samInfo, qualityScores);
		  xDEBUG(DEB_HWIN,fprintf(stderr, "got out line >>%s<<\n", outputLine));
                }
							
								fprintf(outputFilePtr, "%s", outputLine);
							}
							
						} else {
							xDEBUG(DEB_FAIL_SW,
										 fprintf(stderr, "fail sw %s M %d m %d \n",
														 sequenceInfo->sequenceName, anchoringMatches, anchoringMismatches));
							if (skeletonScore<sequenceInfo->bestSkeletonScore*9/10) {
								predSkelScore ++;
							}
						}
					}
				}
			}

			// if yes, figure out vertical match start, vertical match stop, and perform alignment
			// if alignment score is >= best sw score, dump alignment location
		} else {
			xDEBUG(DEB_WEAK_CANDIDATES, fprintf(stderr, "[%d][%s] bAC %d vs %d\n",
																	sequenceId, sequenceInfo->sequenceName, bestMatchScore, sequenceInfo->bestAnchoringScore*3/4));
		}
	}
	xDEBUG(DEB_PERFORM_COLLATION, fprintf(stderr, "stop collation \n"));
}


int bandedSW(int *scoringMatrix, char* verticalSequence, char *horizontalSequence, int sizeVerticalSequence, int band) {
	int i, j, leftVal, upVal, diagVal;
	int *prevLineH, *currentLineH;
	int score, bestScore;
	int matchGain = 1;
	int mismatchPenalty = -2;
	int gapPenalty = -3;
	int *tmpLine;
	bestScore = 0;
	int bandPlus1 = band+1;
	int bestGlobalScore = 0;

	xDEBUG(DEB_BSW, fprintf(stderr, "aligning %s vs %s\n", verticalSequence, horizontalSequence));
	prevLineH = scoringMatrix;
	for (j=0; j<=bandPlus1; j++) {
		prevLineH[j]=0;
	}
	for (i=0; i<sizeVerticalSequence;i++) {
		currentLineH = prevLineH+band+2;
		currentLineH[0]=0;
		xDEBUG(DEB_BSW, fprintf(stderr, "vert base %c\n", verticalSequence[i]));
		for (j=1; j<=band; j++) {
			leftVal = currentLineH[j-1]+gapPenalty;
			upVal = prevLineH[j+1]+gapPenalty;
			if (verticalSequence[i]==horizontalSequence[i+j-1]) {
				diagVal = prevLineH[j]+matchGain;
			} else {
				diagVal = prevLineH[j]+mismatchPenalty;
			}
			if (diagVal>upVal) {
				bestScore = diagVal;
			} else {
				bestScore = upVal;
			}
			if (bestScore<leftVal) {
				bestScore = leftVal;
			}
			if (bestScore>0) {
				currentLineH[j]=bestScore;
				if (bestScore>bestGlobalScore) {
					bestGlobalScore=bestScore;
				}
			} else {
				currentLineH[j]=0;
			}

			xDEBUG(DEB_BSW, fprintf(stderr, "[%d][%d] cmp %c vs %c leftVal=%d upVal=%d  diagH %d, diagVal=%d bestScore=%d\n",
															i, j, verticalSequence[i], horizontalSequence[i+j-1], leftVal, upVal, prevLineH[j], diagVal, currentLineH[j]));
		}
		currentLineH[bandPlus1]=0;
		//tmpLine = prevLineH;
		prevLineH = currentLineH;
		//currentLineH = prevLineH;
	}
	return bestGlobalScore;
}

void filterOutput(char* tmpOutputFileName, FILE *outputFilePtr, SequenceInfo* verticalSequenceInfos, guint32 maxReadMappings, double withinTopPercent) {
	// for now, select only best mappings
	// optimize for best match mapping: don't write additional mappings of the same score during collation step
	FILE *tmpOutputFilePtr= fopen(tmpOutputFileName, "rt");
  if (tmpOutputFilePtr==NULL) {
    fprintf(stderr, "could not open temporary output file %s for reading\n", tmpOutputFileName);
    fflush(stderr);
    exit(2);
  }
  char chromosome[MAX_DEFNAME_SIZE], readName[MAX_DEFNAME_SIZE];
	guint32 chromStart, chromStop, bwScore, sequenceId, readLength, chromLength;
	char strand;
	char tmpLine[20*MAX_LINE_LENGTH];
	while(fgets(tmpLine, 20*MAX_LINE_LENGTH-1, tmpOutputFilePtr) !=NULL) {
		sscanf(tmpLine, "%d %d", &sequenceId, &bwScore);
		SequenceInfo * sequenceInfo = verticalSequenceInfos+sequenceId;
      		xDEBUG(DEB_TOPPERCENT, fprintf(stderr, "2: got mapping for for %d, score %d \n", sequenceId, bwScore));
		if (sequenceInfo->passingMappings<=maxReadMappings&& bwScore >=sequenceInfo->bestSWScore*withinTopPercent) {
			sequenceInfo->passingMappings ++;
      xDEBUG(DEB_TOPPERCENT, fprintf(stderr, "passing mapping for %s: %d vs %d %g %g\n",
                                     sequenceInfo->sequenceName, bwScore, sequenceInfo->bestSWScore, withinTopPercent, sequenceInfo->bestSWScore*withinTopPercent ));
		}
	}
  fclose(tmpOutputFilePtr);
  
  tmpOutputFilePtr= fopen(tmpOutputFileName, "rt");
  if (tmpOutputFilePtr==NULL) {
    fprintf(stderr, "could not open temporary output file %s for reading\n", tmpOutputFileName);
    fflush(stderr);
    exit(2);
  }
  
	while(fgets(tmpLine, 20*MAX_LINE_LENGTH-1, tmpOutputFilePtr) !=NULL) {
		sscanf(tmpLine, "%d %d", &sequenceId, &bwScore);
		SequenceInfo * sequenceInfo = verticalSequenceInfos+sequenceId;
      		xDEBUG(DEB_TOPPERCENT, fprintf(stderr, "2: got mapping for for %d, score %d \n", sequenceId, bwScore));
		xDEBUG(DEB_FILTER_OUTPUT,
					 fprintf(stderr, "got %s\t%d\t%d\t%s\t%c\t%d; bestScore=sequenceInfo->bestSWScore; numMappings=%d vs %d\n",
							chromosome, chromStart, chromStop, readName, strand, bwScore, sequenceInfo->bestScoreMappings, maxReadMappings));
		if (sequenceInfo->passingMappings<=maxReadMappings && bwScore >=sequenceInfo->bestSWScore*withinTopPercent) {
			int idx;
			int lenTmpLine = strlen(tmpLine);
                        tmpLine[lenTmpLine-1]=' ';
			for (idx=0; tmpLine[idx]!='$' && idx<lenTmpLine ; idx++);
			if (idx<lenTmpLine) {
				fprintf(outputFilePtr, "%s\t%d\n", tmpLine+idx+2, sequenceInfo->passingMappings);
			} else {
				fprintf(stderr, "incorrect line %s", tmpLine);
			}
		}
	}
  fclose(tmpOutputFilePtr);
}



int bandedSWAlignmentInfo(int *scoringMatrix,
													char* verticalSequence, char *horizontalSequence,
													int sizeVerticalSequence, int band,
													float targetScore, AlignmentSummary* alignmentSummary) {
	int i, j, leftVal, upVal, diagVal;
	int *prevLineH, *currentLineH;
	int score, bestScore;
	int matchGain = 1;
	int mismatchPenalty = -2;
	int gapPenalty = -3;
	int *tmpLine;
	bestScore = 0;
	int bandPlus1 = band+1;
	int bestGlobalScore = 0;
	int bestScoreRow, bestScoreCol;

	bestScoreRow = 0;
	bestScoreCol = 0;

	prevLineH = scoringMatrix;
	for (j=0; j<=bandPlus1; j++) {
		prevLineH[j]=0;
	}
	for (i=0; i<sizeVerticalSequence;i++) {
		currentLineH = prevLineH+band+2;
		currentLineH[0]=0;
		xDEBUG(DEB_BSW, fprintf(stderr, "vert base %c line %x\n", verticalSequence[i], currentLineH));
		for (j=1; j<=band; j++) {
			leftVal = currentLineH[j-1]+gapPenalty;
			upVal = prevLineH[j+1]+gapPenalty;
			if (verticalSequence[i]==horizontalSequence[i+j-1]) {
				diagVal = prevLineH[j]+matchGain;
			} else {
				diagVal = prevLineH[j]+mismatchPenalty;
			}
			if (diagVal>upVal) {
				bestScore = diagVal;
			} else {
				bestScore = upVal;
			}
			if (bestScore<leftVal) {
				bestScore = leftVal;
			}
			if (bestScore>0) {
				currentLineH[j]=bestScore;
				if (bestScore>bestGlobalScore) {
					bestScoreCol = j;
					bestScoreRow = i;
					bestGlobalScore=bestScore;
				}
			} else {
				currentLineH[j]=0;
			}

			xDEBUG(DEB_BSW, fprintf(stderr, "[%d][%d] cmp %c vs %c leftVal=%d upVal=%d  diagH %d, diagVal=%d bestScore=%d\n",
															i, j, verticalSequence[i], horizontalSequence[i+j-1], leftVal, upVal, prevLineH[j], diagVal, currentLineH[j]));
		}
		currentLineH[bandPlus1]=0;
		//tmpLine = prevLineH;
		prevLineH = currentLineH;
		//currentLineH = prevLineH;
	}

	if (bestGlobalScore<targetScore || bestGlobalScore==0) {
		return bestGlobalScore;
	}


	// backtrace the best alignment and dump it
	int currentRow = bestScoreRow;
	int currentCol = bestScoreCol;
	int done = 0;

	alignmentSummary->numMatches =0;
  alignmentSummary->numMismatches = 0;
  alignmentSummary->numVerticalGaps= 0;
  alignmentSummary->numHorizontalGaps= 0;
  alignmentSummary->numGapBases = 0;
  alignmentSummary->numBlocks = 0;
	alignmentSummary->verticalStop = currentRow;
  alignmentSummary->horizontalStop = currentRow+currentCol-1;
  int horizontalBlockStop = alignmentSummary->horizontalStop;
  int verticalBlockStop = alignmentSummary->verticalStop;
	TraceStatus status = InBlock;

	currentLineH = scoringMatrix+(band+2)*(currentRow+1);
	while (currentRow>=0 && !done) {
		prevLineH = currentLineH - (band+2);
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "currentRow=%d currentCol=%d currentLine=%x prevLine=%x\n",
									 currentRow, currentCol, currentLineH, prevLineH));
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "current score %d\n", currentLineH[currentCol]));
		// determine direction of score OR stop
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "comp with prev[%d]=%d +M %d -m %d currentLineH[%d]=%d +G %d prevLineH[%d]=%d +G %d %c %c\n",
																	currentCol,	prevLineH[currentCol],
																	prevLineH[currentCol]+matchGain, prevLineH[currentCol]+mismatchPenalty,
																	currentCol-1,
																	currentLineH[currentCol-1],
																	currentLineH[currentCol-1]+gapPenalty,
																	currentCol+1, prevLineH[currentCol+1], prevLineH[currentCol+1]+gapPenalty,
																	verticalSequence[currentRow], horizontalSequence[currentRow+currentCol-1]
																	))
		if (prevLineH[currentCol]+matchGain==currentLineH[currentCol] && verticalSequence[currentRow]==horizontalSequence[currentRow+currentCol-1]) {
			alignmentSummary->numMatches ++;
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"match gain on diag, matches=%d \n", alignmentSummary->numMatches));
			if (prevLineH[currentCol]==0) {
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow+1;
				alignmentSummary->verticalStart = currentRow;
				alignmentSummary->horizontalStart = currentRow+currentCol-1;
				alignmentSummary->numBlocks++;
				done=1;
			}
			if(status==InGap) {
				status=InBlock;
				verticalBlockStop=currentRow;
				horizontalBlockStop = currentCol;
			}
			currentRow-=1;
			currentLineH=prevLineH;
		} else if ( prevLineH[currentCol]+mismatchPenalty==currentLineH[currentCol] &&
							 verticalSequence[currentRow]!=horizontalSequence[currentRow+currentCol-1]) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"mismatch gain on diag\n"));
			alignmentSummary->numMismatches ++;
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"match gain on diag, matches=%d \n", alignmentSummary->numMatches));
			if (prevLineH[currentCol]==0) {
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow+1;
				alignmentSummary->verticalStart = currentRow;
				alignmentSummary->horizontalStart = currentRow+currentCol-1;
				alignmentSummary->numBlocks++;
				done=1;
			}
			if(status==InGap) {
				status=InBlock;
				verticalBlockStop=currentRow;
				horizontalBlockStop = currentCol;
			}
			currentRow-=1;
			currentLineH=prevLineH;
		} else if (currentLineH[currentCol-1]+gapPenalty == currentLineH[currentCol]) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"horiz gap\n"));
			alignmentSummary->numGapBases ++;
			if(status==InBlock) {
				status=InGap;
				alignmentSummary->numHorizontalGaps++;
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow+1;
        //alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow;
				alignmentSummary->numBlocks++;
			}
			currentCol-=1;
		} else if (prevLineH[currentCol+1]+gapPenalty == currentLineH[currentCol]) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"vert gap\n"));
			if(status==InBlock) {
				status=InGap;
				alignmentSummary->numVerticalGaps++;
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow+1;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow;
				alignmentSummary->numBlocks++;
			}
			currentRow-=1;
			currentCol+=1;
			currentLineH=prevLineH;
		}
		//break;
	}
	if (!done) {
		fprintf(stderr, "boo %S %S!\n", verticalSequence, horizontalSequence);
		return 0;
	}
	return bestGlobalScore;
}


/** Use the alignment summary to call base pair variants and methylation status.*/
int callVariants(char* readSequence, char* templateSequence, 
								 SAMInfo *samInfo, AlignmentSummary* alignmentSummary, int bisulfiteSequencing) {
	guint32 blockIndex;
	guint32 hStartBlock, hStopBlock;
	guint32 vStartBlock, vStopBlock;
	guint32 hBpIndex, vBpIndex;
	
	guint32 templateSize = strlen(templateSequence);
	guint32 newConvertedBase, newBPVariantBase, newMethylatedBase;
	
	samInfo->numberOfBasePairVariants=0;
	samInfo->numberOfCGMethylatedBases=0;
	samInfo->numberOfCHGMethylatedBases=0;
	samInfo->numberOfCHHMethylatedBases=0;
	samInfo->numberOfConvertedBases=0;
	xDEBUG(DEB_VARIANTS, fprintf(stderr, "checking variants and methylated bases between %s\n%s\n",
															 readSequence, templateSequence));
	xDEBUG(DEB_VARIANTS, fprintf(stderr,"num blocks: %d\n", alignmentSummary->numBlocks));
	// detect basepair variants
	for(blockIndex=alignmentSummary->numBlocks-1; 1; blockIndex--) {
		hStartBlock = alignmentSummary->horizontalBlockStarts[blockIndex];
		hStopBlock = alignmentSummary->horizontalBlockStarts[blockIndex] + alignmentSummary->blockSizes[blockIndex]-1;
		vStartBlock = alignmentSummary->verticalBlockStarts[blockIndex];
		vStopBlock = alignmentSummary->verticalBlockStarts[blockIndex]+alignmentSummary->blockSizes[blockIndex]-1;
		xDEBUG(DEB_VARIANTS, fprintf(stderr, "block %d h [%d, %d], v [%d,%d]\n", blockIndex,
																	 hStartBlock+1, hStopBlock+1, vStartBlock, vStopBlock));
		
		for (hBpIndex = hStartBlock, vBpIndex = vStartBlock;
				 hBpIndex <= hStopBlock;
				 hBpIndex++, vBpIndex++) {	
			if (bisulfiteSequencing) {
				xDEBUG(DEB_CHH_VARIANTS, fprintf(stderr, "[%d]%c vs [%d]%c\n", vBpIndex, readSequence[vBpIndex],
																				 hBpIndex, templateSequence[hBpIndex]));
				if (readSequence[vBpIndex]=='C' && templateSequence[hBpIndex]=='C') {
					if (templateSequence[hBpIndex+1]=='G') {
						// check for g
					samInfo->cgMethylatedBasesPositions[samInfo->numberOfCGMethylatedBases]=hBpIndex;
					samInfo->cgMethylatedBasesPositionsInRead[samInfo->numberOfCGMethylatedBases]=vBpIndex;
					samInfo->numberOfCGMethylatedBases += 1;
					xDEBUG(DEB_CHH_VARIANTS, fprintf(stderr, "CG M %d\n", hBpIndex));
					} else if (templateSequence[hBpIndex+2]=='G') {
						// chg
						samInfo->chgMethylatedBasesPositions[samInfo->numberOfCHGMethylatedBases]=hBpIndex;
						samInfo->chgMethylatedBasesPositionsInRead[samInfo->numberOfCHGMethylatedBases]=vBpIndex;
						samInfo->numberOfCHGMethylatedBases+= 1;
						xDEBUG(DEB_CHH_VARIANTS, fprintf(stderr, "CHG M %d\n", hBpIndex));
					} else {
						// chh
						samInfo->chhMethylatedBasesPositions[samInfo->numberOfCHHMethylatedBases]=hBpIndex;
						samInfo->chhMethylatedBasesPositionsInRead[samInfo->numberOfCHHMethylatedBases]=vBpIndex;
						samInfo->numberOfCHHMethylatedBases+=1;
						xDEBUG(DEB_CHH_VARIANTS, fprintf(stderr, "CHH M %d\n", hBpIndex));
					}
					continue;
				}
				if (readSequence[vBpIndex]=='T' && templateSequence[hBpIndex]=='C') {
					// converted base
					newConvertedBase = samInfo->numberOfConvertedBases;
					samInfo->numberOfConvertedBases += 1;
					samInfo->convertedBasesPositions[newConvertedBase] = hBpIndex;
					samInfo->convertedBasesPositionsInRead[newConvertedBase] = vBpIndex;
					xDEBUG(DEB_VARIANTS, fprintf(stderr, "Converted base at %d: %c vs %c\n",
																				 hBpIndex, readSequence[vBpIndex], templateSequence[hBpIndex]));
					xDEBUG(DEB_CHH_VARIANTS, fprintf(stderr, "T %d\n", hBpIndex));
					continue;
				}
				if (readSequence[vBpIndex] != templateSequence[hBpIndex]) {
					// SNP !!
					newBPVariantBase = samInfo->numberOfBasePairVariants;
					samInfo->numberOfBasePairVariants += 1;
					samInfo->basePairVariantsPositions[newBPVariantBase] = hBpIndex;
					samInfo->basePairVariantsPositionsInRead[newBPVariantBase] = vBpIndex;
					samInfo->basePairVariantAlleles[newBPVariantBase] = readSequence[vBpIndex];
					xDEBUG(DEB_VARIANTS, fprintf(stderr, "Variant base at %d: %c vs %c\n",
																				 hBpIndex, readSequence[vBpIndex], templateSequence[hBpIndex]));

					continue;
				}
			} else {
				if (readSequence[vBpIndex] != templateSequence[hBpIndex]) {
					// SNP !!
					newBPVariantBase = samInfo->numberOfBasePairVariants;
					samInfo->numberOfBasePairVariants += 1;
					samInfo->basePairVariantsPositions[newBPVariantBase] = hBpIndex;
					samInfo->basePairVariantsPositionsInRead[newBPVariantBase] = vBpIndex;
					samInfo->basePairVariantAlleles[newBPVariantBase] = readSequence[vBpIndex];
					xDEBUG(DEB_VARIANTS, fprintf(stderr, "Variant base at %d: %c vs %c\n",
																				 hBpIndex, readSequence[vBpIndex], templateSequence[hBpIndex]));
					continue;
				}
			}
		}
		
		if (blockIndex==0) {
			break;
		}
	};
}

void resetSamInfo(SAMInfo* samInfo) {
	samInfo->numberOfBasePairVariants=0;
	samInfo->numberOfCGMethylatedBases=0;
	samInfo->numberOfCHGMethylatedBases=0;
	samInfo->numberOfCHHMethylatedBases=0;
	samInfo->numberOfConvertedBases=0;
}

int bandedSWAlignmentInfoBisulfiteSeq(int *scoringMatrix,
													char* verticalSequence, char *horizontalSequence,
													int sizeVerticalSequence, int band,
													float targetScore, AlignmentSummary* alignmentSummary) {
	int i, j, leftVal, upVal, diagVal;
	int *prevLineH, *currentLineH;
	int score, bestScore;
	int matchGain = 1;
	int mismatchPenalty = -2;
	int gapPenalty = -3;
	int *tmpLine;
	bestScore = 0;
	int bandPlus1 = band+1;
	int bestGlobalScore = 0;
	int bestScoreRow, bestScoreCol;
	char cH, cV;
	bestScoreRow = 0;
	bestScoreCol = 0;

	prevLineH = scoringMatrix;
	for (j=0; j<=bandPlus1; j++) {
		prevLineH[j]=0;
	}
	for (i=0; i<sizeVerticalSequence;i++) {
		currentLineH = prevLineH+band+2;
		currentLineH[0]=0;
		xDEBUG(DEB_BSW, fprintf(stderr, "vert base %c line %x\n", verticalSequence[i], currentLineH));
		cV = verticalSequence[i];
		for (j=1; j<=band; j++) {
			leftVal = currentLineH[j-1]+gapPenalty;
			upVal = prevLineH[j+1]+gapPenalty;
			cH = horizontalSequence[i+j-1];
			if (cV==cH || (cV=='T' && cH=='C')) {
				diagVal = prevLineH[j]+matchGain;
			} else {
				diagVal = prevLineH[j]+mismatchPenalty;
			}
			if (diagVal>upVal) {
				bestScore = diagVal;
			} else {
				bestScore = upVal;
			}
			if (bestScore<leftVal) {
				bestScore = leftVal;
			}
			if (bestScore>0) {
				currentLineH[j]=bestScore;
				if (bestScore>bestGlobalScore) {
					bestScoreCol = j;
					bestScoreRow = i;
					bestGlobalScore=bestScore;
				}
			} else {
				currentLineH[j]=0;
			}

			xDEBUG(DEB_BSW, fprintf(stderr, "[%d][%d] cmp %c vs %c leftVal=%d upVal=%d  diagH %d, diagVal=%d bestScore=%d\n",
															i, j, verticalSequence[i], horizontalSequence[i+j-1], leftVal, upVal, prevLineH[j], diagVal, currentLineH[j]));
		}
		currentLineH[bandPlus1]=0;
		//tmpLine = prevLineH;
		prevLineH = currentLineH;
		//currentLineH = prevLineH;
	}

	if (bestGlobalScore<targetScore || bestGlobalScore==0) {
		return bestGlobalScore;
	}


	// backtrace the best alignment and dump it
	int currentRow = bestScoreRow;
	int currentCol = bestScoreCol;
	int done = 0;

	alignmentSummary->numMatches =0;
  alignmentSummary->numMismatches = 0;
  alignmentSummary->numVerticalGaps= 0;
  alignmentSummary->numHorizontalGaps= 0;
  alignmentSummary->numGapBases = 0;
  alignmentSummary->numBlocks = 0;
	alignmentSummary->verticalStop = currentRow;
  alignmentSummary->horizontalStop = currentRow+currentCol-1;
  int horizontalBlockStop = alignmentSummary->horizontalStop;
  int verticalBlockStop = alignmentSummary->verticalStop;
	TraceStatus status = InBlock;
	//WSchar cV, cH;
	
	currentLineH = scoringMatrix+(band+2)*(currentRow+1);
	while (currentRow>=0 && !done) {
		prevLineH = currentLineH - (band+2);
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "currentRow=%d currentCol=%d currentLine=%x prevLine=%x\n",
									 currentRow, currentCol, currentLineH, prevLineH));
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "current score %d\n", currentLineH[currentCol]));
		// determine direction of score OR stop
		xDEBUG(DEB_TRACEBACK, fprintf(stderr, "comp with prev[%d]=%d +M %d -m %d currentLineH[%d]=%d +G %d prevLineH[%d]=%d +G %d %c %c\n",
																	currentCol,	prevLineH[currentCol],
																	prevLineH[currentCol]+matchGain, prevLineH[currentCol]+mismatchPenalty,
																	currentCol-1,
																	currentLineH[currentCol-1],
																	currentLineH[currentCol-1]+gapPenalty,
																	currentCol+1, prevLineH[currentCol+1], prevLineH[currentCol+1]+gapPenalty,
																	verticalSequence[currentRow], horizontalSequence[currentRow+currentCol-1]
																	))
		cV = verticalSequence[currentRow];
		cH = horizontalSequence[currentRow+currentCol-1];
		if (prevLineH[currentCol]+matchGain==currentLineH[currentCol] && (cV==cH || (cV=='T' && cH=='C'))) {
			alignmentSummary->numMatches ++;
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"match gain on diag, matches=%d \n", alignmentSummary->numMatches));
			if (prevLineH[currentCol]==0) {
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow+1;
				alignmentSummary->verticalStart = currentRow;
				alignmentSummary->horizontalStart = currentRow+currentCol-1;
				alignmentSummary->numBlocks++;
				done=1;
			}
			if(status==InGap) {
				status=InBlock;
				verticalBlockStop=currentRow;
				horizontalBlockStop = currentCol;
			}
			currentRow-=1;
			currentLineH=prevLineH;
		} else if ( prevLineH[currentCol]+mismatchPenalty==currentLineH[currentCol] &&
							 cV!=cH && (cV !='T' || cH!='C')) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"mismatch gain on diag\n"));
			alignmentSummary->numMismatches ++;
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"match gain on diag, matches=%d \n", alignmentSummary->numMatches));
			if (prevLineH[currentCol]==0) {
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow+1;
				alignmentSummary->verticalStart = currentRow;
				alignmentSummary->horizontalStart = currentRow+currentCol-1;
				alignmentSummary->numBlocks++;
				done=1;
			}
			if(status==InGap) {
				status=InBlock;
				verticalBlockStop=currentRow;
				horizontalBlockStop = currentCol;
			}
			currentRow-=1;
			currentLineH=prevLineH;
		} else if (currentLineH[currentCol-1]+gapPenalty == currentLineH[currentCol]) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"horiz gap\n"));
			alignmentSummary->numGapBases ++;
			if(status==InBlock) {
				status=InGap;
				alignmentSummary->numHorizontalGaps++;
				//alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol-1;
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow+1;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow;

				alignmentSummary->numBlocks++;
			}
			currentCol-=1;
		} else if (prevLineH[currentCol+1]+gapPenalty == currentLineH[currentCol]) {
			xDEBUG(DEB_TRACEBACK, fprintf(stderr,"vert gap\n"));
			if(status==InBlock) {
				status=InGap;
				alignmentSummary->numVerticalGaps++;
				alignmentSummary->verticalBlockStarts[alignmentSummary->numBlocks]=currentRow+1;
        alignmentSummary->horizontalBlockStarts[alignmentSummary->numBlocks]=currentRow+currentCol;
        alignmentSummary->blockSizes[alignmentSummary->numBlocks] = verticalBlockStop - currentRow;
				alignmentSummary->numBlocks++;
			}
			currentRow-=1;
			currentCol+=1;
			currentLineH=prevLineH;
		}
		//break;
	}
	if (!done) {
		fprintf(stderr, "boo %S %S!\n", verticalSequence, horizontalSequence);
		return 0;
	}
	return bestGlobalScore;
}

int outputRegularPashLine(guint32 sequenceId, int swScore, SequenceInfo* sequenceInfo,
                          guint32 chromLength, char* currentSequence,
                          guint32 alignmentHorizontalStart,
                          AlignmentSummary* alignmentSummary,
                          char strand, char *outputLine, SAMInfo* samInfo, const char* qualityScores) {
  guint32 sequenceLength = sequenceInfo->sequenceLength;
  outputLine[0]='\0';
  sprintf(outputLine+strlen(outputLine),
          "%d\t%d\t%d\t%d\t$\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t",
    sequenceId,
    swScore,
    sequenceLength,
    chromLength,
    currentSequence,
    alignmentHorizontalStart+1+alignmentSummary->horizontalStart,
    alignmentHorizontalStart+1+alignmentSummary->horizontalStop,
    sequenceInfo->sequenceName,
    1+alignmentSummary->verticalStart, 1+alignmentSummary->verticalStop,
    strand,
    alignmentSummary->numMatches, alignmentSummary->numMismatches,
    alignmentSummary->numVerticalGaps+alignmentSummary->numHorizontalGaps,
    alignmentSummary->numGapBases,
    alignmentSummary->numBlocks);

  int i;
  
  for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
    sprintf(outputLine+strlen(outputLine), "%d,", alignmentSummary->blockSizes[i]);
  }
  sprintf(outputLine+strlen(outputLine), "\t");
  for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
    sprintf(outputLine+strlen(outputLine), "%d,",
            alignmentHorizontalStart+1+alignmentSummary->horizontalBlockStarts[i]);
  }
  sprintf(outputLine+strlen(outputLine), "\t");
  for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
    sprintf(outputLine+strlen(outputLine), "%d,", 1+alignmentSummary->verticalBlockStarts[i]);
  }
  
  sprintf(outputLine+strlen(outputLine), "\t");
  sprintf(outputLine+strlen(outputLine), "%d\t", samInfo->numberOfBasePairVariants);
  if (samInfo->numberOfBasePairVariants>0) {
    for (i=0; i<samInfo->numberOfBasePairVariants; i++) {
      sprintf(outputLine+strlen(outputLine), "%d/%c/%d," ,
              alignmentHorizontalStart+1+samInfo->basePairVariantsPositions[i],
              samInfo->basePairVariantAlleles[i],
              strand=='+'?qualityScores[samInfo->basePairVariantsPositionsInRead[i]]:
                          qualityScores[sequenceLength-1-samInfo->basePairVariantsPositionsInRead[i]]);
    }
  } else {
    sprintf(outputLine+strlen(outputLine), "-");
  }
    
  sprintf(outputLine+strlen(outputLine), "\n");
  return 0;
}


int outputBisulfiteMappingLine(guint32 sequenceId, int swScore, SequenceInfo* sequenceInfo,
                          guint32 chromLength, char* currentSequence,
                          guint32 alignmentHorizontalStart,
                          AlignmentSummary* alignmentSummary,
                          char *outputLine, SAMInfo* samInfo,
                          const char* qualityScores,
                          PashParameters* pp) {
  
  outputLine[0]='\0';
  guint32 sequenceLength = sequenceInfo->sequenceLength;
  if (pp->reverseStrandDnaMethMapping) {
    sprintf(outputLine+strlen(outputLine),
          "%d\t%d\t%d\t%d\t$\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t",
    sequenceId,
    swScore,
    sequenceLength,
    chromLength,
    pp->actualChromName,
    chromLength-(alignmentHorizontalStart+alignmentSummary->horizontalStop),
    chromLength-(alignmentHorizontalStart+alignmentSummary->horizontalStart),
    sequenceInfo->sequenceName,
    1+alignmentSummary->verticalStart, 1+alignmentSummary->verticalStop,
    '-',
    alignmentSummary->numMatches, alignmentSummary->numMismatches,
    alignmentSummary->numVerticalGaps+alignmentSummary->numHorizontalGaps,
    alignmentSummary->numGapBases,
    alignmentSummary->numBlocks);

    int i;
    for(i=0; i<alignmentSummary->numBlocks; i++) {
      sprintf(outputLine+strlen(outputLine), "%d,", alignmentSummary->blockSizes[i]);
    }
    sprintf(outputLine+strlen(outputLine), "\t");
    for(i=0; i<alignmentSummary->numBlocks; i++) {
      sprintf(outputLine+strlen(outputLine), "%d,",
              chromLength+1-
              (alignmentSummary->blockSizes[i]+alignmentHorizontalStart+alignmentSummary->horizontalBlockStarts[i]));
    }
    sprintf(outputLine+strlen(outputLine), "\t");
    for(i=0; i<alignmentSummary->numBlocks;  i++) {
      sprintf(outputLine+strlen(outputLine), "%d,",
              sequenceLength+1-(alignmentSummary->verticalBlockStarts[i]+alignmentSummary->blockSizes[i]));
    }
    
    // output SAM Info
    sprintf(outputLine+strlen(outputLine), "\t");
    sprintf(outputLine+strlen(outputLine), "%d\t", samInfo->numberOfBasePairVariants);
    if (samInfo->numberOfBasePairVariants>0) {
      for (i=samInfo->numberOfBasePairVariants-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%c/%d," ,
                chromLength-(alignmentHorizontalStart+samInfo->basePairVariantsPositions[i]),
                revComplementQuick(samInfo->basePairVariantAlleles[i]),
                qualityScores[samInfo->basePairVariantsPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }

    sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCGMethylatedBases);
    if (samInfo->numberOfCGMethylatedBases>0) {
      for (i=samInfo->numberOfCGMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                chromLength-(alignmentHorizontalStart+samInfo->cgMethylatedBasesPositions[i]),
                qualityScores[samInfo->cgMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }

		sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCHGMethylatedBases);
    if (samInfo->numberOfCHGMethylatedBases>0) {
      for (i=samInfo->numberOfCHGMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                chromLength-(alignmentHorizontalStart+samInfo->chgMethylatedBasesPositions[i]),
                qualityScores[samInfo->chgMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
		
		sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCHHMethylatedBases);
    if (samInfo->numberOfCHHMethylatedBases>0) {
      for (i=samInfo->numberOfCHHMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                chromLength-(alignmentHorizontalStart+samInfo->chhMethylatedBasesPositions[i]),
                qualityScores[samInfo->chhMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }

    sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfConvertedBases);
    if (samInfo->numberOfConvertedBases>0) {
      for (i=samInfo->numberOfConvertedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                chromLength-(alignmentHorizontalStart+samInfo->convertedBasesPositions[i]),
                qualityScores[samInfo->convertedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
		
    sprintf(outputLine+strlen(outputLine), "\n");
  } else {
    sprintf(outputLine+strlen(outputLine),
          "%d\t%d\t%d\t%d\t$\t%s\t%d\t%d\t%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t",
    sequenceId,
    swScore,
    sequenceLength,
    chromLength,
    currentSequence,
    alignmentHorizontalStart+alignmentSummary->horizontalStart+1,
    alignmentHorizontalStart+1+alignmentSummary->horizontalStop,
    sequenceInfo->sequenceName,
    1+alignmentSummary->verticalStart, 1+alignmentSummary->verticalStop,
    '+',
    alignmentSummary->numMatches, alignmentSummary->numMismatches,
    alignmentSummary->numVerticalGaps+alignmentSummary->numHorizontalGaps,
    alignmentSummary->numGapBases,
    alignmentSummary->numBlocks);

    int i;
    for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
      sprintf(outputLine+strlen(outputLine), "%d,", alignmentSummary->blockSizes[i]);
    }
    sprintf(outputLine+strlen(outputLine), "\t");
    for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
      sprintf(outputLine+strlen(outputLine), "%d,",
              alignmentHorizontalStart+1+alignmentSummary->horizontalBlockStarts[i]);
    }
    sprintf(outputLine+strlen(outputLine), "\t");
    for(i=alignmentSummary->numBlocks-1; i>=0; i--) {
      sprintf(outputLine+strlen(outputLine), "%d,", 1+alignmentSummary->verticalBlockStarts[i]);
    }
    
    // output SAM Info
    sprintf(outputLine+strlen(outputLine), "\t");
    sprintf(outputLine+strlen(outputLine), "%d\t", samInfo->numberOfBasePairVariants);
    if (samInfo->numberOfBasePairVariants>0) {
      for (i=0; i<samInfo->numberOfBasePairVariants; i++) {
        sprintf(outputLine+strlen(outputLine), "%d/%c/%d," ,
                alignmentHorizontalStart+1+samInfo->basePairVariantsPositions[i],
                samInfo->basePairVariantAlleles[i],
                qualityScores[samInfo->basePairVariantsPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
    
    
		sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCGMethylatedBases);
    if (samInfo->numberOfCGMethylatedBases>0) {
      for (i=samInfo->numberOfCGMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                alignmentHorizontalStart+1+samInfo->cgMethylatedBasesPositions[i],
                qualityScores[samInfo->cgMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }

		sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCHGMethylatedBases);
    if (samInfo->numberOfCHGMethylatedBases>0) {
      for (i=samInfo->numberOfCHGMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                alignmentHorizontalStart+1+samInfo->chgMethylatedBasesPositions[i],
                qualityScores[samInfo->chgMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
		
		sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfCHHMethylatedBases);
    if (samInfo->numberOfCHHMethylatedBases>0) {
      for (i=samInfo->numberOfCHHMethylatedBases-1; i>=0; i--) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                alignmentHorizontalStart+1+samInfo->chhMethylatedBasesPositions[i],
                qualityScores[samInfo->chhMethylatedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
		
    sprintf(outputLine+strlen(outputLine), "\t%d\t", samInfo->numberOfConvertedBases);
    if (samInfo->numberOfConvertedBases>0) {
      for (i=0; i<samInfo->numberOfConvertedBases; i++) {
        sprintf(outputLine+strlen(outputLine), "%d/%d," ,
                alignmentHorizontalStart+1+samInfo->convertedBasesPositions[i],
                qualityScores[samInfo->convertedBasesPositionsInRead[i]]);
      }
    } else {
      sprintf(outputLine+strlen(outputLine), "-");
    }
    
    sprintf(outputLine+strlen(outputLine), "\n");
  }
  return 0;
}


