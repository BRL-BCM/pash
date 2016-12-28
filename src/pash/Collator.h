/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef PASH30_COLLATOR_H___
#define PASH30_COLLATOR_H___

#include <glib.h>
#include "HiveHash.h"
#include "someConstants.h"

/** Data structure enabling the traversal of a stream of matches for a vertical kmer.*/
typedef struct {
  /** stream of vertical kmers */
  IntListRunner intListRunner;
	guint32 okey;
	guint32 verticalSeqID;
  // current vertical offset
	//guint16 verticalOffset;
	//guint16 horizontalOffset;
	// current vertical sequence id
	// current horizontal offset
	/// Current kmer.
	//	guint32 kmer;
} MatchStream;

/** Data structure holding the information for a matching pair of positions.
*/
typedef struct {
  /** The sequence ID of the vertical chunk containing the matching position.*/
  guint32 verticalSeqID;
  /** Vertical offset of the matching pair.*/
  guint32 verticalOffset;
  /** Horizontal offset of the matching pair.*/
  guint32 horizontalOffset;
  /** Diagonal of the current match pair.*/
  int diagonal;
  /** Actual kmer.*/
  guint32 kmer;
  /** The position of the previous match pair in the best run for the current match pair.*/
  int previousMatchPairInRun;
  /** Score of the best run ending in current match pair.*/
  guint32 bestRunScore;
  int isMatchRunStart;
} MatchPair;


/** Data structure containing information necessary for the collation.*/
typedef struct {
  /** Current stream of matches.*/
	MatchStream* matchStreams;
	MatchStream** matchStreamPtrs;
  /** Maximum match streams capacity.*/
	guint32 maxMatchStreams;
  /** Number of valid match streams for current horizontal chunk.*/
	guint32 validMatchStreams;
  /** Number of diagonal capacity.*/
  guint32 numberOfDiagonals ;
  /** Match pairs.*/
  MatchPair* matchPairs;
  /** Overall match pairs capacity.*/
  guint32 matchPairsCapacity;
  /** Starting positions in the match pairs of match runs.*/
  guint32 matchRunStarts[MAX_MATCH_STREAMS];
	/** Match run scores.*/
  guint32 matchRunScores[MAX_MATCH_STREAMS];
	char targetTemplate[3*MAX_READ_SIZE+2*DEFAULT_BAND];
	char readTemplate[MAX_READ_SIZE+1];
	long targetTemplateStart;
	int *bswMemory;

} CollatorControl;

/** Construct & initialize a CollatorControl data structure.*/
CollatorControl* initCollatorControl();
/** Add a new match stream to the collator control if there is a matching vertical kmer.*/
void addMatchStreamCollatorControl(CollatorControl *c, guint32 kmer, HiveHash* hh);
/** Reset the collator control for a new collation.*/
void resetCollatorControl(CollatorControl *c);
/** Release collator control resources.*/
void freeCollatorControl(CollatorControl *c);
/** Do the actual collation.*/
void performCollation(FILE* outputFilePtr, CollatorControl* c, SequenceHash* sequenceHash, PashParameters *pp,
											char* currentSequence, guint32 start, guint32 stop, guint32 currentChrom);
/*
#define STREAM_LESS_THAN(matchStream1,matchStream2)  			             \
  ( (matchStream1.verticalSeqID < matchStream2.verticalSeqID) ||       \
    ( (matchStream1.verticalSeqID == matchStream2.verticalSeqID) &&    \
      (matchStream1.verticalOffset < matchStream2.verticalOffset))  || \
    ( (matchStream1.verticalSeqID == matchStream2.verticalSeqID) &&    \
      (matchStream1.verticalOffset == matchStream2.verticalOffset) &&  \
      (matchStream1.horizontalOffset < matchStream2.horizontalOffset)) \
  )
*/

/*
#define STREAM_LESS_THAN(matchStream1,matchStream2)  			             \
  ( (matchStream1->verticalSeqID < matchStream2->verticalSeqID) ?1:           \
      (matchStream1->verticalSeqID > matchStream2->verticalSeqID) ?0:         \
          (matchStream1->verticalOffset < matchStream2->verticalOffset)?1:    \
           (matchStream1->verticalOffset > matchStream2->verticalOffset)? 0:  \
            (matchStream1->horizontalOffset < matchStream2->horizontalOffset))
*/

#define STREAM_LESS_THAN(matchStream1,matchStream2)  			             \
  ((matchStream1->verticalSeqID < matchStream2->verticalSeqID) ?1:           \
      (matchStream1->verticalSeqID > matchStream2->verticalSeqID) ?0:         \
          (matchStream1->okey < matchStream2->okey))

#define STREAM_LESS_THAN2(matchStream1,matchStream2)  			             \
  ( (matchStream1->verticalSeqID < matchStream2->verticalSeqID) ?1:           \
      (matchStream1->verticalSeqID > matchStream2->verticalSeqID) ?0:         \
          (matchStream1->verticalOffset < matchStream2->verticalOffset)?1:    \
           (matchStream1->verticalOffset > matchStream2->verticalOffset)? 0:  \
            (matchStream1->horizontalOffset < matchStream2->horizontalOffset))


/*
#define STREAM_LESS_THAN(matchStream1,matchStream2)  			\
  (matchStream1->key< matchStream2->key)
*/

#define STREAM_PARENT(i) (i/2)
#define LEFT_CHILD(i) (2*i)
#define RIGHT_CHILD(i) (2*i+1)
#define TOP_HORIZONTAL_MATCH(matchStreams) (matchStreams[1].horizontalSeqID)
#define TOP_VERTICAL_MATCH(matchStreams) (matchStreams[1].verticalSeqID)




/** Match information.*/
typedef struct {
	/** number of matching bases.*/
	guint32 numMatchingBases;
	/** number of gaps.*/
	guint32 numGaps;
	/** total gaps size.*/
	guint32 numGapBases;
	/** horizontal match start (within chunk).*/
	guint32 horizontalStart;
	/** horizontal match stop (within chunk).*/
	guint32 horizontalStop;
	/** vertical match start (within chunk).*/
	guint32 verticalStart;
	/** vertical match stop (within chunk).*/
	guint32 verticalStop;
	/** vertical chunk ID.*/
	guint32 verticalChunkID;
	/** horizontal chunk ID.*/
	guint32 horizontalChunkID;
	/** adaptive score, if used.*/
	unsigned long adaptiveScore;
	/** Flag whether the vertical matches belong to the forward or to the reverse complement sequence.*/
	int isFwd;
	int numBlocks;
	int blockSizes[MAX_GAPS];
	int vBlockStarts[MAX_GAPS];
	int hBlockStarts[MAX_GAPS];
} MatchInfo;


/** Inserts a stream in the streams priority queue.*/
int insertMatchStream(MatchStream* matchStreams, int numStreams, MatchStream oneMatchStream);
/** Extract first stream from streams priority queue.*/
MatchStream extractFirstStream(MatchStream* matchStreams, int numStreams);
/** Restore the heap property for the match streams priority queue.*/
void HeapifyStreams(MatchStream* matchStreams, int numStreams, int pos);
/** Setup a match stream and insert it into the priority queue.*/
int setMatchStreamAndInsertInQueue(MatchStream matchStream, MatchStream* matchStreamPriorityQueue, int numStreams,
																		guint32 kmer, guint32 horizontalOffset);
/** Setup the match stream for a horizontal kmer, querying the hive hash.*/
void addMatchStreamCollatorControl(CollatorControl *c, guint32 kmer, HiveHash* hh, guint32 horizontalOffset);
/** Advance the top match stream and restore the priority queue.*/
int advanceTopMatchStream(MatchStream* matchStreamPriorityQueue, int numStreams, int numDiagonals);

int deleteHeapMin(MatchStream* matchStreams, int numStreams);
void filterOutput(char* tmpOutputFileName, FILE *outputFilePtr,
                  SequenceInfo* verticalSequenceInfos, guint32 maxReadMappings, double withinTopPercent) ;

#endif
