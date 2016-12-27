/*
Copyright (c) 2004 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/

/******************************************************
 * NAME: FixedHashKey.h
 *
 * DESCRIPTION:
 * This module provides functionality for converting a sequence of up to
 * 16 letters from the alphabet "AGTCagtc" to a 32-bit int (the "key") and
 * furthermore provides functions for accessing the contents of the key as
 * bases, not bits.
 *
 * REQUIRES:
 *  glib library (GTK+, multiplatform C library)
 *  boolean.h and byte.h for data-type standardization.
 * *
 * AUTHOR: Andrew R. Jackson (andrewj@bcm.tmc.edu)
 *
 */
#ifndef ARJ_FIXED_HASH_KEY
#define ARJ_FIXED_HASH_KEY

/******************************************************
 *INCLUDES
 *****************************************************/
#include <glib.h>
#include "boolean.h"
#include "byte.h"
#include "buffers.h"
#include "PashDebug.h"

/******************************************************
 *DEFINES FOR CONSTANTS
 *****************************************************/

#define BAD_KEY ((keytype) 0)   // placeholder in key cache for a bad key
                                // note that this means that poly-A is always
                                // ignored

/******************************************************
 *DEFINES FOR INLINE FUNCTIONS
 *****************************************************/

/******************************************************
 *TYPEDEFS
 *****************************************************/


// use for weighting particular HBE hits, e.g. by k-mer frequency
typedef byte ScoreFactorType;
typedef bytebuff ScoreFactorList;

#define keytype guint32		// for hash keys


/******************************************************
 *EXTERNS
 *****************************************************/
extern const guint32 A2bits;
extern const guint32 T2bits;
extern const guint32 G2bits;
extern const guint32 C2bits;
extern const char twoBitsToBase[];

#define DEB_GET_KEY 0
/******************************************************
 * SCORE FACTOR FUNCTIONS
 ******************************************************
 ******************************************************/
/** Initialize the score factor list.*/
int initScoreFactorList(ScoreFactorList *scoreFactorList, FILE *scoreFactorFile, guint32 totalKeys);
/** Get score factor for a key.*/
ScoreFactorType getScoreFactor(const ScoreFactorList* const list, guint32 key);
/** Free score factor list memory.*/
void destroyScoreFactorList(ScoreFactorList *scoreFactorList);

#define maxScoreFactor 255
#define minScoreFactor 0

/******************************************************
 * PRIMARY FUNCTION DECLARATIONS
 *****************************************************/

/**   Fills in the preallocated string buffer  with the sequence corresponding to the key provided. */
boolean getSeqForKey(guint32 key, guint32 keyLen, char** seq);
/**   Fills in the preallocated string buffer  with the sequence corresponding to the key provided. */
boolean getSeqForKey1(guint32 key, guint32 keyLen, char*seq);
/******************************************************
 * SUPPLEMENTAL FUNCTION DECLARATIONS (not needed by most programs...gives access to some internals)
 *****************************************************/
 /******************************************************
 * FUNCTION: getNthBaseBits
 *
 * Usage   : getNthBaseBits(key, 5)
 * Function:
 *           Returns a gunint32 with just the two bits corresponding to the Nth
 *           base in the sequence set to 1. All other bits set to 0.
 * Returns : guint32
 * Args    :  guint32, guint32
 *
 * NOW INLINED FOR SPEED (WATCH YOUR *TYPES* THOUGH!): guint32 getNthBaseBits(guint32 key, guint32 n)
 */
#define getNthBaseBits(key, n) (((key)>>(2*(n))) & 0x03)

 /******************************************************
 * FUNCTION: setNthBaseBits
 *
 * Usage   : setNthBaseBits(key, 5, T2bits)
 * Function:
 *           Sets the 2 bits corresponding to the Nth base in the sequence to
 *           the provided 2 bits. Returns the new resultant key.
 * Returns : guint32
 * Args    : guint32, guint32, guint32
 *
 * NOW INLINED FOR SPEED (WATCH YOUR *TYPES* THOUGH!): guint32 setNthBaseBits(guint32 key, guint32 n, guint32 setBits);
 */
#define setNthBaseBits(key, n, setBits) ((setBits)!=0 ? ((key) | ((setBits) << (2*(n)))) : ((key) & ~(0x03 << (2*(n)))))

/******************************************************
 * FUNCTION: getNthBase
 *
 * Usage   : getNthBase(key, 5)
 * Function:
 *           Returns the Nth base (a char) stored in the key.
 * Returns : char
 * Args    : guint32, guint32
 *
 * NOW INLINED FOR SPEED (WATCH YOUR *TYPES* THOUGH!): unsigned int getNthBase(guint32 key, guint32 n);
 */
#define getNthBase(key, n) ((n) > 15 ? '\0' : twoBitsToBase[(getNthBaseBits((key), (n)))])

/** Prints numBitsToPrint from key on STDOUT, but in binary (1's & 0'). */
void printBinary(guint32 x, guint32 nBits);

/** Sets the 2 bits corresponding to the Nth base in the sequence to
 the provided letter, and returns the new resultant key.
 @returns the new key 
 @param key input key
 @param n base index
 @param base  new base
*/
static inline guint32 setNthBase(guint32 key, const guint32 n, const char base) {
	//	if(n > 15u) { return key; } already chexked this case
	
	switch(base)
	{
	case 'A' :
  case 'a' :
		return setNthBaseBits(key, n, A2bits);
	case 'T' :
  case 't' :
		return setNthBaseBits(key, n, T2bits);
	case 'G' :
  case 'g':
		return setNthBaseBits(key, n, G2bits);
	case 'C' :
  case 'c' :
		return setNthBaseBits(key, n, C2bits);
	default :
		return key;
	}
}

/** Fills key in with the 32bit value corresponding to seq.
@return  TRUE on success, FALSE otherwise.
@param seq sequence
@param keyPtr pointer to key
\note Usage   : getKeyForSeq(seq, &key)
\pre seq != NULL
*/
static inline boolean getKeyForSeq(const char* const seq, guint32* const keyPtr) {
	int ii = 0;
	//int len = 0;
	char* seqWalker = (char *)seq;
	guint32 lclKey = 0x00u;
  
	// (seq==NULL || ((len = strlen(seq))>16)) { return FALSE; }
	/*if (seq==NULL) { return FALSE; }
	else
  */
	/* loop through each character and set the key bits for it */
	//for(ii=0; ii<len; ii++, seqWalker++)
	xDEBUG(DEB_GET_KEY, fprintf(stderr, "converting kmer %s\n", seq));
	for(ii=0; *seqWalker != '\0'; ii++, seqWalker++) {
		lclKey = setNthBase(lclKey, ii, *seqWalker);
    xDEBUG(DEB_GET_KEY, fprintf(stderr, "[%d]\t%s\t%d\t%x\n", ii, seqWalker, lclKey, lclKey));
	}
	*keyPtr = lclKey;
	return TRUE;
}

#endif /* FixedHashKey.h */
