/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef KJK_MASK_H
#define KJK_MASK_H

/************************************************************************
 * Mask.h
 * utilities for dealing with sampling patterns
 * See "Pash Documentation for developers.doc" for more information about sampling patterns
 *
 ************************************************************************/

#include <glib.h>
#include "byte.h"
#include "pashtypes.h"

#define MAX_MASK_LEN 50		// allows us to allocate a fixed amount of space
// to store masks in a bindump section header; increase arbitrarily if there is ever a need

// mask record (stored in bindump section header): maskLen (guint32), keyLen
// (guint32), mask (MAX_MASK_LEN bytes)
#define MASK_RECORD_SIZE (sizeof(guint32)*2 + sizeof(byte)*MAX_MASK_LEN)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct masktype
{
	byte mask[MAX_MASK_LEN];  // array of 0s and 1s corresponding to unsampled and sampled positions, respectively
	guint32 maskLen;  // length of pattern (number of bases between the first and last 1)
	guint32 keyLen;  // number of sampled positions in the pattern
	/// Array that describes how many new bases a kmer adds to a previous kmer based on the number of overlapped bases.
	guint32 maskOverlapContribution[MAX_MASK_LEN];
} Mask;

/** Generate the sampling pattern based on the user's specification.*/
int setMask(char const * pattern, Mask *mask);

/** Copy a sampling pattern.*/
static inline int copyMask(Mask *dest, Mask *source) {
	guint32 ii;
	//*check for errors
	if(source->maskLen<1) return 1;	
	for(ii=0;ii<source->maskLen;ii++) {
		dest->mask[ii]=source->mask[ii];
	}
	dest->maskLen=source->maskLen;
	dest->keyLen=source->keyLen;
	return 0;
}

/** Compute the new bases a kmer adds to a previous kmer based on the number of overlapped bases, for
   all range of overlaps.
 */
void computeOverlapContribution(Mask *mask);  

/** Given a sampling mask, list of matching words, and word offset, determines the number and positions of matching bases detected.*/
wordtype countMatches(const guint32 wordOffsetGap, const wordtype numHits,
		wordtype * const hits, Mask *mask, byte *matchingBases,
		const guint32 length, wordtype *start, wordtype *stop);

#ifdef __cplusplus
}
#endif

#endif
