/*
Copyright (c) 2004 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/

/************************************************************************
 * Mask.c
 * utilities for dealing with sampling patterns
 * See "Pash Documentation for developers.doc" for more information about sampling patterns
 *
 ************************************************************************/

#include <stdio.h>
 #include <glib.h>
#include "pashtypes.h"
#include <stdlib.h>
#include "byte.h"
#include "Pattern.h"
#include "Mask.h"
#include "PashDebug.h"

#define DEB_MASK_OVERLAP_CONTRIBUTION 0

/** Generates the sampling pattern of weight 14 and length 22.
@param mask string holding the pattern
111#111##1#11##1#1111
*/
void fourteenOfTwentyOne(byte *mask)
{
	mask[ 0]=1;
	mask[ 1]=1;
	mask[ 2]=1;
	mask[ 3]=0;
	mask[ 4]=1;
	mask[ 5]=1;
	mask[ 6]=1;
	mask[ 7]=0;
	mask[ 8]=0;
	mask[ 9]=1;
	mask[10]=0;
	mask[11]=1;
	mask[12]=1;
	mask[13]=0;
	mask[14]=0;
	mask[15]=1;
	mask[16]=0;
	mask[17]=1;
	mask[18]=1;
	mask[19]=1;
	mask[20]=1;
}
/** Generates the sampling pattern of weight 13 and length 21.
@param mask string holding the pattern
*/
void thirteenOfTwentyone(byte *mask)
{
	mask[ 0]=1;
	mask[ 1]=1;
	mask[ 2]=1;
	mask[ 3]=0;
	mask[ 4]=1;
	mask[ 5]=1;
	mask[ 6]=0;
	mask[ 7]=1;
	mask[ 8]=1;
	mask[ 9]=0;
	mask[10]=0;
	mask[11]=0;
	mask[12]=1;
	mask[13]=1;
	mask[14]=0;
	mask[15]=1;
	mask[16]=0;
	mask[17]=1;
	mask[18]=0;
	mask[19]=1;
	mask[20]=1;
}



/** Generates the sampling pattern of weight 11 and length 18 used by PatternHunter (Ma et al 2002).
@param mask string holding the pattern
*/
void elevenOfEighteen(byte *mask)
{
	mask[0]= 1;
	mask[1]= 1;
	mask[2]= 1;
	mask[3]= 0;
	mask[4]= 1;
	mask[5]= 0;
	mask[6]= 0;
	mask[7]= 1;
	mask[8]= 0;
	mask[9]= 1;
	mask[10]= 0;
	mask[11]= 0;
	mask[12]= 1;
	mask[13]= 1;
	mask[14]= 0;
	mask[15]= 1;
	mask[16]= 1;
	mask[17]= 1;
} // elevenOfEighteen

/** Generate the 12-weight sampling pattern  111*1*11*1**11*111
described in "Good spaced seed for homology search" Choi et al, 20(7): 1053 Bioinformatics
@param mask string holding the pattern
*/
void twelveOfEighteen(byte *mask)
{
	mask[0]= 1;
	mask[1]= 1;
	mask[2]= 1;
	mask[3]= 0;
	mask[4]= 1;
	mask[5]= 0;
	mask[6]= 1;
	mask[7]= 1;
	mask[8]= 0;
	mask[9]= 1;
	mask[10]= 0;
	mask[11]= 0;
	mask[12]= 1;
	mask[13]= 1;
	mask[14]= 0;
	mask[15]= 1;
	mask[16]= 1;
	mask[17]= 1;
} // twelveOfEighteen

/** Generate the 9-weight sampling pattern 11*11*1*1***111
described in "Good spaced seed for homology search" Choi et al, 20(7): 1053 Bioinformatics
@param mask string holding the pattern
*/
void nineOfFifteen(byte *mask)
{
	mask[0]= 1;
	mask[1]= 1;
	mask[2]= 0;
	mask[3]= 1;
	mask[4]= 1;
	mask[5]= 0;
	mask[6]= 1;
	mask[7]= 0;
	mask[8]= 1;
	mask[9]= 0;
	mask[10]= 0;
	mask[11]= 0;
	mask[12]= 1;
	mask[13]= 1;
	mask[14]= 1;
} // nineOfFifteen


/** Generate the 8-weight sampling pattern 11**1**1*1*111
described in "Good spaced seed for homology search" Choi et al, 20(7): 1053 Bioinformatics
@param mask string holding the pattern
*/
void eightOfFourteen(byte *mask)
{
	mask[0]= 1;
	mask[1]= 1;
	mask[2]= 0;
	mask[3]= 0;
	mask[4]= 1;
	mask[5]= 0;
	mask[6]= 0;
	mask[7]= 1;
	mask[8]= 0;
	mask[9]= 1;
	mask[10]= 0;
	mask[11]= 1;
	mask[12]= 1;
	mask[13]= 1;
} // eightOfFourteen

/** Generate the 10-weight spaced seed 11*11***11*1*111
described in "Good spaced seed for homology search" Choi et al, 20(7): 1053 Bioinformatics
@param mask string holding the pattern
*/
void tenOfSixteen(byte *mask)
{
	mask[0]= 1;
	mask[1]= 1;
	mask[2]= 0;
	mask[3]= 1;
	mask[4]= 1;
	mask[5]= 0;
	mask[6]= 0;
	mask[7]= 0;
	mask[8]= 1;
	mask[9]= 1;
	mask[10]= 0;
	mask[11]= 1;
	mask[12]= 0;
	mask[13]= 1;
	mask[14]= 1;
	mask[15]=1;
} // tenOfSixteen


/** Generate the sampling pattern based on the user's specification.
\note If a user specifies pattern length and weight but not the pattern,
   assign a default pattern if available, or pick a random one.
@return  0 (always succeeds)
@param maskLen - size of the sampled word, like the "model" of patternhunter
@param keyLen - number of bases used, like the "weight" of patterhunter
\par SIDE EFFECTS
 <ul> <li>  memory allocated for mask
 <li> values of mask initialized
 <li> for keyLen of 9, 10, 11, 12, and 13 use patterns from PatternHunter or Choi et al.
 <li> if keyLen==0 or maskLen==0 use patternhunter 11 of 18 (keyLen and  maskLen are changed appropriately) </ul>
*/
int setMask(Mask *mask)
{
	if( (mask->keyLen==13&&mask->maskLen==21) ||
	    (mask->keyLen==13&&mask->maskLen==0)  ||
	    (mask->keyLen==0&&mask->maskLen==21)) {
		mask->keyLen=13;
		mask->maskLen=21;
		thirteenOfTwentyone(mask->mask);
        } else if ( mask->keyLen==14 && mask->maskLen!=14) {
		mask->keyLen=14;
		mask->maskLen=21;
		fourteenOfTwentyOne(mask->mask);
	} else if ( mask->keyLen==12 && mask->maskLen!=12) {
		mask->keyLen=12;
		mask->maskLen=18;
		twelveOfEighteen(mask->mask);
	} else if ( mask->keyLen==9 && mask->maskLen!=9) {
		mask->keyLen=9;
		mask->maskLen=15;
		nineOfFifteen(mask->mask);
	} else if ( mask->keyLen==8 && mask->maskLen!=8) {
		mask->keyLen=8;
		mask->maskLen=14;
		eightOfFourteen(mask->mask);
	}  else 	if ( mask->keyLen==10 && mask->maskLen!=10) {
		mask->keyLen=10;
		mask->maskLen=16;
		tenOfSixteen(mask->mask);
	} else 	if( (mask->keyLen==0 || mask->maskLen==0) ||
		        (mask->keyLen==11 && mask->maskLen==18)) {
		mask->keyLen=11;
		mask->maskLen=18;
		elevenOfEighteen(mask->mask);
	} else {
		// generate random pattern to use for the job.  This is untested, and a user has reported some errors when attempting this
		getPattern(mask->mask, mask->keyLen, mask->maskLen);
	}
	computeOverlapContribution(mask);
	return 0;
}

/**Given a sampling mask, list of matching words, and word offset,
  determines the number and positions of matching bases detected
 @return 0 if error, otherwise total number of matching bases detected
 \par ERROR CONDITIONS
 <ul> <li> NULL pointers        (0)
 <li> wordOffsetGap < 1       (0)
 <li> maskLen < 1          (0)
 </ul>
 \pre pointers allocated to hold enough data, in particular *matchingBases
 holds (greatest element of hits[])*wordOffsetGap + maskLen - 1
 \par SIDE EFFECTS
 <ul> <li> list of matching positions is built into the byte array
 <li> position of the first matching position is stored in *start
 <li> position of last matching position is stored in *stop
 @param wordOffsetGap - distance between the start of consecutive words
 \note follows convention that first possible word is #0
*/
wordtype countMatches(const guint32 wordOffsetGap, const wordtype numHits,
                      wordtype * const hits, Mask *mask, byte *matchingBases,
                      const guint32 length, wordtype *start, wordtype *stop)
{
	guint32 ii=0,jj=0;
	wordtype numBases=0;

	// initialize matching bases array, do before error checking because this is
	// desired outcome if numHits==0
	for(ii=0;ii<length;ii++) matchingBases[ii]=0;

	//error checking
	if(wordOffsetGap<1) return 0;
	if(numHits<1) return 0;
	if(hits==NULL) return 0;
	if(mask->maskLen<1) return 0;
	if(matchingBases==NULL) return 0;

	for(ii=0;ii<numHits;ii++)
		for(jj=0;jj<mask->maskLen;jj++)
			if(mask->mask[jj])
				if((hits[ii]*wordOffsetGap + jj) < length)  // don't overrun array
					matchingBases[ hits[ii]*wordOffsetGap + jj ] = 1;

	for(numBases=0,ii=0, *start=length;ii<length;ii++)
		if(matchingBases[ii])
		{
			numBases++;
			*stop=ii;
			if(*start==length) *start=ii;
		}
	return numBases;
} // countMatches

/** Copy a sampling pattern
@param dest destination mask
@param source source mask
@return 0 if OK, -1 if error
*/


/** Compute the new bases a kmer adds to a previous kmer based on the number of overlapped bases, for
   all range of overlaps
   @param mask mask
*/
void computeOverlapContribution(Mask *mask) {
	int overlap, i;
	int matchingBases[2*MAX_MASK_LEN];
	int matchCount;

	xDEBUG(DEB_MASK_OVERLAP_CONTRIBUTION, fprintf(stderr, "computeOverlapContribution START\n"));
	mask->maskOverlapContribution[0] = mask->keyLen;
	for (overlap=1; overlap<mask->maskLen; overlap++) {
		for(i=0; i<2*mask->maskLen; i++) {
			matchingBases[i] = 0;
		}
		for (i=0; i<mask->maskLen; i++) {
			if (mask->mask[i]==1) {
				matchingBases[i] = 1;
			}
		}
		for (i=0; i<mask->maskLen; i++) {
			if (mask->mask[i]==1) {
				matchingBases[i+overlap] = 1;
			}
		}
		// now count total number of matches
		for (matchCount=0, i = 0; i<2*mask->maskLen; i++) {
			if (matchingBases[i]==1) {
				matchCount +=1;
			}
		}
		xDEBUG(DEB_MASK_OVERLAP_CONTRIBUTION, fprintf(stderr, "[%d] matchingBases: %d contribution %d\n",
																									overlap, matchCount, matchCount-mask->keyLen));
		mask->maskOverlapContribution[overlap] = matchCount - mask->keyLen;
	}
	xDEBUG(DEB_MASK_OVERLAP_CONTRIBUTION, fprintf(stderr, "computeOverlapContribution STOP\n"));
}
