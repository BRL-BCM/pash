/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


/******************************************************
 * NAME: FixedHashKey.c
 *
 * DESCRIPTION:
 * This module provides functionality for converting a sequence of up to
 * 16 letters from the alphabet "AGTCagtc" to a 32-bit int (the "key") and
 * furthermore provides functions for accessing the contents of the key as
 * bases, not bits.
 *
 * REQUIRES:
 * �glib library (GTK+, multiplatform C library)
 *  boolean.h and byte.h for data-type standardization.
 * *
 * AUTHOR: Andrew R. Jackson (andrewj@bcm.tmc.edu)
 *
 */
/******************************************************
 *INCLUDES
 *****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <string.h>
#include "boolean.h"
#include "byte.h"
#include "error.h"
#include "FixedHashKey.h"


/******************************************************
 *EXTERN'D CONSTANTS
 *****************************************************/
const guint32 A2bits = 0x00u;
const guint32 T2bits = 0x01u;
const guint32 G2bits = 0x02u;
const guint32 C2bits = 0x03u;
const char twoBitsToBase[] = { 'A', 'T', 'G', 'C' };

/******************************************************
 *EXTERN'D FUNCTIONS
 *****************************************************/


/** Fills in the preallocated string buffer (pointer to a char*)
    with the sequence corresponding to the key provided. Sequence will
    be all uppercase.
  \note It makes sense to allocat 17 bytes for the buffer...16bases, plus one for "\0".
 
  \note You must provide the key length, which is the number of bases
      stored in the key since 32bits can fit 16 bases, but you might be
      interested in 9 mers.
  @return TRUE if succeeds, FALSE otherwise.
  @param key key
  @param keyLen key length
  @param seq sequence
*/
boolean getSeqForKey(const guint32 key, const guint32 keyLen, char** seq) {
	guint32 ii = 0;
	char* seqWalker = NULL;

	if(keyLen==0u || (keyLen>16u)) { return FALSE; }
	else
	{
		/* allocate space to store keyLen chars +1 for the '\0' */
		if((seqWalker = *seq = (char*)calloc(keyLen+1, sizeof(char)))==NULL) { /* FAILED */ MEM_ALLOC_CROAK(); }
		/* loop through each pair of bit bases and get the letter for it */
		for(ii=0; ii<keyLen; ii++, seqWalker++)
		{
			*seqWalker = getNthBase(key, ii);
		}
		return TRUE;
	}
}


/** Fills in the preallocated string buffer 
    with the sequence corresponding to the key provided. Sequence will
    be all uppercase.
  \note It makes sense to allocate 17 bytes for the buffer...16bases, plus one for "\0".
 
  \note You must provide the key length, which is the number of bases
      stored in the key since 32bits can fit 16 bases, but you might be
      interested in 9 mers.
  @return TRUE if succeeds, FALSE otherwise.
  @param key key
  @param keyLen key length
  @param seq sequence
*/
boolean getSeqForKey1(const guint32 key, const guint32 keyLen, char*seq) {
	guint32 ii = 0;
	char* seqWalker = NULL;

	if(keyLen==0u || (keyLen>16u)) { return FALSE; }
	else
	{
		/* allocate space to store keyLen chars +1 for the '\0' */
		seqWalker = seq;
    if(seqWalker == NULL) { /* passed a NULL buffer */ MEM_ALLOC_CROAK(); }
		/* loop through each pair of bit bases and get the letter for it */
		for(ii=0; ii<keyLen; ii++, seqWalker++)
		{
			*seqWalker = getNthBase(key, ii);
		}
		return TRUE;
	}
}

/** Initialize the score factor list
@return 0 if success, 1 otherwise
@param list score factor list
@param scoreFactorFile file handle containing the score factors
*/
int initScoreFactorList(ScoreFactorList *list, FILE *scoreFactorFile,
	guint32 totalKeys) {
	ScoreFactorType test=0;

	// allocate memory, initialize vars
	list->capacity=totalKeys;
	list->buff=NULL;
	list->buff=(ScoreFactorType*) malloc(list->capacity * sizeof(ScoreFactorType));
	if(list->buff==NULL)
	{
		fprintf(stderr,"memory allocation for score factor list (%ld items) failed\n", list->capacity * sizeof(ScoreFactorType));
		return 1;
	}
	list->content=0;
	list->pos=0;

	// read input file
	list->content = fread((void*) list->buff, sizeof(ScoreFactorType), list->capacity, scoreFactorFile);
	if(list->content != list->capacity)
	{
		fprintf(stderr,"fatal: premature end of score factor list\nread %d items, expected %d\n", list->capacity, list->content);
		return 1;
	}

	fread((void*) &test, sizeof(ScoreFactorType), 1, scoreFactorFile);  // try to read one element to make sure file has ended
	if(!feof(scoreFactorFile))
	{
		fprintf(stderr,"fatal: score factor list longer than expected (expected %d items)\n", list->capacity);
		return 1;
	}
	return 0;
}

/** Get score factor for a key
@return score factor for the given key
@param list score factor list
@param key input key
*/
ScoreFactorType getScoreFactor(const ScoreFactorList* const list, guint32 key)
{
	if(key >= list->capacity) return 0;
	return (list->buff)[key];
}


/** Free score factor list memory
@param list score factor list
*/
void destroyScoreFactorList(ScoreFactorList *list)
{
	list->buff=NULL;
	list->capacity=0;
	list->content=0;
	list->pos=0;
}

/** Prints numBitsToPrint from key on STDOUT, but in binary (1's & 0'); 
   bit 0 is the LSb.
@param x  input key
@param nBits number of bits to print
 */
void printBinary(const guint32 x, const guint32 nBits) {
	int ii;
	for(ii=nBits-1u; ii>=0; ii--)
	{
		printf("%i", (x>>ii) & 01);
	}
}

/******************************************************
 *MAIN, FOR BASIC UNIT TESTING
 *****************************************************/
#ifdef TEST_FIXED_HASH_KEY
int main(int argc, char* argv[])
{
	char* s1 = "AGTC";
	char* s2 = "AAAAAAAAAAAAAAAA";
	char* s3 = "GGGGGGGGGGGGGGGGGGG";
	char* s4 = "GGCTCTCTATTTCCAA";
	char* s5 = "TAATCGGCCCATG";
	char* getSeqBuff = NULL;
	guint32 k1 = 0x1B;
	guint32 k2 = 0x00;
	guint32 k3 = 0x00;
	guint32 k4 = 0x00;
	guint32 k5 = 0x00;
	boolean retVal = FALSE;

	printf("1st base bits of 00011011: %X\t", getNthBaseBits(k1, 0)); printBinary(getNthBaseBits(k1, 0), 32); putchar('\n');
	printf("2nd base bits of 00011011: %X\t", getNthBaseBits(k1, 1)); printBinary(getNthBaseBits(k1, 1), 32); putchar('\n');
	printf("3rd base bits of 00011011: %X\t", getNthBaseBits(k1, 2)); printBinary(getNthBaseBits(k1, 2), 32); putchar('\n');
	printf("4th base bits of 00011011: %X\t", getNthBaseBits(k1, 3)); printBinary(getNthBaseBits(k1, 3), 32); putchar('\n');
	printf("\n");
	printf("set 1st base bits of 00011011 to 00 (makes 00011000): %X\t", setNthBaseBits(k1, 0, 0x0)); printBinary(setNthBaseBits(k1, 0, 0x0),32); putchar('\n');
	printf("set 2st base bits of 00011011 to 10 (makes 00011011): %X\t", setNthBaseBits(k1, 1, 0x2)); printBinary(setNthBaseBits(k1, 1, 0x2),32); putchar('\n');
	printf("set 3rd base bits of 00011011 to 11 (makes 00111011): %X\t", setNthBaseBits(k1, 2, 0x3)); printBinary(setNthBaseBits(k1, 2, 0x3),32); putchar('\n');
	printf("set 4th base bits of 00011011 to 01 (makes 01011011): %X\t", setNthBaseBits(k1, 3, 0x1)); printBinary(setNthBaseBits(k1, 3, 0x1),32); putchar('\n');
	printf("\n");
	printf("get 1st base of 00011011 ('C'): %c\t", getNthBase(k1, 0)); putchar('\n');
	printf("get 2nd base of 00011011 ('G'): %c\t", getNthBase(k1, 1)); putchar('\n');
	printf("get 3rd base of 00011011 ('T'): %c\t", getNthBase(k1, 2)); putchar('\n');
	printf("get 4th base of 00011011 ('A'): %c\t", getNthBase(k1, 3)); putchar('\n');
	printf("\n");
	getSeqForKey(k1, 4, &getSeqBuff); printf("translated key (CGAT): %s\n", getSeqBuff);
	printf("set 1st *base* of 00011011 to A (makes 00011000): %X\t", setNthBase(k1, 0, 'A')); printBinary(setNthBase(k1, 0, 'A'),32); putchar('\t');

	getSeqForKey(setNthBase(k1, 0, 'A'), 4, &getSeqBuff); printf("translated key (AGTA): %s\n", getSeqBuff);
	printf("set 2st *base* of 00011011 to G (makes 00011011): %X\t", setNthBase(k1, 1, 'G')); printBinary(setNthBase(k1, 1, 'G'),32); putchar('\t');

	getSeqForKey(setNthBase(k1, 1, 'G'), 4, &getSeqBuff); printf("translated key (CGTA): %s\n", getSeqBuff);
	printf("set 3rd *base* of 00011011 to C (makes 00111011): %X\t", setNthBase(k1, 2, 'C')); printBinary(setNthBase(k1, 2, 'C'),32); putchar('\t');

	getSeqForKey(setNthBase(k1, 2, 'C'), 4, &getSeqBuff); printf("translated key (CGCA): %s\n", getSeqBuff);
	printf("set 4th *base* of 00011011 to T (makes 01011011): %X\t", setNthBase(k1, 3, 'T')); printBinary(setNthBase(k1, 3, 'T'),32); putchar('\t');

	getSeqForKey(setNthBase(k1, 3, 'T'), 4, &getSeqBuff); printf("translated key (CGTT): %s\n", getSeqBuff);
	printf("set 4th *base* of 00011011 to N (makes 00011011): %X\t", setNthBase(k1, 3, 'N')); printBinary(setNthBase(k1, 3, 'N'),32); putchar('\t');

	getSeqForKey(setNthBase(k1, 0, 'N'), 4, &getSeqBuff); printf("translated key (CGAT): %s\n", getSeqBuff);
	printf("\n");

	k1 = 0x0;
	if((retVal = getKeyForSeq(s1, &k1))==TRUE)
	{
		printf("Key for %s: (%X)\t", s1, k1); printBinary(k1, 32); putchar('\t');
		getSeqForKey(k1, 4, &getSeqBuff); printf("translated key (%s): %s\n", s1, getSeqBuff);
	}
	else { printf("Key for %s FAILED\n", s1); }
	if((retVal = getKeyForSeq(s2, &k2))==TRUE)
	{
		printf("Key for %s: (%X)\t", s2, k2); printBinary(k2, 32); putchar('\t');
		getSeqForKey(k2, 16, &getSeqBuff); printf("translated key (%s): %s\n", s2, getSeqBuff);
	}
	else { printf("Key for %s FAILED\n", s2); }
	if((retVal = getKeyForSeq(s3, &k3))==TRUE)
	{
		printf("Key for %s: (%X)\t", s3, k3); printBinary(k3, 32); putchar('\t');
		getSeqForKey(k3, 19, &getSeqBuff); printf("translated key (%s): %s\n", s3, getSeqBuff);
	}
	else { printf("Key for %s FAILED\n", s3); }
	if((retVal = getKeyForSeq(s4, &k4))==TRUE)
	{
		printf("Key for %s: (%X)\t", s4, k4); printBinary(k4, 32); putchar('\t');
		getSeqForKey(k4, 16, &getSeqBuff); printf("translated key (%s): %s\n", s4, getSeqBuff);
	}
	else { printf("Key for %s FAILED\n", s4); }
	if((retVal = getKeyForSeq(s5, &k5))==TRUE)
	{
		printf("Key for %s: (%X)\t", s5, k5); printBinary(k5, 32); putchar('\t');
		getSeqForKey(k5, 13, &getSeqBuff); printf("translated key (%s): %s\n", s5, getSeqBuff);
	}
	else { printf("Key for %s FAILED\n", s5); }

	printf("\n...about to croak due to memory (not really, just a test)...\n\n");
	MEM_ALLOC_CROAK();
}


#endif
