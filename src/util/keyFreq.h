#ifndef KJK_KEYFREQ_H
#define KJK_KEYFREQ_H

#include <glib.h>

#include "IgnoreList.h"
#include "FixedHashKey.h"

// data modes

/** SEQ means currently within sequence.*/ 
#define SEQ 0
/** NEWLN means a newline has just been encountered, so a defline could be starting.*/ 
#define NEWLN 1
/** DEFLINE means currently within a defline, can be cleared by hitting a newln.*/
#define DEFLINE 2
/** START means an input file is just starting (unused).*/
#define START 3

/** Kmer frequency entry.*/
typedef struct kmerfreqentry {
   /** kmer index.*/
	 guint32 kmerIndex;
	 /** kmer frequency.*/
   guint32 kmerFrequency; 
} KmerFreqEntry;


/** Print usage information.*/
void printUsage();

/** Parse arguments.*/
int parse(int argc, char ** argv);

/** Read input file.*/
static inline int readInputFile(FILE *input, int *mode,
	charbuff *readBuff, charbuff *seqBuff, 
	int patternLen);

/** Converts base to a 0-4 number.*/
static inline int baseToNum(char c);

/** Checks if the letter corresponds to a base.*/
static inline int isBase(char c);

/** Print the program arguments*/
void printParseResults(int argc, char ** argv, int firstFile);

/** Process input file.*/
void processInputFile(const char *inputName, guint32* freq, KmerFreqEntry *kfreq);

/** Compute reverse complement.*/
static inline void revComp(const char *seq, char *rev, int len);

/** Compare two kmer frequency entries based on frequency.*/
int compareKmerFrequencyEntry(const void *p1, const void*p2);

#endif // KJK_KEYFREQ_H
