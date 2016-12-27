#ifndef __PAsh__To__SAM__Converter__H__
#define __PAsh__To__SAM__Converter__H__

#include <glib.h>
#include "FastQUtil.h"

struct SequenceInfo {
	guint32 sequenceId;
	guint32 numMappings;
};

class Pash2SAMConverter {
  PashFastqUtil *pashFastQUtil;
  char fastqFile[MAX_FILE_NAME];
  char referenceSequencesFile[MAX_FILE_NAME];
  char mappingFile[MAX_FILE_NAME];
  char outputFile[MAX_FILE_NAME];
  GHashTable* sequenceHashTable;
  GHashTable *referenceSizesHashTable;
	SequenceInfo* sequenceInfos;
	FILE* outPtr;
	char sampleName[MAX_LINE_LENGTH], centerName[MAX_LINE_LENGTH];
	int bisulfiteSeqFlag;
public:
  Pash2SAMConverter();
  ~Pash2SAMConverter();
  int parseParameters(int argc, char* argv[]);
  int loadReferenceSizes();
  int convertPash2SAM();
private:
  void usage();
	/** Performs the actual conversion to SAM, mapping by mapping.*/
	int convertMappings();
	/** Build a hash from read name to read id.*/
	int prepareReadsHash();
	/** Generate extended CIGAR information for a read mapping.*/
	int buildCigarString(char* cigar, guint32 numBlocks,
											 guint32 readStart, guint32 readStop,
											 guint32* blockSizesArray,
											 guint32* horizontalStartsArray,
											 guint32* verticalStartsArray,
											 guint32 readLength);
	/** Output unmapped reads.*/
	int outputUnmappedReads();
  /** Obtain the reverse complement for a read.*/
  int getReverseComplement(const char* fwdRead, char* revComplementRead, int maxReadLength);
  /** Obtain the reverse of a string.*/
  int getReverseString(const char* source, char* destination, int maxReadLength);
  /** Backconvert a bisulfite sequencing read, for the purpose of using it w/ generic tools such as samtools or picard to call genotypes.*/
  int bisulfiteBackConvertRead(char* backConvertedRead, char strand,
															guint32 numT, char* tBases,
															guint32 numBlocks,
											 				guint32* blockSizesArray,
															guint32* horizontalStartsArray,
															guint32* verticalStartsArray, const char* fwdRead);
};

#endif

