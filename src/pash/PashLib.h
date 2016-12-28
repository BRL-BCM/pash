/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef PASHLIB_H
#define PASHLIB_H

#include <glib.h>		// standard data types
#include <stdio.h>		// FILE
#include <string.h>		// strncpy
#include "Mask.h"		// Sampling patterns
#include "FastaUtil.h"
#include "FastQUtil.h"
#include "pashtypes.h"
#include "byte.h"
#include "IgnoreList.h"


/** Updating this as bugs are fixed, features improved, etc.*/
#define MAX_FILE_NAME_SIZE 2048
#define FASTA_FILES 2
#define VERTICAL_FASTA_FILE 0
#define HORIZONTAL_FASTA_FILE 1
#define DEFAULT_HIVE_HASH_MEMORY 4096 
#define DEFAULT_NUMBER_OF_DIAGONALS 500
#define DEFAULT_MIN_SCORE 40
#define DEFAULT_WORD_OFFSET 6

enum SensitivityMode {HighSensitivity, MediumSensitivity, LowSensitivity, FastSensitivity, UserDefinedSensitivity };

typedef struct {
	/// Vertical sequence file (typically reads/chromosomes/genome).
	char verticalFile[MAX_FILE_NAME_SIZE+1];
	/// Fasta Utility for the vertical sequence.
	//FastaUtil* fastaUtilVertical;
	/// FastQ Utility for the vertical sequence.
	PashFastqUtil* verticalFastqUtil;
	guint32 lastVerticalSequenceMapped;
	SequenceInfo *verticalSequencesInfos;
	void* hiveHash;
	// Bisulfite sequencing support
	int bisulfiteSequencingMapping;
	// dna meth support
	int reverseStrandDnaMethMapping;
	char actualChromName[MAX_FILE_NAME_SIZE+1];
	guint32 reverseComplementSequenceLength ;


	/// Horizontal sequence file (typically chromosomes/genome).
	char horizontalFile[MAX_FILE_NAME_SIZE+1];
	/// Fasta Utility for the horizontal sequence.
	FastaUtil* fastaUtilHorizontal;
	/// Location of scratch directory (hopefully soon to become obsolete).
	char scratchDirectory[MAX_FILE_NAME_SIZE+1];
	/// Output file
	char outputFile[MAX_FILE_NAME_SIZE+1];
	/** Whether the sampling pattern has been defined yet --
	  default sampling patterns can be used in certain cases, in which case
	  the sampling pattern is not set until after parameter parsing.*/
	char ignoreListFile[MAX_FILE_NAME_SIZE];
	int useIgnoreList;
	IgnoreList ignoreList;
	/// Number of Pash diagonals.
	guint32 numberOfDiagonals;
	/// Minimum score.
	guint32 minScore;
	/// Vertical sequence word offset.
	guint32 wordOffset;
	/// Pash sensitivity mode
	SensitivityMode sensitivityMode;
	/// percentOfKmers to keep
	guint32 keepHashedKmersPercent;
	/// Sampling pattern used by Pash.
	Mask mask;
	/// Flag whether the sampling pattern was defined.
	gboolean isMaskDefined;
	/// Avoid redundancy when doing self comparison -- any match found below the main diagonal is ignored.
	gboolean selfComparison;
	/// Hive hash memory limit.
	guint32 hiveHashMemoryLimit;
	/// List of kmers to ignore.
	/// Flag whether the user desires gzipped output.
	FILE* outputFilePtr;
	gboolean useGzippedOutput;
	guint32 maxMappings;
	double topPercent;
} PashParameters;

typedef struct {
	/// Hive hash.
	void *hiveHash;
	/// last chunk id
	guint32 lastSequenceId;
	/// offset of sequence buffer in the actual sequence
	guint32 offsetOfSequenceBufferInRealSequence;
	/// number of chunks in the current sequence.
	guint32 numberOfChunksInCurrentSequence;
	/// index of current chunk in the current sequence
	guint32 currentSequenceChunk;
} SequenceHash;

/// Parse command-line options and setup the Pash parameters.
PashParameters* parseCommandLine(int argc, char**argv);
/// Print Pash parameters.
void PashUsage();
/// Hash vertical sequence until the hive hash size reaches a user-specified limit.
int sizeCurrentVerticalSequencesBatch(PashParameters* pashParams);
int hashCurrentVerticalSequencesBatch(PashParameters* pashParams);
/// Scan the horizontal sequence (typically chromosome/genome) agains the hivehash.
int scanHorizontalSequence(PashParameters* pashParams, SequenceHash* sequenceHash);
/// Initialize the sequence hash
SequenceHash *initSequenceHash();

/// Print current time
void printNow();


#endif
