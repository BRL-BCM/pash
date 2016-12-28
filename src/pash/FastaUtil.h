/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef __fasta_Util_h____
#define __fasta_Util_h____

#include "SequencePool.h"
#include "SequenceInfo.h"

#define RAW_BUFFER_SIZE 16000
#define MAX_DEFNAME_SIZE 200
#define SEQUENCE_BUFFER_SIZE 40000
#define MAX_FILE_NAME 2048
#define INIT_SEQUENCE_INFO 50000

typedef enum {DefLine, Sequence, DetermineLineType, Comment, Unknown} ParserStates;

typedef struct {
    SequencePool* seqPool;
    /// sequence information for all FASTA/fof sequences
    SequenceInfo* sequencesInformation;
    /// array of FASTA file(s) names
    char** fileArray;
    /// number of FASTA files
    int numFiles;
    /// index of current FASTA file
    int currentFileIndex;
    /// file pointer for the current FASTA file
    FILE* currentFile;
    /// number of FASTA sequences
    guint32 numberOfSequences;
    /// number of allocated sequences
    guint32 numberOfAllocatedSequences;
    /// index of current FASTA sequence
    guint32 currentSequenceIndex;
    /// index of current line in FASTA file
    guint32 currentLine;
    /// raw read buffer
    char rawBuffer[RAW_BUFFER_SIZE];
    /// characters read in the raw buffer
    guint32 readChars;
    /// characters attempted to read in the raw buffer
    guint32 bufferRead;
    /// position in the read buffer
    guint32 positionInRawBuffer;
    /// did we reach EOF in current file ?
    int isEof;
    /// did we consume the eof for current file ?
    int eofConsumed;
    /// current defline buffer
    char deflineBuffer[MAX_DEFNAME_SIZE+1];
    /// index in the defline buffer
    guint32 currentDeflineBufferPos;
    /// maximum defline size
    guint32 maximDeflineSize;
    /// maximum sequence buffer size
    guint32 maximSequenceSize;
    /// current nucleotide sequence buffer
    char sequenceBuffer[SEQUENCE_BUFFER_SIZE+1];
    /// index in the sequence buffer
    guint32 currentSequenceBufferPos;
		/// position in the actual sequence
		guint32 currentActualSequencePos;
    /// boolean flag indicating whether we finished the current sequence
    int isCurrentSequenceFinished;
    /// Parser state.
    ParserStates parserStatus;
    /// Parsing finished indicator.
    int parsingDone;
} FastaUtil;

int parseFastaFileFirstPassFastaUtil(FastaUtil* fastaUtil);
int addFileFastaUtil(char *fileName, FastaUtil* fastaUtil);
FastaUtil* initFastaUtil(char *fileName);
int rewindFastaUtil(FastaUtil* fastaUtil);
int nextChunkFastaUtil(FastaUtil* fastaUtil);
void fastaUtilKeepPartialBuffer(FastaUtil* fastaUtilKeepPartialBuffer, int basesToKeep);

#endif
