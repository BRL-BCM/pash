/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef _PASH_FASTQ_UTIL___H__
#define _PASH_FASTQ_UTIL___H__

#include <glib.h>
#include "someConstants.h"
#include "SequencePool.h"

typedef enum {FastaOnly, FastaAndQualityScores} ReadsSequenceType;


class PashFastqUtil {
  char fastqFile[MAX_FILE_NAME];
  guint32 numberOfSequences;
  guint32 availableEntries;
  char** sequenceNames;
  char** forwardSequences;
  char** reverseSequences;
  char** qualities;
  SequencePool* sequencePool;
  ReadsSequenceType readsSequenceType;
  int reverseSequencesAvailable;
  
  public:
    PashFastqUtil(char* fileName, ReadsSequenceType sequenceType);
    ~PashFastqUtil();
    int loadSequences(int loadReverseComplement);
    const char* retrieveSequence(guint32 sequenceId);
    const char* retrieveRevComplementSequence(guint32 sequenceId);
    const char* retrieveQualityScores(guint32 sequenceId);
    const char* retrieveDefName(guint32 sequenceId);
    guint32 getNumberOfSequences();
  private:
    char baseComplement(char b);      
};



#endif

