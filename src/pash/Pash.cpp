/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


/*******************************************************
 * NAME: Pash.c
 *
 *
 * AUTHOR: Cristian Coarfa (coarfa@bcm.edu), based on work by
 * Ken Kalafus (kkalafus@bcm.tmc.edu) and Andrew Jackson (andrewj@bcm.tmc.edu).
 *
 * FIXES AND MINOR CHANGES: Piotr Pawliczek
*/


#include <sys/types.h>
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include "PashLib.h"
#include "PashDebug.h"

#define DEB_MAIN 1 

int runPash(int argc, char** argv);

int main(int argc, char** argv) {
  int result;
  result = runPash(argc, argv);
  return result;
}

int runPash(int argc, char** argv) {
  PashParameters *pashParams;
  pashParams = parseCommandLine(argc, argv);

  printNow();
  pashParams->outputFilePtr = fopen(pashParams->outputFile, "wt");
  if (pashParams->outputFilePtr==NULL) {
    fprintf(stderr, "could not open temporary output file %s\n", pashParams->outputFile);
    fflush(stderr);
    exit(2);
  }

  pashParams->verticalFastqUtil = new PashFastqUtil(pashParams->verticalFile, FastaAndQualityScores);
  pashParams->verticalFastqUtil->loadSequences(1);

  fprintf(stderr, "initialized  vertical sequence util\n");
  printNow();
  pashParams->fastaUtilHorizontal = initFastaUtil(pashParams->horizontalFile);
  fprintf(stderr, "initialized  horizontal sequence util\n");
  printNow();

  // print SAM header
  fprintf(pashParams->outputFilePtr, "@HD\tVN:1.0\n");
  fprintf(pashParams->outputFilePtr, "@PG\tID:pash3\tPN:Pash\tVN:3.01.03\n");
  //fprintf(pashParams->outputFilePtr, "@RG\tID:--\tCN:BRL\n");
  for ( unsigned i = 1;  i <= pashParams->fastaUtilHorizontal->numberOfSequences;  ++i ) {
    char const * name = pashParams->fastaUtilHorizontal->sequencesInformation[i].sequenceName;
    if ( strncmp(name, "#RC.pash.", 9) == 0 ) continue;
    fprintf(pashParams->outputFilePtr, "@SQ\tSN:%s\tLN:%u\n", name, pashParams->fastaUtilHorizontal->sequencesInformation[i].sequenceLength);
  }
  
  SequenceHash * horizontalSequenceHash = initSequenceHash();
  pashParams->lastVerticalSequenceMapped = 0;
  guint32 numberOfVerticalReads = pashParams->verticalFastqUtil->getNumberOfSequences();
  pashParams->verticalSequencesInfos = (SequenceInfo*) malloc((numberOfVerticalReads+1)*sizeof(SequenceInfo));
  if (pashParams->verticalSequencesInfos==NULL) {
    fprintf(stderr,"could not allocate memory for the reads\n");
    exit(2);
  }
  // load as much stuff into memory as specified by user
  xDEBUG(DEB_MAIN, fprintf(stderr, "starting hashed the vertical sequence\n"));
  sizeCurrentVerticalSequencesBatch(pashParams);
  hashCurrentVerticalSequencesBatch(pashParams);
  xDEBUG(DEB_MAIN, fprintf(stderr, "done hashed the vertical sequence\n"));
  fprintf(stderr, "hashed the vertical sequence\n");
  printNow();
  fflush(stderr);
  // pass the hive hash to the horizontal scan sequence
  horizontalSequenceHash->hiveHash = pashParams->hiveHash;
  fprintf(stderr, "starting horizontal sequence\n");
  printNow();
  fflush(stderr);
  // run genome against it, perform anchoring/alignment, write output
  scanHorizontalSequence(pashParams, horizontalSequenceHash);
  // free memory
  fprintf(stderr, "finished traversing horizontal sequence\n");
  printNow();
  fprintf(stderr, "freed resources\n");
  printNow();
  fflush(stderr);

  fclose(pashParams->outputFilePtr);
  return 0;
}
