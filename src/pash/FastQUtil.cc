#include <stdlib.h>
#include <string.h>

#include "FastQUtil.h"
#include "BRLGenericUtils.h"
#include "generic_debug.h"

#define DEB_LOAD_SEQUENCES 0
#define DEB_INIT  0
#define DEB_TRAV 0

PashFastqUtil::PashFastqUtil(char* fileName, ReadsSequenceType sequenceType) {
  strcpy(fastqFile, fileName);
  numberOfSequences=0;
  availableEntries = 1;
  sequenceNames = (char**) malloc(sizeof(char*)*availableEntries);
  forwardSequences= (char**) malloc(sizeof(char*)*availableEntries);
  reverseSequences = (char**) malloc(sizeof(char*)*availableEntries);
  qualities = (char**) malloc(sizeof(char*)*availableEntries);
  if (qualities==NULL || sequenceNames==NULL || forwardSequences==NULL || reverseSequences==NULL) {
    fprintf(stderr, "Insufficient memory, exiting ...\n");
    exit(2);
  }
  forwardSequences[0]=NULL;
  reverseSequences[0]=NULL;
  qualities[0]=NULL;
  sequenceNames[0]=NULL;
  sequencePool = new SequencePool();
  readsSequenceType = sequenceType;
  xDEBUG(DEB_INIT, fprintf(stderr, "init finished\n"));
}

PashFastqUtil::~PashFastqUtil() {
  delete sequencePool;
  free(forwardSequences);
  free(reverseSequences);
  free(qualities);
  free(sequenceNames);
}

int PashFastqUtil::loadSequences(int loadReverseComplement) {
  FILE* fastqPtr = BRLGenericUtils::openTextGzipBzipFile(fastqFile);
  if (fastqPtr==NULL) {
    fprintf(stderr, "could not open fastq file %s\n", fastqFile);
    exit(2);
  }
  reverseSequencesAvailable = loadReverseComplement;
  char currentLine[MAX_LINE_LENGTH], currentDef[MAX_LINE_LENGTH];
  char *buffer = (char*) malloc(sizeof(char)*(DEFAULT_BUFFER_SIZE+5*MAX_LINE_LENGTH));;
  guint32 bufferReadPos = 0; 
  guint32 numCharsRead = 0;
  guint32 newLinePos = 0;
  guint32 targetChars;
  int endOfFile;
  guint32 bufferIndex;
  guint32 newLineIdx1, newLineIdx2, newLineIdx3, newLineIdx4;
  guint32 actualBufferSize;
  guint32 currentSequenceId = 0;
  guint32 sequenceBufferIndexStart;
  
  for (endOfFile =0; !endOfFile; ) {
    targetChars = DEFAULT_BUFFER_SIZE;
    numCharsRead = fread(buffer+bufferReadPos, sizeof(char), targetChars, fastqPtr);
    actualBufferSize = bufferReadPos+numCharsRead;
    xDEBUG(DEB_LOAD_SEQUENCES,
           fprintf(stderr, "bufferReadPos=%d targetChars=%d actualBufferSize=%d\n",
                   bufferReadPos, targetChars, actualBufferSize));
    if (numCharsRead<targetChars) {
      endOfFile = 1;
      if (buffer[actualBufferSize-1]!='\n') {
        buffer[actualBufferSize]='\n';
        actualBufferSize+=1;
      }
    }
    sequenceBufferIndexStart = 0;
    int done;
    for(done=0; !done; ) {
      // search exactly 4 new lines
      xDEBUG(DEB_LOAD_SEQUENCES,
           fprintf(stderr, "sequenceBufferIndexStart=%d \n",
                   sequenceBufferIndexStart));
      for (bufferIndex  = sequenceBufferIndexStart;
           bufferIndex<actualBufferSize && buffer[bufferIndex]!='\n';
           bufferIndex++);
      if (bufferIndex==actualBufferSize) {
        done = 1;
        continue;
      } else {
        newLineIdx1 = bufferIndex;
      }
      
      for (bufferIndex++;
           bufferIndex<actualBufferSize && buffer[bufferIndex]!='\n';
           bufferIndex++);
      if (bufferIndex==actualBufferSize) {
        done = 1;
        continue;
      } else {
        newLineIdx2 = bufferIndex;
      }
      
      for (bufferIndex++;
           bufferIndex<actualBufferSize && buffer[bufferIndex]!='\n';
           bufferIndex++);
      if (bufferIndex==actualBufferSize) {
        done = 1;
        continue;
      } else {
        newLineIdx3 = bufferIndex;
      }
      
      for (bufferIndex++;
           bufferIndex<actualBufferSize && buffer[bufferIndex]!='\n';
           bufferIndex++);
      if (bufferIndex==actualBufferSize) {
        xDEBUG(DEB_LOAD_SEQUENCES, fprintf(stderr, "oops: %d %s ", sequenceBufferIndexStart, &buffer[sequenceBufferIndexStart])); 
        done = 1;
        continue;
      } else {
        newLineIdx4 = bufferIndex;
      }
      currentSequenceId += 1;
      buffer[newLineIdx1]='\0';
      buffer[newLineIdx2]='\0';
      buffer[newLineIdx3]='\0';
      buffer[newLineIdx4]='\0';
      xDEBUG(DEB_LOAD_SEQUENCES, fprintf(stderr, "[%d] 4 lines:\n%s\n%s\n%s\n%s\n", currentSequenceId,
                                         &buffer[sequenceBufferIndexStart],
                                         &buffer[newLineIdx1+1],
                                         &buffer[newLineIdx2+1],
                                         &buffer[newLineIdx3+1]));
      // add sequence and quality to appropriate repository
      if (strlen(&buffer[newLineIdx1+1]) <= MAX_READ_SIZE) {
        numberOfSequences++;
        if (availableEntries == numberOfSequences) {
          guint32 oldAvailableEntries = availableEntries;
          availableEntries = 5*availableEntries/4+1;
          sequenceNames = (char**) realloc(sequenceNames, sizeof(char*)*availableEntries);
          forwardSequences= (char**) realloc(forwardSequences, sizeof(char*)*availableEntries);
          if (loadReverseComplement) {
            reverseSequences = (char**) realloc(reverseSequences, sizeof(char*)*availableEntries);
          }
          if (readsSequenceType==FastaAndQualityScores) {
            qualities = (char**) realloc(qualities, sizeof(char*)*availableEntries);
          }
          if (qualities==NULL || sequenceNames==NULL || forwardSequences==NULL || reverseSequences==NULL) {
            fprintf(stderr, "Insufficient memory, exiting ...\n");
            exit(2);
          }       
        }
        char fixedName[MAX_LINE_LENGTH]; 
        sscanf(&buffer[sequenceBufferIndexStart+1], " %s ", fixedName );
        sequenceNames[numberOfSequences] = sequencePool->addSequence(fixedName);
        // convert sequence to uppercase
        guint32 fwdSequenceIndex = 0;
        char *fwdSequence = &buffer[newLineIdx1+1];
        guint32 sequenceLength =   strlen(fwdSequence);
        
        for (fwdSequenceIndex=0; fwdSequenceIndex<sequenceLength; fwdSequenceIndex++) {
          if (fwdSequence[fwdSequenceIndex]>='a') {
            fwdSequence[fwdSequenceIndex] = fwdSequence[fwdSequenceIndex] - ('a'-'A');
          }
        }
        forwardSequences[numberOfSequences] = sequencePool->addSequence(fwdSequence);
        xDEBUG(DEB_LOAD_SEQUENCES, fprintf(stderr,"added fwd seq %s\n", fwdSequence));
  
        if (loadReverseComplement) {
          char reverseComplementSequence [MAX_LINE_LENGTH];
          int revComplementIndex;
          for (revComplementIndex=0; revComplementIndex<sequenceLength; revComplementIndex++) {
            reverseComplementSequence[revComplementIndex]= baseComplement(fwdSequence[sequenceLength-1-revComplementIndex]);
          }
          reverseComplementSequence[sequenceLength]='\0';
          reverseSequences[numberOfSequences] = sequencePool->addSequence(reverseComplementSequence);
          xDEBUG(DEB_LOAD_SEQUENCES, fprintf(stderr,"added rev seq %s\n", reverseComplementSequence));
        }
        
        if (readsSequenceType==FastaAndQualityScores) {
          qualities[numberOfSequences] = sequencePool->addSequence(&buffer[newLineIdx3+1]);
        }
      
      }
      
      sequenceBufferIndexStart=newLineIdx4+1;
    }
    xDEBUG(DEB_LOAD_SEQUENCES,
           fprintf(stderr, ">> bufferIndex=%d, sequenceBufferIndexStart=%d \n",
                   bufferIndex, sequenceBufferIndexStart));
    if (bufferIndex == actualBufferSize && sequenceBufferIndexStart<actualBufferSize) {
      if (endOfFile) {
        // should not get here
        fprintf(stderr, "incorrect input; could not find 4x number of sequences lines\n");
      } else {
        guint32 copyIdx;
        bufferReadPos = actualBufferSize-sequenceBufferIndexStart;
        for (copyIdx=0; copyIdx<bufferReadPos; copyIdx++) {
          buffer[copyIdx] = buffer[sequenceBufferIndexStart+copyIdx];
        }
      }
      // roll over leftover buffer
    } else {
      bufferReadPos=0;
    }
  }
  
  free(buffer);
  fclose(fastqPtr);
  return 0;
}

const char* PashFastqUtil::retrieveSequence(guint32 sequenceId) {
  if (sequenceId<=numberOfSequences) {
    return forwardSequences[sequenceId];
  } else {
    return NULL;
  }
}

const char* PashFastqUtil::retrieveRevComplementSequence(guint32 sequenceId) {
  if (reverseSequencesAvailable && sequenceId<=numberOfSequences) {
    return reverseSequences[sequenceId];
  } else {
    return NULL;
  }
}

const char* PashFastqUtil::retrieveQualityScores(guint32 sequenceId) {
  if (readsSequenceType==FastaAndQualityScores && sequenceId<=numberOfSequences) {
    return qualities[sequenceId];
  } else {
    return NULL;
  }
}


char PashFastqUtil::baseComplement(char b) {
  switch(b) {
    case 'a':
    case 'A':
      return 'T';
      break;
    case 't':
    case 'T':
      return 'A';
      break;
    case 'c':
    case 'C':
      return 'G';
      break;
    case 'g':
    case 'G':
      return 'C';
      break;
    default:
      return 'N';
      break;
  }
}


const char* PashFastqUtil::retrieveDefName(guint32 sequenceId) {
  if (sequenceId<=numberOfSequences) {
    return sequenceNames[sequenceId];
  } else {
    return NULL;
  }
}

guint32 PashFastqUtil::getNumberOfSequences() {
  return numberOfSequences;
}

