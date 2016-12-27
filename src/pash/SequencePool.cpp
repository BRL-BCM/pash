#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SequencePool.h"
#include "generic_debug.h"

#define DEB_FREE 0

/*
struct ASeqPool {
  char *sequence;
  guint32 used;
  guint32 capacity;
};

class SequencePool {
  char** poolSkeleton;
  guint32 poolSize;
  guint32 poolCapacity;
public:
  StringPool();
  ~StringPool();
  char* addSequences(char* fwd, char* revComplement);
};

*/

SequencePool::SequencePool() {
  poolCapacity = 10;
  poolSize = 1;
  poolSkeleton = (ASeqPool*) malloc(sizeof(ASeqPool)*poolCapacity);
  if (poolSkeleton==NULL) {
    fprintf(stderr, "could not allocate pool skeleton");
    exit(2);
  }
  latestPool = 0;
  int i;
  poolSkeleton[0].used = 0;
  poolSkeleton[0].capacity = POOL_INCREMENT;
  poolSkeleton[0].sequence = (char*) malloc(sizeof(char)*POOL_INCREMENT);
  if (poolSkeleton[0].sequence==NULL) {
    fprintf(stderr, "could not allocate seq pool");
    exit(2);
  }
}


SequencePool::~SequencePool() {
  int i;
  xDEBUG(DEB_FREE, fprintf(stderr, "freeing %d pools\n", poolCapacity));
  for (i=0; i<poolCapacity; i++) {
    xDEBUG(DEB_FREE, fprintf(stderr, "about to free pool %d %x\n", i, poolSkeleton[i].sequence));
    free(poolSkeleton[i].sequence);
  }
  free(poolSkeleton);
}


char* SequencePool::addSequence(char* fwd) {
  int needSize = strlen(fwd)+1;
  if ((poolSkeleton[latestPool].used + needSize)>poolSkeleton[latestPool].capacity) {
    latestPool ++;
    if (latestPool>=poolSize) {
      poolCapacity = poolCapacity+1;
      poolSkeleton = (ASeqPool*) realloc(poolSkeleton, sizeof(ASeqPool)*poolCapacity);
      if (poolSkeleton==NULL) {
        fprintf(stderr, "could not reallocate pool skeleton\n");
        exit(2);
      }
    }
    poolSkeleton[latestPool].used = 0;
    poolSkeleton[latestPool].capacity = POOL_INCREMENT;
    poolSkeleton[latestPool].sequence = (char*) malloc(sizeof(char)*POOL_INCREMENT);
    if (poolSkeleton[latestPool].sequence==NULL) {
      fprintf(stderr, "could not allocate seq pool");
      exit(2);
    }
  }
  guint32 used = poolSkeleton[latestPool].used;
  guint32 readLen = strlen(fwd);
  char* returnPtr = &poolSkeleton[latestPool].sequence[used];
  strcpy(&poolSkeleton[latestPool].sequence[used], fwd);
  poolSkeleton[latestPool].used =used+readLen+1;
  return returnPtr;
}
