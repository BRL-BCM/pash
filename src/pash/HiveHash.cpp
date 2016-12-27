#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define xDEBUG(flag, code) if (flag) {code; fflush(stdout); fflush(stderr);}

#include "HiveHash.h"
#define DEB_CONSTR 0
#define DEB_ADD_ENTRY 0
#define DEB_XALLOC 0
#define DEB_MARK_ENTRY 0 
#define DEB_GET_RUNNER 0
#define DEB_HASH_TRUALLOC 0
#define DEB_ADD_ENTRYXX 0
#define DEB_CHECKHASHXX 0 

/**
  Add a pair (offset, value) to hash entry key, extending the current bin accordingly.
  Check invariant that offsets are monotonically increasing.
  @param key current hash key
  @param offset current offset
  @param value current value
  @return 0 for success, 1 for failure
*/
int
HiveHash::markEntry(guint32 key) {
  unsigned long currentBin;
  guint32 lastValue;
  int binSize;
  xDEBUG(DEB_MARK_ENTRY, fprintf(stderr, "B markEntry %d %x\n", key, key));

  currentBin = (unsigned long) hashSkeleton[key];
  hashSkeleton[key]= (guint32*)(currentBin+1); 
  xDEBUG(DEB_MARK_ENTRY, fprintf(stderr, "E markEntry %d %d\n", key, (unsigned long)hashSkeleton[key]));
  
  return 0;
}

/**
  Add a pair (offset, value) to hash entry key, extending the current bin accordingly.
  Check invariant that offsets are monotonically increasing.
  @param key current hash key
  @param offset current offset
  @param value current value
  @return 0 for success, 1 for failure
*/
int
HiveHash::addEntryX(guint32 key, guint32 value, guint32 offset) {
  guint32* currentBin;
  guint32 lastValue;
  int binSize;
  xDEBUG(DEB_XALLOC, fprintf(stderr, "B addEntryX %d %x --> %d %d\n", key, key, value, offset));
  if (value==0) {
    fprintf(stderr, "adding value of 0\n");
    exit(0);
  }

  currentBin = hashSkeleton[key];
  if (currentBin==NULL) {
    fprintf(stderr, "bin should not be empty\n");
    exit(0);
  } else {
    if (currentBin[1]==0) {
      binSize = currentBin[0];
      xDEBUG(DEB_XALLOC, fprintf(stderr, "about to alloc final %d --> %d entries\n", key, 2*binSize+1));
      currentBin = (guint32*) realloc(currentBin, (2*binSize+1)*sizeof(guint32));
      xDEBUG(DEB_XALLOC, fprintf(stderr, "done alloc final %d --> %d entries\n", key, 2*binSize+1));
      currentBin[0]=1;
      hashSkeleton[key]=currentBin;
    } else {
      currentBin[0] +=1;  
    }
    binSize = currentBin[0];
    currentBin[2*binSize-1] = value;
    currentBin[2*binSize] = offset;
  }
  if (DEB_ADD_ENTRY) {
        int numberOfBinEntries;
        int pos;
        numberOfBinEntries = currentBin[0];
        fprintf(stderr, "key %d has data: %d (value,offset) pairs: [\n",
                key, numberOfBinEntries/2);
        for (pos=1; pos<numberOfBinEntries;pos+=2) {
          fprintf(stderr, "(%d, %d)\t", key, currentBin[pos], currentBin[pos+1]);
        }
        fprintf(stderr,"]\n");
  }

  xDEBUG(DEB_XALLOC, fprintf(stderr, "E addEntryX \n"));
  return 0;
}

/**
  Add a pair (offset, value) to hash entry key, extending the current bin accordingly.
  Check invariant that offsets are monotonically increasing.
  @param key current hash key
  @param offset current offset
  @param value current value
  @return 0 for success, 1 for failure
*/
int
HiveHash::addEntry(guint32 key, guint32 value, guint32 offset) {
  guint32* currentBin;
  guint32 lastValue;
  int binSize;
  xDEBUG(DEB_XALLOC, fprintf(stderr, "B addEntry \n"));

  currentBin = hashSkeleton[key];
  if (currentBin==NULL) {
    numberOfKeys ++;
    currentBin = (guint32*) malloc(3*sizeof(guint32));
    hashSkeleton[key]= currentBin;
    currentBin[0]=3;
    currentBin[1]=value;
    currentBin[2]=offset;
    memoryFootprint += 3.0;
  } else {
    // check that values are monotonically increasing
    binSize = currentBin[0];
    lastValue = currentBin[binSize-2];
    if (lastValue>value) {
      fprintf(stderr, "for key %d current value %d is smaller than last value in bin %d, at position %d\n",
              key, value, lastValue, binSize-2);
      return 1;
    } else {
      binSize = currentBin[0];
      // increase bin size to make room for current value
      xDEBUG(DEB_XALLOC, fprintf(stderr, "about to reallocate %d --> %d\n", binSize, binSize+2));
      currentBin = (guint32*) realloc(currentBin, (binSize+2)*sizeof(guint32));
      xDEBUG(DEB_XALLOC, fprintf(stderr, "done reallocating %d --> %d\n", binSize, binSize+2));
      if (currentBin == NULL) {
        fprintf(stderr, "could not reallocate HiveHash bin from %d to %d bytes\n",
                binSize*sizeof(guint32), (binSize+2)*sizeof(guint32));
        return 1;
      } else {
        hashSkeleton[key]=currentBin;
        currentBin[0] = binSize+2;
        currentBin[binSize] = value;
        currentBin[binSize+1] = offset;
        memoryFootprint+=2.0;
      }
    }
  }
  if (DEB_ADD_ENTRY) {
        int numberOfBinEntries;
        int pos;
        numberOfBinEntries = currentBin[0];
        fprintf(stderr, "key %d has data: %d (value,offset) pairs: [\n",
                key, numberOfBinEntries/2);
        for (pos=1; pos<numberOfBinEntries;pos+=2) {
          fprintf(stderr, "(%d, %d)\t", key, currentBin[pos], currentBin[pos+1]);
        }
        fprintf(stderr,"]\n");
  }

  numberOfHashValues++;
  xDEBUG(DEB_XALLOC, fprintf(stderr, "E addEntry \n"));
  return 0;
}

/** Checks the content of the hive hash; useful for debugging.
*   @param filePtr  output stream
*/
void
HiveHash::checkHash() {
  guint32* currentBin;
  int numberOfBinEntries, numberOfOffsets, pos;
  guint32 key;
  int offsetIndex, valueIndex;
  for (key=0; key<hashSize; key++) {
    if (hashSkeleton[key] != NULL) {
      currentBin = hashSkeleton[key];
      numberOfBinEntries = currentBin[0];
      for (pos=1; pos<numberOfBinEntries;pos+=2) {
        if (currentBin[pos]==0) {
          fprintf(stderr, "hash with seq values 0\n");
          exit(0);
        }
      }
    }
  }
}


void
HiveHash::checkHashXX() {
  guint32* currentBin;
  int numberOfBinEntries, numberOfOffsets, pos;
  guint32 key;
  int offsetIndex, valueIndex;
  for (key=0; key<hashSize; key++) {
    if (hashSkeleton[key] != NULL) {
      currentBin = hashSkeleton[key];
      numberOfBinEntries = currentBin[0];
      for (pos=2; pos<2+numberOfBinEntries*2;pos+=2) {
        xDEBUG(DEB_CHECKHASHXX, fprintf(stderr, "[%d][%d]: %d\n", key, pos, currentBin[pos]));
        if (currentBin[pos]==0) {
          fprintf(stderr, "hash with seq values 0\n");
          exit(0);
        }
      }
    }
  }
}

/** Dumps the content of the hive hash; useful for debugging.
*   @param filePtr  output stream
*/
void
HiveHash::dumpHash(FILE *filePtr) {
  guint32* currentBin;
  int numberOfBinEntries, numberOfOffsets, pos;
  guint32 key;
  int offsetIndex, valueIndex;
  fprintf(filePtr, "dumping hiveHash %x %p\n", this);
  for (key=0; key<hashSize; key++) {
    if (hashSkeleton[key] != NULL) {
      currentBin = hashSkeleton[key];
      numberOfBinEntries = currentBin[0];
      fprintf(filePtr, "key %d has data: %d (value,offset) pairs: [\n",
              key, numberOfBinEntries/2);
      for (pos=1; pos<numberOfBinEntries;pos+=2) {
        fprintf(filePtr, "(%d, %d)\t", currentBin[pos], currentBin[pos+1]);
      }
      fprintf(filePtr,"]\n");
    }
  }
  printStatistics(filePtr);
}
/** Sets up a (value,offset) running pointer
*   @param key kmer
*   @param listRunner  list runner data structure
*   @return 0 for success, 1 otherwise
*   */
int
HiveHash::getIntListRunner(guint32 key, IntListRunner* listRunner) {
  guint32 * currentBin = hashSkeleton[key];
  guint32 pos;
  if (currentBin != NULL) {
    listRunner->left = 2*currentBin[0];
    listRunner->list = currentBin+2;
    xDEBUG(DEB_GET_RUNNER, fprintf(stderr, "key %d %x runner %p %x size %d list %p\n", 
         key, key, listRunner, listRunner, currentBin[0], currentBin+2));
    xDEBUG(DEB_GET_RUNNER, 
      for (pos=2; pos<=2*currentBin[0];pos+=2) {
        fprintf(stderr, "(%d, %d)\t", currentBin[pos], currentBin[pos+1]);
      }
      ; fprintf(stderr, "\n");
    )
  } else {
    listRunner->list = NULL;
    listRunner->left = 0;
  }
  return 0;
}

/** Initializes a HiveHash object
*    @param size number of entries in the hive hash
*/
HiveHash::HiveHash(int size, guint32 keepKmerPercent) {
  int key;
  hashSize = size;
  if (size < 0) {
    fprintf(stderr, "HiveHash size should be greater than 0!\nExiting ...\n");
    exit(1);
  }
  hashSkeleton = (guint32**) malloc(hashSize*sizeof(guint32*));
  if (hashSkeleton == NULL) {
    fprintf(stderr, "could not allocate HiveHash skeleton\n");
    exit(1);
  }
  for (key=0; key<size; key++) {
    hashSkeleton[key]=(guint32*)0;
  }
  memoryFootprint = 0.0;
  numberOfHashValues = 0;
  numberOfKeys = 0;
  kmerPercent = (double)keepKmerPercent*1.0/100.0;
}

/** Destroys a HiveHash object.*/
HiveHash::~HiveHash() {
  int key;
  /*for (key=0; key<hashSize; key++) {
    if (hashSkeleton[key] != NULL) {
      free(hashSkeleton[key]);
    }
  }
  */
  free(hashSkeleton);
}

/** Returns the size in bytes of the hive hash.
  @return size in bytes of hive hash
*/
double
HiveHash::getMemoryFootprint() {
  return (memoryFootprint*sizeof(guint32)+hashSize*sizeof(guint32*))/(1024.0*1024.0);
}

/** Print hive hash statistics; useful for debugging.
*   @param filePtr  output stream
*/
void
HiveHash::printStatistics(FILE* filePtr) {
  fprintf(filePtr, "Hive hash %p %x; memory footprint %g elements, %gMB\n", this, this,
          memoryFootprint, getMemoryFootprint());
  fprintf(filePtr, "number of values: %d number of keys %d\n",
          numberOfHashValues, numberOfKeys);
}


int HiveHash::allocateHashMemory() {
  int key;
  xDEBUG(DEB_HASH_TRUALLOC, fprintf(stderr, "starting truAlloc\n"));
  
  poolList = NULL;
  individualPoolSize = 1*1024*1024;
  //individualPoolSize = 10;
  
  currentPool = (guint32*)malloc( individualPoolSize* sizeof(guint32));
  poolList = g_slist_append (poolList, (gpointer) currentPool);
  guint32 currentPoolUsed = 0;
  guint32 currentPoolLeft = individualPoolSize-currentPoolUsed;
  // get statistics about the has occupancy
  double nBins, sumKmers, sumX2Kmers;
  double maxKmers = 0;
  guint32 binSize;
  
  /*for (nBins=0, sumKmers=0, key=0; key<hashSize; key++) {
    guint32 binSize = (unsigned long)hashSkeleton[key];
    sumKmers+= (double)binSize;
    sumX2Kmers+= ((double)binSize)*((double)binSize);
    if (binSize>0) {
      if (maxKmers<binSize) {
         maxKmers = binSize;
      }
      nBins += 1.0;
    }
  }*/
  double mean ; //= sumKmers*1.0/(nBins*1.0);
  double stddev ;//= sqrt( sumX2Kmers/nBins - mean*mean);
  
  
  
  guint32 *kmerFreqHist = (guint32*) malloc(250000*sizeof(guint32));
  guint32 *kmerOccurencesHist = (guint32*) malloc(250000*sizeof(guint32));
  guint32 ii;
  double nIndividualKmers=0, nKmerOccurences=0;
  
  for (ii=0; ii<250000; ii++) {
    kmerFreqHist[ii] = 0;
    kmerOccurencesHist[ii]=0;
  }
  guint32 r;
  for (key=0; key<hashSize; key++) {
    binSize = (unsigned long)hashSkeleton[key];
    if (binSize>0) {
      if (maxKmers<binSize) {
        maxKmers = binSize;
      }

      r=binSize/2;
      if (r>=250000) {
        r=250000-1;
      }
      kmerFreqHist[r]+=1;
      kmerOccurencesHist[r]+= binSize;
      nIndividualKmers += 1;
      nKmerOccurences += binSize;
    }
  }
  
  double threshold;
  double sumIndividualKmers=0;
  double sumKmerOccurence=0; 
  double sumIndividKmersNorm=0;
  double sumKmerOccurenceNorm=0;
  fprintf(stderr, "Total kmer %g total kmer occurances %g\n", nIndividualKmers, nKmerOccurences);
  for(r=0; r<250000; r++) {
    if (kmerFreqHist[r]>0) {
      fprintf(stderr, "[%d] %d outof %g %d outof %g ====> ",
              r, kmerFreqHist[r], nIndividualKmers, kmerOccurencesHist[r], nKmerOccurences);
      sumIndividualKmers += kmerFreqHist[r];
      sumKmerOccurence += kmerOccurencesHist[r];
      fprintf(stderr, "%g %g ---> ", sumIndividualKmers, sumKmerOccurence);
      sumIndividKmersNorm = sumIndividualKmers/nIndividualKmers;
      sumKmerOccurenceNorm = sumKmerOccurence/nKmerOccurences;
      fprintf(stderr, "[%d]"
                       "\t%d\t%g\t%g\t||\t"
                       "\t%d\t%g\t%g\n",
              r,
              kmerFreqHist[r], sumIndividualKmers, sumIndividKmersNorm,
              kmerOccurencesHist[r], sumKmerOccurence, sumKmerOccurenceNorm);
      if (sumKmerOccurenceNorm>=kmerPercent) {
        threshold = r*2+1;
        break;
      }
    }
    //break when sumCummNorm>=0.9
  }
  fflush(stderr);
  //exit(0);
  
  //double threshold = mean; 
  //double threshold = mean+stddev; 
  //double threshold = mean+10*stddev;
  
  xDEBUG(1, fprintf(stderr, "%g nonempty bins sum %g mean %g max %g threshold %g\n", 
     nIndividualKmers, sumKmerOccurence, nIndividualKmers, maxKmers, threshold)); 

  for (key=0; key<hashSize; key++) {
    guint32 binSize = (unsigned long)hashSkeleton[key];
    xDEBUG(DEB_HASH_TRUALLOC,  if (binSize>0) { fprintf(stderr, "individual bin size %d\n", binSize);} );
    if (binSize==0 || binSize>threshold) {
      hashSkeleton[key]=NULL;
    } else {
      guint32 neededSize = 2+binSize*2;
      if (neededSize > individualPoolSize) {
        individualPoolSize = neededSize;
        xDEBUG(1, fprintf(stderr, "pool size increased to %d\n", neededSize));
      }
      if (neededSize>currentPoolLeft) {
        currentPool = (guint32*)malloc( individualPoolSize* sizeof(guint32));
        poolList = g_slist_append (poolList, (gpointer) currentPool);
        currentPoolUsed = 0;
        currentPoolLeft = individualPoolSize-currentPoolUsed;
        xDEBUG(DEB_HASH_TRUALLOC, fprintf(stderr, "reallocate mempool: used %d needed %d left %d\n",
                                          currentPoolUsed, neededSize, currentPoolLeft));
      }
      hashSkeleton[key] = currentPool + currentPoolUsed;
      guint32* currentBin = hashSkeleton[key];
      currentBin[0]=binSize;
      currentBin[1]=0;
      currentPoolUsed += neededSize;
      currentPoolLeft -= neededSize;
      xDEBUG(DEB_HASH_TRUALLOC, fprintf(stderr, "after key %d: used %d needed %d left %d\n",
                                          key, currentPoolUsed, neededSize, currentPoolLeft));
    }
  }
  return 0;
}


/**
  Add a pair (offset, value) to hash entry key, extending the current bin accordingly.
  Check invariant that offsets are monotonically increasing.
  @param key current hash key
  @param offset current offset
  @param value current value
  @return 0 for success, 1 for failure
*/
int
HiveHash::addEntryXX(guint32 key, guint32 value, guint32 offset) {
  guint32* currentBin;
  guint32 lastValue;
  guint32 binSize, binEntries;;
  xDEBUG(DEB_ADD_ENTRYXX, fprintf(stderr, "B addEntryXX %d %x --> %d %d\n", key, key, value, offset));
  if (value==0) {
    fprintf(stderr, "adding value of 0\n");
    exit(0);
  }

  currentBin = hashSkeleton[key];
  if (currentBin==NULL) {
    //  fprintf(stderr, "bin should not be empty\n");
    return 0;
  } else {
    binSize = currentBin[0];
    binEntries = currentBin[1];
    xDEBUG(DEB_ADD_ENTRYXX, fprintf(stderr, "key %d binSize %d binEntries %d\n", key, binSize, binEntries));
    if (binEntries==binSize) {
      fprintf(stderr, "trying to add more keys to key %d than capacity %d\n", key, binSize);
      exit(0);
    }
    currentBin[1] +=1;  
    binEntries = currentBin[1];
    currentBin[2*binEntries] = value;
    currentBin[2*binEntries+1] = offset;
  }
  
  if (DEB_ADD_ENTRYXX) {
        int numberOfBinEntries;
        int pos;
        numberOfBinEntries = currentBin[0];
        fprintf(stderr, "key %d has data: %d (value,offset) pairs: [\n",
                key, numberOfBinEntries);
        for (pos=2; pos<2+numberOfBinEntries*2;pos+=2) {
          fprintf(stderr, "(%d, %d)\t",  currentBin[pos], currentBin[pos+1]);
        }
        fprintf(stderr,"]\n");
  }

  xDEBUG(DEB_ADD_ENTRYXX, fprintf(stderr, "E addEntryXX \n"));
  return 0;
}
