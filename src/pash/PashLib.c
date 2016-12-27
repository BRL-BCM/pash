/*
Copyright (c) 2004-2008 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/

/***********************************************************************
* PashLib.c
* various support functions for Pash
***********************************************************************/

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <getopt.h>
#include "PashLib.h"
#include "PashDebug.h"
#include "HiveHash.h"
#include "FastaUtil.h"
#include "FastQUtil.h"
#include "FixedHashKey.h"
#include "IgnoreList.h"
#include "BisulfiteKmerGenerator.h"

#define DEB_HIVE_HASH 0
#define DEB_HASH_VERTICAL_SEQ 0
#define DEB_DUMP_HIVE_HASH 0
#define DEB_DUMP_HIVE_HASH_1 0
#define DEB_MASK_LEN 0
#define DEB_LOAD_VERT_HASH 0
#define DEB_SIZE_VERT_HASH 0
#define DEB_LOAD_PER_READ 0

/**Get the compliment of a base
  \note Used to get reverse compliment sequences, for hash key determination
  \pre upper case
@param base base
@return complementa of a base
*/

static inline char complement(char base)
{
  switch(base) {
    case 'a':
    case 'A':
      return 'T';
      break;
    case 'c':
    case 'C':
      return 'G';
      break;
    case 'g':
    case 'G':
      return 'C';
      break;
    case 'T':
    case 't':
      return 'A';
    case 'N':
    case 'n':
    default:
      return 'N';
      break;
  }
}

/** Interpret  command line args and initialize control variables
@param argc number of command line arguments +1
@param argv array of command line arguments
@param p pash control parameters
@return  pashParameters structure, if OK, NULL if error
*/
PashParameters* parseCommandLine(int argc, char**argv) {
	PashParameters* pp;
	char maskStr [MAX_MASK_LEN+1];
	int opt=0;
	int ii;
	int keyCount;
	// command line long options
  static struct option long_options[] = {
    {"text", no_argument, 0, 't'},
    {"scratch", required_argument, 0, 'S'},
    {"diagonals", required_argument, 0, 'd'},
    {"verticalInput", required_argument, 0, 'v'},
    {"horizontalInput", required_argument, 0, 'h'},
    {"verticalWordOffset", required_argument, 0,'G'},
    {"patternWeight", required_argument, 0, 'k'},
    {"patternLength", required_argument, 0, 'n'},
    {"samplingPattern",required_argument,0,'m'},
    {"outputFile",required_argument,0,'o'},
    {"ignoreList",required_argument,0,'L'},
    {"maxMappings",required_argument,0,'N'},
    {"topPercent",required_argument,0,'P'},
    {"score", required_argument, 0,'s'},
    {"indexMemory", required_argument, 0, 'M'},
    {"bisulfiteSequencingMapping", no_argument, 0, 'B'},
    {"highSensitivity", no_argument, 0, '0'},
    {"mediumSensitivity", no_argument, 0, '1'},
    {"lowSensitivity", no_argument, 0, '2'},
    {"fastMode", no_argument, 0, '3'},
    {"keepHashedKmersPercent",required_argument,0,'K'},
    {"self", no_argument, 0, 'A'},
    {0, 0, 0, 0}
  };
	if (argc==1) {
		PashUsage();
		exit(0);
	}
	if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-help") ) {
		PashUsage();
		exit(0);
	}
	pp = (PashParameters*) malloc(sizeof(PashParameters));
	xDieIfNULL(pp, fprintf(stderr, "could not allocate memory for Pash parameters\n"), 1);
	// set up sensible defaults
	strcpy(pp->verticalFile, "");
	strcpy(pp->horizontalFile, "");
	pp->verticalFastqUtil = NULL;
	pp->fastaUtilHorizontal = NULL;
	strcpy(pp->scratchDirectory, "");
	strcpy(pp->outputFile, "");
	pp->numberOfDiagonals = DEFAULT_NUMBER_OF_DIAGONALS;
	pp->minScore = DEFAULT_MIN_SCORE;
	pp->wordOffset = DEFAULT_WORD_OFFSET;
	pp->isMaskDefined = FALSE;
	pp->hiveHashMemoryLimit = DEFAULT_HIVE_HASH_MEMORY;
	pp->useIgnoreList = FALSE;
	pp->useGzippedOutput = FALSE;
	int option_index = 0;
	pp->mask.maskLen = 0;
	pp->mask.keyLen = 0;
	pp->topPercent = 0.01;
        pp->useIgnoreList = 0;
	pp->maxMappings = 1;
        pp->bisulfiteSequencingMapping=0;
	pp->sensitivityMode = MediumSensitivity;
        pp->keepHashedKmersPercent=99;
	while((opt=getopt_long(argc,argv,
                         ":S:M:d:v:h:L:g:G:k:n:m:o:s:tBA:N:P:0123K:",
                         long_options, &option_index))!=-1) {
		switch(opt) {
		case 'S':  // scratch directory location
			strncpy(pp->scratchDirectory,optarg,MAX_FILE_NAME_SIZE);
			break;
		case 'M':  // maximum amount of RAM available for hive hash
			pp->hiveHashMemoryLimit=atoi(optarg);
			if(pp->hiveHashMemoryLimit < 0) {
				xDie(fprintf(stderr, "Vertical sequence hash limit should be positive\n"), 1);
			}
			break;
		case 'd':  // shortcut for setting 'paralellogram' geometry of this diagonal length
			pp->numberOfDiagonals=atoi(optarg);
			if(pp->hiveHashMemoryLimit < 0) {
				xDie(fprintf(stderr, "Number of diagonals should be positive\n"), 1);
			}
			break;
		case 'v':  // name of vertical input file
			strncpy(pp->verticalFile, optarg, MAX_FILE_NAME_SIZE);
			break;
    case 'L':  // name of vertical input file
			strncpy(pp->ignoreListFile, optarg, MAX_FILE_NAME_SIZE);
      pp->useIgnoreList = 1;
			break;
		case 'h':  // name of horizontal input file
			strncpy(pp->horizontalFile, optarg, MAX_FILE_NAME_SIZE);
		break;
		case 'G':  // offset between starts of adjacent words to be processed
			pp->wordOffset=atoi(optarg);
			pp->sensitivityMode = UserDefinedSensitivity;
			break;
      case '0':
      fprintf(stderr, "High sensitivity mode\n");
      pp->sensitivityMode=HighSensitivity;
      break;
		case '1':
      fprintf(stderr, "Medium sensitivity mode\n");
      pp->sensitivityMode=MediumSensitivity;
      break;
		case '2':
      fprintf(stderr, "Low sensitivity mode\n");
      pp->sensitivityMode=LowSensitivity;
      break;
		case '3':
                fprintf(stderr, "Fast mode\n");
                pp->sensitivityMode=FastSensitivity;
                pp->keepHashedKmersPercent=93;
                break;
		case 'k':  // sampling pattern 'weight'
			pp->mask.keyLen=atoi(optarg);
			break;
		case 'n':  // sampling pattern length
			pp->mask.maskLen=atoi(optarg);
			break;  // sampling pattern
		case 'm':
			strncpy(maskStr,optarg, MAX_MASK_LEN);
			pp->mask.maskLen = strlen(maskStr);
			for(ii=0; ii<pp->mask.maskLen; ii++) {
				pp->mask.mask[ii]=0;
			}
			for(ii=0; ii<pp->mask.maskLen; ii++) {
				if(maskStr[ii]=='1') {
					pp->mask.mask[ii]=1;
					keyCount++;
				}
			}
			pp->mask.keyLen = keyCount;
			pp->isMaskDefined = TRUE;
			computeOverlapContribution(&pp->mask);
			xDEBUG(DEB_MASK_LEN, fprintf(stderr, "sampling pattern  %d %d \n", pp->mask.keyLen, pp->mask.maskLen));
			break;
		case 'o':  // output file name
			strncpy(pp->outputFile, optarg, MAX_FILE_NAME_SIZE);
			break;
		case 's':  // score cutoff (base score)
			pp->minScore=atoi(optarg);
			break;
		case 'N':  // offset between starts of adjacent words to be processed
			pp->maxMappings=atoi(optarg);
			break;
    case 'P':  // offset between starts of adjacent words to be processed
			pp->topPercent=atof(optarg)/100.0;
			break;
		case 'z': // GZIP-ed output (default is text)
			pp->useGzippedOutput=TRUE;
			break;
                case 'K': // GZIP-ed output (default is text)
			pp->keepHashedKmersPercent=atoi(optarg);
                        if (pp->keepHashedKmersPercent<90) {
			  pp->keepHashedKmersPercent=90;
                        }
                        if (pp->keepHashedKmersPercent>100) {
			  pp->keepHashedKmersPercent=100;
                        }
			break;
		case 'A':  // self-comparison; ignore any matches below the main diagonal
			pp->selfComparison=TRUE;
			break;
    case 'B':
      fprintf(stderr, "Performing bisulfite sequencing mapping\n");
      pp->bisulfiteSequencingMapping=1;
      break;
		case ':':
			xDie(fprintf(stderr,"Warning: missing argument for -%c\n",optopt),1);
			break;
		case '?':
			xDie(fprintf(stderr,"Warning: -%c is not a defined option\n",optopt),1);
			break;
		}
	}

	// validate parameters
	if (!strcmp(pp->verticalFile, "")) {
		xDie(fprintf(stderr, "no vertical sequence file provided\n"), 1);
	}
	if(!strcmp(pp->horizontalFile, "")) {
		xDie(fprintf(stderr, "no horizontal sequence file provided\n"), 1);
	}
	if (!strcmp(pp->scratchDirectory, "")) {
		xDie(fprintf(stderr, "no scratch directory provided\n"), 1);
	}
	if (!strcmp(pp->outputFile, "")) {
		xDie(fprintf(stderr, "no output file provided\n"), 1);
	}
	// generate sampling pattern
  if(pp->isMaskDefined==FALSE) {
		if (pp->mask.maskLen == 0 || pp->mask.keyLen == 0) {
			xDie(fprintf(stderr, "no sampling pattern provided\n"), 1);
		} else  if(pp->mask.keyLen > pp->mask.maskLen) {
			xDie(fprintf(stderr, "sampling pattern weight is greater than sampling pattern length\n"), 1);
		} else {
			setMask(&pp->mask);
			pp->isMaskDefined= TRUE;
		}
  }
	return pp;
}

/** print Pash usage information.*/
void PashUsage() {
	fprintf(stdout,"Pash version " PASH_BUILD_VERSION "\nUsage:\nPash.exe\n"
  " --verticalInput         | -v <verticalFile>  Vertical sequence as a fasta input file with full path; if file ends\n"
  "                              in '.fof', it is assumed the named file contains a list of Fasta files\n"
  " --horizontalInput       | -h <horizontalFile> Horizontal sequence: as a fasta input file with full path; if file ends\n"
  "                              in '.fof', it is assumed the named file contains a list of Fasta files\n"
  " --diagonals             | -d <number of diagonals> \n"
  " --patternWeight         | -k <pattern weight> Number of sampled positions in the sampling pattern\n"
  " --patternLength         | -n total length of sampling pattern, including unsampled positions\n"
" --samplingPattern       | -m sampling pattern (e.g. 11011 would sample the two positions, skip\n"
"                              one position, then sample the next two\n"
" --verticalWordOffset    | -G <vertical word offset gap - must be a multiple of diagonal offset gap>\n"
" --outputFile            | -o <output file name>\n"
" --score                 | -s <scoreCutoff>\n"
" --gzip                  | -z  request gzip-ed output (default is text)\n"
" --scratch               | -S Scratch directory location \n"
" --indexMemory           | -M index of the vertical sequence hash in MB(default 1024)\n"
" --ignoreList            | -L ignore the kmers present in the ignore list file\n"
" --maxMappings           | -N maximum number of mappings per read\n"
" --topPercent            | -P top percent from the best alignment score to be reported for each read; use numbers in the interval 0-100; default 1\n"
" --bisulfiteSeq          | -B perform mapping of bisulfite sequencing reads\n"
" --highSensitivity       | -0 run pash in high-sensitivity mode \n"
" --mediumSensitivity     | -1 run pash in medium-sensitivity mode (default setting)\n"
" --lowSensitivity        | -2 run pash in low-sensitivity mode \n"
" --fastMode              | -3 run pash in fast mode \n"
" --keepHashedKmersPercent| -K percent amount of hashed kmers to keep; this value should be between 90 and 100; default is 99 \n"
" --help print usage info and exit\n"
  "\n");
}

/// Determine the memory requirements of the current vertical batch.
int sizeCurrentVerticalSequencesBatch(PashParameters* pp) {
	guint32 hiveHashSize;
	//FastaUtil* fastaUtilVertical = pp->fastaUtilVertical;
  PashFastqUtil* verticalFastqUtil = pp->verticalFastqUtil;
	int i;
	guint32 currentHashSize;
  guint32 currentForwardChunkStart, currentForwardChunkStop;
  guint32 currentReverseChunkStop, currentReverseChunkStart;
  guint32 sequenceLength;
  char currentKmer [MAX_MASK_LEN+1];
  char currentReverseKmer[MAX_MASK_LEN+1];
  int maskLen, maskWeight, kmerPos, maskPos;
  int startOffset, maxOffset, minOffset, offsetGap;
  Mask mask = pp->mask;
  int currentSequencePos;
  int numberOfDiagonals;
  int minPosition;
  int sequenceToKeep;
  char fwdSeq [2000], revCSeq[2000];
  guint32 crtGap;
  BisulfiteKmerGenerator *bisulfiteKmerGenerator;
  int bisulfiteSequencingMapping = pp->bisulfiteSequencingMapping;
  int kmersPerRead = 0; 

  int maxSeeds=512;
	guint32* guintSeeds;
	int actualSeeds;
  if (bisulfiteSequencingMapping) {
    bisulfiteKmerGenerator= new BisulfiteKmerGenerator();
    guintSeeds = (guint32*) malloc(maxSeeds*sizeof(guint32));
  }
  
  guint32 forwardKey, reverseKey;
  guint32 memoryTarget;
  maskLen =pp->mask.maskLen;
  maskWeight = pp->mask.keyLen;
  currentKmer[maskWeight] = '\0';
  currentReverseKmer[maskWeight] = '\0';
  offsetGap = pp->wordOffset;
  numberOfDiagonals = pp->numberOfDiagonals;
  memoryTarget = (pp->hiveHashMemoryLimit*95)/100;
  HiveHash *hiveHash = (HiveHash*)pp->hiveHash;
  double totalKmers;
	xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "START sizeCurrentVerticalSequencesBatch\n"));
  IgnoreList ignoreList; // see IgnoreList.h
  FILE *ignoreListFile;  // file handle for the ignore list
  // read ignore list, if applicable.  depends on sampling pattern (i.e. must call setMask before this point)
  int useIgnoreList = pp->useIgnoreList;
  if(useIgnoreList) {
    ignoreListFile=fopen(pp->ignoreListFile,"r");
    if(ignoreListFile==NULL) {
      fprintf(stdout,"failed to open ignore list file %s, bailing out\n", pp->ignoreListFile);
    }
    readIgnoreList(ignoreListFile, &ignoreList, pp->mask.keyLen); // ignoreList is allocated, initialized and read
    fclose(ignoreListFile);
    ignoreListFile=NULL;
    pp->ignoreList=ignoreList;
  }


  if (hiveHash == NULL) {
		hiveHashSize = 1;
		for ( i=0; i< pp->mask.keyLen; i++) {
			hiveHashSize *= 4;
		}
		xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "hive hash size = %d\n", hiveHashSize));
		hiveHash = new HiveHash(hiveHashSize, pp->keepHashedKmersPercent);
		pp->hiveHash = hiveHash;
	}
 
 
  // if at limit of memory, then stop, because we have enough info to resume the vert hash filling
  guint32 numberOfVerticalSequences = verticalFastqUtil->getNumberOfSequences();
  guint32 lastVerticalSequenceMapped = pp->lastVerticalSequenceMapped;
  guint32 currentVerticalSequence;
  kmersPerRead = 0;
  guint32 maxReadLength = 1;
  for (currentVerticalSequence=1;
       currentVerticalSequence<=numberOfVerticalSequences;
       currentVerticalSequence++) {
     kmersPerRead = 0;
     // print current sequence length
    const char* currentSequence = verticalFastqUtil->retrieveSequence(currentVerticalSequence);
    const char* currentDefName = verticalFastqUtil->retrieveDefName(currentVerticalSequence);
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
           fprintf(stderr, ">> got current sequence [%d] %s -- %s: %d characters\n",
                  currentVerticalSequence,
                  currentDefName,
                  currentSequence,
                  strlen(currentSequence)
                  ));
    sequenceLength=strlen(currentSequence); 
    if (maxReadLength < sequenceLength) {
			maxReadLength = sequenceLength ;
		} 
    switch(pp->sensitivityMode) {
			case HighSensitivity:
				offsetGap = (sequenceLength<=50?(sequenceLength<=36?2:3):(sequenceLength<76?4:6));
				break;		
			case MediumSensitivity:
				offsetGap = (sequenceLength<76?(sequenceLength<=36?2:(sequenceLength<=50?3:4)):(sequenceLength<100?8:12));
				break;
			case LowSensitivity:
				offsetGap = (sequenceLength<76?(sequenceLength<=36?2:(sequenceLength<=50?4:4)):(sequenceLength<100?10:(sequenceLength<150?14:18)));
				break;
                        case FastSensitivity:
				offsetGap = (sequenceLength<=100?(sequenceLength<=36?2:(sequenceLength<=50?3:4)):(sequenceLength<=200?8:10));
				break;

			case UserDefinedSensitivity:
				offsetGap=pp->wordOffset;
				break;
			default:
				fprintf(stderr, "offset gap not determined !\n");
				exit(0);
		} 
    int fwdIdx;
    const char* fwdSeq = currentSequence;
    sequenceLength=strlen(fwdSeq);
    if (!bisulfiteSequencingMapping) {
      for (fwdIdx =0; fwdIdx<sequenceLength; fwdIdx++) {
        revCSeq[sequenceLength-1-fwdIdx] = complement(fwdSeq[fwdIdx]);
      }  
      revCSeq[sequenceLength] = '\0';
    }
    
    xDEBUG(DEB_SIZE_VERT_HASH,
           fprintf(stderr, "S got current sequence [%d] %s >>%s<< >>%s<< ~~%s~~ of length %d: %d characters\n",
                  currentVerticalSequence,
                  currentDefName,
                  fwdSeq, revCSeq,
                  currentSequence,
                  sequenceLength, strlen(fwdSeq)));

    //xDEBUG(DEB_DUMP_HIVE_HASH, hiveHash->dumpHash(stderr));
    //offsetGap = (sequenceLength<=50?(sequenceLength<=36?2:4):(sequenceLength<=76?6:12));
    // start filling the vertical hash table
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
           fprintf(stderr, "S processing short sequence (maybe read)\n"));
    maxOffset = sequenceLength-maskLen;
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
           fprintf(stderr, "S %d Gap: %d \n", maxOffset, offsetGap));
    for (startOffset = 0; startOffset<=maxOffset; startOffset+= offsetGap) {
      for (currentSequencePos=startOffset,maskPos = 0,kmerPos = 0;
           kmerPos < maskWeight;
           currentSequencePos++, maskPos++) {
        if (mask.mask[maskPos]) {
          currentKmer [kmerPos] = currentSequence[currentSequencePos];
          kmerPos++;
        }
      }
      getKeyForSeq(currentKmer, &forwardKey);
      xDEBUG(DEB_HASH_VERTICAL_SEQ,
             fprintf(stderr, "VERT SEQ HASH: found forward kmer %s %d %x, v seq %d at v offset %d\n",
                      currentKmer, forwardKey, forwardKey, currentVerticalSequence,
                      startOffset));
      if(!useIgnoreList || !isIgnored(forwardKey, ignoreList)) {
        if (bisulfiteSequencingMapping) {
          actualSeeds=bisulfiteKmerGenerator->generateKmerList(currentKmer, guintSeeds, maxSeeds);
          for (i=0; i<actualSeeds; i++) {
                xDEBUG(DEB_SIZE_VERT_HASH, fprintf(stderr, "mark entry %d of %d\n", i+1, actualSeeds));
          	hiveHash->markEntry(guintSeeds[i]);
                kmersPerRead += 1;
          }  
        } else {
          xDEBUG(DEB_SIZE_VERT_HASH, fprintf(stderr, "mark entry %d \n", forwardKey ));
          hiveHash->markEntry(forwardKey);
          kmersPerRead += 1;
        }
        xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "done adding to hive hash\n")); 
      } else {
        // fprintf(stderr, "ignore %s %u\n", currentKmer, forwardKey);
      }
      xDEBUG(DEB_DUMP_HIVE_HASH, hiveHash->dumpHash(stderr));
    }
    
    if (!bisulfiteSequencingMapping) {
      minOffset = maskLen -1;
      for (startOffset = sequenceLength-1; startOffset>=minOffset; startOffset-= offsetGap) {
        for (currentSequencePos=startOffset,maskPos = 0,kmerPos = 0;
             kmerPos < maskWeight;
             currentSequencePos--, maskPos++) {
          if (mask.mask[maskPos]) {
            currentReverseKmer[kmerPos] = complement(currentSequence[currentSequencePos]);
            kmerPos++;
          }
        }
        getKeyForSeq(currentReverseKmer, &reverseKey);
        xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "found reverse kmer %s %d %x, adding value %d at offset %d\n",
                                              currentReverseKmer, reverseKey, reverseKey,
                                              currentVerticalSequence, sequenceLength-1-startOffset));
        if(!useIgnoreList || !isIgnored(reverseKey, ignoreList)) {
          hiveHash->markEntry(reverseKey);
          kmersPerRead += 1;
        } else {
          // fprintf(stderr, "ignore %s %u\n", currentReverseKmer, reverseKey);
        }
      }
    }
    totalKmers += kmersPerRead;
    xDEBUG(DEB_LOAD_PER_READ, fprintf(stderr, "rl: %s\t%d\n",
       currentDefName, kmersPerRead)); 
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
      fprintf(stderr, "++ sequence status %d\n",
              currentVerticalSequence));
	}
	xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "end of parsing or fill hash capacity\n"));
  xDEBUG(DEB_DUMP_HIVE_HASH_1, hiveHash->dumpHash(stderr));
	xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "STOP hashCurrentVerticalSequencesBatch\n"));
  
  if (pp->bisulfiteSequencingMapping) {
    //delete bisulfiteKmerGenerator;
    free(guintSeeds);
  }
  xDEBUG(DEB_LOAD_PER_READ, fprintf(stderr, "Total kmers: %g\n", totalKmers));
	hiveHash->allocateHashMemory();
  pp->numberOfDiagonals = maxReadLength ;
  if (pp->numberOfDiagonals<100) {
	pp->numberOfDiagonals = 100;
  }
  fprintf(stderr, "Number of diagonals set to %d\n", pp->numberOfDiagonals);

  return 0;
}

/// Hash vertical sequence until the hive hash size reaches a user-specified limit.
int hashCurrentVerticalSequencesBatch(PashParameters* pp) {
	guint32 hiveHashSize;
	//FastaUtil* fastaUtilVertical = pp->fastaUtilVertical;
  PashFastqUtil* verticalFastqUtil = pp->verticalFastqUtil;
	int i;
	guint32 currentHashSize;
  guint32 currentForwardChunkStart, currentForwardChunkStop;
  guint32 currentReverseChunkStop, currentReverseChunkStart;
  guint32 sequenceLength;
  char currentKmer [MAX_MASK_LEN+1];
  char currentReverseKmer[MAX_MASK_LEN+1];
  int maskLen, maskWeight, kmerPos, maskPos;
  int startOffset, maxOffset, minOffset, offsetGap;
  Mask mask = pp->mask;
  int currentSequencePos;
  int numberOfDiagonals;
  int minPosition;
  int sequenceToKeep;
  char fwdSeq [2000], revCSeq[2000];
  guint32 crtGap;
  BisulfiteKmerGenerator *bisulfiteKmerGenerator;
  int bisulfiteSequencingMapping = pp->bisulfiteSequencingMapping;
  
  int maxSeeds=512;
	guint32* guintSeeds;
	int actualSeeds;
  if (bisulfiteSequencingMapping) {
    bisulfiteKmerGenerator= new BisulfiteKmerGenerator();
    guintSeeds = (guint32*) malloc(maxSeeds*sizeof(guint32));
  }
  
  guint32 forwardKey, reverseKey;
  guint32 memoryTarget;
  maskLen =pp->mask.maskLen;
  maskWeight = pp->mask.keyLen;
  currentKmer[maskWeight] = '\0';
  currentReverseKmer[maskWeight] = '\0';
  offsetGap = pp->wordOffset;
  numberOfDiagonals = pp->numberOfDiagonals;
  memoryTarget = (pp->hiveHashMemoryLimit*95)/100;
  HiveHash *hiveHash = (HiveHash*)pp->hiveHash;
	xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "START hashCurrentVerticalSequencesBatch\n"));
  IgnoreList ignoreList; // see IgnoreList.h
  FILE *ignoreListFile;  // file handle for the ignore list
  // read ignore list, if applicable.  depends on sampling pattern (i.e. must call setMask before this point)
  int useIgnoreList = pp->useIgnoreList;
  if(useIgnoreList) {
    ignoreListFile=fopen(pp->ignoreListFile,"r");
    if(ignoreListFile==NULL) {
      fprintf(stdout,"failed to open ignore list file %s, bailing out\n", pp->ignoreListFile);
    }
    readIgnoreList(ignoreListFile, &ignoreList, pp->mask.keyLen); // ignoreList is allocated, initialized and read
    fclose(ignoreListFile);
    ignoreListFile=NULL;
    pp->ignoreList=ignoreList;
  }


  if (hiveHash == NULL) {
		hiveHashSize = 1;
		for ( i=0; i< pp->mask.keyLen; i++) {
			hiveHashSize *= 4;
		}
		xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "hive hash size = %d\n", hiveHashSize));
		hiveHash = new HiveHash(hiveHashSize, pp->keepHashedKmersPercent);
		pp->hiveHash = hiveHash;
	}
 
 
  // if at limit of memory, then stop, because we have enough info to resume the vert hash filling
  guint32 numberOfVerticalSequences = verticalFastqUtil->getNumberOfSequences();
  guint32 lastVerticalSequenceMapped = pp->lastVerticalSequenceMapped;
  guint32 currentVerticalSequence;
  for (currentVerticalSequence=1;
       currentVerticalSequence<=numberOfVerticalSequences;
       currentVerticalSequence++) {
     // print current sequence length
    const char* currentSequence = verticalFastqUtil->retrieveSequence(currentVerticalSequence);
    const char* currentDefName = verticalFastqUtil->retrieveDefName(currentVerticalSequence);
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
           fprintf(stderr, ">> got current sequence [%d] %s -- %s: %d characters\n",
                  currentVerticalSequence,
                  currentDefName,
                  currentSequence,
                  strlen(currentSequence)
                  ));

    SequenceInfo* currentSequenceInfo = &pp->verticalSequencesInfos[currentVerticalSequence];
    sequenceLength = strlen(currentSequence);
    currentSequenceInfo->sequenceName = currentDefName;
    currentSequenceInfo->sequenceLength = sequenceLength;
    currentSequenceInfo->bestAnchoringScore = 0;
    currentSequenceInfo->bestSWScore = 0;
		currentSequenceInfo->bestSkeletonScore = 0;
    currentSequenceInfo->bestScoreMappings = 0;
    currentSequenceInfo->passingMappings = 0;
   
   
    // TODO: in the first pass stage, find a suitable # of diagonals
    if (sequenceLength>numberOfDiagonals) {
      fprintf(stderr, "number of diagonals %d should be greater than sequence length %d\n",
              numberOfDiagonals, sequenceLength);
      exit(1);
    }
   
    int fwdIdx;
    const char* fwdSeq = currentSequence;
    if (!bisulfiteSequencingMapping) {
      for (fwdIdx =0; fwdIdx<sequenceLength; fwdIdx++) {
        revCSeq[sequenceLength-1-fwdIdx] = complement(fwdSeq[fwdIdx]);
      }  
      revCSeq[sequenceLength] = '\0';
    }
    
    xDEBUG(DEB_LOAD_VERT_HASH,
           fprintf(stderr, "got current sequence [%d] %s >>%s<< >>%s<< ~~%s~~ of length %d: %d characters\n",
                  currentVerticalSequence,
                  currentDefName,
                  fwdSeq, revCSeq,
                  currentSequence,
                  sequenceLength, strlen(fwdSeq)));
    switch(pp->sensitivityMode) {
			case HighSensitivity:
				offsetGap = (sequenceLength<=50?(sequenceLength<=36?2:3):(sequenceLength<76?4:6));
				break;		
			case MediumSensitivity:
				offsetGap = (sequenceLength<76?(sequenceLength<=36?2:(sequenceLength<=50?3:4)):(sequenceLength<100?8:12));
				break;
			case LowSensitivity:
				offsetGap = (sequenceLength<76?(sequenceLength<=36?2:(sequenceLength<=50?4:4)):(sequenceLength<100?10:(sequenceLength<150?14:18)));
				break;
			case FastSensitivity:
				offsetGap = (sequenceLength<=100?(sequenceLength<=36?2:(sequenceLength<=50?3:4)):(sequenceLength<=200?8:10));
				break;
			case UserDefinedSensitivity:
				offsetGap=pp->wordOffset;
				break;
			default:
				fprintf(stderr, "offset gap not determined !\n");
				exit(0);
		}
    //xDEBUG(DEB_DUMP_HIVE_HASH, hiveHash->dumpHash(stderr));
    //offsetGap = (sequenceLength<=50?(sequenceLength<=36?2:4):(sequenceLength<=76?6:12));
    // start filling the vertical hash table
    xDEBUG(DEB_HASH_VERTICAL_SEQ,
           fprintf(stderr, "processing short sequence (maybe read)\n"));
    maxOffset = sequenceLength-maskLen;
    for (startOffset = 0; startOffset<=maxOffset; startOffset+= offsetGap) {
      for (currentSequencePos=startOffset,maskPos = 0,kmerPos = 0;
           kmerPos < maskWeight;
           currentSequencePos++, maskPos++) {
        if (mask.mask[maskPos]) {
          currentKmer [kmerPos] = currentSequence[currentSequencePos];
          kmerPos++;
        }
      }
      getKeyForSeq(currentKmer, &forwardKey);
      xDEBUG(DEB_HASH_VERTICAL_SEQ,
             fprintf(stderr, "VERT SEQ HASH: found forward kmer %s %d %x, v seq %d at v offset %d\n",
                      currentKmer, forwardKey, forwardKey, currentVerticalSequence,
                      startOffset));
      if(!useIgnoreList || !isIgnored(forwardKey, ignoreList)) {
        if (bisulfiteSequencingMapping) {
          actualSeeds=bisulfiteKmerGenerator->generateKmerList(currentKmer, guintSeeds, maxSeeds);
          for (i=0; i<actualSeeds; i++) {
                xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "add entry %d of %d\n", i+1, actualSeeds));
          	hiveHash->addEntryXX(guintSeeds[i], currentVerticalSequence, startOffset);
          }  
        } else {
          hiveHash->addEntryXX(forwardKey, 2*currentVerticalSequence, startOffset);
        }
        xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "done adding to hive hash\n")); 
      } else {
        // fprintf(stderr, "ignore %s %u\n", currentKmer, forwardKey);
      }
      xDEBUG(DEB_DUMP_HIVE_HASH, hiveHash->dumpHash(stderr));
    }
    
    if (!bisulfiteSequencingMapping) {
      minOffset = maskLen -1;
      for (startOffset = sequenceLength-1; startOffset>=minOffset; startOffset-= offsetGap) {
        for (currentSequencePos=startOffset,maskPos = 0,kmerPos = 0;
             kmerPos < maskWeight;
             currentSequencePos--, maskPos++) {
          if (mask.mask[maskPos]) {
            currentReverseKmer[kmerPos] = complement(currentSequence[currentSequencePos]);
            kmerPos++;
          }
        }
        getKeyForSeq(currentReverseKmer, &reverseKey);
        xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "found reverse kmer %s %d %x, adding value %d at offset %d\n",
                                              currentReverseKmer, reverseKey, reverseKey,
                                              currentVerticalSequence, sequenceLength-1-startOffset));
        if(!useIgnoreList || !isIgnored(reverseKey, ignoreList)) {
          hiveHash->addEntryXX(reverseKey, 2*currentVerticalSequence+1, sequenceLength-1-startOffset);
        } else {
          // fprintf(stderr, "ignore %s %u\n", currentReverseKmer, reverseKey);
        }
      }
    }

    xDEBUG(DEB_HASH_VERTICAL_SEQ,
      fprintf(stderr, "++ sequence status %d\n",
              currentVerticalSequence));
	}
  pp->lastVerticalSequenceMapped = currentVerticalSequence;
	xDEBUG(DEB_HASH_VERTICAL_SEQ, fprintf(stderr, "end of parsing or fill hash capacity\n"));
  xDEBUG(DEB_DUMP_HIVE_HASH_1, hiveHash->dumpHash(stderr));
	xDEBUG(DEB_HIVE_HASH, fprintf(stderr, "STOP hashCurrentVerticalSequencesBatch\n"));
  
  if (pp->bisulfiteSequencingMapping) {
    //delete bisulfiteKmerGenerator;
    free(guintSeeds);
  }

  // set number of diagonals based on maximum read length
  hiveHash->checkHashXX();
    return 0;
}


/// Initialize the sequence hash
SequenceHash *initSequenceHash() {
  SequenceHash* sequenceHash = (SequenceHash*) malloc(sizeof(SequenceHash));
  sequenceHash->lastSequenceId = 0;
  sequenceHash->numberOfChunksInCurrentSequence = 0;
  sequenceHash->offsetOfSequenceBufferInRealSequence=0;
  sequenceHash->currentSequenceChunk = 0;
	sequenceHash->hiveHash = NULL;
}

void printNow() {
  time_t now;
  time(&now);
  fprintf(stderr, "%.24s\n", ctime(&now));
}
