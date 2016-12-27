#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>

#include "Pash2SAMConverter.h"
#include "BRLGenericUtils.h"
#include "generic_debug.h"

#define DEB_PARSE_INTLIST 0
#define DEB_HASH_SEQS     0
#define DEB_UNMAPPED      0
#define DEB_BACKCONVERT   0
#define DEB_FIX_BS_READSTARTBUG 0

Pash2SAMConverter::Pash2SAMConverter() {
	pashFastQUtil = NULL;
	sequenceHashTable = NULL;
	referenceSizesHashTable = NULL;
}

Pash2SAMConverter::~Pash2SAMConverter() {
}

int Pash2SAMConverter::parseParameters(int argc, char*argv[]) {
	static struct option long_options[] = {
		{"pashMappings", required_argument, 0, 'p'},
		{"fastqFile", required_argument, 0, 'f'},
		{"referenceSequnces", required_argument, 0, 'r'},
		{"SAMFile", required_argument, 0, 'S'},
		{"sample", required_argument, 0, 's'},
		{"center", required_argument, 0, 'C'},
		{"help", no_argument, 0, 'h'},
		{"bisulfiteSeq", no_argument, 0, 'B'},
		{0, 0, 0, 0}
	};
 // setup params
		if (argc==1) {
	usage();
	exit(0);
	}
	if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-help") ||
      !strcmp(argv[1],"-h") || !strcmp(argv[1],"--h")) {
		usage();
		exit(0);
	}

  int option_index = 0;
  char opt;

  strcpy(fastqFile, "");
  strcpy(mappingFile, "");
  strcpy(referenceSequencesFile, "");
  strcpy(outputFile, "");
  strcpy(centerName, "");
  strcpy(sampleName, "");
  bisulfiteSeqFlag = 0;
	
  while((opt=getopt_long(argc,argv,
                         "f:p:r:S:h:s:C:B",
                         long_options, &option_index))!=-1) {
    switch(opt) {
    case 'p':
      strcpy(mappingFile, optarg);
      break;
    case 'S':
      strcpy(outputFile, optarg);
      break;
    case 'f':
      strcpy(fastqFile, optarg);
      break;
    case 'r':
      strcpy(referenceSequencesFile, optarg);
      break;
    case 's':
      strcpy(sampleName, optarg);
      break;
    case 'C':
      strcpy(centerName, optarg);
      break;
    case 'h':
      usage();
      exit(2);
		case 'B':
			bisulfiteSeqFlag=1;
			break;
	  default:
      fprintf(stderr, "unknown option %s \n", optarg);
      usage();
      exit(2);
    }
  }

  // validate the user inputs
	if (!strcmp(centerName,"")) {
		fprintf(stderr, "Center name was not specified\n");
		exit(2);
	}
	if (!strcmp(sampleName,"")) {
		fprintf(stderr, "The sample name was not specified\n");
		exit(2);
	}
  return 0;
} // Pash2SAMConverter::parseParameters

int Pash2SAMConverter::loadReferenceSizes() {
	FILE* refseqPtr = BRLGenericUtils::openTextGzipBzipFile(referenceSequencesFile);
	char chromosome[MAX_LINE_LENGTH];
	guint32 chromSize;
	char tmpLine[MAX_LINE_LENGTH];
	while(fgets(tmpLine, MAX_LINE_LENGTH-1, refseqPtr) !=NULL) {
		sscanf(tmpLine, "%s %d", chromosome, &chromSize);
		fprintf(outPtr, "@SQ	SN:%s	LN:%d\n", chromosome, chromSize);
	}
	fclose(refseqPtr);
}


int Pash2SAMConverter::convertPash2SAM() {
  outPtr = fopen(outputFile, "w");
  // print header
	fprintf(outPtr,
					"@HD	VN:1.0\n"
					"@PG	ID:pash\n"
					"@RG  ID:-- SM:%s CN:%s\n",
					sampleName, centerName)		;			
						//@RG	ID:--	LB:--	PU:--	CN:HGSC	DT:--PL:--

// load references
	loadReferenceSizes();
// load fastq
  pashFastQUtil = new PashFastqUtil(fastqFile, FastaAndQualityScores);
// load reads
  pashFastQUtil->loadSequences(0);
// setup sequence hash
	prepareReadsHash();
//
	convertMappings();
  fclose(outPtr);
  return 0;
}

/** Performs the actual conversion to SAM, mapping by mapping
 * @return 0 for success, 1 otherwise
 */
int Pash2SAMConverter::convertMappings() {
	int matchGain = 1;
  int gapPenalty = -3;
  int mismatchPenalty=-2;
  int numMatches, numMismatches, numGaps, numGapBases, numBlocks;
  char chromosomeName[MAX_LINE_LENGTH], readName[MAX_LINE_LENGTH];
  char blockSizes[MAX_LINE_LENGTH], startReadBlocks[MAX_LINE_LENGTH], startChromosomeBlocks[MAX_LINE_LENGTH];
  int chromosomeSize, chromosomeStart, chromosomeStop;
  int readSize, readStart, readStop;
  char strand;
  int matchScore;
	char samLine[4*MAX_LINE_LENGTH];
	int numBPVariations;
	char bpVariations[MAX_LINE_LENGTH];
	int maxAlignmentBlocks = 1000;
	guint32 *blockSizesArray       = (guint32*) malloc(maxAlignmentBlocks*sizeof(guint32));
	guint32 *horizontalStartsArray = (guint32*) malloc(maxAlignmentBlocks*sizeof(guint32));
	guint32 *verticalStartsArray   = (guint32*) malloc(maxAlignmentBlocks*sizeof(guint32));
	guint32 alignmentScore;
	char tmpLine[MAX_LINE_LENGTH];
	const char* fwdRead, *qualityScores;
  char revComplementRead[MAX_LINE_LENGTH+1], revQualityScores[MAX_LINE_LENGTH];
	char cigar[MAX_LINE_LENGTH];
	guint32 readId;
	guint32 numCG, numCHG, numCHH, numT;
  char cgBases[MAX_LINE_LENGTH];
  char chgBases[MAX_LINE_LENGTH];
  char chhBases[MAX_LINE_LENGTH];
  char tBases[MAX_LINE_LENGTH];
  guint32 methylatedBases[MAX_READ_SIZE];
  guint32 nonMethylatedBases[MAX_READ_SIZE];
  guint32 numMethylatedBases;
	
	FILE* mapPtr = BRLGenericUtils::openTextGzipBzipFile(mappingFile);
	if (mapPtr==NULL) {
		fprintf(stderr, "cannot open mapping file %s\n", mappingFile);
		perror("");
		exit(2);
	}
	// determine the number of matches per read
	while(fgets(tmpLine, MAX_LINE_LENGTH-1, mapPtr) !=NULL) {
		sscanf(tmpLine, "%s %d %d %s %d %d %c %d %d %d %d %d %s %s %s %d %s %d %s %d %s %d %s %d %s",
                     chromosomeName, &chromosomeStart, &chromosomeStop,
                     readName, &readStart, &readStop,
                     &strand,
                     &numMatches, &numMismatches, &numGaps, &numGapBases, &numBlocks,
                     blockSizes, startChromosomeBlocks, startReadBlocks,
                     &numBPVariations, bpVariations,
                     &numCG, cgBases, &numCHG, chgBases, &numCHH, chhBases, &numT, tBases);
              xDEBUG(DEB_PARSE_INTLIST,
                fprintf(stderr, "read %d matches, %d mismatches,"
                  "%d Gaps, %d GapBases, strand %c, read %s: %d-%d, chr %s %d-%d, numBlocks: %d\n",
                  numMatches, numMismatches, numGaps, numGapBases, strand,
                  readName, readStart, readStop,
                  chromosomeName, chromosomeStart, chromosomeStop, numBlocks));
		SequenceInfo* readInfo = (SequenceInfo*) g_hash_table_lookup(sequenceHashTable, readName);
		if (readInfo ==0) {
			fprintf(stderr, "read %s not found in fastq file\n", readName);
			
			continue;
		}
		readInfo->numMappings++;
	}
	
	fclose(mapPtr);
	mapPtr = BRLGenericUtils::openTextGzipBzipFile(mappingFile);
	
	// perform SAM conversion of mapped reads
	while(fgets(tmpLine, MAX_LINE_LENGTH-1, mapPtr) !=NULL) {
		sscanf(tmpLine, "%s %d %d %s %d %d %c %d %d %d %d %d %s %s %s %d %s %d %s %d %s %d %s %d %s",
                     chromosomeName, &chromosomeStart, &chromosomeStop,
                     readName, &readStart, &readStop,
                     &strand,
                     &numMatches, &numMismatches, &numGaps, &numGapBases, &numBlocks,
                     blockSizes, startChromosomeBlocks, startReadBlocks,
                     &numBPVariations, bpVariations,
                     &numCG, cgBases, &numCHG, chgBases, &numCHH, chhBases, &numT, tBases);
              xDEBUG(DEB_PARSE_INTLIST,
                fprintf(stderr, "read %d matches, %d mismatches,"
                  "%d Gaps, %d GapBases, strand %c, read %s: %d-%d, chr %s %d-%d, numBlocks: %d\n",
                  numMatches, numMismatches, numGaps, numGapBases, strand,
                  readName, readStart, readStop,
                  chromosomeName, chromosomeStart, chromosomeStop, numBlocks));
		xDEBUG(DEB_PARSE_INTLIST,
					 fprintf(stderr, "numBlocks: %d blockSizes %s readStarts %s genomeStarts %s\n",
							numBlocks, blockSizes, startReadBlocks, startChromosomeBlocks));
		alignmentScore = matchGain*numMatches + mismatchPenalty*numMismatches + gapPenalty*numGapBases;
		BRLGenericUtils::parseListOfIntegers(blockSizesArray, blockSizes, numBlocks);
		BRLGenericUtils::parseListOfIntegers(horizontalStartsArray, startChromosomeBlocks, numBlocks);
		BRLGenericUtils::parseListOfIntegers(verticalStartsArray, startReadBlocks, numBlocks);
		
		SequenceInfo* readInfo = (SequenceInfo*) g_hash_table_lookup(sequenceHashTable, readName);
		if (readInfo ==0) {
			fprintf(stderr, "read %s not found in fastq file\n", readName);
			continue;
		}
		readId = readInfo->sequenceId;
		fwdRead = pashFastQUtil->retrieveSequence(readId);
		qualityScores = pashFastQUtil->retrieveQualityScores(readId);
		samLine [0]='\0';
		guint32 Hi0=0;
		if (numGapBases==0 && numMismatches==0) {
			Hi0 = readInfo->numMappings;
		}
    if (strand!='+') {
      getReverseComplement(fwdRead, revComplementRead, MAX_LINE_LENGTH);
      getReverseString(qualityScores, revQualityScores, MAX_LINE_LENGTH);
    }
    if (strand=='-' && bisulfiteSeqFlag) {
      guint32 tmpStart = strlen(fwdRead)-readStop+1;
      guint32 tmpStop = strlen(fwdRead)-readStart+1;
      readStart = tmpStart;
      readStop = tmpStop;
      xDEBUG(DEB_FIX_BS_READSTARTBUG, fprintf(stderr, "new readStart %d new readStop %d\n", readStart, readStop));
    }
    buildCigarString(cigar, numBlocks,
				readStart, readStop,
				blockSizesArray, horizontalStartsArray, verticalStartsArray, strlen(fwdRead));
		
		if (bisulfiteSeqFlag) {
			char bisulfiteBackConvertedRead[MAX_READ_SIZE+1];
			if (strand=='+') {
				strcpy(bisulfiteBackConvertedRead, fwdRead);
			} else {
				strcpy(bisulfiteBackConvertedRead, revComplementRead);
			}
			int flag = 0;
			if (numT==0) {
				flag = 0;
			} else {
			  flag=bisulfiteBackConvertRead(bisulfiteBackConvertedRead, strand, numT, tBases,
						numBlocks, blockSizesArray, horizontalStartsArray, verticalStartsArray, fwdRead);
			}
			if (!flag){
				sprintf(samLine+strlen(samLine), "%s\t%d\t%s\t%d\t99\t%s\t*\t0\t0\t%s\t%s\tAS:i:%d\tNM:i:%d\tH0:i:%d\tBR:Z:%s\n",
					readName, strand=='+'?0:16, chromosomeName, chromosomeStart, cigar,
					bisulfiteBackConvertedRead,
					strand=='+'?qualityScores:revQualityScores,
					alignmentScore, numMismatches+numGapBases, Hi0, fwdRead);
			}
		} else {
			sprintf(samLine+strlen(samLine), "%s\t%d\t%s\t%d\t99\t%s\t*\t0\t0\t%s\t%s\tAS:i:%d\tNM:i:%d\tH0:i:%d\tNH:i:%d",
				readName, strand=='+'?0:16, chromosomeName, chromosomeStart, cigar,
				strand=='+'?fwdRead:revComplementRead,
				strand=='+'?qualityScores:revQualityScores,
				alignmentScore, numMismatches+numGapBases, Hi0, readInfo->numMappings);
		}
		fprintf(outPtr, "%s\n", samLine);
	}
	
	// SAM output of unmapped reads
	outputUnmappedReads();
	fclose(mapPtr);
}


void Pash2SAMConverter::usage() {
  fprintf(stderr, "pash2SAM.exe\n"
"    Utility that converts Pash mappings into the SAM format. \n"
"  --pashMappings      | -p    pash mappings file\n"
"  --fastqFile         | -f    FASTQ file for the mapped reads. When using .fa and .qual file, use our utility faqualToFastq.rb\n"
"  --referenceSequnces | -r    file containing the reference sequences and their lengths\n"
"  --bisulfiteSeqFlag  | -B    converting the mappings of a bisulfite sequencing run\n"
"  --SAMFile           | -S    output SAM file\n"
"  --sample            | -s    sample name\n"
"  --center            | -C    center name\n"
"\nExample usage\n"
"pash2SAM.exe -f h1.h3k36me3.fastq.bz2 -p h1.h3k36me3.pash.gz -r hg18.chromosomes -S h1.h3k36me3.sam\n");
}

/** Build a hash from read name to read id
 * @return 0 for succes, 1 otherwise
 */
int Pash2SAMConverter::prepareReadsHash() {
	sequenceHashTable = g_hash_table_new(g_str_hash, g_str_equal);
  xDEBUG(DEB_HASH_SEQS, fprintf(stderr, "build hash table at address %p\n", sequenceHashTable)); 
	guint32 numSequences = pashFastQUtil->getNumberOfSequences();
	guint32 i;
	sequenceInfos = (SequenceInfo*) malloc((numSequences+1)*sizeof(SequenceInfo));
	for (i=1; i<=numSequences; i++) {
		sequenceInfos[i].sequenceId = i;
		sequenceInfos[i].numMappings  = 0;
		const char* seqName = pashFastQUtil->retrieveDefName(i);
		g_hash_table_insert(sequenceHashTable, g_strdup(seqName), (gpointer) &sequenceInfos[i]);
	}
	return 0;
}

/** Generate extended CIGAR information for a read mapping.*/
int Pash2SAMConverter::buildCigarString(char* cigar, guint32 numBlocks,
																				guint32 readStart, guint32 readStop,
																				guint32* blockSizesArray,
																				guint32* horizontalStartsArray,
																				guint32* verticalStartsArray,
																				guint32 readLength) {
	guint32 blockIdx;
	cigar[0]='\0';
	guint32 crtVBlockStop, crtHBlockStop, nextVBlockStart, nextHBlockStart, Vdist, Hdist;
	if (readStart>1) {
		sprintf(cigar+strlen(cigar), "%dS", readStart-1);
	}
	for (blockIdx=1; blockIdx<numBlocks; blockIdx++) {
		sprintf(cigar+strlen(cigar), "%dM", blockSizesArray[blockIdx-1]);
		// call an insert or a deletion from the read
		crtVBlockStop = verticalStartsArray[blockIdx-1]+blockSizesArray[blockIdx-1]-1;
		nextVBlockStart = verticalStartsArray[blockIdx];
		crtHBlockStop = horizontalStartsArray[blockIdx-1]+blockSizesArray[blockIdx-1]-1;
		nextHBlockStart = horizontalStartsArray[blockIdx];
		Vdist = nextVBlockStart-crtVBlockStop-1;
		Hdist = nextHBlockStart-crtHBlockStop-1;
    if (Vdist<0 || Hdist<0) {
      fprintf(stderr, "erroneous line\n");
      exit(2);
    }
		if (Vdist>Hdist) {
			sprintf(cigar+strlen(cigar), "%dI", Vdist);
		} else {
			sprintf(cigar+strlen(cigar), "%dD", Hdist);
		}
	}
	sprintf(cigar+strlen(cigar), "%dM", blockSizesArray[numBlocks-1]);
	if (readStop<readLength) {
		sprintf(cigar+strlen(cigar), "%dS", readLength-readStop);
	}
}


/** Output unmapped reads.*/
int Pash2SAMConverter::outputUnmappedReads() {
	guint32 readId, numSequences;
	char samLine[MAX_LINE_LENGTH];
	const char* fwdRead, *qualityScores;
	const char* readName;
	numSequences = pashFastQUtil->getNumberOfSequences();
	for (readId=1; readId<numSequences; readId++) {
		const char* tmpReadName= pashFastQUtil->retrieveDefName(readId);
		xDEBUG(DEB_UNMAPPED, fprintf(stderr, "evaluating [%d] %s\n", readId, tmpReadName));
		if (sequenceInfos[readId].numMappings==0) {
			xDEBUG(DEB_UNMAPPED, fprintf(stderr, "unmapped [%d] %s\n", readId,  tmpReadName));
			fwdRead = pashFastQUtil->retrieveSequence(readId);
			qualityScores = pashFastQUtil->retrieveQualityScores(readId);
			readName = pashFastQUtil->retrieveDefName(readId);
			samLine [0]='\0';
			
			sprintf(samLine+strlen(samLine), "%s\t0\t*\t0\t99\t*\t*\t0\t0\t%s\t%s",
							readName, fwdRead, qualityScores);
			fprintf(outPtr, "%s\n", samLine);		
		}
	}
}


/** Obtain the reverse complement for a read.
@param fwdRead string containing the forward sequence
@param revComplementRead  location for the reverse complement sequence
@param maxReadLength maximum read length
@return 0 if success, 1 otherwise
*/
int Pash2SAMConverter::getReverseComplement(const char* fwdRead, char* revComplementRead, int maxReadLength) {
  int fwdIdx, revComplementIdx;
  int seqLength;
  seqLength = strlen(fwdRead);
  if (seqLength>maxReadLength) {
    return 1;
  }
  for (fwdIdx=0, revComplementIdx=seqLength-1; fwdIdx<seqLength; fwdIdx++, revComplementIdx--) {
    switch (fwdRead[fwdIdx]) {
      case 'a':
      case 'A':
        revComplementRead[revComplementIdx] = 'T';
        break;
      case 't':
      case 'T':
        revComplementRead[revComplementIdx] = 'A';
        break;
      case 'c':
      case 'C':
        revComplementRead[revComplementIdx] = 'G';
        break;
      case 'g':
      case 'G':
        revComplementRead[revComplementIdx] = 'C';
        break;
      default:
        revComplementRead[revComplementIdx] = 'N';
    }
  }
  revComplementRead[seqLength] = '\0';
  return 0;
}

/** Obtain the reverse of a string

/** Obtain the reverse complement for a read.
@param source string to be reversed
@param destination location for the reverse sequence
@param maxReadLength maximum read length
@return 0 if success, 1 otherwise
*/
int Pash2SAMConverter::getReverseString(const char* source, char* destination, int maxReadLength) {
  int fwdIdx, reverseIdx;
  int seqLength;
  seqLength = strlen(source);
  if (seqLength>maxReadLength) {
    return 1;
  }
  for (fwdIdx=0, reverseIdx=seqLength-1; fwdIdx<seqLength; fwdIdx++, reverseIdx--) {
    destination[reverseIdx] = source[fwdIdx];
  }
  destination[seqLength] = '\0';
  return 0;
}

/** Backconvert a bisulfite sequencing read, for the purpose of using it w/ generic tools such as samtools or picard to call genotypes*/
int Pash2SAMConverter::bisulfiteBackConvertRead(char* backConvertedRead, char strand,
															guint32 numT, char* tBases,
															guint32 numBlocks,
											 				guint32* blockSizesArray,
															guint32* horizontalStartsArray,
															guint32* verticalStartsArray,
															const char* fwdRead) {
	guint32 nonMethylatedBases[MAX_READ_SIZE];
															 
	int flag=0;
	if (numT>0) {
	 flag=BRLGenericUtils::parseListOfIntegers(&nonMethylatedBases[0], tBases, numT);
	}
	if (flag) {
		fprintf(stderr, "wrong line %s for read %s\n", tBases , fwdRead);
		return 1;
	}
	int idx1, idx2;
	guint32 tmp;
	// sort converted bases
	if (numT>1){ 
		// sort methylated bases
		for (idx1=0; idx1<numT-1; idx1++) {
			for (idx2=idx1+1; idx2<numT; idx2++) {
				if (nonMethylatedBases[idx1]>nonMethylatedBases[idx2]) {
					tmp = nonMethylatedBases[idx1];
					nonMethylatedBases[idx1]=nonMethylatedBases[idx2];
					nonMethylatedBases[idx2] = tmp;
				}
			}
		}
	}
	guint32 nonMethylatedBasesIndex = 0;
	
	guint32 readSize = strlen(backConvertedRead);
	guint32 blockIdx;
	guint32 crtVBlockStop, crtHBlockStop;
  guint32 crtVBlockStart, crtHBlockStart;
	int inBlockIdx;
	int inReadIdx;
	
	for (blockIdx=0; blockIdx<numBlocks && nonMethylatedBasesIndex<numT; blockIdx++) {
		crtVBlockStop  = verticalStartsArray[blockIdx]   + blockSizesArray[blockIdx]-1;
		crtHBlockStop  = horizontalStartsArray[blockIdx] + blockSizesArray[blockIdx]-1;
		crtHBlockStart = horizontalStartsArray[blockIdx];
		crtVBlockStart = verticalStartsArray[blockIdx];
	  for (inBlockIdx=crtHBlockStart; inBlockIdx<=crtHBlockStop && nonMethylatedBasesIndex<numT ; inBlockIdx++) {
			if(inBlockIdx==nonMethylatedBases[nonMethylatedBasesIndex]) {
				nonMethylatedBasesIndex++;
				if (strand=='+') {
					inReadIdx = inBlockIdx-crtHBlockStart+crtVBlockStart-1;
					backConvertedRead[inReadIdx]='C'; // back convert T to C
					xDEBUG(DEB_BACKCONVERT, fprintf(stderr, "back convert base %d --> %d %c %c\n", inBlockIdx, inReadIdx+1, fwdRead[inReadIdx], backConvertedRead[inReadIdx]));
					if (fwdRead[inReadIdx]!='T') {
						fprintf(stderr, "incorrect back conversion base %d --> %d %c %c", inBlockIdx, inReadIdx+1, fwdRead[inReadIdx], backConvertedRead[inReadIdx]);
						return 1;
					}
				} else {
					inReadIdx=inBlockIdx-crtHBlockStart+crtVBlockStart-1;
					backConvertedRead[inReadIdx]='G';
					xDEBUG(DEB_BACKCONVERT, fprintf(stderr, "back convert RC base %d --> %d %c %c\n", inBlockIdx, inReadIdx+1, fwdRead[readSize-1-inReadIdx], backConvertedRead[inReadIdx] ));
				  if (fwdRead[readSize-1-inReadIdx]!='T') {
						fprintf(stderr, "incorrect back conversion base %d --> %d %c %c", inBlockIdx, inReadIdx+1, fwdRead[readSize-1-inReadIdx], backConvertedRead[inReadIdx]);
						return 1;
					}
				}
			}
		}
	}
	
	return 0;
}


