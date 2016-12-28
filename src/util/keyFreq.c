/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

/* keyFreq.c
see usage information below in function printUsage()
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <ctype.h>
#include <assert.h>
#include "keyFreq.h"
#include "FixedHashKey.h"
#include "err.h"
#include "PashDebug.h"
#include "../pash/Mask.h"

#define KEYFREQ_VERSION "0.9"

#define DEB_SORT_KMERS 0
#define DEB_FOF 0

#define MAX_FOF_LINE 2048
#define DEB_REVCOMP 0

//***** GLOBAL VARIABLES
/** Whether output should be human-readable.*/
gint32 readable=0;
/** Array giving sampling pattern.*/
gint32 *pattern=NULL;
Mask mask;
/** Location in argv of output file name.*/
gint32 outputID=0;
/** Number of bases in sampling pattern.*/
gint32 patternLen=0;
/** Number of sampled positions in pattern.*/
gint32 patternWeight=0;

guint32 debug=0;

/** Print usage information*/
void printUsage() {

	fprintf(stderr,"keyFreq - tool distributed with Pash version 3.01.03\n"
			"For a specified sampling pattern, generates global list of k-mer frequencies\n"
			"from the specified set of sequence files in  .fa format.\n"
			"Usage: \n"
			"pash3_keyFreq -o <outputFile> -p <101001101> <-h> myInput1.fa ...\n"
			"  '-h' if human readable output is desired (otherwise lists in a flat binary format, 4 bytes per key)\n"
			"  '-p' <sampling pattern> (e.g. 11011 would sample the two positions, skip one position, then\n"
			"       sample the next two), to use predefined pattern choose one of the following: 8from14, 9from15,\n"
			"       10from16, 11from18, 12from18, 13from21, 14from21 (default is 12from18)\n"
			"\n"
			"Predefined sampling patterns:\n"
			"  8from14  = 11001001010111\n"
			"  9from15  = 110110101000111\n"
			"  10from16 = 1101100011010111\n"
			"  11from18 = 111010010100110111\n"
			"  12from18 = 111010110100110111\n"
			"  13from21 = 111011011000110101011\n"
			"  14from21 = 111011100101100101111\n"
			"\n"
			"\n");
	return;
}



/** Converts base to a 0-4 number
@param c base letter */
static inline int baseToNum(char c)
{
	if(c=='a' || c=='A') return 0;
	if(c=='t' || c=='T') return 1;
	if(c=='g' || c=='G') return 2;
	if(c=='c' || c=='C') return 3;
	if(c=='n' || c=='N') return 4;
	return 5;
}

/** Checks if the letter corresponds to a base
@param c alleged base letter*/
static inline int isBase(char c)
{
	if(c=='a' || c=='A') return 1;
	if(c=='t' || c=='T') return 1;
	if(c=='g' || c=='G') return 1;
	if(c=='c' || c=='C') return 1;
	return 0;
}

static inline void revComp(const char* const seq, char *rev, int len)
{
	int ii=0;
	char c=0;
	xDEBUG(DEB_REVCOMP, fprintf(stderr, "revcomp for %s\n",seq););
	for(ii=0;ii<len;ii++)
	{
		c=toupper(seq[ii]);
		if(c=='A') rev[len-ii-1]='T';
		else if(c=='T') rev[len-ii-1]='A';
		else if(c=='G') rev[len-ii-1]='C';
		else if(c=='C') rev[len-ii-1]='G';
		else            rev[len-ii-1]='N';
		xDEBUG(DEB_REVCOMP, fprintf(stderr, "[%c]->[%c]\n",seq[ii], rev[len-ii-1]););
	}
}

int main(int argc, char ** argv)
{
	guint32 ii=0,   // input file (location in argv of name)
			firstFile=0;    // first input file (location in argv of name)
	FILE *output;
	guint32* freq = NULL;
	KmerFreqEntry *kfreq=NULL;
	char *seq=(char*) malloc(17);
	long freqSum, freqMin, freqMax, freqMean, freqMedian;
	int nonZeroFreqKmers, firstNonZeroIndex;
	FILE* fofFilePtr;
	char fofLine[MAX_FOF_LINE];
	char fileName[MAX_FOF_LINE], singleFastaFileName[MAX_FOF_LINE];
	int len;

	firstFile=parse(argc, argv);
	if(debug) printParseResults(argc, argv, firstFile);

	if(strcmp("-",argv[outputID])==0) output=stdout;
	else output=fopen(argv[outputID],"w");
	if(output==NULL) die("unable to open output file");

	if (readable) {
		kfreq=(KmerFreqEntry*)calloc(power_int(4,patternWeight),sizeof(KmerFreqEntry));
		if (kfreq == NULL) {
			die("could not allocate memory for kfreq\n");
		}

		for (ii =0; ii < power_int(4, patternWeight); ii++) {
			kfreq[ii].kmerIndex = ii;
			// kfreq[ii].kmerFrequency = 0;
		}
	}

	freq=(guint32*)calloc(power_int(4,patternWeight),sizeof(guint32));

	/** fprintf(stderr, "kFreq: %ud bytes, freq: %ud bytes\n",
	power_int(4,patternWeight)*sizeof(KmerFreqEntry), 	power_int(4,patternWeight)*sizeof(guint32));
	 */
	if (freq == NULL) {
		die("could not allocate memory for freq\n");
	}


	fprintf(stderr, "starting reading input file(s)\n");
	for(ii=firstFile; ii<argc; ii++) { // loop ii over input files
		// check if .fof file
		strcpy(fileName, argv[ii]);
		len=strlen(fileName);
		if( len>4 && toupper(fileName[len-4]) == '.' && toupper(fileName[len-3]) == 'F' &&
				toupper(fileName[len-2]) == 'O' &&	toupper(fileName[len-1]) == 'F' ) {
			fofFilePtr = fopen(fileName, "rt");
			while (fgets(fofLine, MAX_FOF_LINE, fofFilePtr)) {
				sscanf(fofLine, " %s", singleFastaFileName);
				processInputFile(singleFastaFileName, freq, kfreq);
			}
			fclose(fofFilePtr);
		} else {
			processInputFile(fileName, freq, kfreq);
		}
	}
	fprintf(stderr, "finished reading input file(s)\n");

	// sort kmers by frequency
	if (readable) {
		nonZeroFreqKmers = 0;
		freqSum = freqMean = freqMedian = 0;
		fprintf(stderr, "starting sorting\n");
		qsort(kfreq, power_int(4, patternWeight), sizeof(KmerFreqEntry),
				compareKmerFrequencyEntry);
		fprintf(stderr, "finished sorting\n");
		// print the sorted list
		fprintf(stderr, "starting writing human readable output\n");

		freqMin = freqMax = kfreq[ power_int(4, patternWeight)-1].kmerFrequency;
		for(ii=0; ii<power_int(4,patternWeight); ii++)
		{
			if (kfreq[ii].kmerFrequency>0) {
				nonZeroFreqKmers ++;
				freqSum += kfreq[ii].kmerFrequency;
				if (freqMin>kfreq[ii].kmerFrequency) {
					freqMin = kfreq[ii].kmerFrequency;
				}
				getSeqForKey(kfreq[ii].kmerIndex, patternWeight, &seq);
				fprintf(output,"%s\t%i\t%i\n", seq ,  kfreq[ii].kmerFrequency, kfreq[ii].kmerIndex);
			}
		}
		// determine the median

		for(ii=0; ii<power_int(4,patternWeight); ii++) {
			if (kfreq[ii].kmerFrequency>0) {
				firstNonZeroIndex = ii;
				break;
			}
		}
		//printf("The first nonzero index is %d\n", ii);
		if (nonZeroFreqKmers%2 == 1) {
			freqMedian = kfreq[ii+(nonZeroFreqKmers/2)].kmerFrequency;
		} else {
			freqMedian = (kfreq[ii-1+(nonZeroFreqKmers/2)].kmerFrequency+kfreq[ii+(nonZeroFreqKmers/2)].kmerFrequency)/2;
		}

		freqMean = freqSum / nonZeroFreqKmers;
		fprintf(stderr, "Statistics: min=%ld max=%ld mean=%ld median=%ld\n",
				freqMin, freqMax, freqMean, freqMedian);
		fprintf(stderr, "finished writing human readable output\n");
	}
	if(!readable)  {
		fwrite( (void*) freq,  sizeof(guint32), power_int(4,patternWeight), output);
	}
	return 0;
}  // main


/** Parse command line parameters
@return index of first input file in argv*/
int parse(int argc, char ** argv)
{
	int ii=0, jj=0;

	if (argc==1) {
		printUsage();
		exit(0);
	}
	if (!strcmp(argv[1],"--help")) {
		printUsage();
		exit(0);
	}

	for(ii=1;ii<argc;ii++)
	{
		if(strcmp("-h",argv[ii])==0) readable=1;
		else if(strcmp("-p",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify pattern after -p");
			setMask(argv[ii], &mask);
			patternLen=mask.maskLen;
			patternWeight = mask.keyLen;
			if(patternLen==0) die("bad pattern");
			if(pattern!=NULL) die("attempt to specify pattern twice");
			pattern=(int*)malloc(sizeof(gint32)*patternLen);
			if(pattern==NULL) die("pattern memory allocation failed");
			for(jj=0;jj<patternLen;jj++)
			{
				pattern[jj]=mask.mask[jj];
			}
			if(patternWeight<=0 || patternWeight > MAXPATTERN) die("Pattern must have between 1 and 15 sampled positions");
		}
		else if(strcmp("-o",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify output file after -o");
			outputID=ii;
		}
		else if(strcmp("-d",argv[ii])==0) {debug++; fprintf(stderr,"\npash3_keyFreq version 3.01.03.  Debugging activated.\n");}
		else
		{
			if(patternLen==0) {
				/* set default pattern */
				setMask("12from18", &mask);
				patternLen=mask.maskLen;
				patternWeight = mask.keyLen;
				pattern=(int*)malloc(sizeof(gint32)*patternLen);
				if(pattern==NULL) die("pattern memory allocation failed");
				for(jj=0;jj<patternLen;jj++)
				{
					pattern[jj]=mask.mask[jj];
				}
			}
			if(outputID==0) die("must specify output file");
			if(patternLen > BUFFSIZE) die("pattern length greater than size of input buffer; use a shorter pattern or recompile with increased BUFFSIZE");
			return ii;
		}
	}
	die("must specify at least one input file");
	return argc;
}

/** Reads input file. Fill seqBuff with bases, wrapping as needed
Deflines are stripped, and seqBuff only contains sequences from one record at a time
@param input input file
@param mode data mode
@param charbuff read buffer
@param seqBuff  sequence buffer
@param patternLen pattern length
@return 0 if input buffer is at EOF and readBuff is empty, otherwise 1*/
static inline int readInputFile(FILE *input, int *mode, charbuff *readBuff, charbuff *seqBuff, int patternLen)
{
	char nextChar;
	seqBuff->pos=0;
	if(debug > 1) fprintf(stderr,"starting readInputFile");
	if(feof(input) && readBuff->content==0) return 0;

	if(debug) fprintf(stderr,"mode is %i\n", *mode);
	if(*mode==DEFLINE) seqBuff->content=0;
	if(*mode==SEQ)   // if looking at sequence, wrap end of buffer back to beginning
	{
		if(seqBuff->content >= patternLen)
		{
			memmove((void*)seqBuff->buff, (void*)(seqBuff->buff +seqBuff->content -patternLen +1), patternLen-1);
			seqBuff->content=patternLen-1;
		}
	}
	while(1)  // loop once per byte of readBuff processed
	{
		// refill read buffer, if necessary
		if(debug > 2) fprintf(stderr,"%i\n",readBuff->pos);
		if(readBuff->content==0 || readBuff->pos >= readBuff->content)
		{
			readBuff->content=fread((void*)(readBuff->buff), 1, readBuff->capacity, input);
			if(ferror(input)) die("error reading from input file");
			readBuff->pos=0;
		}
		if(readBuff->content==0) return 1;	// no more sequence in file
		nextChar=readBuff->buff[readBuff->pos];
		if(*mode==NEWLN && nextChar=='>') {*mode=DEFLINE; readBuff->pos++; return 1;}
		if(*mode==NEWLN && nextChar!='>') {*mode=SEQ;}
		if(nextChar=='\n') {*mode=NEWLN;}
		if(*mode==SEQ && nextChar!='\0') {seqBuff->buff[seqBuff->content]=nextChar; seqBuff->content++;}
		readBuff->pos++;
		if(seqBuff->content == seqBuff->capacity) return 1;	// sequence buffer is full
	}
	assert(0);  // should never get here
	return 0;
}

/** Get sequences from input file, convert to keys, store frequencies in freq
@param inputName name of input file
@param freq frequency array
@param kfreq sortable kmer frequency array */
void processInputFile(const char * inputName, guint32 *freq, KmerFreqEntry* kfreq)
{
	int mode=NEWLN;	// type of data encountered at end of last read
	guint32 jj=0,   // position in input buffer
			kk=0,   // position in sampling pattern
			pos=0,  // base position in current word
			badWord=0,      // set to 1 if current word contains bad characters
			badReverseWord=0;      // set to 1 if current reverse word contains bad characters
	charbuff readBuff, seqBuff;
	FILE *input=NULL;
	char seq[17], revseq[17], rev[17];
	char currentBase;
	char currentRevBase;

	guint32 key;  // hash key for current word

	readBuff.buff=(char*) calloc(BUFFSIZE,sizeof(char));
	if(readBuff.buff==NULL) die("failed to allocate memory for read buffer");
	readBuff.capacity=BUFFSIZE;
	readBuff.content=0;
	readBuff.pos=0;

	seqBuff.buff=(char*) calloc(BUFFSIZE,sizeof(char));
	if(seqBuff.buff==NULL) die("failed to allocate memory for sequence buffer");
	seqBuff.capacity=BUFFSIZE;
	seqBuff.content=0;
	seqBuff.pos=0;

	fprintf(stderr,"reading file: %s\n",inputName);
	if(strcmp("-",inputName)==0) input=stdin;
	else input=fopen(inputName,"r");
	if(input==NULL) die("failed to open input file");

	while(readInputFile(input, &mode, &readBuff, &seqBuff, patternLen)>0)
	{
		if(debug) fprintf(stderr,"%i bytes in input buffer\n%i bytes in sequence buffer\n", readBuff.content, seqBuff.content);
		if(debug) {for(jj=0; jj<seqBuff.content; jj++) {
			if(seqBuff.buff[jj]==0) fprintf(stderr,"."); else
				fprintf(stderr,"%c",seqBuff.buff[jj]);} fprintf(stderr,"\n");}

		if( seqBuff.content + 1 > patternLen )
			for(jj=0;jj<seqBuff.content-patternLen+1;jj++)
			{
				badWord=0;
				badReverseWord=0;
				pos=0;

				for(kk=0;kk<patternLen;kk++)
				{
					currentBase = seqBuff.buff[jj+kk];
					currentRevBase = seqBuff.buff[jj+patternLen-1-kk];
					xDEBUG(DEB_REVCOMP, fprintf(stderr, "base %c revBase %c\n", currentBase, currentRevBase););
					if (currentBase > 'a' && currentBase < 'z')  {
						currentBase += 'A'-'a';
					}
					if (currentRevBase > 'a' && currentRevBase < 'z')  {
						currentRevBase += 'A'-'a';
					}
					if(pattern[kk]) {
						if (!(currentBase == 'A' || currentBase == 'C' ||
								currentBase == 'G' || currentBase=='T') ) {
							badWord=1;
						}
						seq[pos]=currentBase;
						if (!(currentRevBase == 'A' || currentRevBase == 'C' ||
								currentRevBase == 'G' || currentRevBase=='T') ) {
							badReverseWord=1;
						}
						rev[patternWeight-pos-1]=currentRevBase;
						pos ++;
					}
				} // build word from current sequence
				seq[pos]='\0';
				rev[pos]='\0';
				revseq[pos]='\0';
				xDEBUG(DEB_REVCOMP, fprintf(stderr, "pos=%d patternLen=%d\n", pos, patternLen));
				if(!badWord) {
					getKeyForSeq(seq, &key);
					if(debug) fprintf(stderr,"getting key for %s : %u\n",seq, key);
					if((freq[key] + 1) > 0 ) {
						freq[key]++;  // don't increment if it would overflow
						if (readable) {
							kfreq[key].kmerFrequency ++;
						}
					}
				}  // good word, make and store
				if(!badReverseWord){ 	// add reverse compliment
					revComp(rev, revseq, pos);
					getKeyForSeq(revseq, &key);
					if(debug) fprintf(stderr,"getting key for %s->%s : %u\n",rev, revseq, key);
					if((freq[key] + 1) > 0 ) {
						freq[key]++;  // don't increment if it would overflow
						if (readable) {
							kfreq[key].kmerFrequency ++;
						}
					}
				}
			}  // loop jj over sequence buffer
		if(debug) fprintf(stderr,"\n");
	}  // while read from input file
	fclose(input); input=NULL;
	free(readBuff.buff);
	free(seqBuff.buff);
} // processInputFile

/** Print the program arguments*/
void printParseResults(int argc, char ** argv, int firstFile)
{
	guint32 ii=0;
	fprintf(stderr,"patternLen=%i\npattern=",patternLen);
	for(ii=0;ii<patternLen;ii++) fprintf(stderr,"%i",pattern[ii]);
	fprintf(stderr,"\n");
	fprintf(stderr,"output file: %s\n",argv[outputID]);
	fprintf(stderr,"%i input file(s):\n",argc-firstFile);
	for(ii=firstFile;ii<argc;ii++) fprintf(stderr,"%s\n",argv[ii]);
	fprintf(stderr,"\n");
}




/** Compare two kmer frequency entries based on frequency.*/
int compareKmerFrequencyEntry(const void *p1, const void*p2 ) {
	if ( ((KmerFreqEntry*)p1)->kmerFrequency < ((KmerFreqEntry*)p2)->kmerFrequency) {
		return -1;
	} else 	if ( ((KmerFreqEntry*)p1)->kmerFrequency > ((KmerFreqEntry*)p2)->kmerFrequency) {
		return 1;
	} else {
		return 0;
	}
}

