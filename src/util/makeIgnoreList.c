/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

/* makeIgnoreList.c
   Reads key frequency table generated by keyFreq.exe (in binary format) and outputs
   'ignore list' file for Pash.  Ignore list is a bit vector, where '1' means ignore and
   '0' means don't ignore
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <glib.h>
#include "FixedHashKey.h"
#include "IgnoreList.h"
#include "err.h"

//***** GLOBAL VARIABLES
gint32	maxCutoff=-1,	// max frequency cutoff to apply
minCutoff=-1,	// min frequency cutoff to apply
inputID=0,	// location in argv of input file name
outputID=0;	// location in argv of output file name

int debug=0;	// activate debug mode
//***** GLOBAL VARIABLES


//***** printUsage
void printUsage() {
	fprintf(stderr,"makeIgnoreList - tool distributed with Pash version 3.01.03\n"
			"Reads key frequency table generated by pash3_keyFreq (in binary format) and\n"
			"outputs 'ignore list' file for Pash.  Keys with occurance > specified\n"
			"maxCutoff or < specified minCutoff are ignored.  May only specify a maximum or\n"
			"minimum, in which case the other is ignored.  Ignore list output is a bit vector, \n"
			"where '1' means ignore and '0' means don't ignore.\n"
			"Usage: \n"
			"pash3_makeIgnoreList -i <inputFile> -o <outputFile> -c <max frequency cutoff> -m <min frequency cutoff>\n"
			"\n");
}
//***** printUsage


//***** parse_makeIgnoreList
// parse input parameters
void parse_makeIgnoreList(int argc, char ** argv)
{
	int ii=0;

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
		if(strcmp("-i",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify input file after -i");
			inputID=ii;
		}
		else if(strcmp("-o",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify output file after -o");
			outputID=ii;
		}
		else if(strcmp("-c",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify max cutoff after -c");
			maxCutoff=atoi(argv[ii]);
			if(maxCutoff<0) die("maximum cutoff must be non-negative");
		}
		else if(strcmp("-m",argv[ii])==0)
		{
			ii++;
			if(ii>=argc) die("must specify min cutoff after -c");
			minCutoff=atoi(argv[ii]);
			if(minCutoff<0) die("minimum cutoff must be non-negative");
		}
		else if(strcmp("-d",argv[ii])==0) debug++;
		else die("parameter error");
	}
	if(inputID==0 || outputID==0 || (maxCutoff==-1 && minCutoff==-1)) die("missing parameter");
	if(minCutoff > maxCutoff && minCutoff!=-1 && maxCutoff!=-1) die("minCutoff must not be larger than maxCutoff");
}
//***** parse_makeIgnoreList


//***** processFile
void processFile(FILE *inputFile, FILE *outputFile)
{
	guint32buff readBuff;
	IgnoreList outputBuff;
	int ii=0, jj=0;
	long int total=0, rejected=0;
	int percent1, percent2;

	guint8 *testBuff=(guint8*) calloc(8,sizeof(guint8));

	readBuff.capacity = BUFFSIZE;
	readBuff.buff = (guint32*) malloc(sizeof(guint32) * readBuff.capacity);
	if(readBuff.buff==NULL) die("allocation of file read buffer failed");
	readBuff.pos = 0;
	readBuff.content = 0;

	outputBuff.capacity = BUFFSIZE/8;
	outputBuff.pos = 0;
	outputBuff.buff = (guint8 *) malloc(sizeof(guint8) * outputBuff.capacity);
	if(outputBuff.buff==NULL) die("allocation of output buffer failed");
	outputBuff.content = 0;

	while(!feof(inputFile))
	{
		readBuff.content=fread((void*) readBuff.buff, sizeof(guint32), readBuff.capacity, inputFile);
		for(ii=0; ii<readBuff.content; ii++)
		{
			total++;
			if((minCutoff!=-1 && readBuff.buff[ii] < minCutoff) || (maxCutoff!=-1 && readBuff.buff[ii] > maxCutoff)) {
				readBuff.buff[ii]=1;
				++rejected;
			} else {
				readBuff.buff[ii]=0;
			}
			if(debug) {
				if( readBuff.buff[ii]) {
					fprintf(stderr,"key %ld ignored\n", total-1);
				} else {
					fprintf(stderr,"key %ld considered\n", total-1);
				}
			}
			if((( ii + 1) % 8) == 0)
			{
				outputBuff.buff[outputBuff.content++]=setBits(readBuff.buff+ii-7);
				if(debug && outputBuff.buff[outputBuff.content-1])
				{
					fprintf(stderr,"writing %i, bits are ",outputBuff.buff[outputBuff.content-1]);
					getBits(outputBuff.buff[outputBuff.content-1], testBuff);
					for(jj=0;jj<8;jj++) fprintf(stderr,"%i",testBuff[jj]);
					fprintf(stderr," (least to most significant)\n");
				}
			} // store 8th byte
		} // process input buffer content
		fwrite((void*) outputBuff.buff, sizeof(guint8), outputBuff.content, outputFile);
		if(debug) fprintf(stderr,"%i bytes written to output\n", outputBuff.content);
		outputBuff.content=0;
	}
	ii=0;

	percent1 = rejected * 100 / total;
	percent2 = (rejected * 10000 / total) % 100;
	printf("Procent of rejected kmers: %d.%02d\n", percent1, percent2);

	if(inferPatternWeight(total, &ii))
		fprintf(stderr,"Warning: %li total keys processed -- not a power of 4.  Maybe input is wrong?\n", total);
	else if(debug) fprintf(stderr,"%li total keys processed.  Inferred key length is %i\n",total, ii);

	free(readBuff.buff);
	free(outputBuff.buff);
}
//***** processFile


//***** main
int main(int argc, char ** argv)
{
	FILE *inputFile=NULL, *outputFile=NULL;
	parse_makeIgnoreList(argc, argv);
	if(strcmp("-",argv[inputID])==0) inputFile=stdin;
	else inputFile=fopen(argv[inputID],"r");
	if(inputFile==NULL) die("failed to open input file");

	if(strcmp("-",argv[outputID])==0) outputFile=stdout;
	else outputFile=fopen(argv[outputID],"w");
	if(outputFile==NULL) die("failed to open output file");

	processFile(inputFile, outputFile);
	fclose(inputFile);
	fclose(outputFile);
	return 0;
}
//**** main
