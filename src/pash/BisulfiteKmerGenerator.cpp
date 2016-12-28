/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <stdlib.h>

#include "PashDebug.h"
#include "FixedHashKey.h"
#include "BisulfiteKmerGenerator.h"

#define DEB_SETUP 0
#define DEB_GENERATE 0

BisulfiteKmerGenerator::BisulfiteKmerGenerator() {
    
}

int BisulfiteKmerGenerator::generateKmerList(char* kmer, guint32* kmerList, int maxKmerNumber) {
    //maxKmerNumber = 300;
    xDEBUG(DEB_GENERATE, fprintf(stderr, "start generating kmer list\n"));
    int kmerSize = strlen(kmer);
    int generatedKmers = 0;
    int crtIndex, direction;
    int validKmerIndex = 0;
    bases[kmerSize]='\0';
    // fill in bisulfite base
    int allC, allT;  
    for (crtIndex=0;  crtIndex<kmerSize;
         crtIndex++) {
        if (kmer[crtIndex]=='T') {
            bisulfiteBases[crtIndex].numberOfValues=2;
            bisulfiteBases[crtIndex].baseValues[0]='C';
            bisulfiteBases[crtIndex].baseValues[1]='T';
            bisulfiteBases[crtIndex].index = 0;
        } else {
            bisulfiteBases[crtIndex].numberOfValues=1;
            bisulfiteBases[crtIndex].baseValues[0]=kmer[crtIndex];
            bisulfiteBases[crtIndex].index = 0;
        }
        xDEBUG(DEB_SETUP, fprintf(stderr, "setup crtIndex %d\n", crtIndex));
    }
    
    for (crtIndex=0, direction=1;crtIndex>=0;) {
        xDEBUG(DEB_GENERATE, fprintf(stderr, "index %d direction %d kmerSize %d\n",
                                     crtIndex, direction, kmerSize));
        if (crtIndex==kmerSize) {
            validKmerIndex += 1;
            int j;
            for (j=0, allC=1, allT=1; j<kmerSize; j++) {
               if (bases[j]=='T' && kmer[j]=='T') {
                allC = 0;
               }
               if (bases[j]=='C' && kmer[j]=='T') {
                allT = 0;
               }
            }
            xDEBUG(DEB_GENERATE, fprintf(stderr," currentKmer %s allC %d allT %d\n", bases, allC, allT)); 
            //if(allC || allT || (generatedKmers<maxKmerNumber && (validKmerIndex%2==0))) {
            if(allC || allT || (generatedKmers<maxKmerNumber) ) {
            xDEBUG(DEB_GENERATE, fprintf(stderr, "current kmer=%s\n", bases));
						getKeyForSeq(bases, &kmerList[generatedKmers]);
            
            generatedKmers += 1;
            //if (generatedKmers==maxKmerNumber) {
            //   break;
            //}
            }
            
            crtIndex -=1;
            direction = 0;
        } else if (direction==1) {
            bisulfiteBases[crtIndex].index = 0;
            bases[crtIndex] = bisulfiteBases[crtIndex].baseValues[0];
            crtIndex +=1;
        } else { // direction = 0
            bisulfiteBases[crtIndex].index +=1;
            if (bisulfiteBases[crtIndex].index==bisulfiteBases[crtIndex].numberOfValues) {
                crtIndex-=1;
                direction = 0;
            }   else {
                bases[crtIndex] = bisulfiteBases[crtIndex].baseValues[bisulfiteBases[crtIndex].index];
                crtIndex +=1;
                direction = 1;
            }
        }
    }
    xDEBUG(DEB_GENERATE, fprintf(stderr, "done generating %d kmer list\n", generatedKmers));
    return generatedKmers;
}

#define DEBUG_BKG 0
#if DEBUG_BKG
int main(void) {
    BisulfiteKmerGenerator *bisulfiteKmerGenerator = new BisulfiteKmerGenerator();
    char kmerSample[MAX_KMER_SIZE+1];
		guint32 guintSeeds[40];
    strcpy(kmerSample, "ACTGTTT");
		int i;
		int numSeeds, actualSeeds;
		numSeeds=4;
    actualSeeds=bisulfiteKmerGenerator->generateKmerList(kmerSample, guintSeeds, numSeeds);
    for (i=0; i<actualSeeds; i++) {
				fprintf(stderr, "ssed [%d]=%u\n", i, guintSeeds[i]);
		}
		fprintf(stderr, "----------\n");
		numSeeds=7;
    actualSeeds=bisulfiteKmerGenerator->generateKmerList(kmerSample, guintSeeds, numSeeds);
    for (i=0; i<actualSeeds; i++) {
				fprintf(stderr, "ssed [%d]=%u\n", i, guintSeeds[i]);
		}
		fprintf(stderr, "----------\n");
		numSeeds=32;
    actualSeeds=bisulfiteKmerGenerator->generateKmerList(kmerSample, guintSeeds, numSeeds);
    for (i=0; i<actualSeeds; i++) {
				fprintf(stderr, "ssed [%d]=%u\n", i, guintSeeds[i]);
		}
		fprintf(stderr, "----------\n");
		
    return 1;
}
#endif


