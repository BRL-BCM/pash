#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "PashDebug.h"

#define DEB_BSW 1

int bandedSW(int *scoringMatrix, char* verticalSequence, char *horizontalSequence, int sizeVerticalSequence, int band) {
	int i, j, leftVal, upVal, diagVal;
	int *prevLineH, *currentLineH;
	int score, bestScore;
	int matchGain = 1;
	int mismatchPenalty = -2;
	int gapPenalty = -3;
	int *tmpLine;
	bestScore = 0;
	int bandPlus1 = band+1;
	int bestGlobalScore = 0;

	prevLineH = scoringMatrix;
	for (j=0; j<=bandPlus1; j++) {
		prevLineH[j]=0;
	}
	for (i=0; i<sizeVerticalSequence;i++) {
		currentLineH = prevLineH+band+2;
		currentLineH[0]=0;
		fprintf(stderr, "vert base %c\n", verticalSequence[i]);
		for (j=1; j<=band; j++) {
			leftVal = currentLineH[j-1]+gapPenalty;
			upVal = prevLineH[j+1]+gapPenalty;
			if (verticalSequence[i]==horizontalSequence[i+j-1]) {
				diagVal = prevLineH[j]+matchGain;
			} else {
				diagVal = prevLineH[j]+mismatchPenalty;
			}
			if (diagVal>upVal) {
				bestScore = diagVal;
			} else {
				bestScore = upVal;
			}
			if (bestScore<leftVal) {
				bestScore = leftVal;
			}
			if (bestScore>0) {
				currentLineH[j]=bestScore;
				if (bestScore>bestGlobalScore) {
					bestGlobalScore=bestScore;
				}
			} else {
				currentLineH[j]=0;
			}

			xDEBUG(DEB_BSW, fprintf(stderr, "[%d][%d] cmp %c vs %c leftVal=%d upVal=%d  diagH %d, diagVal=%d bestScore=%d\n",
															i, j, verticalSequence[i], horizontalSequence[i+j-1], leftVal, upVal, prevLineH[j], diagVal, currentLineH[j]));
		}
		currentLineH[bandPlus1]=0;
		//tmpLine = prevLineH;
		prevLineH = currentLineH;
		//currentLineH = prevLineH;
	}
	return bestGlobalScore;
}

int main(int argc, char**argvc) {


	char *verticalSequence = "ACTGAACTG";


	char *verticalSequence2 = "ACTGACTGAG";
	char *horizSequence     = "ACTGAACTGAG@@@@@";
/*
	char *verticalSequence2 =    "ACTGAATCTGAG";
	char *horizSequence =     "@@@ACTGAACTGAG@@@@@";
*/
	int *memoryMatrix = malloc(200*sizeof(int));

	int score = bandedSW(memoryMatrix, verticalSequence2, horizSequence, strlen(verticalSequence2), atoi(argvc[1]));
	fprintf(stderr, "%d\n", score);
}
