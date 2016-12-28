/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <glib.h>


#include "generic_debug.h"
#include "someConstants.h"
#include "BRLGenericUtils.h"

#define DEB_PARSE_COMMA_LIST  0
#define DEB_PARSE_INTLIST     0
#define DEB_ZIP 0

/** Open a gzip-ed or a bzip-ed file by setting up a pipe between (b)zcat and the current process.
 * @return a FILE pointer to current file
 * @param fileName name of the file to open
 */
FILE* BRLGenericUtils::openTextGzipBzipFile(char* fileName) {
  char catCommand[2*MAX_FILE_NAME];
  FILE *inputFilePtr;
  
  if (!strcmp( &fileName[strlen(fileName)-3], ".gz")) {
    sprintf(catCommand, "gzip -d -c %s", fileName);
    xDEBUG(DEB_ZIP, fprintf(stderr, "about to open gzip-ed file %s\n", catCommand));
    inputFilePtr = popen(catCommand, "r");
  } else {
    if (!strcmp(&fileName[strlen(fileName)-4], ".bz2")) {
      sprintf(catCommand, "bunzip2 -c  %s", fileName);
      xDEBUG(DEB_ZIP, fprintf(stderr, "about to open bzip-ed file %s\n", catCommand));
      inputFilePtr = popen(catCommand, "r");
    } else {
      xDEBUG(DEB_ZIP, fprintf(stderr, "opening text file %s\n", fileName));
      inputFilePtr = fopen(fileName, "rt");
    }
  }

  return inputFilePtr; 
}

/** Print current time
 * @param outPtr file stream to print to
 */
void BRLGenericUtils::printNow(FILE* outPtr=stdout) {
  time_t now;
	time(&now);
	fprintf(outPtr, "%.24s\n", ctime(&now));
}

/** Parse a comma separated list, and set an array of strings with the result.
@return 0 if success, 1 otherwise
*/
int BRLGenericUtils::parseCommaSeparatedList(char* commaSeparatedList, char*** stringArray, guint32 *numberOfStrings) {
  xDEBUG(DEB_PARSE_COMMA_LIST, fprintf(stderr, "starting out with comma separated list %s\n", commaSeparatedList));
  guint32 stringCount;
  guint32 stringIndex;
  guint32 inputStringLen = strlen(commaSeparatedList);
  
  for(stringIndex=0, stringCount=1; stringIndex<inputStringLen; stringIndex++) {
    if (commaSeparatedList[stringIndex]==',') {
      stringCount +=1;
    }
  }
  
  *stringArray = (char**) malloc(stringCount*sizeof(char*));
  if (stringArray==NULL) {
    fprintf(stderr, "could not allocate string array\n");
    return 1;
  }
  
  guint32 startIdx, endIdx, idx;
  *numberOfStrings = 0;
  for (startIdx = 0, stringIndex=0; startIdx<inputStringLen;) {
    for (idx=startIdx; idx<inputStringLen && commaSeparatedList[idx]!=','; idx++) {
    }
    endIdx = idx;
    if (endIdx<inputStringLen) {
      commaSeparatedList[endIdx]='\0';
    }
    if (strlen(&commaSeparatedList[startIdx])>0) {
      (*stringArray)[stringIndex]= g_strdup(&commaSeparatedList[startIdx]);
      stringIndex++;
      *numberOfStrings +=1;
      xDEBUG(DEB_PARSE_COMMA_LIST, fprintf(stderr, "added string %s\n", (*stringArray)[stringIndex-1]));
    } else {
      fprintf(stderr, "empry input string !!!\n");
      return 1;
    }
    startIdx = endIdx+1;
  }
  xDEBUG(DEB_PARSE_COMMA_LIST, fprintf(stderr, "filled in %d strings at address %p\n",
                                       *numberOfStrings, numberOfStrings));
  return 0; 
}


int BRLGenericUtils::parseListOfIntegers(guint32* intArray, char* listOfIntegers, int numIntegers) {
  int i;
  int value, index;
  int startPos;
  int len;

  xDEBUG(DEB_PARSE_INTLIST, fprintf(stderr,
                                    "parsing >> %s << for %d integers\n", listOfIntegers, numIntegers));
  index = 0;
  len = strlen(listOfIntegers);
  startPos = 0;
  for (i=0; i<len; i++) {
     xDEBUG(DEB_PARSE_INTLIST,
         fprintf(stderr,
           "string[%d]=%c\n", i, listOfIntegers[i]));
    if (listOfIntegers[i]==',') {
       xDEBUG(DEB_PARSE_INTLIST,
         fprintf(stderr,
           "found a comma at %d\n", i));
      if (startPos >= i) {
        fprintf(stderr, "incorrect start position string\n");
        exit(2);
      } else {
        value = strtol(&listOfIntegers[startPos], NULL, 10);
        intArray[index]=value;
        index++;
        startPos = i+1;
         xDEBUG(DEB_PARSE_INTLIST,
         fprintf(stderr,
           "found integer value %d\n", value));
      }
    }
  }
  if (startPos<len) {
       value = strtol(&listOfIntegers[startPos], NULL, 10);
        intArray[index]=value;
        index++;
         xDEBUG(DEB_PARSE_INTLIST,
         fprintf(stderr,
           "found integer value %d\n", value));
  }
  if (index!=numIntegers) {
    fprintf(stderr, "Parse error: expected %d integers; read %d integers\n", numIntegers, index);
    return 1;
  }
  return 0;
}
