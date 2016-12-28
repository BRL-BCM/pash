/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef __Generic_Utils___H_
#define __Generic_Utils___H_

#include <stdio.h>
#include <glib.h>

class BRLGenericUtils {
public:
  static FILE* openTextGzipBzipFile(char* fileName);
  static void printNow(FILE* outPtr);
  static int parseCommaSeparatedList(char* commaSeparatedList,
                                     char*** stringArray, guint32 *numberOfStrings);
	static int parseListOfIntegers(guint32* intArray, char* listOfIntegers, int numIntegers);

};


#endif
