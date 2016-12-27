/*
Copyright (c) 2004 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/


/******************************************************
 * NAME: Pattern.c
 *
 * DESCRIPTION:
 * 
 * Function to generate random sampling patterns 
 * Using this is not recommended -- intelligently selected patterns will perform better than random ones
 * Furthermore, the implementation can generate patterns that do not start/end with 1, thus the effective pattern length is shorter than desired
 * (furthermore, a user has reported that patterns not starting with 1 result in a program crash)
 *
 * These have not been used/tested since Ruben authored them.
 *
 * AUTHOR: Ruben Valas (rvalas@andrew.cmu.edu)
 */

/******************************************************
 *INCLUDES
 *****************************************************/

#include <stdlib.h>	// random number generation
#include <glib.h>
#include "Pattern.h"


/*****************************************************
 *Functions
 *****************************************************/
/******************************************************
 * FUNCTION: getPattern
 *
 * Usage   : getPattern(kmerSize, rangeSize) 
 *
 * generate a random sampling mask
 *
 */

guint8 getPattern(guint8* pattern, guint16 kmerSize, guint16 rangeSize)
{
	
	int s, k;
	float r, limit;
	for(k=s=0; k<rangeSize; k++)/*loop through sequence*/
		{
		
			r=((float)rand())/((float)RAND_MAX);
			limit=((float)(kmerSize-s))/((float)(rangeSize-k));
			if(r<limit)
				{
					s++;
					pattern[k]=1;
				}
			else
				pattern[k]=0;
			if(s==kmerSize)/*made enough 1's*/
				break;

		}
	return 0;
}
