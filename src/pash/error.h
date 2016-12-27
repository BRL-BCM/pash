/*
Copyright (c) 2004 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/

/******************************************************
 * NAME: error.h
 *
 * DESCRIPTION:
 * Defines basic error-printing messages that also report the fileName and
 * line number when printing the error. Defines a specific macro for
 * reporting out-of-memory problems.
 *
 * NOTE: these macros cause program termination. Don't call exit() on your own,
 * 			 call one of these or make some new wrappers/functions in here.
 *
 * REQUIRES:
 *  glib library (GTK++, multiplatform C library)
 *  boolean.h and byte.h for data-type standardization.
 *
 * AUTHOR: Andrew R. Jackson (andrewj@bcm.tmc.edu)
 *
 */

#ifndef ARJ_ERROR
#define ARJ_ERROR


/******************************************************
 *INCLUDES
 *****************************************************/
#include <stdio.h>

/******************************************************
 * DEFINED CONSTANTS
 *****************************************************/
#define ARJ_NOT_OK 8

/******************************************************
 * FUNCTION: CROAK
 *
 * Usage   : CROAK(msgStr, fileStr, lineNum)
 * Function:
 * Returns : [none]
 * Args    : string constant, char*, int
 */
#define CROAK(MSG, FILE, LINE) fprintf(stderr, MSG " in '%s', line %d\n", (FILE), (LINE)) ; exit(ARJ_NOT_OK)

/******************************************************
 * FUNCTION: MEM_ALLOC_CROAK
 *
 * Usage   : MEM_ALLOC_CROAK()
 * Function:
 * Returns :  [none]
 * Args    :  [none]
 */
#define MEM_ALLOC_CROAK() CROAK("ERROR: Memory allocation failed", __FILE__, __LINE__)


#endif /* error.h */
