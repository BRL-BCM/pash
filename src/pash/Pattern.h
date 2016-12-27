/*
Copyright (c) 2004 Baylor College of Medicine.
Use of this software is governed by a license.  See the included file
LICENSE.TXT for details.
*/


/******************************************************
 * NAME: Pattern.h
 *
 * DESCRIPTION:
 * Header file for Pattern.c.
 * 
 * Function to generate random sampling patterns 
 * Using this is not recommended -- intelligently selected patterns will perform better than random ones
 * Furthermore, the implementation can generate patterns that do not start/end with 1, thus the effective pattern length is shorter than desired
 * (furthermore, a user has reported that patterns not starting with 1 result in a program crash)
 *
 * These have not been used/tested since Ruben authored them.
 * 
 * AUTHOR: Ruben Valas (rvalas@andrew.cmu.edu)
 *
 */

#ifndef ARJ_PATTERN
#define ARJ_PATTERN

guint8 getPattern(guint8* pattern, guint16 kmerSize, guint16 rangeSize);

#endif
