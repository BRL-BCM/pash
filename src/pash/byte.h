/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef ARJ_BYTE
#define ARJ_BYTE

/******************************************************
 * NAME: byte.h
 *
 * DESCRIPTION:
 * Defines the 'byte' keyword. We do this in case we want to change the
 * "implementation" of a byte without changing other code.
 *
 * Right now, a 'byte' is a 'guint8' from the GTK++ library. Thus, it is
 * unsigned.
 *
 * REQUIRES:
 * �glib library (GTK++, multiplatform C library)
 *
 * AUTHOR: Andrew R. Jackson (andrewj@bcm.tmc.edu)
 *
 */


/******************************************************
 *INCLUDES
 *****************************************************/
#include <glib.h>

/******************************************************
 *DEFINES FOR TYPES
 *****************************************************/
#define byte guint8


#endif /* byte.h */
