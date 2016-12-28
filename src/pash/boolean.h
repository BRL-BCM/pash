/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef ARJ_BOOLEAN
#define ARJ_BOOLEAN

/******************************************************
 * NAME: boolean.h
 *
 * DESCRIPTION:
 * Defines the 'boolean' keyword. We do this in case we want to change the
 * "implementation" of a boolean without changing other code.
 *
 * Right now, a 'boolean' is a 'gboolean' from the GTK++ library.
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
#define boolean gboolean

/******************************************************
 *DEFINES FOR INLINE FUNCTIONS
 *****************************************************/

/******************************************************
 * FUNCTION: BOOLSTR
 *
 * Usage   :�BOOLSTR(var)
 * Function:
 *         ��Returns string "TRUE" or "FALSE" depending on whether value passed
 *           is true or false, in terms of C's view on this.
 * Returns : �"TRUE" or "FALSE"
 * Args    : �variable
 */
#define BOOLSTR(x) (x) ? "TRUE" : "FALSE"

#endif /* boolean.h */
