/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef KJK_ERR_H
#define KJK_ERR_H

#include <stdlib.h>
#include <assert.h>

#define die(msg)   {fprintf(stderr,"*****\nFatal error\n%s\n*****\n",(msg)); printUsage(); assert(0); exit(1);}

#define dieNoUsage(msg)   {fprintf(stderr,"*****\nFatal error\n%s\n*****\n",(msg)); assert(0); exit(1);}

#define warning(msg) {fprintf(stderr,"*****\nWarning:\n%s\n*****\n", (msg));}

#endif

