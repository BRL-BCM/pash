#ifndef KJK_ERR_H
#define KJK_ERR_H

#include <stdlib.h>
#include <assert.h>

#define die(msg)   {fprintf(stderr,"*****\nFatal error\n%s\n*****\n",(msg)); printUsage(); assert(0); exit(1);}

#define dieNoUsage(msg)   {fprintf(stderr,"*****\nFatal error\n%s\n*****\n",(msg)); assert(0); exit(1);}

#define warning(msg) {fprintf(stderr,"*****\nWarning:\n%s\n*****\n", (msg));}

#endif

