/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef GENERIC_xDEBUG_UTILS__H
#define GENERIC_xDEBUG_UTILS__H

#define xDEBUG(flag, code) if (flag) {code;}
/** Execute a block of code and exit.*/
#define xDie(code,exitCode) {code; fflush(stdout); fflush(stderr); exit(exitCode);}
#define xDieIfNULL(var, code, exitCode) if (var==NULL) { \
    code; fflush(stdout); fflush(stderr); exit(exitCode);}
    
#endif
