/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef PASH_DEBUG____H
#define PASH_DEBUG____H


#include "generic_debug.h"

#define SIMPLE_FASTA_CC_OK 0
#define SIMPLE_FASTA_CC_ERROR -1


// PashLib
#define DEB_VERTICAL_KEY_CACHE 0
#define DEB_CACHE_KMER 0
#define DEB_HASH_CACHED 0
#define DEB_NEXT_CHUNK 0
#define DEB_DUMP_COLLATED_LIST 0
#define DEB_GET_COORDINATE 0
#define DEB_GET_CHUNK_ERROR 0
#define DEB_COUNT_MATCHES 0
#define DEB_COLLATION 0
#define DEB_FIRST_FASTA_PASS 0
#define DEB_FILL_REVERSE_CHUNK 0
#define DEB_OPEN_BINDUMP 0
#define DEB_ZLIB 0
#define DEB_COLLATOR1 0
#define DEB_MULTI_DIAG_COLLATION_HASH 0
#define DEB_MULTI_DIAG_COLLATION 			0
#define DEB_MULTI_DIAG_SCORING 				0
#define DEB_MULTI_DIAG_OUTPUT 				0
#define DEB_MULTI_DIAG_HSEQ_ADVANCE 	0
#define DEB_MBLOCK_OUTPUT							0
#define DEB_MULTI_RUNS 								0
#define DEB_ADAPTIVE_SCORES 					0
#define DEB_INDEXING1                 0
#define DEB_LARGEBINS									0

// various optimization steps
#define DEB_MULT_HASH_GROW 	0
#define DEB_HASH_INVERT_BINDUMP 0
#define DEB_HASH_FASTFREE 0
#define DEB_SPLIT_OUTPUT 0

// reciprocal.exe flags
#define DEB_POPULATE_FASTA_METADATA 0
#define DEB_INIT_SFFM 0
#define DEB_FIND_CHUNK 0
#define DEB_POPULATE_TABLE 0
#define DEB_READ_INPUT 0
#define DEB_WRITE_OUTPUT 0
#define DEB_REC_ZLIB 0

extern char***bindumps;
#endif
