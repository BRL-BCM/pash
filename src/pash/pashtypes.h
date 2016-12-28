/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/


#ifndef KJK_PASHTYPES_H 
#define KJK_PASHTYPES_H
/*******************************
 * pashtypes.h
 *
 * defines datatypes used for diagonals, words, and anything else
 * based on some type of int that we may want to change some day
 * Also some constants that are needed by both hashing/inversion and collation
 *******************************/
#define NUM_FASTA_FILES 2       // number of input files processed per job
				// for arrays relating to input files:
#define VERTICAL_FILE_NUM 0	// item 0 is the vertical file
#define HORIZONTAL_FILE_NUM 1	// item 1 is the horizontal file
// in principle one might be able to simultaneously compare a third sequence but it is not clear whether there would be any substantial advantages to this, nor how difficult it would be to modify pash in this manner

#define diagtype gint16	     // diagonals are signed, reflecting that negative numbered diagonals are supported
#define wordtype guint16     // word position identifiers

typedef struct scoreType {   // currently Pash can compute 3 different scores for a match
  wordtype bases;	// count of base identities detected in a match
  gdouble bitScore;	// uses mutual information to compute a score, influenced by both the number of base/word identities, and also their spatial organization.  See Pash manuscript for details
  guint32 scoreFactorScore;	// when using per-base scoring, this score is computed by adding the scores of all the individual k-mers in the match
} scoretype;

#define binsizetype guint16   // Used in the hashbin datatype, variables of this type are used to store the number of HashBinEntries in a hashbin
			     // NOTE if this changes, must update MAX_HASHBIN_CAPACITY in HashBin.h

#define off_tt guint32	     // datatype used for file pointers, in place of off_t.  This is left over from the days of supporting RedHat 6.2 linux machines that had some weird file pointer behavior.
			     // as these systems are no longer supported, it would most likely be okay to start using the standard off_t datatype
			     
#endif
