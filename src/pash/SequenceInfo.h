/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef __SEQUENCE_INFO_____H___
#define __SEQUENCE_INFO_____H___


typedef struct {
    const char* sequenceName;
    guint32 sequenceLength;
    guint32 bestAnchoringScore;
    guint32 bestSWScore;
		guint32 bestSkeletonScore;
    guint32 bestScoreMappings;
    guint32 bestChrom;
    guint32 bestStart;
    guint32 passingMappings;
} SequenceInfo;


#endif

