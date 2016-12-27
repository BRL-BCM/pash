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

