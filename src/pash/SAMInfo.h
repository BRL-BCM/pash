#ifndef _SAM__INFO___H__
#define _SAM__INFO___H__

struct SAMInfo {
	guint32 basePairVariantsPositions[MAX_READ_SIZE];
  guint32 basePairVariantsPositionsInRead[MAX_READ_SIZE];
	char basePairVariantAlleles[MAX_READ_SIZE];
	guint32 numberOfBasePairVariants;
	
	guint32 cgMethylatedBasesPositions[MAX_READ_SIZE];
	guint32 cgMethylatedBasesPositionsInRead[MAX_READ_SIZE];
	guint32 numberOfCGMethylatedBases;
	guint32 chgMethylatedBasesPositions[MAX_READ_SIZE];
	guint32 chgMethylatedBasesPositionsInRead[MAX_READ_SIZE];
	guint32 numberOfCHGMethylatedBases;
	guint32 chhMethylatedBasesPositions[MAX_READ_SIZE];
	guint32 chhMethylatedBasesPositionsInRead[MAX_READ_SIZE];
	guint32 numberOfCHHMethylatedBases;
	guint32 convertedBasesPositions[MAX_READ_SIZE];
	guint32 convertedBasesPositionsInRead[MAX_READ_SIZE];
	guint32 numberOfConvertedBases;
};


#endif



