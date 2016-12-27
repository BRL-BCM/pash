#include <stdio.h>
#include <string.h>
#include <glib.h>
#include "PashDebug.h"
#include "FixedHashKey.h"

#define MAX_KMER_SIZE 16

typedef struct {
    char baseValues[4];
    int numberOfValues;
    int index;
} BisulfiteBase ;
    
class BisulfiteKmerGenerator {
    BisulfiteBase bisulfiteBases[MAX_KMER_SIZE];
    char bases[MAX_KMER_SIZE];
public:
    BisulfiteKmerGenerator();
    int generateKmerList(char* kmer, guint32* kmerList, int maxKmerNumber);
};


