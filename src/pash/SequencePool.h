#ifndef _SequencePool_____H_____
#define _SequencePool_____H_____

#include <glib.h>

#define POOL_INCREMENT 10000

struct ASeqPool {
  char *sequence;
  guint32 used;
  guint32 capacity;
};

class SequencePool {
  ASeqPool* poolSkeleton;
  guint32 poolSize;
  guint32 poolCapacity;
  guint32 latestPool;
public:
  SequencePool();
  ~SequencePool();
  char* addSequence(char* fwd);
};

#endif
