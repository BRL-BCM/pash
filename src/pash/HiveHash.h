/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#ifndef HiveHash__H__
#define HiveHash__H__

#include <glib.h>
#include <stdio.h>

/* Collapsed hash class
*/

class IntListRunner {
public:
  guint32 *list;
  guint32 left;
};

/**
* The main difference between Pash 2.0 and Pash 3.0: a hive kmer hash that
* collapses hashtables for multiple offsets.
*/
class HiveHash {
protected:
  int hashSize;

  guint32 numberOfEntriesIndex;
  double memoryFootprint;
  int numberOfHashValues;
  int numberOfKeys;
  
  GSList* poolList;
  guint32 *currentPool;
  guint32 individualPoolSize;
  double kmerPercent;
public:
  HiveHash(int size, guint32 keepKmerPercent);
  int addEntry(guint32 key, guint32 value, guint32 offset);
  int addEntryX(guint32 key, guint32 value, guint32 offset);
  int markEntry(guint32 key);
  double getMemoryFootprint();
  void dumpHash(FILE *filePtr);
  void checkHash();
  int getIntListRunner(guint32 key, IntListRunner* listRunner);
  void printStatistics(FILE* filePtr);
  ~HiveHash();
  guint32** hashSkeleton;
  int addEntryXX(guint32 key, guint32 value, guint32 offset);
  int allocateHashMemory();
  void checkHashXX();
};

#endif
