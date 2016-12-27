#ifndef KJK_IGNORELIST_H
#define KJK_IGNORELIST_H

#include <stdio.h>
#include <glib.h>
#include "buffers.h"

#define BUFFSIZE (1024*64)  // read buffer size (number of members).  must be a multiple of 8, so an integral number of ignore list bytes can be derived from a buffer full of of key counts
#define MAXPATTERN 15

typedef bytebuff IgnoreList;

void initIgnoreList(IgnoreList* buff, guint32 capacity);

guint32 power_int(guint32 base, guint32 exp);
guint8 setBits(guint32 *bits);
guint8 setBitsFromByteArray(guint8 *bits);

void getBits(guint8 b, guint8 *bits);
void getBitsToGuint32Array(guint8 b, guint32 *bits);

int inferPatternWeight(guint32 keysRead, int *weight);

int readIgnoreList(FILE *input, IgnoreList *buff, int weight);
int isIgnored(guint32 key, const IgnoreList ignoreList);
guint32 numKeysInIgnoreList(const IgnoreList ignoreList);

#endif // KJK_IGNORELIST_H
