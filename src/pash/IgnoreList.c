/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "IgnoreList.h"
#include "err.h"
#include "buffers.h"
#include "PashDebug.h"

void initIgnoreList(IgnoreList* buff, guint32 capacity)
  {
  initByteBuff(buff, capacity);
  }

//***** power_int
// integer exponent
guint32 power_int(guint32 base, guint32 exp)
  {
  guint32 ii=0, res=1;
  for(ii=1;ii<=exp;ii++)
    res=res*base;
  return res;
  }
//***** power_int

//***** setBits
// return byte set from array.  bits[0] is least significant, bits[7] is most significant
guint8 setBits(guint32 *bits)
  {
  guint8 r=0;
  if(bits[0]) r+=1;
  if(bits[1]) r+=2;
  if(bits[2]) r+=4;
  if(bits[3]) r+=8;
  if(bits[4]) r+=16;
  if(bits[5]) r+=32;
  if(bits[6]) r+=64;
  if(bits[7]) r+=128;
  return r;
  }
//***** setBits

//***** setBitsFromByteArray
guint8 setBitsFromByteArray(guint8 *bits)
  {
  guint8 r=0;
  if(bits[0]) r+=1;
  if(bits[1]) r+=2;
  if(bits[2]) r+=4;
  if(bits[3]) r+=8;
  if(bits[4]) r+=16;
  if(bits[5]) r+=32;
  if(bits[6]) r+=64;
  if(bits[7]) r+=128;
  return r;
  }
//***** setBitsFromByteArray

//***** getBits
// set array values from bits of b.  bits[0] is least significant, bits[7] is most significant
// array must be pre-allocated
void getBits(guint8 b, guint8 *bits)
  {
  int ii=0;
  for(ii=0; ii<8; ii++)
    {
    bits[ii]=0;
    if((b % 2)==1) bits[ii]=1;
    b=b/2;
    }
  }
//***** getBits


//***** getBitsToGuint32Array
void getBitsToGuint32Array(guint8 b, guint32 *bits)
  {
  int ii=0;
  for(ii=0; ii<8; ii++)
    {
    bits[ii]=0;
    if((b % 2)==1) bits[ii]=1;
    b=b/2;
    }
  }
//***** getBitsToGuint32Array


//***** inferPatternWeight
// determine weight of pattern from total number of keys
// i.e. (int) log(4,totalKeys)
// return 1 if keysRead is not a power of 4 (suggests likely error)
int inferPatternWeight(guint32 totalKeys, int *weight)
  {
  *weight=0;
  if(totalKeys==0) return 0;
  if(totalKeys<4) {*weight=1; return 1;}
  while(totalKeys>1)
    {
    (*weight)++;
    if(totalKeys % 4 > 0 && totalKeys != 1) return 1;
    if(totalKeys==4) return 0;
    totalKeys=totalKeys/4;
    }
  assert(0); // should never get here
  return 0;
  }
//***** inferPatternWeight


//***** readIgnoreList
// read ignore list file into byte buffer.  returns inferred key length
// if weight > 0 start with the assumption that this is the pattern weight
// (but adjust buffer size if it's wrong)
int readIgnoreList(FILE *input, IgnoreList *buff, int weight)
  {
  bytebuff readBuff;
  guint8 *temp=NULL;
  guint32 ii=0;
  gint32 oldCapacity = 0;
  int jj;

  if(input==NULL) dieNoUsage("readIgnoreList(): failed to read input file");
  if(weight==0) weight=11;

  readBuff.content=0;
  readBuff.pos=0;
  readBuff.capacity=power_int(4,weight)/8 + 1;  // make buffer slightly larger than expectation so we hit EOF on the first read, thus avoid growing the buffer unnecessarily
  if(readBuff.capacity<2) readBuff.capacity = 2;
  readBuff.buff=(guint8*) malloc( sizeof(guint8) * readBuff.capacity);
  if(readBuff.buff==NULL) dieNoUsage("failed to allocate ignore list read buffer (increasing size)");

  // read entire input file into readbuff
  readBuff.content=fread((void*) readBuff.buff, sizeof(guint8), readBuff.capacity, input);

  while(!feof(input))
    {
    oldCapacity = readBuff.capacity;
    readBuff.capacity = 4 * readBuff.capacity + 1;
    temp=(guint8*) realloc((void*) readBuff.buff, readBuff.capacity);
    if(temp==NULL) dieNoUsage("failed to reallocate memory for read buffer (decreasing size)");
    readBuff.buff=temp;
    readBuff.content+=fread( (void*) (readBuff.buff+oldCapacity), sizeof(guint8), readBuff.capacity - oldCapacity, input);

    }

  if(readBuff.capacity != readBuff.content)
    {
    temp=(guint8*) realloc((void*) readBuff.buff, readBuff.content);
    if(temp==NULL) dieNoUsage("failed to reallocate memory for read buffer");
    readBuff.buff=temp;
    readBuff.capacity=readBuff.content;
    }
  buff->capacity=8*readBuff.content;
  buff->content=0;
  buff->pos=0;
  buff->buff= (guint8*)malloc(sizeof(guint8)*buff->capacity);
  if(buff->buff==NULL) dieNoUsage("failed to allocate memory for ignore list");

  for(ii=0; ii<readBuff.content; ii++) {
    getBits(readBuff.buff[ii], buff->buff+ii*8);
    xDEBUG(0, for ( jj=0; jj<7; jj++) {
       if (!buff->buff[8*ii+jj]) {
 	fprintf(stderr, "key %d considered\n", ii*8+jj);
       }
    });
  }

  buff->content=8*readBuff.content;
  if(inferPatternWeight(buff->content, &weight)) dieNoUsage("number of keys is not a power of 4.  Maybe ignore list is bad.  Aborting.");

  free(readBuff.buff);

  return weight;
  }
//***** readIgnoreList

int isIgnored(guint32 key, const IgnoreList ignoreList)
  {
  if(key > ignoreList.content) dieNoUsage("requested sequence key exceeds size of ignore list");
  return ignoreList.buff[key];
  }

guint32 numKeysInIgnoreList(const IgnoreList ignoreList) {return ignoreList.content;}
