/*
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
*/

#include <stdio.h>
#include <glib.h>
#include "buffers.h"
#include "err.h"

void initGuint32Buff(guint32buff * buff, guint32 capacity)
  {
  buff->pos=0;
  buff->capacity=capacity;
  buff->content=0;
  buff->buff=(guint32*)malloc(sizeof(guint32) * buff->capacity);
  if(buff->buff==NULL) dieNoUsage("failed to allocate buffer memory");
  }

void initByteBuff(bytebuff * buff, guint32 capacity)
  {
  buff->pos=0;
  buff->capacity=capacity;
  buff->content=0;
  buff->buff=(guint8*)malloc(sizeof(guint8) * buff->capacity);
  if(buff->buff==NULL) dieNoUsage("failed to allocate buffer memory");
  }

void initCharBuff(charbuff * buff, guint32 capacity)
  {
  buff->pos=0;
  buff->capacity=capacity;
  buff->content=0;
  buff->buff=(char*)malloc(sizeof(char) * buff->capacity);
  if(buff->buff==NULL) dieNoUsage("failed to allocate buffer memory");
  }

gint32 bufferedRead(bytebuff *readBuff, FILE* input)
  {
  readBuff->content=fread((void*)readBuff->buff, sizeof(guint8), readBuff->capacity, input);
  if(ferror(input)) dieNoUsage("error reading input file");
  readBuff->pos=0;
  return readBuff->content;
  }

int bufferedWrite(const bytebuff* const writeBuff, FILE * output)
  {
  gint32 written=fwrite((void*)writeBuff->buff, sizeof(guint8), writeBuff->content, output);
  if(ferror(output) || written!=writeBuff->content) dieNoUsage("error writing output file");
  return writeBuff->content;
  }

void freeByteBuff(bytebuff * buff)
  {
  buff->pos=0;
  buff->capacity=0;
  buff->content=0;
  if(buff->buff!=NULL) free(buff->buff); buff->buff=NULL;
  }

void freeCharBuff(charbuff * buff)
  {
  buff->pos=0;
  buff->capacity=0;
  buff->content=0;
  if(buff->buff!=NULL) free(buff->buff); buff->buff=NULL;
  }

void freeGuint32Buff(guint32buff * buff)
  {
  buff->pos=0;
  buff->capacity=0;
  buff->content=0;
  if(buff->buff!=NULL) free(buff->buff); buff->buff=NULL;
  }
