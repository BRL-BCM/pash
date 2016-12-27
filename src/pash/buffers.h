#ifndef KJK_BUFFERS_H
#define KJK_BUFFERS_H

#include <glib.h>

typedef struct Guint32Buff {
  guint32 *buff;
  gint32 capacity; // total bytes that can be held
  gint32 content;  // total bytes currently held
  gint32 pos;      // subscript of buff currently under consideration
} guint32buff;

typedef struct ByteBuff {
  guint8 *buff;
  gint32 capacity; // total bytes that can be held
  gint32 content;  // total bytes currently held
  gint32 pos;      // subscript of buff currently under consideration
} bytebuff;

typedef struct CharBuff {
  char *buff;
  gint32 capacity; // total bytes that can be held
  gint32 content;  // total bytes currently held
  gint32 pos;      // subscript of buff currently under consideration
} charbuff;

void initByteBuff(bytebuff * buff, guint32 capacity);
void initCharBuff(charbuff * buff, guint32 capacity);
void initGuint32Buff(guint32buff * buff, guint32 capacity);

void freeByteBuff(bytebuff * buff);
void freeCharBuff(charbuff * buff);
void freeGuint32Buff(guint32buff * buff);

int bufferedRead(bytebuff *readbuff, FILE* input);
int bufferedWrite(const bytebuff * const writebuff, FILE* output);


#endif  // KJK_BUFFERS_H

