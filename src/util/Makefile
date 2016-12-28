include ../Makefile.include

CFLAGS+=$(COMMON_COMPILE_FLAGS) $(GLIB_INCLUDE) -I. -I../pash
TARGETS=pash3_keyFreq pash3_makeIgnoreList

all: $(TARGETS)
VPATH=../pash

keyFreq_OBJECTS=IgnoreList.o buffers.o FixedHashKey.o keyFreq.o ../pash/Mask.o
makeIgnoreList_OBJECTS= makeIgnoreList.o IgnoreList.o buffers.o

pash3_keyFreq: $(keyFreq_OBJECTS)
	$(CC) -o $@ $+ $(GLIB_LIB)

pash3_makeIgnoreList: $(makeIgnoreList_OBJECTS) 
	$(CC) -o $@ $+  $(GLIB_LIB)

clean:
	rm -f *.o $(TARGETS)
