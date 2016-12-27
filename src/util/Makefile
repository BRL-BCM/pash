include ../../Makefile.include
PASH_INCLUDES=$(GLIB_INCLUDE) -I. -I../pash
COMPILE_FLAGS= $(PASH_INCLUDES) $(COMMON_COMPILE_FLAGS) -DKJK_NO_GSL -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE
TARGETS=keyFreq.exe makeIgnoreList.exe

all: $(TARGETS)
VPATH=../pash

keyFreq_OBJECTS=IgnoreList.o buffers.o FixedHashKey.o keyFreq.o 
makeIgnoreList_OBJECTS= makeIgnoreList.o IgnoreList.o buffers.o

keyFreq.exe: $(keyFreq_OBJECTS)
	$(CC) $(COMPILE_FLAGS) -o $@ $+  $(GLIB_LIB) -lm 
	cp -f $@ ../../bin/	

makeIgnoreList.exe: $(makeIgnoreList_OBJECTS) 
	$(CC) $(COMPILE_FLAGS) -o $@ $+ $(GLIB_LIB)
	cp -f $@ ../../bin/	

%.o: %.c
	$(CC) $(COMPILE_FLAGS) -o $@ -c $< 

clean:
	rm -f *.o $(TARGETS)
