include ../../Makefile.include

PASH_INCLUDES=$(GLIB_INCLUDE)  -I../utils -I. -I../pash
COMPILE_FLAGS= $(PASH_INCLUDES) $(COMMON_COMPILE_FLAGS) -DKJK_NO_GSL -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

PASH=pash-3.0lx.exe
TARGETS=$(PASH) 
#reciprocal.exe

all: $(TARGETS)

Pash_OBJECTS=Pash.o FastaUtil.o PashLib.o Mask.o Pattern.o HiveHash.o FixedHashKey.o Collator.o SequencePool.o IgnoreList.o buffers.o FastQUtil.o BRLGenericUtils.o BisulfiteKmerGenerator.o

reciprocal_OBJECTS= reciprocal.o simpleFasta.o

$(PASH): $(Pash_OBJECTS)
	$(CC) $(COMPILE_FLAGS) -o $@ $+  $(GLIB_LIB)  $(ZLIB_LIB) -lm  # -L$(HOME)/pkgs/lib/ -lprofiler  # -ltcmalloc
	mkdir -p ../../bin
	cp -f $@ ../../bin/	

reciprocal.exe: $(reciprocal_OBJECTS)	
	$(CC) $(COMPILE_FLAGS) -o $@ $+ $(GLIB_LIB) $(ZLIB_LIB) -lm
	cp -f $@ ../../bin/	

%.o: %.c
	$(CC) $(COMPILE_FLAGS) -o $@ -c $<

%.o: %.cpp
	$(CC) $(COMPILE_FLAGS) -o $@ -c $<

%.o: %.cc
	$(CC) $(COMPILE_FLAGS) -o $@ -c $<

reciprocal.o: reciprocal.c
	$(CC) $(COMPILE_FLAGS) -DRECIPROCAL_MAIN -o $@ -c $<

clean:
	rm -f *.o $(TARGETS)
