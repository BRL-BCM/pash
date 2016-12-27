include ../../Makefile.include

PASH_INCLUDES=$(GLIB_INCLUDE) -I. -I../pash
COMPILE_FLAGS= $(PASH_INCLUDES) $(COMMON_COMPILE_FLAGS) -DKJK_NO_GSL -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE

PASH=pash2SAM.exe
TARGETS=$(PASH) 
VPATH=../pash

all: $(TARGETS)

Pash2SAM_OBJECTS= SequencePool.o FastQUtil.o BRLGenericUtils.o Pash2SAM.o Pash2SAMConverter.o driver.cc

$(PASH): $(Pash2SAM_OBJECTS)
	$(CC) $(COMPILE_FLAGS) -o $@ $+  $(GLIB_LIB)  $(ZLIB_LIB) -lm  # -L$(HOME)/pkgs/lib/ -lprofiler  # -ltcmalloc
	mkdir -p ../../bin
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
