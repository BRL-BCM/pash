include ../Makefile.include

CFLAGS+=$(COMMON_COMPILE_FLAGS) $(GLIB_INCLUDE) -I../utils -I.
CXXFLAGS+=$(COMMON_COMPILE_FLAGS) $(GLIB_INCLUDE) -I../utils -I.

TARGETS=pash3 

all: $(TARGETS)

Pash_OBJECTS=Pash.o FastaUtil.o PashLib.o Mask.o Pattern.o HiveHash.o FixedHashKey.o Collator.o SequencePool.o 
Pash_OBJECTS+=IgnoreList.o buffers.o FastQUtil.o BRLGenericUtils.o BisulfiteKmerGenerator.o


pash3: $(Pash_OBJECTS)
	$(CXX) -o $@ $+ -static-libstdc++ -static-libgcc $(GLIB_LIB) 

clean:
	rm -f *.o $(TARGETS)
