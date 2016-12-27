#include <stdio.h>
#include <stdlib.h>

#include "Pash2SAMConverter.h"

int main(int argc, char* argv[]) {
  Pash2SAMConverter *pash2SAMConverter  = new Pash2SAMConverter();
  if (pash2SAMConverter->parseParameters(argc,  argv)) {
    return 1;
  }
  pash2SAMConverter->convertPash2SAM();
  return 0;
}


