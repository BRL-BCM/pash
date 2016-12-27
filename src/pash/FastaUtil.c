#include <stdio.h>
#include <string.h>
#include <glib.h>
#include <stdlib.h>
#include "FastaUtil.h"
#include "BRLGenericUtils.h"

/** Check a flag and execute a block of code.*/
#define xDEBUG(flag, code) if (flag) {code; fflush(stdout); fflush(stderr);}
/** Check if a variable is NULL; if yes, execute a block of code and exit.*/
#define xDieIfNULL(var, code) if (var==NULL) {code; fflush(stdout); fflush(stderr); exit(1);}

#define DEB_FIRST_PASS 0
#define DEB_FIRST_PASS_CHAR 0
#define DEB_RAW_READ 0
#define DEB_NEXT_CHUNK 0
#define DEB_DRIVER 0
#define DEB_KEEP_PARTIAL_BUFFER 0
#define DEB_PROGRESS_HORIZ 1

int parseFastaFileFirstPassFastaUtil(FastaUtil* fastaUtil) {
	int done;
  int i,j,k;
  ParserStates status;
  int readChars, bufferRead;
  int positionInRawBuffer;
  int positionInSequenceBuffer;
  int positionInDeflineBuffer;
  int isEof;
  char *rawBuffer, *sequenceBuffer, *deflineBuffer;
  SequenceInfo* sequencesInfo;
  char currentChar;
  int haveSequence;
  int currentSequenceLength;
  int nextSequenceIndex;
  int numberOfSequences, numberOfAllocatedSequences;
  int parsingDone;

  guint32 currentLine ;
  rawBuffer = fastaUtil->rawBuffer;
  sequenceBuffer = fastaUtil->sequenceBuffer;
  deflineBuffer = fastaUtil->deflineBuffer;
  sequencesInfo = fastaUtil->sequencesInformation;
  numberOfSequences = fastaUtil->numberOfSequences;
  numberOfAllocatedSequences = fastaUtil->numberOfAllocatedSequences;
  xDEBUG(DEB_FIRST_PASS,
        fprintf(stderr, "performing first pass for file >>%s<<\n",
          fastaUtil->fileArray[fastaUtil->numFiles-1]));
  //fastaUtil->currentFile =  fopen(fastaUtil->fileArray[fastaUtil->numFiles-1], "rt");
  fastaUtil->currentFile =  BRLGenericUtils::openTextGzipBzipFile(fastaUtil->fileArray[fastaUtil->numFiles-1]);

  perror("opening fasta file");
  xDieIfNULL(fastaUtil->currentFile,
           fprintf(stderr, "could not open fasta file >>%s<<\n",
                   fastaUtil->fileArray[fastaUtil->numFiles-1]));
  parsingDone = 0;
  bufferRead = RAW_BUFFER_SIZE;
  readChars = fread(rawBuffer, sizeof(char),
                  bufferRead, fastaUtil->currentFile);
  xDEBUG(DEB_RAW_READ, fprintf(stderr, "read %d chars , requested %d\n", readChars, bufferRead));
  positionInRawBuffer = 0;
  status = DetermineLineType;
  if (readChars<bufferRead) {
      isEof = 1;
  }
  if (readChars ==0) {
      parsingDone = 1;
  }
  fastaUtil->currentLine = 0;
  currentSequenceLength = 0;
  haveSequence = 0;
  currentLine = 0;

  while (!parsingDone) {
    switch(status) {
      case DetermineLineType:
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "DetermineLineType at %d\n", currentLine));
        while(!parsingDone) {
          if (positionInRawBuffer == readChars) { // buffer full, read again
            readChars = fread(rawBuffer, sizeof(char), bufferRead, fastaUtil->currentFile);
            xDEBUG(DEB_FIRST_PASS,
              fprintf(stderr, "refresh raw buffer, read %d out of %d requested\n", readChars, bufferRead));
            positionInRawBuffer = 0;
            if (readChars<bufferRead) {
                isEof = 1;
            }
            if (readChars ==0) {
                parsingDone = 1;
                continue;
            }
          }
          currentChar = rawBuffer[positionInRawBuffer];
          xDEBUG(DEB_FIRST_PASS_CHAR, fprintf(stderr, "[%c]\n", currentChar));
          if (currentChar!=' ' && currentChar!='\t') {
              break;
          } else {
              positionInRawBuffer ++;
          }
        }
        if (parsingDone) {
            continue;
        }
        switch(currentChar) {
            case ';':
                status = Comment;
                break;
            case '>':
                status = DefLine;
                break;
            case '\n':
                currentLine ++;
                positionInRawBuffer++;
                status = DetermineLineType;
                break;
            default:
                status = Sequence;
                break;
        }
        continue;
        break;
      case DefLine:
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "DefLine at %d\n", currentLine));
        // update info for previous sequence
        if (haveSequence) {
          if (currentSequenceLength == 0) {
            fprintf(stderr, "zero -sized sequence at line %d\n", currentLine-1);
            exit(1);
          } else {
            fastaUtil->sequencesInformation[numberOfSequences].sequenceLength = currentSequenceLength;
            xDEBUG(DEB_PROGRESS_HORIZ, fprintf(stderr, "++ sequence %d has size %d\n", numberOfSequences, currentSequenceLength));
          }
        }
        // increment sequence info storage, if necessary
        numberOfSequences++;
        if (numberOfSequences > numberOfAllocatedSequences) {
          numberOfAllocatedSequences = numberOfAllocatedSequences+numberOfAllocatedSequences/5;
          fastaUtil->sequencesInformation = (SequenceInfo*)
            realloc(fastaUtil->sequencesInformation, numberOfAllocatedSequences*sizeof(SequenceInfo));
          xDieIfNULL(fastaUtil->sequencesInformation,
                     fprintf(stderr, "could not reallocate sequences information array\n"));
        }
        // search for new line
        haveSequence =1;
        while(!parsingDone) {
          if (positionInRawBuffer == readChars) {  // buffer full, read again
            readChars = fread(rawBuffer, sizeof(char),  bufferRead, fastaUtil->currentFile);
            xDEBUG(DEB_FIRST_PASS,
              fprintf(stderr, "refresh raw buffer, read %d out of %d requested\n", readChars, bufferRead));
            positionInRawBuffer = 0;
            if (readChars<bufferRead) {
                isEof = 1;
            }
            if (readChars==0) {
                parsingDone = 1;
                continue;
            }
          }
          currentChar = rawBuffer[positionInRawBuffer];
          xDEBUG(DEB_FIRST_PASS_CHAR, fprintf(stderr, "[%c]\n", currentChar));
          if (currentChar=='\n') {
              break;
          } else {
              positionInRawBuffer ++;
          }
        }
        if (parsingDone) {
            continue;
        } else if (currentChar!='\n') {
          fprintf(stderr, "empty fasta sequence at line %d\n", currentLine);
          exit(1);
        } else {
          currentLine ++;
          currentSequenceLength = 0;
          positionInRawBuffer ++;
          status = DetermineLineType;
        }
        break;
      case Sequence:
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "Sequence at line %d\n", currentLine));
        // no blanks on sequence lines
        while(!parsingDone) {
          if (positionInRawBuffer == readChars) {  // buffer full, read again
            readChars = fread(rawBuffer, sizeof(char),  bufferRead, fastaUtil->currentFile);
            xDEBUG(DEB_FIRST_PASS,
              fprintf(stderr, "refresh raw buffer, read %d out of %d requested\n", readChars, bufferRead));
            positionInRawBuffer = 0;
            if (readChars<bufferRead) {
                isEof = 1;
            }
            if (readChars==0) {
                parsingDone = 1;
                continue;
            }
          }
          currentChar = rawBuffer[positionInRawBuffer];
          xDEBUG(DEB_FIRST_PASS_CHAR, fprintf(stderr, "[%c]\n", currentChar));
          if (currentChar=='\n') {
            positionInRawBuffer ++;
            currentLine ++;
            status = DetermineLineType;
            break;
          } else {
            currentSequenceLength ++;
            positionInRawBuffer ++;
          }
        }
        break;
      case Comment:
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "Comment at %d\n", currentLine));
        while(!parsingDone) {
          if (positionInRawBuffer == readChars) {  // buffer full, read again
            readChars = fread(rawBuffer, sizeof(char),  bufferRead, fastaUtil->currentFile);
            xDEBUG(DEB_FIRST_PASS,
              fprintf(stderr, "refresh raw buffer, read %d out of %d requested\n", readChars, bufferRead));
            positionInRawBuffer = 0;
            if (readChars<bufferRead) {
                isEof = 1;
            }
            if (readChars==0) {
                parsingDone = 1;
                continue;
            }
          }
          currentChar = rawBuffer[positionInRawBuffer];
          xDEBUG(DEB_FIRST_PASS_CHAR, fprintf(stderr, "[%c]\n", currentChar));
          if (currentChar=='\n') {
              status = DetermineLineType;
              currentLine++;
              positionInRawBuffer++;
              break;
          } else {
              positionInRawBuffer ++;
          }
        }
        break;
      case Unknown:
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "Unknown at %d\n"));
      default:
          xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "Default\n"));
          fprintf(stderr, "unknown parser type\n");
          exit(1);
          break;
    }
  }
  if (haveSequence) {
      if (currentSequenceLength==0) {
          fprintf(stderr, "0-size sequence at line %d\n", currentLine);
          exit(1);
      } else {
          fastaUtil->sequencesInformation[numberOfSequences].sequenceLength = currentSequenceLength;
          xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "-- sequence %d has size %d\n",
                                         numberOfSequences, currentSequenceLength));
      }
  }
  fastaUtil->numberOfAllocatedSequences = numberOfAllocatedSequences;
  fastaUtil->numberOfSequences = numberOfSequences;
  fastaUtil->currentSequenceIndex = numberOfSequences-1;
  fastaUtil->currentLine = currentLine;
  fclose(fastaUtil->currentFile);
}

int addFileFastaUtil(char *fileName, FastaUtil* fastaUtil) {
  int numFiles;
	numFiles = fastaUtil->numFiles +1;
  fastaUtil->fileArray = (char**)realloc(fastaUtil->fileArray, numFiles*sizeof(char));
  xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "adding FASTA file %s\n", fileName));
	if (fastaUtil->fileArray==NULL) {
		fprintf(stderr, "couldn't realloc fileArray for file %s\n", fileName);
		exit(1);
	}
  fastaUtil->fileArray[numFiles-1]=(char*) malloc(sizeof(char)*(strlen(fileName)+1));
	strcpy(fastaUtil->fileArray[numFiles-1], fileName);
  fastaUtil->numFiles = numFiles;
	parseFastaFileFirstPassFastaUtil(fastaUtil);
}

/** Initializes a  FastaUtil from a FASTA/FOF (list of FASTA) file.
 * @param fastaUtil fasta utility object
 * @return new FastaUtil object
 */
FastaUtil* initFastaUtil(char *fileName) {
  FastaUtil* fastaUtil;
  FILE *fofFile;
	char currentLine[MAX_FILE_NAME+2];
  char currentFastaFile[MAX_FILE_NAME+1];
	fastaUtil = (FastaUtil*) malloc(sizeof(FastaUtil));
	fastaUtil->seqPool = new SequencePool();
	fastaUtil->fileArray = NULL;
	fastaUtil->numFiles = 0;
	fastaUtil->currentFile = NULL;
  fastaUtil->bufferRead = RAW_BUFFER_SIZE;
	fastaUtil->deflineBuffer[0]='\0';
  fastaUtil->maximDeflineSize = MAX_DEFNAME_SIZE;
  fastaUtil->maximSequenceSize = SEQUENCE_BUFFER_SIZE;
	fastaUtil->sequencesInformation =
    (SequenceInfo*) malloc(sizeof(SequenceInfo)*INIT_SEQUENCE_INFO);
  xDieIfNULL(fastaUtil->sequencesInformation,
             fprintf(stderr, "could not allocate sequence information array\n"));
  fastaUtil->numberOfAllocatedSequences = INIT_SEQUENCE_INFO;
  fastaUtil->numberOfSequences = 0;
	if (fastaUtil == NULL){
		fprintf(stderr, "could not allocate memory\n");
		exit(1);
	} else {
		if (strstr(fileName+strlen(fileName)-4, ".fof")!=NULL) {
			// have fof file
			fofFile = BRLGenericUtils::openTextGzipBzipFile(fileName);
			while(fgets(currentLine, MAX_FILE_NAME+2, fofFile)!=NULL) {
        sscanf(currentLine, "%s", currentFastaFile);
				if (strlen(currentFastaFile)>MAX_FILE_NAME) {
					fprintf(stderr, "file name longer than %d at line %d in %s\n",
							MAX_FILE_NAME, fastaUtil->numFiles+1, fileName);
					exit(1);
				}
        xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "adding current file: %s\n", currentLine));
				addFileFastaUtil(currentFastaFile, fastaUtil);
			}
			fclose(fofFile);
		} else {
			// have single file
      xDEBUG(DEB_FIRST_PASS, fprintf(stderr, "adding a single file: %s\n", fileName));
			addFileFastaUtil(fileName, fastaUtil);
		}
	}
  return fastaUtil;
}


void fillRawBuffer(FastaUtil* fastaUtil) {
  fastaUtil->readChars = fread(fastaUtil->rawBuffer, sizeof(char),
                               fastaUtil->bufferRead,
                               fastaUtil->currentFile);
  xDEBUG(DEB_RAW_READ, fprintf(stderr, "read %d chars , requested %d\n",
                               fastaUtil->readChars, fastaUtil->bufferRead));
  fastaUtil->positionInRawBuffer = 0;
  if (fastaUtil->readChars<fastaUtil->bufferRead) {
    if (!fastaUtil->isEof) {
      xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "[fill buffer] setting EOF\n"));
      fastaUtil->isEof = 1;
    } else {
      // assertion would be more appropriate
      if (fastaUtil->readChars ==0 && fastaUtil->isEof) {
        // attempt to go to next file
        fclose(fastaUtil->currentFile);
        fastaUtil->currentFileIndex++;
        if (fastaUtil->currentFileIndex == fastaUtil->numFiles) {
          fastaUtil->parsingDone = 1;
          fastaUtil->currentFileIndex--;
        } else {
          fastaUtil->currentFile = BRLGenericUtils::openTextGzipBzipFile(fastaUtil->fileArray[fastaUtil->currentFileIndex]);
          xDieIfNULL(fastaUtil->currentFile,
            fprintf(stderr, "could not open file %s\n", fastaUtil->fileArray[fastaUtil->currentFileIndex]));
          fastaUtil->positionInRawBuffer = 0;
          fastaUtil->currentSequenceIndex = 0;
          fastaUtil->currentDeflineBufferPos = 0;
          fastaUtil->currentSequenceBufferPos = 0;
          fastaUtil->currentLine = 0;
          fastaUtil->isEof = 0;
          fillRawBuffer(fastaUtil);
        }
      } else {
        fprintf(stderr, "wrong assertion on eof\n");
        exit(1);
      }
    }
  }
}

/** Rewinds a FastaUtil to the first file.
 * @param fastaUtil fasta utility object
 */
int rewindFastaUtil(FastaUtil* fastaUtil) {
	xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "rewinding fasta util..."));
	fastaUtil->currentFileIndex = 0;

	fastaUtil->currentFile = BRLGenericUtils::openTextGzipBzipFile(fastaUtil->fileArray[fastaUtil->currentFileIndex]);
	xDieIfNULL(fastaUtil->currentFile,
		fprintf(stderr, "could not open file %s\n", fastaUtil->fileArray[fastaUtil->currentFileIndex]));
	fastaUtil->positionInRawBuffer= 0;
	fastaUtil->currentSequenceIndex = 0;
	fastaUtil->currentDeflineBufferPos = 0;
	fastaUtil->currentSequenceBufferPos = 0;
  fastaUtil->currentSequenceIndex = 0;
  fastaUtil->isEof = 0;
  fastaUtil->eofConsumed = 0;
  fastaUtil->parsingDone = 0;
  fastaUtil->currentLine= 0;
  fastaUtil->parserStatus = DetermineLineType;
  fillRawBuffer(fastaUtil);
	xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "done \n"));
  return 1;
}

inline char getNextCharacter(FastaUtil* fastaUtil) {
  char result;
	int pos = fastaUtil->positionInRawBuffer;
  if (pos<fastaUtil->readChars) {
    fastaUtil->positionInRawBuffer++;
    xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "[%d]\t[%c]\t[%x]\n", fastaUtil->positionInRawBuffer-1,
                                   fastaUtil->rawBuffer[pos], fastaUtil->rawBuffer[pos]));
		result =  fastaUtil->rawBuffer[pos];
		return result;
  } else {
    xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "consumed current buffer\n"));
    if (fastaUtil->isEof && !fastaUtil->eofConsumed) {
      xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "[%d]\t[EOF]\n", fastaUtil->positionInRawBuffer));
      fastaUtil->eofConsumed = 1;
      return -1;
    } else {
      xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "attempt to refill\n"));
      fillRawBuffer(fastaUtil);
      if (fastaUtil->readChars>0) {
        pos = 0;
        fastaUtil->positionInRawBuffer = 1;
        xDEBUG(DEB_NEXT_CHUNK,
              fprintf(stderr, "buffer refill successful\n[%d]\t[%c]\t[%x]\n",
                fastaUtil->positionInRawBuffer-1, fastaUtil->rawBuffer[pos], fastaUtil->rawBuffer[pos]));
				result = fastaUtil->rawBuffer[pos];
				return result;
      } else {
        if (fastaUtil->parsingDone) {
          xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "returning EndOfParsing\n"));
          return -2;
        } else if (fastaUtil->isEof) {
          xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "returning EOF 2nd time\n"));
          return -1;
        }
      }
    }
  }
}

inline void putCharacterBackInRawBuffer(FastaUtil* fastaUtil) {
  if (fastaUtil->positionInRawBuffer==0) {
    fprintf(stderr, "incorrect put char back operation\n");
    exit(1);
  } else {
    fastaUtil->positionInRawBuffer --;
  }
}

/** Reads EITHER the next read from the fasta file OR the
 * next chunk of sequence of size up to SEQUENCE_BUFFER_SIZE.
 * @param fastaUtil fasta utility structure
 * @return the number of DNA nucleotideds read, or -1 in case of error
 */
int nextChunkFastaUtil(FastaUtil* fastaUtil) {
  char *rawBuffer, *sequenceBuffer, *deflineBuffer;
  SequenceInfo* sequencesInfo;
  char currentChar;
  int haveSequence;
  char* sequenceName;

  rawBuffer = fastaUtil->rawBuffer;
  sequenceBuffer = fastaUtil->sequenceBuffer;
  deflineBuffer = fastaUtil->deflineBuffer;
  sequencesInfo = fastaUtil->sequencesInformation;

  haveSequence = fastaUtil->currentSequenceBufferPos>0;
	xDEBUG(DEB_NEXT_CHUNK,
        fprintf(stderr, "getting next chunk from file %s; haveSequence=%d ->%d\n",
          fastaUtil->fileArray[fastaUtil->currentFileIndex], haveSequence,
					fastaUtil->currentSequenceBufferPos));
  currentChar = 0;
  //while (!fastaUtil->isEof && !fastaUtil->parsingDone) {
  while (currentChar>=0) {
    switch(fastaUtil->parserStatus) {
      case DetermineLineType:
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "DetermineLineType at %d\n", fastaUtil->currentLine));
        while((currentChar = getNextCharacter(fastaUtil))>=0) {
          if (currentChar!=' ' && currentChar!='\t') {
              break;
          }
        }
        if (currentChar<0) {
            continue;
        } else {
          switch(currentChar) {
              case ';':
                  fastaUtil->parserStatus = Comment;
                  break;
              case '>':
                  fastaUtil->parserStatus = DefLine;
                  break;
              case '\n':
                  fastaUtil->currentLine ++;
                  fastaUtil->parserStatus = DetermineLineType;
                  break;
              case ' ':
                  fprintf(stderr, "shouldn't get there\n");
                  exit(1);
                  break;
              default:
                  fastaUtil->parserStatus = Sequence;
                  putCharacterBackInRawBuffer(fastaUtil);
                  break;
          }
        }
        break;
      case DefLine:
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "DefLine at %d\n", fastaUtil->currentLine));

        if (haveSequence) {
          // return with previous sequence information
					return 1;
        }
        // search for first non-blank
        while((currentChar = getNextCharacter(fastaUtil)) >=0) {
          if (currentChar!=' ' && currentChar !='\t' && currentChar !='\n') {
            break;
          }
        }
        if (currentChar == '\n') {
            fprintf(stderr, "incorrect FASTA defline in file %s at line %d\n",
                    fastaUtil->fileArray[fastaUtil->currentFileIndex],
                    fastaUtil->currentLine);
            exit(1);
          }
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "found first defname char at rawbufferpos %d: %c\n",
                                       fastaUtil->positionInRawBuffer, currentChar));
        // start populating the defline
        fastaUtil->currentSequenceIndex ++;
        fastaUtil->currentDeflineBufferPos = 0;
        fastaUtil->deflineBuffer[fastaUtil->currentDeflineBufferPos]=currentChar;
        fastaUtil->currentDeflineBufferPos ++;
        // search for space new line
        while( (currentChar = getNextCharacter(fastaUtil))>=0) {
          if (currentChar=='\n' || currentChar==' ' || currentChar=='\t') {
              break;
          } else {
            if (fastaUtil->currentDeflineBufferPos>=fastaUtil->maximDeflineSize) {
              fprintf(stderr, "long-named sequence at line %d\n", fastaUtil->currentLine);
              exit(1);
            }
            fastaUtil->deflineBuffer[fastaUtil->currentDeflineBufferPos]=currentChar;
            fastaUtil->currentDeflineBufferPos ++;
          }
        }
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "[Defline] ended search for whitespace/newline at rawbufferpos %d: %c\n",
                                       fastaUtil->positionInRawBuffer, currentChar));

        if (currentChar<0) {
          fprintf(stderr, "empty fasta sequence at line %d\n", fastaUtil->currentLine);
          exit(1);
        }
        if (currentChar!='\n') {
          xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "[Defline] searching for newline at rawbufferpos %d: %c\n",
                                       fastaUtil->positionInRawBuffer, currentChar));
          while( (currentChar = getNextCharacter(fastaUtil))>=0) {
            if (currentChar=='\n') {
                break;
            }
          }
          if (currentChar<0) {
            fprintf(stderr, "empty fasta sequence at line %d\n", fastaUtil->currentLine);
            exit(1);
          }
        }
        fastaUtil->deflineBuffer[fastaUtil->currentDeflineBufferPos]='\0';
        // replace by using a more efficient glib function
        // or a list of large character arrays containing sequence names
        sequenceName = (char*) malloc(sizeof(char)*(fastaUtil->currentDeflineBufferPos+1));
        strcpy(sequenceName, fastaUtil->deflineBuffer);
        fastaUtil->sequencesInformation[fastaUtil->currentSequenceIndex].sequenceName = sequenceName;
        fastaUtil->currentLine ++;
        fastaUtil->currentSequenceBufferPos = 0;
				fastaUtil->currentActualSequencePos = 0;
        fastaUtil->parserStatus = DetermineLineType;
        xDEBUG(DEB_NEXT_CHUNK,
               fprintf(stderr, "found definition [%d] %s %s %s\n",
                fastaUtil->currentSequenceIndex, fastaUtil->deflineBuffer,
                sequenceName, fastaUtil->sequencesInformation[fastaUtil->currentSequenceIndex].sequenceName));
        break;
      case Sequence:
				if (fastaUtil->currentActualSequencePos < fastaUtil->sequencesInformation[fastaUtil->currentSequenceIndex].sequenceLength) {
					haveSequence = 1;
				}
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "Sequence at line %d\n", fastaUtil->currentLine));
        // no blanks on sequence lines
        while( (currentChar = getNextCharacter(fastaUtil))>=0) {
          if (currentChar=='\n') {
              break;
          } else {
						xDEBUG(DEB_NEXT_CHUNK,
               fprintf(stderr, "added char %d %c\n", currentChar, currentChar));
						if (currentChar>='a' && currentChar<='z') {
							currentChar = currentChar-'a';
							currentChar = currentChar+'A';
						}
						fastaUtil->sequenceBuffer[fastaUtil->currentSequenceBufferPos]=currentChar;
            fastaUtil->currentSequenceBufferPos ++;
						fastaUtil->currentActualSequencePos ++;
            // check for sequence buffer max capacity
            if (fastaUtil->currentSequenceBufferPos >= fastaUtil->maximSequenceSize) {
              return 1;
            }
          }
        }
        if (currentChar < 0) {
            continue;
        }
        if (currentChar=='\n') {
          fastaUtil->currentLine ++;
          fastaUtil->parserStatus = DetermineLineType;
        }
        break;
      case Comment:
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "Comment at %d\n", fastaUtil->currentLine));
        while( (currentChar = getNextCharacter(fastaUtil))>=0) {
          if (currentChar=='\n') {
              break;
          }
        }
        if (currentChar<0) {
            continue;
        }
        if (currentChar=='\n') {
          fastaUtil->currentLine ++;
          fastaUtil->parserStatus = DetermineLineType;
        }
        break;
      case Unknown:
        xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "Unknown at %d\n"));
      default:
          xDEBUG(DEB_NEXT_CHUNK, fprintf(stderr, "Default\n"));
          fprintf(stderr, "unknown parser type\n");
          exit(1);
          break;
    }
  }

	xDEBUG(DEB_NEXT_CHUNK,
				 fprintf(stderr, "parsing loop stopped, current char %d, EndOfParsing %d\n",
								 currentChar, fastaUtil->parsingDone));
  if (haveSequence) {
      if (fastaUtil->currentSequenceBufferPos==0) {
          fprintf(stderr, "%s:%d 0-size sequence at line %d\n",
									__FILE__, __LINE__, fastaUtil->currentLine);
          exit(1);
      } else {
        fastaUtil->parserStatus = DetermineLineType;
        return 1;
      }
  }
  return 0;
}

void fastaUtilKeepPartialBuffer(FastaUtil* fastaUtil, int sequenceBasesToKeep) {
  xDEBUG(DEB_KEEP_PARTIAL_BUFFER, fprintf(stderr, "keeping partial buffer of %d bases\n",
                                          sequenceBasesToKeep));
  if (sequenceBasesToKeep==0) {
    fastaUtil->currentSequenceBufferPos = 0;
  } else if (sequenceBasesToKeep>fastaUtil->currentSequenceBufferPos) {
    fprintf(stderr, "asked to keep %d characters from a buffer with only %d characters\n",
            sequenceBasesToKeep, fastaUtil->currentSequenceBufferPos);
  } else {
    int destIndex, targetIndex;
    for (destIndex = 0, targetIndex = fastaUtil->currentSequenceBufferPos - sequenceBasesToKeep;
         targetIndex<fastaUtil->currentSequenceBufferPos; destIndex++, targetIndex++) {
      fastaUtil->sequenceBuffer[destIndex] = fastaUtil->sequenceBuffer[targetIndex];
    }
    fastaUtil->currentSequenceBufferPos = destIndex;
    xDEBUG(DEB_KEEP_PARTIAL_BUFFER, fprintf(stderr, "reset current buffer pos to %d\n",
                                          fastaUtil->currentSequenceBufferPos));
  }
}

#define UNIT_TEST_FASTA_UTIL_1______  1
#undef UNIT_TEST_FASTA_UTIL_1______
#ifdef UNIT_TEST_FASTA_UTIL_1______

int main(int argc, char**argv) {
  int i;
  FastaUtil* fastaUtil;
  if (argc == 1) {
    fprintf(stderr, "testFastaUtil <fasta/fof file>\n");
    exit(0);
  }
  fprintf(stderr, "performing first pass for file %s\n", argv[1]);
  fastaUtil = initFastaUtil(argv[1]);
	rewindFastaUtil(fastaUtil);
	while(!fastaUtil->parsingDone) {
		nextChunkFastaUtil(fastaUtil);
    xDEBUG(DEB_DRIVER, {
			fastaUtil->deflineBuffer[fastaUtil->currentDeflineBufferPos]='\0';
			fastaUtil->sequenceBuffer[fastaUtil->currentSequenceBufferPos]='\0';
			fprintf(stderr, "current fasta util state\n");
					 fprintf(stderr, "current fasta file [%d]->%s\n",
									 fastaUtil->currentFileIndex,
									 fastaUtil->fileArray[fastaUtil->currentFileIndex]);
					 fprintf(stderr, "current line %d\n", fastaUtil->currentLine);
					 fprintf(stderr, "def: %s \n seq: %s\n",
									 fastaUtil->deflineBuffer, fastaUtil->sequenceBuffer);
		}	);
		fastaUtilKeepPartialBuffer(fastaUtil, 0);
	}
  return 0;
}
#endif
