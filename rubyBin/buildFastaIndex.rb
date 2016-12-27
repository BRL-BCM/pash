#!/usr/bin/env ruby
# Author: Cristian Coarfa (coarfa@bcm.edu)
# This utility creates an index (offset + sequence files) for a FASTA file OR a list of fasta files

require 'brl/util/textFileUtil'
require 'brl/util/util'

class ReadIndexBuilder
  DEBUG=false
	def initialize (optsHash)
		@optsHash = optsHash
		setParameters
	end
	
	def setParameters()
		@readsFile = @optsHash['--fastaFile']
		@fofFile = @optsHash['--fofFile']
		if (@readsFile !=nil && @fofFile!=nil) then
			$stderr.puts "You need to specify a reads file OR a fof file, not both"
			exit(2)
		end
		if (@readsFile==nil && @fofFile == nil) then
      	$stderr.puts "You need to specify a reads file OR a fof file"
			exit(2)
		end
		@readsOffsetFile = @optsHash['--offsetFile']
		@readsSequenceFile = @optsHash['--sequenceFile']
    $stderr.puts "write sequence in #{@readsSequenceFile} and offsets in #{@readsOffsetFile}" if (DEBUG)
	end
	
	def ReadIndexBuilder.processArguments()
		optsArray =	[	['--fastaFile', '-r', GetoptLong::OPTIONAL_ARGUMENT],
									['--offsetFile', '-o', GetoptLong::REQUIRED_ARGUMENT],
									['--sequenceFile', '-s', GetoptLong::REQUIRED_ARGUMENT],
									['--fofFile', '-f', GetoptLong::OPTIONAL_ARGUMENT],
									['--help', '-h', GetoptLong::NO_ARGUMENT]
								]
		
		progOpts = GetoptLong.new(*optsArray)
		optsHash = progOpts.to_hash
		ReadIndexBuilder.usage() if(optsHash.key?('--help'));
		
		unless(progOpts.getMissingOptions().empty?)
			ReadIndexBuilder.usage("USAGE ERROR: some required arguments are missing") 
		end
		ReadIndexBuilder.usage() if(optsHash.empty?);
		return optsHash
	end
	
	def indexFastaFile(reader,writer,offsetWriter)
    line = nil
		crtDef = ""
		crtSeqSize = 0
		ff = nil
		$stdout.sync = true
		offsetWriter.sync = true
		reader.each { |line|
      line.strip!
			if (line =~/^>/) then
				unless (crtDef == "") 
					$stderr.puts "adding #{crtDef}->#{crtSeqSize}" if (DEBUG)
					offsetWriter.puts "#{crtDef}\t#{@currentOffset}\t#{crtSeqSize}"
					@currentReadIndex += 1
					@currentOffset += crtSeqSize
				  if (@currentReadIndex %100000 ==0) then
            $stdout.print "."
				  end
				end
				line =~/^\s*>\s*([^\s]+)(\s|$)/
				crtDef = $1.to_s
				$stderr.puts "found def #{crtDef} #{$2} #{$3}" if (DEBUG)
				crtSeqSize = 0
			else
				writer.print "#{line}"
        crtSeqSize+= line.size
			end
		}
		unless (crtDef == "")
      $stderr.puts "adding #{crtDef}->#{crtSeqSize}" if (DEBUG)
			offsetWriter.puts "#{crtDef}\t#{@currentOffset}\t#{crtSeqSize}"
      @currentReadIndex += 1
      @currentOffset += crtSeqSize
      if (@currentReadIndex %100000 ==0) then
        $stdout.print "."
      end
		end
	end
	
	
	def buildIndex
    @currentReadIndex = 0
    @currentOffset = 0
    offsetFile = BRL::Util::TextWriter.new(@readsOffsetFile)
    @offsetStruct = Struct.new("ReadOffsetStruct", :readIndex, :offset, :size)
    writer = BRL::Util::TextWriter.new(@readsSequenceFile,"w")
    if (@readsFile !=nil) then
      fastaFileReader = BRL::Util::TextReader.new(@readsFile)
      $stderr.puts "indexing #{@readsFile}" if (DEBUG)
      indexFastaFile(fastaFileReader, writer, offsetFile)
    else
      fofReader = BRL::Util::TextReader.new(@fofFile)
      fastaFileIndex = 0  
      line = ""
      fofReader.each { |line|
        fastaFileName = line.strip
        $stderr.puts "indexing #{fastaFileName}" if (DEBUG)
        fastaFileReader = BRL::Util::TextReader.new(fastaFileName)
        indexFastaFile(fastaFileReader,writer, offsetFile)
        fastaFileIndex+=1
        fastaFileReader.close()
      }
    end
    writer.close()
    offsetFile.close()
	end
	
	
	def ReadIndexBuilder.usage(msg='')
			unless(msg.empty?)
				puts "\n#{msg}\n"
			end
			puts "
PROGRAM DESCRIPTION:
  This utility creates an index for a fasta file OR a fof list of fasta files.
COMMAND LINE ARGUMENTS:
  --fastaFile        |-r   => FASTA sequence file
  --fofFile          |-f   => list of fasta files
  --offsetFile       |-o   => FASTA index files
  --sequenceFile     |-s   => file containing concatenated nucleotide sequence
  --help             |-h   => [optional flag] Output this usage info and exit

USAGE:
  buildFastaIndex.rb  -f allReads.fof -o offsets.txt -s sequence.txt 
";
			exit(2);
	end
end

###########################################################################
# MAIN
###########################################################################

# Process command line options
optsHash = ReadIndexBuilder.processArguments()
# Instantiate read index builder using the program arguments
fastaIndexer = ReadIndexBuilder.new(optsHash)
# Index
fastaIndexer.buildIndex()
exit(0);

