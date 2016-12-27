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
		@outputFile = @optsHash['--outputFile']
	end
	
	def ReadIndexBuilder.processArguments()
		optsArray =	[	['--fastaFile', '-r', GetoptLong::OPTIONAL_ARGUMENT],
					['--outputFile', '-o', GetoptLong::REQUIRED_ARGUMENT],
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
	
  def indexFastaFile(reader, offsetWriter)
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
					offsetWriter.puts "#{crtDef}\t#{crtSeqSize}"
				end
				line =~/^\s*>\s*([^\s]+)(\s|$)/
				crtDef = $1.to_s
				$stderr.puts "found def #{crtDef} #{$2} #{$3}" if (DEBUG)
				crtSeqSize = 0
			else
                          crtSeqSize+= line.size
			end
		}
		unless (crtDef == "")
                  $stderr.puts "adding #{crtDef}->#{crtSeqSize}" if (DEBUG)
		  offsetWriter.puts "#{crtDef}\t#{crtSeqSize}"
		end
	end
	
	
	def buildIndex
          outputFileWriter = BRL::Util::TextWriter.new(@outputFile)
          fastaFileReader = BRL::Util::TextReader.new(@readsFile)
          $stderr.puts "indexing #{@readsFile}" if (DEBUG)
          indexFastaFile(fastaFileReader, outputFileWriter)
          outputFileWriter.close()
	end
	
	
	def ReadIndexBuilder.usage(msg='')
			unless(msg.empty?)
				puts "\n#{msg}\n"
			end
			puts "
PROGRAM DESCRIPTION:
  This utility extracts the chromosomes names and lengths from a reference fasta file.

COMMAND LINE ARGUMENTS:
  --fastaFile        |-r   => FASTA sequence file
  --outputFile       |-o   => FASTA index files
  --help             |-h   => [optional flag] Output this usage info and exit

USAGE:
  extractChromLengths  -r ref.fa -o ref.chroms.txt 
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

