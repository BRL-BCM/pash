#!/usr/bin/env ruby
# Author: Cristian Coarfa (coarfa@bcm.edu)
# This utility creates an index (offset + sequence files) for a FASTA file OR a list of fasta files



class FastaQualToFastqConverter
  DEBUG=false
	def initialize (optsHash)
		@optsHash = optsHash
		setParameters
	end
	
	def setParameters()
		@fastaFile    = @optsHash['--fastaFile']
		@qualityFile  = @optsHash['--qualityFile']
		@injectQualityScore = @optsHash['--injectQualityScore'].to_i
		@fastqFile    = @optsHash['--fastqFile']
		$stderr.puts "convert sequence #{@fastaFile} and quality scores  #{@qualityFile} or #{@injectQualityScore} into fastq file #{@fastqFile}" if (DEBUG)
	end
	
	def FastaQualToFastqConverter.processArguments()
		optsArray =	[	['--qualityFile', 				'-q', GetoptLong::OPTIONAL_ARGUMENT],
								  ['--injectQualityScore', 	'-i', GetoptLong::OPTIONAL_ARGUMENT],
									['--fastaFile', 					'-f', GetoptLong::REQUIRED_ARGUMENT],
									['--fastqFile', 					'-o', GetoptLong::REQUIRED_ARGUMENT],
									['--help', 								'-h', GetoptLong::NO_ARGUMENT]
								]
		optsHash={}
		progOpts = GetoptLong.new(*optsArray)
		progOpts.each do |opt, arg|
      case opt
        when '--help'
          FastaQualToFastqConverter.usage()
        when '--fastaFile'
          optsHash['--fastaFile']=arg
        when '--qualityFile'
          optsHash['--qualityFile']=arg
        when '--fastqFile'
          optsHash['--fastqFile']=arg
				when '--injectQualityScore'
					optsHash['--injectQualityScore']=arg
      end
    end

		FastaQualToFastqConverter.usage() if(optsHash.empty?);
		FastaQualToFastqConverter.usage() if(optsHash.key?('--help'));
		if (!optsHash.key?('--fastqFile') || !optsHash.key?('--fastaFile')) then
      FastaQualToFastqConverter.usage("Some argument are missing")
      return nil
		end
		if (!optsHash.key?('--qualityFile') && !optsHash.key?('--injectQualityScore')) then
      FastaQualToFastqConverter.usage("You need to provide either a quality scores file or a quality score")
			return nil
		end
		if (optsHash.key?('--qualityFile') && optsHash.key?('--injectQualityScore')) then
      FastaQualToFastqConverter.usage("You need to provide either a quality scores file or a quality score")
      return nil
		end

				
		return optsHash
	end
	
	def myFileOpen(fileName)
    $stderr.puts "Trying to open #{fileName}" if (DEBUG)
    reader = nil
    if (fileName=~/\.gz$/) then
      reader = File.popen("zcat #{fileName}", "r")
    elsif (fileName=~/\.bz2$/) then
      reader = File.popen("bzcat #{fileName}", "r")
    elsif
      reader = File.open(fileName)
    end
    return reader
	end
	
	def loadSequences()
	  line = nil
		crtDef = ""
		crtSeq=nil
		ff = nil
		reader = myFileOpen(@fastaFile)
		reader.each { |line|
      line.strip!
			if (line =~/^>/) then
				unless (crtDef == "") 
					$stderr.puts "adding #{crtDef}->#{crtSeq}" if (DEBUG)
          s = @fastqStruct.new(crtDef, crtSeq, nil)
          @fastqHash[crtDef]=s
				end
				line =~/^\s*>\s*([^\s]+)(\s|$)/
				crtDef = $1.to_s
				$stderr.puts "found def #{crtDef} #{$2} #{$3}" if (DEBUG)
				crtSeq = ""
			else
        crtSeq << line
			end
		}
		unless (crtDef == "")
      $stderr.puts "adding #{crtDef}->#{crtSeq}" if (DEBUG)
      s = @fastqStruct.new(crtDef, crtSeq, nil)
      @fastqHash[crtDef]=s
		end
		reader.close()
	end
	
	def assignQualityScore()
		@injectQualityScore += 33
		if (@injectQualityScore<93) then
			@injectQualityScore = 93
		end
		@fastqHash.keys.each {|k|
      s = @fastqHash[k]
      if (!s.seq.nil?) then
				s.qual=@injectQualityScore.chr*s.seq.size
      end
    }
	end

  def loadQualityScores()
    line = nil
    crtDef = ""
    crtSeq=nil
    ff = nil
    reader = myFileOpen(@qualityFile)
    adjustedQual = 0
    reader.each { |line|
      line.strip!
      if (line =~/^>/) then
        unless (crtDef == "") 
          $stderr.puts "adding #{crtDef}->#{crtSeq}" if (DEBUG)
          if (!@fastqHash.key?(crtDef)) then
            s = @fastqStruct.new(crtDef, nil, crtSeq)
            @fastqHash[crtDef]=s    
          else
            @fastqHash[crtDef].qual = crtSeq
          end
        end
        line =~/^\s*>\s*([^\s]+)(\s|$)/
        crtDef = $1.to_s
        $stderr.puts "found def #{crtDef} #{$2} #{$3}" if (DEBUG)
        crtSeq = ""
      else
        f = line.split(/\s+/)
        f.each {|n|
          q = n.to_i
          if (q>93) then
            q = 93
          end
          adjustedQual = q + 33;
          crtSeq<< adjustedQual.chr
        }
      end
    }
    unless (crtDef == "")
      $stderr.puts "adding #{crtDef}->#{crtSeq}" if (DEBUG)
      if (!@fastqHash.key?(crtDef)) then
        s = @fastqStruct.new(crtDef, nil, crtSeq)
        @fastqHash[crtDef]=s    
      else
        @fastqHash[crtDef].qual = crtSeq
      end
    end
    reader.close()
  end
	
	# Traverse the hash of sequences and qualities
	# and output the resulting FASTQ file
	def outputFastqFile()
    writer = File.open(@fastqFile, "w")
    @fastqHash.keys.each {|k|
      s = @fastqHash[k]
      if (s.seq.nil?) then
        $stderr.puts "Tag #{k} does not have fasta sequence"
      elsif(s.qual.nil?) then
        $stderr.puts "Tag #{k} does not have quality scores"
      elsif(s.qual.size != s.seq.size) then
        $stderr.puts "Tag #{k} has sequence of size #{s.seq.size} and quality scores of size #{s.qual.size}"
      else
        writer.puts "@#{s.def}\n#{s.seq}\n+\n#{s.qual}"
      end
    }
    writer.close()
	end
	
	# Drive the conversion from fasta+qual files into fastq
	# The steps are
	#* load the sequences
	#* load the quality scores
	#* traverse the sequence and quality scores and output the fastq files
	#** check for errors in this step
	def convert()
    @fastqStruct = Struct.new("FastQ", :def, :seq, :qual)
    @fastqHash ={}
    loadSequences()
    if (@qualityFile!=nil) then
			loadQualityScores()
    else
			assignQualityScore()
    end
    outputFastqFile()
	end
	
	# Print usage information and exit
	def FastaQualToFastqConverter.usage(msg='')
			unless(msg.empty?)
				puts "\n#{msg}\n"
			end
			puts "
PROGRAM DESCRIPTION:
  This utility converts a pair of FASTA and QUAL files into a FASTQ file.
COMMAND LINE ARGUMENTS:
  --qualityFile        |-q   => [optional] quality scores file
  --fastaFile          |-f   => fasta file
  --fastqFile          |-o   => list of quality scores files
  --injectQualityScore |-i   => [optional] quality scores are not available; user chooses a quality score
  --help               |-h   => [optional flag] Output this usage info and exit

USAGE:
  convertFastaQualToFastQ.rb  -f file.fa -q file.qual -o file.fastq
";
	exit(2);
	end
end

###########################################################################
# MAIN
###########################################################################

# Process command line options
optsHash = FastaQualToFastqConverter.processArguments()
if (optsHash==nil) then
	FastaQualToFastqConverter.usage()
end
# Instantiate read index builder using the program arguments
fastqConverter= FastaQualToFastqConverter.new(optsHash)
# Index
fastqConverter.convert()
exit(0);

