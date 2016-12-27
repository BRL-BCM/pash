#!/usr/bin/env ruby

# require "socket"

class SubmitHitMapper
  DEBUG = false
	def initialize (optsHash)
		@optsHash = optsHash
		setParameters
	end

	def setParameters()
    @targetGenome = @optsHash['--targetGenome']
    @query      = @optsHash['--query']
    if (@optsHash.key?('--topPercent')) then
	    @topPercent = @optsHash['--topPercent']
    else
			@topPercent = 1
    end
    if (@optsHash.key?('--maxMappings')) then
	    @maxMappings= @optsHash['--maxMappings']
    else
			@maxMappings= 1
    end
    @outputDirectory = @optsHash['--outputDirectory']
    @ignoreList = @optsHash['--ignoreList']
    if (@optsHash.key?('--diagonals')) then
			@diagonals = @optsHash['--diagonals']
    else
			@diagonals = 500
    end
    if (@optsHash.key?('--gap')) then
			@gap = @optsHash['--gap']
		else
			@gap = 6
    end
		@scratchDir = @optsHash['--scratch']
		@ignorePercent = @optsHash['--ignorePercent']
		@nodes = 1
		if (@optsHash.key?('--nodes')) then
			@nodes = @optsHash['--nodes']
		end
    @kw = 12
    @ks = 18
    if(@optsHash.key?('--kweight') && @optsHash.key?('--kspan')) then
      @kw = @optsHash['--kweight'].to_i
      @ks = @optsHash['--kspan'].to_i
    end
    $stderr.puts @optsHash['--bisulfiteSeq'] if (DEBUG)
		if (@optsHash['--bisulfiteSeq']!=nil) then
			@bisulfiteFlag =" -B "
			@splitSize = 200000
		else
			@bisulfiteFlag =" "
			@splitSize = 1000000
		end
		$stderr.puts "hash >>#{@optsHash['--bisulfiteSeq']}<< flag >>#{@bisulfiteFlag}<<" if (DEBUG)
		#exit(2)
	end

	def submitThis
		dName = File.basename(@outputDirectory)
		system("mkdir -p #{@outputDirectory}")
		sleep 2
		thisHost = Socket.gethostname
		if (thisHost =~ /(woodstock)|(ned)|(lewis)|(clark)|(technetium)/) then
			batchType = "lsf"
			scriptName = "#{@outputDirectory}/pashMap.#{dName}.#{Process.pid}.lsf"
			scriptFile = File.open(scriptName, "w")
			scriptFile.puts "#!/bin/bash";
			scriptFile.puts "#BSUB -L /bin/bash";
			scriptFile.puts "#BSUB -q normal"
			scriptFile.puts "#BSUB -J pash.#{dName}.#{Process.pid}";
			scriptFile.puts "#BSUB -o #{@outputDirectory}/#{dName}.pash.%J.o";
			scriptFile.puts "#BSUB -e #{@outputDirectory}/#{dName}.pash.%J.e";
		elsif (thisHost =~ /(sphere)|(gauss)/) then
			batchType = "pbs"
			scriptName = "#{@outputDirectory}/pashMap.#{dName}.#{Process.pid}.pbs"
			scriptFile = File.open(scriptName, "w")
			scriptFile.puts "#!/bin/bash";
			scriptFile.puts "#PBS -q dque";
			scriptFile.puts "#PBS -l nodes=1:ppn=#{@nodes}\n";
			scriptFile.puts "#PBS -l walltime=48:00:00\n";
			scriptFile.puts "#PBS -l cput=48:00:00\n";
			scriptFile.puts "#PBS -M #{ENV["USER"]}\@bcm.tmc.edu\n";
			scriptFile.puts "#PBS -m ea\n";
			scriptFile.puts "#PBS -N pash3.#{dName}.#{Process.pid}\n";
			scriptFile.puts "#PBS -o #{@outputDirectory}/#{dName}.pash3.#{Process.pid}.o";
			scriptFile.puts "#PBS -e #{@outputDirectory}/#{dName}.pash3.#{Process.pid}.e";
		else
			puts "unsupported host"
			exit(2)
		end

    scriptFile.puts "/bin/hostname"
    scriptFile.puts "ulimit -a"
    #scriptFile.puts "mkdir -p #{@outputDirectory}"

    scriptFile.puts "cd #{@outputDirectory}"
    scriptFile.puts "sleep 1"
    ignorePattern = "#{File.dirname(@targetGenome)}/*#{@kw}*#{@ks}*.#{@ignorePercent}p.*il"
    ignoreList = Dir[ignorePattern][0]
    $stderr.puts "Using pattern #{ignorePattern} ---> ignore list #{ignoreList}"
    if (ignoreList.nil?) then
      ignoreList="" 
    else 
      ignoreList = "-L #{ignoreList}"
    end
    Dir.chdir(@outputDirectory)
    # split the input file in small (200k) pieces
    
    scriptFile.puts "splitFastq.rb #{@query} #{@splitSize} #{File.basename(@query)}.split.#{@splitSize}" 
    scriptFile.print "for f in #{File.basename(@query)}.split.#{@splitSize}.*.fastq; "
    outName = "map.`basename $f`.on.#{File.basename(@targetGenome)}.P#{@topPercent}.n#{@maxMappings}.d#{@diagonals}.G#{@gap}.ignore.#{@ignorePercent}"
		scriptFile.print "do echo $f; time pash-3.0d.exe -v $f -h #{@targetGenome} -k #{@kw} -n #{@ks} -S . -s 22 -G #{@gap} -d #{@diagonals} "
		#scriptFile.print " -P #{@topPercent} -n #{@maxMappings} "
		scriptFile.print " -o #{@scratchDir}/#{outName} -N #{@maxMappings} #{@bisulfiteFlag}"
		scriptFile.print  " #{ignoreList}> log.map.`basename $f`.on.#{File.basename(@targetGenome)}.P#{@topPercent}.n#{@maxMappings} 2>&1 ; "
		scriptFile.puts "sleep 2; mv #{@scratchDir}/#{outName} . ; sleep 2; rm $f; done"
		
		scriptFile.close()
		if (batchType == "lsf") then
			command = "bsub < #{scriptName}"
		else
			command = "qsub #{scriptName}"
		end
		system(command)
	end

	def SubmitHitMapper.processArguments()
		# We want to add all the prop_keys as potential command line options

		optsArray =	[ ['--targetGenome', 	   '-t', GetoptLong::REQUIRED_ARGUMENT],
                  ['--query',            '-q', GetoptLong::REQUIRED_ARGUMENT],
                  ['--outputDirectory',  '-o', GetoptLong::REQUIRED_ARGUMENT],
                  ['--scratch',          '-s', GetoptLong::REQUIRED_ARGUMENT],
                  ['--diagonals',        '-d', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--gap',              '-G', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--topPercent',       '-P', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--maxMappings',      '-n', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--nodes',            '-N', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--ignorePercent',    '-i', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--kweight',          '-k', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--kspan',            '-l', GetoptLong::OPTIONAL_ARGUMENT],
                  ['--bisulfiteSeq',     '-B', GetoptLong::NO_ARGUMENT],
	                ['--help', 				     '-h', GetoptLong::NO_ARGUMENT]
								]
		#PROP_KEYS.each { |propName|
		#	argPropName = "--#{propName}"
		#	optsArray << [argPropName, GetoptLong::OPTIONAL_ARGUMENT]
		#}
		progOpts = GetoptLong.new(*optsArray)
		optsHash = progOpts.to_hash
		SubmitHitMapper.usage() if(optsHash.key?('--help'));

		unless(progOpts.getMissingOptions().empty?)
			SubmitHitMapper.usage("USAGE ERROR: some required arguments are missing")
		end

		SubmitHitMapper.usage() if(optsHash.empty?);
		return optsHash
	end

	def SubmitHitMapper.usage(msg='')
		unless(msg.empty?)
			puts "\n#{msg}\n"
		end
		puts "
PROGRAM DESCRIPTION:
  Submits a cluster job that maps fastq files using Pash 3. It can invoke pash for
  - regular mapping
  - bisulfite sequencing mapping, for reads obtained through bisulfite sequencing treatment


COMMAND LINE ARGUMENTS:
  --targetGenome    | -t   => fasta file containing the reference genome
  --scratch         | -s   => scratch directory
  --query           | -q   => file to map
  --outputDirectory | -o   => output directory where results should
                              be generated
  --topPercent      | -P   => [optional] top percent of mappings to be kept
                              (default 1)
  --maxMappings     | -n   => [optional] maximum number of mappings within
                              top percent of best score  (default 1). Reads
                              with a larger number of mappings that this
                              value are discarded from mapping results
  --diagonals       | -d   => number of diagonals, default 500
  --gap             | -G   => gap, default 6
  --nodes           | -N   => processors per node, default 1
  --kweight         | -k   => kmer weight
  --kspan           | -l   => kmer length
  --bisulfiteSeq    | -B   => perform bisulfite sequencing mapping
  --help            | -h   => [optional flag] Output this usage info and exit

USAGE:
  submitPashJob.rb -t human.build.36.fa -q reads.fa -o `pwd` -k 13 -l 21 -G 2 -d 100 -S . 
";
		exit(2);
	end
end


##########################
# MAIN                   #
##########################

# Process command line options
optsHash = SubmitHitMapper.processArguments()
# Instantiate analyzer using the program arguments
submitHitMapper = SubmitHitMapper.new(optsHash)
# Submit this !
submitHitMapper.submitThis()
exit(0);
