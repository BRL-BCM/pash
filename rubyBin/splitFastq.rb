#!/usr/bin/env ruby
DEBUG=false 
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

def usage()
  $stderr.puts "splitFastq.rb <inputFastqFile> <numberOfSequencesPerSplitFile> <outputFileRootName>
   Example:
   splitFastq.rb test.fastq 1000000 test.split.1m
"
  exit(2)
end

#MAIN

if (ARGV.size==0 || ARGV.member?("-h") || ARGV.member?("--help")) then
 usage()
end

bound = ARGV[1].to_i
fbase=ARGV[2].to_s
lCount=0
inFile = myFileOpen(ARGV[0])
outFofIndex=0.to_i
outFile=File.open("#{fbase}.#{outFofIndex}.fastq", "w")
inFile.each {|line|
 if (inFile.lineno%4==1) then
   lCount +=1
   if (lCount>bound) then 
     outFile.close
  #   puts "wrote #{lCount} lines in #{fbase}.#{outFofIndex}"
     outFofIndex = outFofIndex+1
      outFile=File.open("#{fbase}.#{outFofIndex}.fastq", "w")
     lCount = 0 
   end 
 end
 outFile.print line 
}

inFile.close
outFile.close

