#!/usr/bin/env ruby
=begin
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
=end

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
  $stderr.puts "splitFastq - tool distributed with Pash version 3.01.03"
  $stderr.puts "Usage:"
  $stderr.puts "pash3_splitFastq.rb <inputFastqFile> <numberOfSequencesPerSplitFile> <outputFileRootName>"
  $stderr.puts "Example:"
  $stderr.puts "pash3_splitFastq.rb test.fastq 1000000 test.split.1m"
  $stderr.puts ""
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

