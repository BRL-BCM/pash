#!/usr/bin/env ruby
=begin
Copyright (c) 2004-2016 Baylor College of Medicine.
Use of this software is governed by a license.
See the included file LICENSE for details.
=end


require 'zlib'
require 'delegate'
$VERBOSE = false

# ============================= brl/util/util
module BRL ; module Util
module Gzip

  GZIP_EXE = 'gzip -t '
  GZIP_INFO = 'gzip -l'
  GZIP_INFO_RE = /\s*compressed\s*uncompressed\s*ratio\s*uncomp\S+\s*\r?\n\s*(\d+)\s*(\d+)\s*([0-9\.]+)\%\s*(\S+)/

  def Gzip.isGzippedFile?(fileStr)
    fullFilePath = File.expand_path(fileStr)
    return false if((!fullFilePath.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
    cmdStr1 = "#{BRL::Util::Gzip::GZIP_EXE} #{fullFilePath} > /dev/null 2> /dev/null "
    begin
      system(cmdStr1)
      exitCode = $?
      return (exitCode == 0) ? true : false
    rescue
      raise "\n\nERROR: can't execute sub-processes??? (BRL::Util.Util#Gzip.isGzippedFile?)\n\n"
    end
  end

  def Gzip.getInfo(fileStr)
    fullFilePath = File.expand_path(fileStr)
    return nil if((!fileStr.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
    cmdStr = "#{BRL::Util::Gzip::GZIP_INFO} #{fullFilePath} "
    begin
      cmdOutput = `#{cmdStr}`
    rescue => err
      raise(SystemCallError, "\nERROR: 'gzip -l' doesn't work! do you have it installed and reachable from sh? Is your gzip up to date?\nThis command returned:]\n  #{cmdOutput}\nError details: #{err.message}\n" + err.backtrace.join("\n")) if(cmdOutput.nil? or cmdOutput.empty?)
      raise(SystemCallError, "\nERROR: '#{fileStr}' is not a valid gzipped file!") if(cmdOutput =~ /not in gzip format/)
    end
    if(cmdOutput =~ BRL::Util::Gzip::GZIP_INFO_RE)
      return [$1, $2, $3.to_f/100.0, $4]
    else
      raise(SystemCallError, "\nERROR: 'gzip -l' doesn't work! do you have it installed and reachable from sh? Is your gzip up to date?\nThis command returned:\n>>\n#{cmdOutput}\n<<")
    end
  end

  def Gzip.empty?(fileStr)
    fullFilePath = File.expand_path(fileStr)
    if((!fileStr.kind_of?(String)) or !FileTest.exists?(fullFilePath) or !FileTest.readable?(fullFilePath))
      raise(SystemCallError, "\nERROR: the file isn't the name of a gzippe file or the file doesn't exist or you don't have permission to read it!\n    Bad file: #{fileStr}\n")
    end
    testByte = nil
    begin
      # Try to read just 1 byte from gzip file
      # If fails, that's because the compressed file was empty!
      gzIo = Zlib::GzipReader.open(fullFilePath)
      testByte = gzIo.readchar
    rescue
    end
    return testByte.nil? # yes, compressed file was empty!
  end
  
  def Gzip.getCompressedSize(fileStr)
    return BRL::Util::Gzip::getGzipInfo(fileStr)[0]
  end

  def Gzip.getUncompressedSize(fileStr)
    return BRL::Util::Gzip::getGzipInfo(fileStr)[1]
  end

  def Gzip.getRatio(fileStr)
    return BRL::Util::Gzip::getGzipInfo(fileStr)[2]
  end

  def Gzip.getFileName(fileStr)
    return BRL::Util::Gzip::getGzipInfo(fileStr)[3]
  end
end # module Gzip

# ============================= brl/util/textFileUtil
# Author: Andrew R Jackson (andrewj@bcm.tmc.edu)
# Date: 3/31/2004 4:38PM
# Purpose:
# Convenience classes for reading/writing plain text that *may* be gzipped.
# We let the instance figure out whether the text is gzipped or not when
# reading. Uses gzip itself to test the file (not the 'file' command).

  class TextReader < SimpleDelegator
    include BRL::Util
    def initialize(fileStr)
      unless(fileStr.kind_of?(String))
        raise(TypeError, "The file to read from must be provided as a String.");
      end
      unless(FileTest.exists?(fileStr))
        raise(IOError, "The file '#{fileStr}' doesn't exist.");
      end
      unless(FileTest.readable?(fileStr))
        raise(IOError, "The file '#{fileStr}' isn't readable.");
      end
      if(Gzip.isGzippedFile?(fileStr))
        @ioObj = Zlib::GzipReader.open(fileStr);
      else
        @ioObj = File.open(fileStr, "r");
      end
      super(@ioObj)
    end
  end

  class TextWriter < SimpleDelegator
    include BRL::Util
    def initialize(fileStr, modeStr="w+", doGzipOutput=false)
      unless(fileStr.kind_of?(String))
        raise(TypeError, "The file to write to must be provided as a String.");
      end
      if(FileTest.exists?(fileStr) and !FileTest.writable?(fileStr))
        raise(IOError, "The file '#{fileStr}' exists but isn't writable.");
      end
      file = File.open(fileStr, modeStr)
      if(doGzipOutput)
        @ioObj = Zlib::GzipWriter.new(file)
      else
        @ioObj = file
      end
      super(@ioObj)
    end
  end
end ; end # module BRL ; module Util

# ====================================================

if ARGV.size != 2
  $stderr.puts "getRCChrom - tool distributed with Pash version 3.01.03"
  $stderr.puts "Usage:"
  $stderr.puts "pash3_getRCChrom  <fasta_file_with_reference_genome>  <output_fasta_file>"
  $stderr.puts ""
  exit 1
end

DEBUG=true
r = BRL::Util::TextReader.new(ARGV[0])
w = BRL::Util::TextWriter.new(ARGV[1])
line = nil
crtDef = ""
crtSeq = ""
r.each { |line|
  line.strip!
  if (line =~/^>/) then
    if (crtDef != "")  then
      $stderr.puts "adding #{crtDef}" if (DEBUG)
      w.puts ">#{crtDef}\n#{crtSeq}"
      w.puts ">#RC.pash.#{crtDef}\n#{crtSeq.reverse.tr("AaCcGgTtn", "TTGGCCAAN")}"
    end
      line =~/^\s*>\s*([^\s]+)(\s|$)/
      crtDef = $1.to_s
      $stderr.puts "found def #{crtDef} #{$2} #{$3}" if (DEBUG)
      crtSeq = ""
    else
      crtSeq << line.strip 
    end
}
if (crtDef != "")  then
      $stderr.puts "adding #{crtDef}" if (DEBUG)
      w.puts ">#{crtDef}\n#{crtSeq}"
      w.puts ">#RC.pash.#{crtDef}\n#{crtSeq.reverse.tr("AaCcGgTtn", "TTGGCCAAN")}"
end
r.close()
w.close()
