#!/usr/bin/env ruby
require 'brl/util/textFileUtil'
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
