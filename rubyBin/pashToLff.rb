#!/usr/bin/env ruby


def usage()
  $stderr.puts "pashToLff.rb <pash output> <lff output> <new track> <new class>"
end

if (ARGV.size==0 || ARGV[0]=='-h' || ARGV[0]=='--help') then
  usage()
  exit(2)
end

if (ARGV.size != 4) then
  usage()
  exit(2)
end

pashFile = ARGV[0]
lffFile = ARGV[1]
track=ARGV[2]
track =~ /([^:]+):([^d]+)/
nType = $1
nSubtype = $2
nClass = ARGV[3] 


r = File.open(pashFile, "r")
w = File.open(lffFile, "w")
l = nil
cStart = 0
cStop = 0
cAdjust = 0
chrom =  nil
strand = nil
r.each {|l|
  next if (l !~ /chr/)
  ff = l.strip.split(/\s+/)
  chrom = ff[0]
  cStart = ff[1].to_i
  cStop = ff[2].to_i
  strand = ff[6] 
  w.puts "#{nClass}\t#{ff[3]}\t#{nType}\t#{nSubtype}\t#{chrom}\t#{cStart}\t#{cStop}\t#{strand}\t.\t0\t0\t0"
}
r.close()
w.close()

