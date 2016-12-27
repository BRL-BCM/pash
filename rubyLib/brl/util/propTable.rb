#!/usr/bin/env ruby
# ##############################################################################
# Copyright (c) 2004 Baylor College of Medicine.
#
# Use of this software is governed by a license. 
# See the included file LICENSE.TXT for details.
# ##############################################################################

# Author: Andrew R Jackson (andrewj@bcm.tmc.edu)
# Date: 3/31/2004 4:38PM
# Purpose:
# Class that reads properties files (similar to Java ones, perhaps with a few
# more "extras" or assumptions about appropriate data-types, and presents as
# ~Hashes.

# ##############################################################################
# REQUIRED LIBRARIES
# ##############################################################################
require 'brl/util/util'
require 'brl/util/textFileUtil'
$VERBOSE = true

module BRL ; module Util
	include BRL::Util
	class ParseError < Exception; end
	class MissingKeyError < Exception; end

	class PropTable < Hash
		attr_accessor(:splitCSV)

		def initialize(inputSrc=nil, splitCSV=true)
			@splitCSV = splitCSV
			return if(inputSrc.nil?)
			load(inputSrc)
		end

		def load(inputSrc)
			lineCounter = 0
			line = nil
			begin
				if(inputSrc.kind_of?(String) or inputSrc.kind_of?(IO) or inputSrc.kind_of?(BRL::Util::TextReader)) # then init from string
					inputSrc.each { |line|
						lineCounter+=1
						self.parseLine(line)
					}
				elsif(inputSrc.kind_of?(Hash)) # then init from hash
					inputSrc.each {	|key, val|
						valueArray = (@splitCSV and val.kind_of?(String) ? val.split(',') : [val])
						self[key] = (valueArray.length == 1 ? valueArray[0] : valueArray)
					}
				else
					raise(TypeError, "Cannot initialize a PropTable from #{inputSrc.type}. Try a String, Hash, or IO object.");
				end
			rescue BRL::Util::ParseError => err
				raise(BRL::Util::ParseError, "Bad line in properties file/string:\n'#{line}'\nat line #{lineCounter}")
			ensure
				inputSrc.close() if(inputSrc.kind_of?(IO) or inputSrc.kind_of?(BRL::Util::TextReader))
			end
		end # load(inputSrc)

		def store(outputDest, doGzip=false)
			unless(outputDest.kind_of?(String))
				raise(TypeError, "BRL::Util::PropTable#store(outputDest) requires a String that is the path to the output file and optionally a boolean doGzip flag whose default is false");
			end
			writer = BRL::Util::TextWriter.new(outputDest, doGzip)
			self.each { |key, val|
				writer.write("#{key} = #{val}\n")
			}
			writer.close
		end # store(outputDest)

		def verify(propStringArray)
			unless(propStringArray.kind_of?(Array))
				raise(TypeError, "BRL::Util::PropTable#verify(propStringArray) requires an array of Strings to check against")
			end
			propStringArray.each {
				|elem|
				unless(self.key?(elem) and self[elem] !~ /^\s*$/)
					raise(BRL::Util::MissingKeyError, "A required key (#{elem}) is missing in the properties file/string")
				end
			}
		end

		# ##########################################################################
		protected # PROTECTED METHODS
		# ##########################################################################
		def parseLine(line)
			propName, valueStr = nil
			return if(line =~ /^\s*#/) # skip comment lines
			return if(line =~ /^\s*$/) # skip blank lines
			# Isolate the name and the value
			if(line =~ /^\s*(\S+)\s*=\s*(.+)$/)
				propName = $1.strip
				valueStr = $2.strip
				# Is the valueStr quoted? If so, extract what is inside as literal value
				if(valueStr =~ /^(?:(?:\"([^\n\r\"]*)\")|(?:\'([^\n\r\']*)\'))$/)
					if(!$1.nil? and !$1.empty?) # then double-quoted string
						valueArray = [$1]
					elsif(!$2.nil? and !$2.empty?) # then singble-quoted string
						valueArray = [$2]
					else
						valueArray = [ nil ]
					end
				else # not quoted. Attempt to split by commas
					valueArray = (@splitCSV ? valueStr.split(/\s*,\s*/) : [ valueStr ])
				end
				# Return either the array, or if only 1 value just return it not array.
				self[propName] = (valueArray.length == 1 ? valueArray[0] : valueArray)
			else
				raise BRL::Util::ParseError
			end
		end
	end # class PropTable < Hash
end ; end # module BRL ; module Util
