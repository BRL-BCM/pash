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
# Convenience classes for reading/writing plain text that *may* be gzipped.
# We let the instance figure out whether the text is gzipped or not when
# reading. Uses gzip itself to test the file (not the 'file' command).

# ##############################################################################
# REQUIRED LIBRARIES
# ##############################################################################
require 'zlib'
require 'delegate'
require 'brl/util/util'
$VERBOSE = true

module BRL ; module Util

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
