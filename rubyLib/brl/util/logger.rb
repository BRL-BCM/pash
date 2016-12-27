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
# Basic message accumulation logger. Not too sophisticated. Used mostly by
# web-based versions of lffMerger and other tools.

# ##############################################################################
# REQUIRED CLASSES/MODULES
# ##############################################################################
require 'getoptlong'
require 'cgi'
require 'brl/util/util'					
require 'brl/util/propTable'
require 'brl/util/textFileUtil'
$VERBOSE = true

module BRL ; module Util
	
	class ErrorRecord
		attr_accessor :errorTitle, :messageLines
		
		def initialize(errorTitle, messageLines=nil)
			@errorTitle = errorTitle
			@messageLines = messageLines.nil?() ? [] : messageLines
		end
	
		def <<(msgLine)
			@messageLines << msgLine
			return
		end
	
		def to_s()
			errorStr = "#{@errorTitle}\n"
			@messageLines.each { |msg|
				errorStr += "  - #{msg}\n"
			}
			return errorStr
		end
	end
	
	class Logger
		attr_accessor :errorList, :activeError, :logFileName, :logFile, :nextToWriteIdx
	
		def initialize(logFileName=nil)
			@errorList = []
			@activeError = nil
			@nextToWriteIdx = 0
			unless(logFileName.nil?)
				@logFileName = logFileName
				@logFile = BRL::Util::TextWriter.new(logFileName)
			else
				@logFile = nil
			end
		end
	
		def size
			return @errorList.size
		end
	
		def <<(errorRecord)
			@errorList << errorRecord
			@activeError = @errorList.last
		end
	
		def addNewError(mainMsg, msgLines=nil)
			@errorList << ErrorRecord.new(mainMsg, msgLines)
			@activeError = @errorList.last
			return
		end
	
		def addToActiveError(msgLine)
			@activeError << msgLine
			return
		end
	
		def to_s(limit=nil)
			logStr = ''
			errCount = 0
			@errorList.each { |err|
				errCount += 1
				logStr += (errCount.to_s + '.  ' + err.to_s + "\n")
				if(!limit.nil? and errCount > limit)
					logStr += "\n...\n...\n [ There were more errors, but only reporting top #{limit}.\n"
					break
				end
			}
			return logStr
		end
	
		def writeToFile()
			return false if(@logFileName.nil?)
			begin
				logFile = BRL::Util::TextWriter.new(@logFileName)
				(@nextToWriteIdx...@errorList.size).each { |ii|
					logStr += ((ii+1).to_s + '.  ' + @errorList[ii].to_s + "\n\n")
					logFile.print logStr
					@nextToWriteIdx += 1
				}
			rescue => err
				$stderr.puts "\n\nFATAL ERROR: couldn't write to log file. Error details:\n\n#{err.message}\n\n#{err.backtrace}\n\n"
			end
			return true
		end
	
		def close()
			unless(@logFile.nil?)
				begin
					@logFile.close()
				rescue
				end
			end
			return
		end
	
		def clear()
			self.close()
			@errorList.clear()
			return
		end
	end
end ; end
