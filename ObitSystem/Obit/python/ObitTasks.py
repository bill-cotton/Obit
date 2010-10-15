"""Functions providing easy access to Obit task descriptions."""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

class EnvVarError(Exception):
    """Raised when there is a problem with an environment variable."""
    def __init__(self, msg):
        self.msg = msg

def obitTaskDict():
    """Return a dictionary with Obit task name keys and task description
       values."""
    import glob, re, os
    
    # List of all dirs containing TDF files
    try:
        TDFDirs = [ os.environ['OBIT']   + '/TDF',
                    os.environ['OBITSD'] + '/TDF' ] 
    except KeyError:
        msg = "Error: Environment variable OBIT or OBITSD " +\
              "not defined. These variables can be defined by sourcing the " +\
              "Obit setup script before starting ObitTalk."
        raise EnvVarError(msg)
    
    # Create a list of all .TDF files in directories TDFDirs
    fileList = []
    for dir in TDFDirs:
        fileList =  fileList + glob.glob(dir + "/*.TDF") # Full paths to TDFs
    
    taskDict = {}
    
    # Step through each TDF file in list
    for i,file in enumerate(fileList):
        # Extract the task name from the file name
        task = os.path.basename(file).rsplit('.',1)[0] 
    
        # Extract 2nd task description line from the TDF file. 
        # There are multiple description lines. We want the 2nd one.
        TDFFile = open( file )
        taskLine = []
        for line in TDFFile:
            if line.startswith(task): # description lines start w/ task name
                taskLine.append(line)
        taskLine = taskLine[1] # Keep only 2nd line
            
        # Get the task description
        # pattern = taskname, any whitespace, 0 or 1 colon, any whitespace, 
        #           all the rest
        pattern = re.compile( task + r"\s*:?\s*(.*)" )
        matchObj = pattern.match(taskLine) 
        description = matchObj.group(1)
    
        # Add task to dict, with description as the value
        taskDict[task] = description
        
    return taskDict

def obitTaskList():        
    """Print obit task names and descriptions."""
    try:
        taskDict = obitTaskDict()
    except EnvVarError, e:
        print e.msg
        return
    for key in taskDict.keys(): 
        print "%-10s  %s" % (key, taskDict[key])
