#!python
# Script to run gcc preprocessor on files whose name ends in "Def.h"
# or echo the file contents otherwise
# File name is only argument
import os, sys

# Routine to recursively read include
def parseInclude(lline):
    fname = lline.split()[1].replace('"','')
    infil = os.getenv("OBIT")+"include/"+fname
    jnput = open(infil)
    line = " "
    while (line):
        line = jnput.readline() # read next line
        if not line:  # EOF?
            break
        # Def include to be expanded?
        if line.startswith('#include "Obit') and (line.find('Def.h"')>=0):
            parseInclude(line)
            continue
        x=sys.stdout.write(line)
    # end parseInclude

# Get filename argument
infile=sys.argv[1]

# read file recursively including
input = open(infile)
line = " "
while (line):
    line = input.readline() # read next line
    if not line:  # EOF?
        break
    # Def include to be expanded?
    if line.startswith('#include "Obit') and (line.find('Def.h"')>=0):
        parseInclude(line)
        continue
    x=sys.stdout.write(line)
