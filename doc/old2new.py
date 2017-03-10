#!/bin/env python

import sys
import os

if(len(sys.argv) == 1):
  print sys.argv[0], "converts old wordom input files (pre-0.21) to the new format"
  print "Usage:", sys.argv[0], "oldinput"
  sys.exit(0)

oldfile = open(sys.argv[1])
oldinput = oldfile.readlines()

for ii in range(len(oldinput)):
  if oldinput[ii][0] == '$' and oldinput[ii][1] == '$' and oldinput[ii][2] == '$' :
    if oldinput[ii].split()[1] == "BEGIN" :
      print "BEGIN "+oldinput[ii].split()[2]
    elif oldinput[ii].split()[1] == "END" :
      print "END"
  elif oldinput[ii][0] == '$' and oldinput[ii][1] == ' ' and oldinput[ii].split()[1] == "TITLE" :
    print "--TITLE "+oldinput[ii].split()[2]
  elif oldinput[ii][0] == '$' :
    print "--"+oldinput[ii][1:-1]
  elif oldinput[ii][0] == '/' :
    print "--SELE "+oldinput[ii][0:-1]

sys.exit(0)
