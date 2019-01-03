
#!/usr/bin/env python
# coding: utf-8

#***************************************************************************
# Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# Web: https://github.com/charlesrocabert/Evo2Sim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***************************************************************************

import sys
import os
import numpy as np

### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE THE ENVIRONMENT HEATMAP ==="
  print "Usage: python environment_heatmap.py [parameters]"
  print "Parameters are:"
  print "-h, --help:"
  print "    Print this help, then exit."
  print "-p, --path:"
  print "    Give the path of the simulation folder."
  print ""

### READ COMMAND LINE ARGUMENTS ###
# param  string-array args : list of command line arguments
def read_args( args ):
  if len(args) < 2:
    print("Lack of parameters, see help (-h --help)")
    sys.exit()
  else:
    PATH    = ""
    for i in range(len(args)):
      if args[i] == "-h" or args[i] == "--help":
        print_usage()
        sys.exit()
      elif args[i] == "-p" or args[i] == "--path":
        if i+1 >= len(args):
          print("Lack of parameters, see help (-h --help)")
          sys.exit()
        else:
          PATH = args[i+1]
    if PATH == "":
      print("Lack of parameters, see help (-h --help)")
      sys.exit()
    else:
      return PATH

############
#   MAIN   #
############
if __name__ == '__main__':
  
  PATH = read_args(sys.argv)
  
  #-------------------------------------------#
  # 1) open environment file and get measures #
  #-------------------------------------------#
  f = open(PATH+"/statistics/environment_metabolic_amounts.txt", "r")
  maxlength = 0
  nblines   = 0
  l = f.readline()
  while l:
    nblines += 1
    l = l.strip("\n")
    l = l.split(" ")
    if maxlength < len(l):
      maxlength = len(l)
    l = f.readline()
  f.close()

  #-------------------------------------------#
  # 2) build the heatmap                      #
  #-------------------------------------------#
  f = open(PATH+"/statistics/environment_metabolic_amounts.txt", "r")
  heatmap = np.zeros((maxlength, nblines-1))
  l = f.readline()
  count = 0
  while l and count < nblines-1:
    l = l.strip("\n")
    l = l.split(" ")
    for i in range(maxlength):
      if len(l) > i:
        heatmap[i, count] = float(l[i])
      else:
        heatmap[i, count] = 0.0
    count += 1
    l = f.readline()
  f.close()

  #-------------------------------------------#
  # 3) write the heatmap in a file            #
  #-------------------------------------------#
  f = open(PATH+"/statistics/environment_heatmap.txt", "w")
  for i in range(nblines-1):
    line = ""
    for j in range(maxlength):
      line += str(heatmap[j, i])+" "
    line = line.strip(" ")
    line += "\n"
    f.write(line)
  f.close()
