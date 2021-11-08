
#!/usr/bin/env python
# coding: utf-8

#****************************************************************************
# Evo2Sim (Evolution of Evolution Simulator)
# -------------------------------------------
# Digital evolution model dedicated to
# bacterial in silico experimental evolution.
#
# Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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


### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM TROPHIC GROUPS FILES ==="
  print "Usage: python2 trophic_groups_to_js.py [parameters]"
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

### Generate the javascript to render the trophic network ###
# param PATH            : the path of the simulation folder
# param count_threshold : defines the minimum number of cells at which the group is considered
def generate_trophic_groups_map( PATH, count_threshold ):
  groups = {}
  current_groups = {}
  f = open(PATH+"statistics/trophic_network_groups.txt", "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    



############
#   MAIN   #
############
if __name__ == '__main__':

  # START #
  PATH = read_args(sys.argv)

  # GENERATE TROPHIC NETWORK JS #
  generate_trophic_network_js(PATH, 10)
