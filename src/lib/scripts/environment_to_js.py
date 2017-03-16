
#!/usr/bin/env python
# coding: utf-8

#***************************************************************************
# Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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
  print "=== WRITE JAVASCRIPT FILES FROM ENVIRONMENT FILES ==="
  print "Usage: python networks_to_js.py [parameters]"
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

### Generate a list of prime numbers ###
# param maximum : maximum number used to build the list
def build_prime_numbers_list( maximum ):
  prime_numbers = [0]*maximum
  for i in range(1,maximum):
    prime_numbers[i] = 1
  for i in range(2, maximum+1):
    multiple = 2*i
    while multiple <= maximum:
      prime_numbers[multiple-1] = 0
      multiple += i
  liste = []
  for i in range(maximum):
    if prime_numbers[i] == 1:
      liste.append(i+1)
  return liste

### Generate the javascript to render the global environment concentrations ###
# param PATH : the path of the simulation folder
def global_environment_js( PATH ):
  prime_numbers = build_prime_numbers_list(500)

  tags = []
  concs = []
  f = open(PATH+"/statistics/last_environment_metabolic_amounts.txt", "r")
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split()
    tag = int(l[0])
    conc = float(l[1])
    if conc > 0.01:
      tags.append(tag)
      concs.append(conc)
    l = f.readline()
  f.close()

  f = open(PATH+"/viewer/src/js/environment.js", "w")

  f.write("/****************************************************************************\n")
  f.write(" * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon\n")
  f.write(" * E-mail: charles.rocabert@gmail.com\n")
  f.write(" * Web: http://www.evoevo.eu/\n")
  f.write(" *\n")
  f.write(" * This program is free software: you can redistribute it and/or modify\n")
  f.write(" * it under the terms of the GNU General Public License as published by\n")
  f.write(" * the Free Software Foundation, either version 3 of the License, or\n")
  f.write(" * (at your option) any later version.\n")
  f.write(" *\n")
  f.write(" * This program is distributed in the hope that it will be useful,\n")
  f.write(" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
  f.write(" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
  f.write(" * GNU General Public License for more details.\n")
  f.write(" *\n")
  f.write(" * You should have received a copy of the GNU General Public License\n")
  f.write(" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n")
  f.write(" ****************************************************************************/\n\n")

  f.write("var ctx = document.getElementById('environment_chart').getContext('2d');\n")
  f.write("Chart.defaults.global.responsive = true;\n")
  f.write("var data = {\n")

  #--------------#
  # WRITE LABELS #
  #--------------#
  labels_line = ""
  labels_line += "  labels: ["
  for tag in tags:
    if tag in prime_numbers:
      labels_line += "'|"+str(tag)+"|',"
    else:
      labels_line += "'"+str(tag)+"',"
  labels_line = labels_line.strip(",")
  labels_line += "],\n"
  f.write(labels_line)

  f.write("  datasets: [\n")
  f.write("  {\n")
  f.write("    label: 'Global environment state',\n")
  f.write("    fillColor: 'rgba(220,120,0,0.5)',\n")
  f.write("    strokeColor: 'rgba(220,120,0,0.8)',\n")
  f.write("    highlightFill: 'rgba(220,120,0,0.75)',\n")
  f.write("    highlightStroke: 'rgba(220,120,0,1)',\n")

  #------------#
  # WRITE DATA #
  #------------#
  data_line = ""
  data_line += "    data: ["
  for conc in concs:
    data_line += str(conc)+","
  data_line = data_line.strip(",")
  data_line += "],\n"
  f.write(data_line)

  f.write("  }\n")
  f.write("  ]\n")
  f.write("};\n")
  f.write("var myBarChart = new Chart(ctx).Bar(data);\n\n")

  f.close()

### Generate the javascript to render the best individual local environment ###
# param PATH : the path of the simulation folder
def best_individual_local_environment_js( PATH ):
  prime_numbers = build_prime_numbers_list(500)

  tags = []
  concs = []
  f = open(PATH+"/statistics/last_local_environment_metabolic_amounts.txt", "r")
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split()
    tag = int(l[0])
    conc = float(l[1])
    if conc > 0.01:
      tags.append(tag)
      concs.append(conc)
    l = f.readline()
  f.close()

  f = open(PATH+"/viewer/src/js/local_environment.js", "w")

  f.write("/****************************************************************************\n")
  f.write(" * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon\n")
  f.write(" * E-mail: charles.rocabert@gmail.com\n")
  f.write(" * Web: http://www.evoevo.eu/\n")
  f.write(" *\n")
  f.write(" * This program is free software: you can redistribute it and/or modify\n")
  f.write(" * it under the terms of the GNU General Public License as published by\n")
  f.write(" * the Free Software Foundation, either version 3 of the License, or\n")
  f.write(" * (at your option) any later version.\n")
  f.write(" *\n")
  f.write(" * This program is distributed in the hope that it will be useful,\n")
  f.write(" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
  f.write(" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
  f.write(" * GNU General Public License for more details.\n")
  f.write(" *\n")
  f.write(" * You should have received a copy of the GNU General Public License\n")
  f.write(" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n")
  f.write(" ****************************************************************************/\n\n")

  f.write("var ctx = document.getElementById('local_environment_chart').getContext('2d');\n")
  f.write("Chart.defaults.global.responsive = true;\n")
  f.write("var data = {\n")

  #--------------#
  # WRITE LABELS #
  #--------------#
  labels_line = ""
  labels_line += "  labels: ["
  for tag in tags:
    if tag in prime_numbers:
      labels_line += "'|"+str(tag)+"|',"
    else:
      labels_line += "'"+str(tag)+"',"
  labels_line = labels_line.strip(",")
  labels_line += "],\n"
  f.write(labels_line)

  f.write("  datasets: [\n")
  f.write("  {\n")
  f.write("    label: 'Local environment state',\n")
  f.write("    fillColor: 'rgba(0,120,220,0.5)',\n")
  f.write("    strokeColor: 'rgba(0,120,220,0.8)',\n")
  f.write("    highlightFill: 'rgba(0,120,220,0.75)',\n")
  f.write("    highlightStroke: 'rgba(0,120,220,1)',\n")

  #------------#
  # WRITE DATA #
  #------------#
  data_line = ""
  data_line += "    data: ["
  for conc in concs:
    data_line += str(conc)+","
  data_line = data_line.strip(",")
  data_line += "],\n"
  f.write(data_line)

  f.write("  }\n")
  f.write("  ]\n")
  f.write("};\n")
  f.write("var myBarChart = new Chart(ctx).Bar(data);\n\n")

  f.close()


############
#   MAIN   #
############
if __name__ == '__main__':

  # START #
  PATH = read_args(sys.argv)

  # GENERATE GLOBAL ENVIRONMENT JS #
  global_environment_js(PATH)

  # GENERATE BEST INDIVIDUAL LOCAL ENVIRONMENT JS #
  best_individual_local_environment_js(PATH)

