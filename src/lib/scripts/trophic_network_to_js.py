
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

# TROPHIC_LEVEL #
LEVEL_0  = 0
LEVEL_1  = 1
LEVEL_2  = 2
NO_LEVEL = 3

### node class ###
class Node:
  def __init__( self, identifier, production_profile, uptake_profile, release_profile, level, count, appearance, lifespan ):
    self.identifier           = identifier
    self.production_profile   = production_profile
    self.uptake_profile       = uptake_profile
    self.release_profile      = release_profile
    self.level                = level
    self.count                = count
    self.appearance           = appearance
    self.lifespan             = lifespan
    self.necrophagy_links     = []
    self.active_release_links = []

### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM TROPHIC NETWORK FILES ==="
  print "Usage: python2 trophic_network_to_js.py [parameters]"
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
def generate_trophic_network_js( PATH, count_threshold ):

  trophic_network_nodes_filename = PATH+"/statistics/trophic_network_nodes.txt"
  trophic_network_edges_filename = PATH+"/statistics/trophic_network_edges.txt"

  nodes = {}

  #~~~~~~~~~~~~~~~~~~#
  # A) Load nodes    #
  #~~~~~~~~~~~~~~~~~~#

  # id production uptake release level count
  f = open(trophic_network_nodes_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    identifier         = int(l[0])
    production_profile = l[1]
    uptake_profile     = l[2]
    release_profile    = l[3]
    level              = int(l[4])
    count              = int(l[5])
    appearance         = int(l[6])
    lifespan           = int(l[7])
    if count > count_threshold or identifier == 0:
      node = Node(identifier, production_profile, uptake_profile, release_profile, level, count, appearance, lifespan)
      nodes[identifier] = node
    l = f.readline()
  f.close()

  #~~~~~~~~~~~~~~~~~~#
  # B) Load edges    #
  #~~~~~~~~~~~~~~~~~~#
  
  # id1 id2 link_type
  f = open(trophic_network_edges_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    source = int(l[0])
    target = int(l[1])
    ltype  = int(l[2])
    if ltype == 0 and source in nodes.keys() and target in nodes.keys():
      nodes[source].necrophagy_links.append(target)
    elif ltype == 1 and source in nodes.keys() and target in nodes.keys():
      nodes[source].active_release_links.append(target)
    l = f.readline()
  f.close()

  #~~~~~~~~~~~~~~~~~~#
  # C) Write JS file #
  #~~~~~~~~~~~~~~~~~~#

  f = open("./viewer/src/js/trophic_network_graph.js", "w")

  #-----------------------------#
  # C.1) Write header           #
  #-----------------------------#
  f.write("/****************************************************************************\n")
  f.write(" * Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon\n")
  f.write(" * Web: https://github.com/charlesrocabert/Evo2Sim")
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

  #-----------------------------#
  # C.2) Create cytoscape class #
  #-----------------------------#
  f.write("$(function(){ // on dom ready\n")
  f.write("\n")
  f.write("var trophic_network = cytoscape({\n")
  f.write("container: document.getElementById('trophic_network_graph'),\n")
  f.write("\n")

  #-----------------------------#
  # C.3) Declare graph style    #
  #-----------------------------#
  f.write("style: cytoscape.stylesheet()\n")
  f.write("  .selector('node').css({\n")
  f.write("    'content': 'data(name)',\n")
  f.write("    'font-size': 10,\n")
  f.write("    'text-valign': 'center',\n")
  f.write("    'text-outline-width': 2,\n")
  f.write("    'shape': 'data(faveShape)',\n")
  f.write("    'text-outline-color': 'data(faveColor)',\n")
  f.write("    'background-color': 'data(faveColor)',\n")
  f.write("    'color': '#fff',\n")
  f.write("    'width': 30\n")
  f.write("  })\n")
  f.write("  .selector('edge').css({\n")
  f.write("    'opacity': 0.666,\n")
  f.write("    'target-arrow-shape': 'triangle',\n")
  f.write("    'width': 2,\n")
  f.write("    'line-color': 'data(faveColor)',\n")
  f.write("    'source-arrow-color': 'data(faveColor)',\n")
  f.write("    'target-arrow-color': 'data(faveColor)'\n")
  f.write("  })\n")
  f.write("  .selector('.highlighted').css({\n")
  f.write("    'background-color': '#61bffc',\n")
  f.write("    'line-color': '#61bffc',\n")
  f.write("    'target-arrow-color': '#61bffc',\n")
  f.write("    'transition-property': 'background-color, line-color, target-arrow-color',\n")
  f.write("    'transition-duration': '0.5s'\n")
  f.write("  }),\n")
  f.write("\n")

  f.write("elements: {\n")
  f.write("\n")

  #-----------------------------#
  # C.4) Declare the nodes      #
  #-----------------------------#
  f.write("nodes: [\n")

  for i in range(len(nodes.keys())):
    identifier = nodes.keys()[i]
    node       = nodes[identifier]
    name       = ""

    #####################
    # Create group name #
    #####################
    for met in range(len(node.uptake_profile)):
      if node.uptake_profile[met] == '1':
        name += str(met+1)+"-"
    name = name.strip("-")
    name += " ("+str(node.count)+", "+str(node.lifespan)+")"

    line = ""

    if identifier == 0:
      env_name = "ENV ("
      for met in range(len(node.production_profile)):
        if node.production_profile[met] == '1':
          env_name += str(met+1)+"-"
      env_name = env_name.strip("-")
      env_name += ")"
      line += "    { data: { id: '"+str(identifier)+"', name: '"+env_name+"', faveColor: 'black', faveShape: 'roundrectangle' } }"

    elif identifier != 0 and node.level == LEVEL_0:
      line += "    { data: { id: '"+str(identifier)+"', name: '"+name+"', faveColor: 'purple', faveShape: 'roundrectangle' } }"

    elif identifier != 0 and node.level == LEVEL_1:
      line += "    { data: { id: '"+str(identifier)+"', name: '"+name+"', faveColor: '#6FB1FC', faveShape: 'roundrectangle' } }"

    elif identifier != 0 and node.level == LEVEL_2:
      line += "    { data: { id: '"+str(identifier)+"', name: '"+name+"', faveColor: '#86B342', faveShape: 'roundrectangle' } }"

    elif identifier != 0 and node.level == NO_LEVEL:
      line += "    { data: { id: '"+str(identifier)+"', name: '"+name+"', faveColor: 'grey', faveShape: 'roundrectangle' } }"

    if i < len(nodes.keys())-1:
      line += ",\n"
    else:
      line += "\n"

    f.write(line)

  f.write("  ], \n")
  f.write("  \n")
  
  #-----------------------------#
  # C.5) Declare the edges      #
  #-----------------------------#
  f.write("  edges: [\n")

  # For each node #
  for i in range(len(nodes.keys())):
    identifier = nodes.keys()[i]
    node       = nodes[identifier]
    line       = ""

    # Explore necrophagy links #
    for j in range(len(node.necrophagy_links)):
      target = node.necrophagy_links[j]
      line  += "    { data: { id: '"+str(identifier)+"-"+str(target)+"', weight: 1, source: '"+str(identifier)+"', target: '"+str(target)+"', faveColor: 'black' } }"
      line  += ",\n"

    # Explore active release links #
    for j in range(len(node.active_release_links)):
      target = node.active_release_links[j]
      line  += "    { data: { id: '"+str(identifier)+"-"+str(target)+"', weight: 1, source: '"+str(identifier)+"', target: '"+str(target)+"', faveColor: '#6FB1FC' } }"
      if i == len(nodes.keys())-1 and j == len(node.active_release_links)-1:
        line += "\n"
      else:
        line += ",\n"

    f.write(line)

  f.write("  ]\n")

  f.write("  },\n")
  f.write("  \n")

  f.write("  layout: {\n")
  #f.write("    name: 'breadthfirst',\n")
  #f.write("    animate: false,\n")
  #f.write("    circle: true,\n")
  #f.write("    spacingFactor: 4,\n")
  #f.write("    directed: true,\n")
  #f.write("    padding: 10\n")
  f.write("    name: 'cose',\n")
  f.write("    animate: false,\n")
  f.write("    directed: true,\n")
  f.write("    padding: 10,\n")
  f.write("  }\n")
  f.write("});\n")
  f.write("\n")
  f.write("}); // on dom ready\n")

  f.close()

############
#   MAIN   #
############
if __name__ == '__main__':

  # START #
  PATH = read_args(sys.argv)

  # GENERATE TROPHIC NETWORK JS #
  generate_trophic_network_js(PATH, 10)
