
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

# P_TYPE #
E  = 0
TF = 1

### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM INHERITED PROTEINS FILE ==="
  print "Usage: python inherited_proteins_to_js.py [parameters]"
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

### Generate the javascript to render the inherited proteins ###
# param PATH : the path of the simulation folder
def generate_inherited_proteins_js( PATH ):

  inherited_proteins_filename = PATH+"/statistics/best_inherited_proteins.txt"

  nodes = {}
  edges = []

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # A) Load nodes and edges         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # type id parent_id s p km kcat BStag CoEtag free_activity bound_activity window TFtag beta
  node_index = 1
  f = open(inherited_proteins_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")

    # EXTRACT GENETIC UNIT DATA #
    ptype = 0
    if l[0] == "ENZYME":
      ptype = E
    elif l[0] == "TRANSCRIPTION_FACTOR":
      ptype = TF
    functional = int(l[1])
    ident      = int(l[2])
    parent_id  = int(l[3])
    s          = int(l[4])
    p          = int(l[5])
    km         = float(l[6])
    kcat       = float(l[7])
    BStag      = int(l[8])
    CoEtag     = int(l[9])
    free_activity = "false"
    if l[10] == "1":
      free_activity = "true"
    bound_activity = "false"
    if l[11] == "1":
      bound_activity = "true"
    window = int(l[12])
    TFtag  = int(l[13])
    beta   = float(l[14])
  
    # ADD NODE #
    nodes[node_index] = {}
    nodes[node_index]["ptype"] = ptype
    nodes[node_index]["functional"] = functional
    nodes[node_index]["ident"] = ident
    nodes[node_index]["parent_id"] = parent_id
    nodes[node_index]["s"] = s
    nodes[node_index]["p"] = p
    nodes[node_index]["km"] = km
    nodes[node_index]["kcat"] = kcat
    nodes[node_index]["BStag"] = BStag 
    nodes[node_index]["CoEtag"] = CoEtag
    nodes[node_index]["free_activity"] = free_activity
    nodes[node_index]["bound_activity"] = bound_activity
    nodes[node_index]["window"] = window
    nodes[node_index]["TFtag"] = TFtag
    nodes[node_index]["beta"] = beta
    node_index += 1

    l = f.readline()
  f.close()

  # ADD EDGES #
  node_index -= 1
  for i in range(1, node_index):
    edges.append([i, i+1])

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # B) Write JS file                #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  f = open("./viewer/src/js/inherited_proteins_graph.js", "w")

  #-----------------------------#
  # B.1) Write header           #
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
  # B.2) Create cytoscape class #
  #-----------------------------#
  f.write("$(function(){ // on dom ready\n")
  f.write("\n")
  f.write("var inherited_proteins = cytoscape({\n")
  f.write("container: document.getElementById('inherited_proteins_graph'),\n")
  f.write("\n")

  #-----------------------------#
  # B.3) Declare graph style    #
  #-----------------------------#
  f.write("style: cytoscape.stylesheet()\n")
  f.write("  .selector('node').css({\n")
  f.write("    'content': 'data(name)',\n")
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
  f.write("    'target-arrow-shape': 'none',\n")
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
  # B.4) Declare the nodes      #
  #-----------------------------#
  f.write("nodes: [\n")

  for i in range(len(nodes.keys())):
    tag  = nodes.keys()[i]
    node = nodes[tag]
    line  = ""

    #-----------------------#
    # If the node is ENZYME #
    #-----------------------#
    if node["ptype"] == E:
      name = "E"
      if node["functional"] == 1:
        line += "    { data: { id: '"+str(tag)+"', name: '"+name+"', faveColor: 'rgb(188, 84, 77)', faveShape: 'ellipse' } }"
      else:
        line += "    { data: { id: '"+str(tag)+"', name: '"+name+"', faveColor: 'rgb(128, 128, 128)', faveShape: 'ellipse' } }"
    #-------------------------------------#
    # If the node is TRANSCRIPTION FACTOR #
    #-------------------------------------#
    elif node["ptype"] == TF:
      name ="TF"
      if node["functional"] == 1:
        line += "    { data: { id: '"+str(tag)+"', name: '"+name+"', faveColor: 'rgb(126, 103, 162)', faveShape: 'rectangle' } }"
      else:
        line += "    { data: { id: '"+str(tag)+"', name: '"+name+"', faveColor: 'rgb(128, 128, 128)', faveShape: 'rectangle' } }"

    if i < len(nodes.keys())-1:
      line += ",\n"
    else:
      line += "\n"

    f.write(line)

  f.write("  ], \n")
  f.write("  \n")
  
  #-----------------------------#
  # B.5) Declare the edges      #
  #-----------------------------#
  f.write("  edges: [\n")

  for i in range(len(edges)):
    edge = edges[i]
    line = ""
    line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: 'black' } }"
    if i < len(edges)-1:
      line += ",\n"
    else:
      line += "\n"

    f.write(line)

  f.write("  ]\n")

  f.write("  },\n")
  f.write("  \n")

  f.write("  layout: {\n")
  f.write("    name: 'cose',\n")
  f.write("    animate: false,\n")
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

  # GENERATE INHERITED PROTEINS JS #
  generate_inherited_proteins_js(PATH)
