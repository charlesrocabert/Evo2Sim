
#!/usr/bin/env python
# coding: utf-8

#***************************************************************************
# Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# E-mail: charles.rocabert@gmail.com
# Web: http://www.evoevo.eu/
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
E   = 0
TF  = 1
COE = 2

# E_TYPE #
INNER_ENZYME = 0
UPTAKE_PUMP  = 1
RELEASE_PUMP = 2

# TF_TYPE #
REPRESSOR        = 0
ACTIVATOR        = 1
ALWAYS_ACTIVE    = 2
ALWAYS_REPRESSED = 3 

### node class ###
class Node:
  def __init__( self, ident, conc, ptype, s, p, km, kcat, e_type, inherited ):
    self.ident     = ident
    self.conc      = conc
    self.type      = ptype
    self.s         = s
    self.p         = p
    self.km        = km
    self.kcat      = kcat
    self.e_type    = e_type
    self.inherited = inherited
    self.tf_type   = 0
    self.pred      = []
    self.succ      = []

### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM NETWORK FILES ==="
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

### Generate the javascript to render the GRN ###
# param PATH : the path of the simulation folder
def generate_grn_js( PATH ):

  edges_filename = PATH+"/statistics/best_grn_edges.txt"
  nodes_filename = PATH+"/statistics/best_grn_nodes.txt"

  #~~~~~~~~~~~~~~~~~~#
  # A) Load nodes    #
  #~~~~~~~~~~~~~~~~~~#
  nodes = {}
  edges = []

  f = open(nodes_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    ident     = int(l[0])
    conc      = float(l[1])
    ptype     = 0
    if l[2] == "E":
      ptype = E
    elif l[2] == "TF":
      ptype = TF
    elif l[2] == "COE":
      ptype = COE
    s         = int(l[3])
    p         = int(l[4])
    km        = float(l[5])
    kcat      = float(l[6])
    inherited = int(l[7])
    e_type    = 0
    if s == p and kcat > 0.0:
      e_type = UPTAKE_PUMP
    elif s == p and kcat < 0.0:
      kcat = abs(kcat)
      e_type = RELEASE_PUMP
    elif s != p and kcat > 0.0:
      e_type = INNER_ENZYME
    elif s != p and kcat < 0.0:
      tmp = s
      s = p
      p = tmp
      kcat = abs(kcat)
      e_typ = INNER_ENZYME
    node = Node(ident, conc, ptype, s, p, km, kcat, e_type, inherited)
    nodes[ident] = node
    l = f.readline()
  f.close()

  #~~~~~~~~~~~~~~~~~~#
  # B) Load edges    #
  #~~~~~~~~~~~~~~~~~~#
  f = open(edges_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    print l
    l = l.strip("\n")
    l = l.split(" ")
    source = int(l[0])
    target = int(l[1])
    aff    = float(l[2])
    present = False
    for i in range(len(edges)):
      if edges[i][0] == source and edges[i][1] == target:
        edges[i][2] += aff
        present = True
        break
    if not present:
      edges.append([source, target, aff])
    l = f.readline()
  f.close()

  #~~~~~~~~~~~~~~~~~~#
  # C) Complete data #
  #~~~~~~~~~~~~~~~~~~#
  for edge in edges:
    if nodes[edge[0]].type == COE and edge[2] == REPRESSOR:
      nodes[edge[1]].tf_type = REPRESSOR
    elif nodes[edge[0]].type == COE and edge[2] == ACTIVATOR:
      nodes[edge[1]].tf_type = ACTIVATOR
    elif nodes[edge[0]].type == COE and edge[2] == ALWAYS_ACTIVE:
      nodes[edge[1]].tf_type = ALWAYS_ACTIVE
    elif nodes[edge[0]].type == COE and edge[2] == ALWAYS_REPRESSED:
      nodes[edge[1]].tf_type = ALWAYS_REPRESSED

  #~~~~~~~~~~~~~~~~~~#
  # D) Write JS file #
  #~~~~~~~~~~~~~~~~~~#
  f = open("./viewer/src/js/grn_graph.js", "w")

  #-----------------------------#
  # D.1) Write header           #
  #-----------------------------#
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

  #-----------------------------#
  # D.2) Create cytoscape class #
  #-----------------------------#
  f.write("$(function(){ // on dom ready\n")
  f.write("\n")
  f.write("var grn = cytoscape({\n")
  f.write("container: document.getElementById('grn_graph'),\n")
  f.write("\n")
  
  #-----------------------------#
  # D.3) Declare graph style    #
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
  f.write("  .selector('$node > node').css({\n")
  f.write("    'text-outline-color': 'data(faveColor)',\n")
  f.write("    'border-color': 'data(faveColor)',\n")
  f.write("    'border-width': 2,\n")
  f.write("    'background-color': 'white',\n")
  f.write("    'color': 'white',\n")
  f.write("    'padding-top': '10px',\n")
  f.write("    'padding-left': '10px',\n")
  f.write("    'padding-bottom': '10px',\n")
  f.write("    'padding-right': '10px',\n")
  f.write("    'text-valign': 'top',\n")
  f.write("    'text-halign': 'center'\n")
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
  f.write("  })\n")
  f.write("  .selector('edge.negative_regulation').css({\n")
  f.write("    'line-style': 'dashed',\n")
  f.write("    'target-arrow-shape': 'diamond'\n")
  f.write("  }),\n")
  f.write("\n")

  f.write("elements: {\n")
  f.write("\n")

  #-----------------------------#
  # D.4) Declare the nodes      #
  #-----------------------------#
  f.write("nodes: [\n")

  f.write("    { data: { id: 'COE', name: 'Co-enzymes', faveColor: '#86B342', faveShape: 'rectangle' } },\n")
  f.write("    { data: { id: 'TF', name: 'Transcription factors', faveColor: '#F5A45D', faveShape: 'rectangle' } },\n")
  f.write("    { data: { id: 'E', name: 'Enzymes', faveColor: '#6FB1FC', faveShape: 'rectangle' } },\n")

  for i in range(len(nodes.keys())):
    curnode = nodes[nodes.keys()[i]]
    name = ""
    line = ""

    # If the node is an ENZYME #
    if curnode.type == E:
      name = ""
      if curnode.e_type == INNER_ENZYME:
        name = "E("+str(curnode.s)+" -> "+str(curnode.p)+")"
      elif curnode.e_type == UPTAKE_PUMP:
        name = "E("+str(-curnode.s)+" -> "+str(curnode.p)+")"
      elif curnode.e_type == RELEASE_PUMP:
        name = "E("+str(curnode.s)+" -> "+str(-curnode.p)+")"
      line += "    { data: { id: '"+str(curnode.ident)+"', parent: 'E', name: '"+name+"', faveColor: '#6FB1FC', faveShape: 'ellipse' } }"

    # Else if the node is a TRANSCRIPTION FACTOR #
    elif curnode.type == TF:
      name = ""
      if curnode.type == ACTIVATOR and curnode.inherited == 0:
        name = "'TF("+str(curnode.ident)+", A, NI)'"
      elif curnode.type == ACTIVATOR and curnode.inherited == 1:
        name = "'TF("+str(curnode.ident)+", A, I)'"
      elif curnode.type == REPRESSOR and curnode.inherited == 0:
        name = "'TF("+str(curnode.ident)+", R, NI)'"
      elif curnode.type == REPRESSOR and curnode.inherited == 1:
        name = "'TF("+str(curnode.ident)+", R, I)'"
      elif curnode.type == ALWAYS_ACTIVE and curnode.inherited == 0:
        name = "'TF("+str(curnode.ident)+", AA, NI)'"
      elif curnode.type == ALWAYS_ACTIVE and curnode.inherited == 1:
        name = "'TF("+str(curnode.ident)+", AA, I)'"
      elif curnode.type == ALWAYS_REPRESSED and curnode.inherited == 0:
        name = "'TF("+str(curnode.ident)+", AI, NI)'"
      elif curnode.type == ALWAYS_REPRESSED and curnode.inherited == 1:
        name = "'TF("+str(curnode.ident)+", AI, I)'"
      line += "    { data: { id: '"+str(curnode.ident)+"', parent: 'TF', name: "+name+", faveColor: '#F5A45D', faveShape: 'rectangle' } }"

    # Else if the node is a CO-ENZYME #
    elif curnode.type == COE:
      name = "CoE("+str(curnode.s+1)+")"
      line += "    { data: { id: '"+str(curnode.ident)+"', parent: 'COE', name: '"+name+"', faveColor: '#86B342', faveShape: 'triangle' } }"

    if i < len(nodes.keys())-1:
      line += ",\n"
    else:
      line += "\n"

    f.write(line)

  f.write("  ], \n")
  f.write("  \n")
  
  #-----------------------------#
  # D.5) Declare the edges      #
  #-----------------------------#
  f.write("  edges: [\n")

  for i in range(len(edges)):
    edge  = edges[i]
    snode = nodes[edge[0]]
    tnode = nodes[edge[1]]
    line  = ""

    # If the source is an ENZYME #
    if snode.type == E:
      line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#6FB1FC' } }"

    # Else if the source is a TRANSCRIPTION FACTOR #
    elif snode.type == TF:
      if edge[2] < 0.0:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#F5A45D' }, classes: 'negative_regulation' }"
      else:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#F5A45D' }, classes: 'positive_regulation' }"

    # Else if the source is a CO-ENZYME #
    elif snode.type == COE:
      if tnode.tf_type == ACTIVATOR:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#86B342' }, classes: 'positive_regulation' }"
      elif tnode.tf_type == REPRESSOR:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#86B342' }, classes: 'negative_regulation' }"
      elif tnode.tf_type == ALWAYS_REPRESSED:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#86B342' }, classes: 'negative_regulation' }"
      elif tnode.tf_type == ALWAYS_ACTIVE:
        line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(edge[0])+"', target: '"+str(edge[1])+"', faveColor: '#86B342' }, classes: 'positive_regulation' }"

    if i < len(edges)-1:
      line += ",\n"
    else:
      line += "\n"

    f.write(line)

  f.write("  ]\n")

  f.write("  },\n")
  f.write("  \n")

  f.write("  layout: {\n")
  f.write("    name: 'breadthfirst',\n")
  f.write("    animate: false,\n")
  f.write("    circle: true,\n")
  f.write("    spacingFactor: 4,\n")
  f.write("    directed: true,\n")
  f.write("    padding: 10\n")
  f.write("  }\n")
  f.write("});\n")
  f.write("\n")
  f.write("}); // on dom ready\n")

  f.close()

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

### Generate the javascript to render the metabolic network ###
# param PATH : the path of the simulation folder
def generate_metabolic_network_js( PATH ):

  edges_filename = PATH+"/statistics/best_metabolic_edges.txt"
  nodes_filename = PATH+"/statistics/best_metabolic_nodes.txt"
  
  nodes = {}
  edges = []

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # A) Load nodes and edges         #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # s p km kcat delta_g e flux

  f = open(edges_filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    source  = int(l[0])
    product = int(l[1])
    flux    = float(l[6])
    if source not in nodes.keys():
      nodes[source] = 1
    if product not in nodes.keys():
      nodes[product] = 1
    present = False
    for i in range(len(edges)):
      if edges[i][0] == source and edges[i][1] == product:
        edges[i][2] += flux
        present = True
        break
    if not present:
      edges.append([source, product, flux])
    l = f.readline()
  f.close()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # B) Build the prime numbers list #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  prime_numbers = build_prime_numbers_list(1000)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # C) Write JS file                #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  f = open("./viewer/src/js/metabolic_network_graph.js", "w")

  #-----------------------------#
  # C.1) Write header           #
  #-----------------------------#
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

  #-----------------------------#
  # C.2) Create cytoscape class #
  #-----------------------------#
  f.write("$(function(){ // on dom ready\n")
  f.write("\n")
  f.write("var metabolic_network = cytoscape({\n")
  f.write("container: document.getElementById('metabolic_network_graph'),\n")
  f.write("\n")

  #-----------------------------#
  # C.3) Declare graph style    #
  #-----------------------------#
  f.write("style: cytoscape.stylesheet()\n")
  f.write("  .selector('node').css({\n")
  f.write("    'content': 'data(id)',\n")
  f.write("    'text-valign': 'center',\n")
  f.write("    'text-outline-width': 2,\n")
  f.write("    'shape': 'data(faveShape)',\n")
  f.write("    'text-outline-color': 'data(faveColor)',\n")
  f.write("    'background-color': 'data(faveColor)',\n")
  f.write("    'color': '#fff',\n")
  f.write("    'width': 30\n")
  f.write("  })\n")
  f.write("  .selector('$node > node').css({\n")
  f.write("    'text-outline-color': '#6FB1FC',\n")
  f.write("    'border-color': '#6FB1FC',\n")
  f.write("    'border-width': 2,\n")
  f.write("    'background-color': 'white',\n")
  f.write("    'color': 'white',\n")
  f.write("    'padding-top': '10px',\n")
  f.write("    'padding-left': '10px',\n")
  f.write("    'padding-bottom': '10px',\n")
  f.write("    'padding-right': '10px',\n")
  f.write("    'text-valign': 'top',\n")
  f.write("    'text-halign': 'center'\n")
  f.write("  })\n")
  f.write("  .selector('edge').css({\n")
  f.write("    'opacity': 0.666,\n")
  f.write("    'target-arrow-shape': 'triangle',\n")
  f.write("    'width': 2,\n")
  f.write("    'line-color': 'data(faveColor)',\n")
  f.write("    'source-arrow-color': 'data(faveColor)',\n")
  f.write("    'target-arrow-color': 'data(faveColor)'\n")
  f.write("  })\n")
  f.write("  .selector('edge.no_flux').css({\n")
  f.write("    'line-style': 'dashed',\n")
  f.write("    'line-color': 'grey',\n")
  f.write("    'source-arrow-color': 'grey',\n")
  f.write("    'target-arrow-color': 'grey'\n")
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

  f.write("    { data: { id: 'cell', faveColor: 'rgb(200, 200, 220)', faveShape: 'rectangle' } },\n")

  for i in range(len(nodes.keys())):
    tag  = nodes.keys()[i]
    conc = nodes[tag]
    line  = ""

    # If the tag is an inner metabolite #
    if tag > 0:
      if tag in prime_numbers:
        line += "    { data: { id: '"+str(tag)+"', parent: 'cell', faveColor: '#86B342', faveShape: 'rectangle' } }"
      else:
        line += "    { data: { id: '"+str(tag)+"', parent: 'cell', faveColor: '#6FB1FC', faveShape: 'ellipse' } }"

    # Else if the tag is an external metabolite #
    else:
      if tag in prime_numbers:
        line += "    { data: { id: '"+str(tag)+"', faveColor: '#86B342', faveShape: 'rectangle' } }"
      else:
        line += "    { data: { id: '"+str(tag)+"', faveColor: '#6FB1FC', faveShape: 'ellipse' } }"

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

  for i in range(len(edges)):
    edge    = edges[i]
    source  = edge[0]
    product = edge[1]
    flux    = edge[2]
    line  = ""
    if flux > 0.0:
      line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(source)+"', target: '"+str(product)+"', faveColor: 'black' } }"
    else:
      line += "    { data: { id: '"+str(edge[0])+"-"+str(edge[1])+"', weight: 1, source: '"+str(source)+"', target: '"+str(product)+"', faveColor: 'black' }, classes: 'no_flux' }"
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
  f.write("    directed: true,\n")
  f.write("    padding: 10,\n")
  f.write("    \n")
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

  # GENERATE GRN JS #
  generate_grn_js(PATH)

  # GENERATE METABOLIC NETWORK JS #
  generate_metabolic_network_js(PATH)



