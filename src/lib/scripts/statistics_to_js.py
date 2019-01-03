
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

### PRINT USAGE ####
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM STATISTICS ==="
  print "Usage: python statistics_to_js.py [parameters]"
  print "Parameters are:"
  print "-h, --help:"
  print "    Print this help, then exit."
  print "-p, --path:"
  print "    Give the path of the simulation folder."
  print ""

### READ COMMAND LINE ARGUMENTS ###
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

### RESET MEAN DICTIONARY ###
def reset_means( means ):
  for key in means:
    means[key] = 0.0

### COMPUTE MEAN ###
def compute_mean( means, period ):
  for key in means:
    means[key] /= period

### LOAD MEANS IN DATA ###
def load_means_in_data( means, data ):
  for key in means:
    data[key].append(means[key])

### GET FILE DATA ###
def get_file_data( filename, vector_length, requirements ):
  ##############################
  # 1) Get variable names      #
  ##############################
  DATA = {}
  RANK = {}
  MEAN = {}
  f = open(filename, "r")
  l = f.readline()
  l = l.strip("\n")
  l = l.split(" ")
  for i in range(len(l)):
    DATA[l[i]] = []
    RANK[i]    = l[i]
    MEAN[l[i]] = 0.0
  f.close()
  ##############################
  # 2) Get the number of lines #
  ##############################
  NB_LINES = 0
  f = open(filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    NB_LINES += 1
    l = f.readline()
  f.close()
  ##############################
  # 3) Compute the mean period #
  ##############################
  PERIOD = 1
  if NB_LINES > vector_length:
    PERIOD = NB_LINES/vector_length
  ##############################
  # 4) Load data               #
  ##############################
  COUNT = 0
  reset_means(MEAN)
  f = open(filename, "r")
  l = f.readline()
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    if len(l) == len(DATA):
      if COUNT < PERIOD:
        for i in range(len(l)):
          try:
            MEAN[RANK[i]] += float(l[i])
          except:
            sys.exit()
        COUNT += 1
      else:
        compute_mean(MEAN, PERIOD)
        load_means_in_data(MEAN, DATA)
        COUNT = 0
        reset_means(MEAN)
      l = f.readline()
    else:
      break
  f.close()
  DATA_LENGTH = len(DATA[DATA.keys()[0]])
  ##############################
  # 6) Compute additional data #
  ##############################
  if requirements == "DATA":
    return DATA
  elif requirements == "CUMULATED":
    SUM = {}
    CUMULATED = {}
    for key in DATA:
      SUM[key] = 0.0
      CUMULATED[key] = []
    for i in range(DATA_LENGTH):
      for key in DATA:
        SUM[key] += DATA[key][i]
        CUMULATED[key].append(SUM[key])
    return CUMULATED
  elif requirements == "REVERTED":
    REVERTED = {}
    for key in DATA:
      REVERTED[key] = []
    for i in range(DATA_LENGTH):
      for key in DATA:
        REVERTED[key].append(DATA[key][DATA_LENGTH-i-1])
    return REVERTED
  elif requirements == "CUMULATED_REVERTED":
    print "COUCOU COUCOU"
    SUM = {}
    CUMULATED = {}
    for key in DATA:
      SUM[key] = 0.0
      CUMULATED[key] = []
    for i in range(DATA_LENGTH):
      for key in DATA:
        SUM[key] += DATA[key][i]
        CUMULATED[key].append(SUM[key])
    REVERTED_CUMULATED = {}
    for key in DATA:
      REVERTED_CUMULATED[key] = []
    for i in range(DATA_LENGTH):
      for key in DATA:
        REVERTED_CUMULATED[key].append(CUMULATED[key][DATA_LENGTH-i-1])
    return REVERTED_CUMULATED
  else:
    print "ERROR: required option wrong ("+requirements+"). Exit."
    sys.exit()

### WRITE DYGRAPH JS FILE WITH SPECIFIED DATA ###
def write_dygraph_js( variables_list, js_filename, div_id, title, xlab, ylab, dygraph_options, dygraph_name, data ):
  DATA_LENGTH = len(data[data.keys()[0]])
  f = open(js_filename, "w")
  ############################
  # 1) Write header          #
  ############################
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
  ############################
  # 2) write data in js file #
  ############################
  f.write("var data = [];\n")
  for i in range(DATA_LENGTH):
    line = "data.push(["
    for key in variables_list:
      line += str(data[key][i])+", "
    line = line.strip(", ")
    line += "]);"
    f.write(line)
  f.write("\n")
  ############################
  # 3) write dygraph class   #
  ############################
  f.write(dygraph_name+" = new Dygraph(\n")
  f.write("document.getElementById(\""+div_id+"\"),\n")
  f.write("data,\n")
  f.write("{\n")
  f.write("title: '"+title+"',\n")
  f.write("xlabel: '"+xlab+"',\n")
  f.write("ylabel: '"+ylab+"',\n")
  line = "labels: ["
  for key in variables_list:
    line += "'"+key +"', "
  line = line.strip(", ")
  line += "],\n"
  f.write(line)
  ############################
  # 4) Custom options        #
  ############################
  for option in dygraph_options:
    f.write(option+",\n")
  ############################
  # 5) Common options        #
  ############################
  f.write("legend: 'always',\n")
  #f.write("rollPeriod: "+str(period)+",\n")
  #f.write("rollPeriod: 1,\n")
  #f.write("showRoller: true,\n")
  f.write("labelsDivStyles: { 'textAlign': 'right' },\n")
  f.write("animatedZooms: true}\n")
  f.write(");\n\n");
  f.close()

### WRITE POPULATION FIGURES ###
def write_population_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/phenotype_mean.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter  = 1
  ###################
  # POPULATION SIZE #
  ###################
  variables_list = ["t", "population_size"]
  js_filename    = path+"viewer/src/js/population_size.js"
  div_id         = "population_size_div"
  title          = "Population size"
  xlab           = "Simulation time"
  ylab           = "Population size"
  dygraph_name   = "population_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ###############
  # GROWTH RATE #
  ###############
  variables_list = ["t", "growth_rate"]
  js_filename    = path+"viewer/src/js/growth_rate.js"
  div_id         = "growth_rate_div"
  title          = "Population growth rate"
  xlab           = "Simulation time"
  ylab           = "Growth rate"
  dygraph_name   = "population_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1

### WRITE MEAN PHENOTYPE FIGURES ###
def write_mean_phenotype_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/phenotype_mean.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter = 1
  ###################
  # PROTEIN AMOUNTS #
  ###################
  variables_list = ["t", "inherited_TF_amount", "inherited_E_amount", "TF_amount", "E_amount"]
  js_filename    = path+"viewer/src/js/mean_protein_amounts.js"
  div_id         = "mean_protein_amounts_div"
  title          = "Protein amounts"
  xlab           = "Simulation time"
  ylab           = "Concentration"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  #####################
  # METABOLIC AMOUNTS #
  #####################
  variables_list = ["t", "inherited_metabolic_amount", "metabolic_amount"]
  js_filename    = path+"viewer/src/js/mean_metabolic_amounts.js"
  div_id         = "mean_metabolic_amounts_div"
  title          = "Metabolic amounts"
  xlab           = "Simulation time"
  ylab           = "Concentration"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  #########
  # SCORE #
  #########
  variables_list = ["t", "score"]
  js_filename    = path+"viewer/src/js/mean_score.js"
  div_id         = "mean_score_div"
  title          = "Score"
  xlab           = "Simulation time"
  ylab           = "Score"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ##########
  # ENERGY #
  ##########
  variables_list = ["t", "energy"]
  js_filename    = path+"viewer/src/js/mean_energy.js"
  div_id         = "mean_energy_div"
  title          = "Energy"
  xlab           = "Simulation time"
  ylab           = "Energy"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ############
  # LIFESPAN #
  ############
  variables_list = ["t", "lifespan"]
  js_filename    = path+"viewer/src/js/mean_lifespan.js"
  div_id         = "mean_lifespan_div"
  title          = "Lifespan"
  xlab           = "Simulation time"
  ylab           = "Lifespan (simulation steps)"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  #######################
  # NUMBER OF DIVISIONS #
  #######################
  variables_list = ["t", "number_of_divisions"]
  js_filename    = path+"viewer/src/js/mean_nb_divisions.js"
  div_id         = "mean_nb_divisions_div"
  title          = "Number of divisions"
  xlab           = "Simulation time"
  ylab           = "Number of divisions"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  #########################
  # TOXICITY ACCUMULATION #
  #########################
  variables_list = ["t", "toxicity"]
  js_filename    = path+"viewer/src/js/mean_toxicity.js"
  div_id         = "mean_toxicity_div"
  title          = "Toxicity accumulation"
  xlab           = "Simulation time"
  ylab           = "Toxicity accumulation"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ####################
  # CELL'S EXCHANGES #
  ####################
  variables_list = ["t", "metabolic_uptake", "metabolic_release"]
  js_filename    = path+"viewer/src/js/mean_exchanges.js"
  div_id         = "mean_exchanges_div"
  title          = "Exchanges with environment"
  xlab           = "Simulation time"
  ylab           = "Concentrations"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ################################
  # CELL'S METABOLIC GROWTH RATE #
  ################################
  variables_list = ["t", "metabolic_growth_rate", "Dmetabolic_growth_rate"]
  js_filename    = path+"viewer/src/js/mean_metabolic_growth_rate.js"
  div_id         = "mean_metabolic_growth_rate_div"
  title          = "Metabolic growth rate"
  xlab           = "Simulation time"
  ylab           = "Growth rate"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ##############################
  # CELL'S GRN NODES AND EDGES #
  ##############################
  variables_list = ["t", "grn_nb_nodes", "grn_nb_edges"]
  js_filename    = path+"viewer/src/js/mean_grn_nodes_edges.js"
  div_id         = "mean_grn_nodes_edges_div"
  title          = "Genetic regulation network"
  xlab           = "Simulation time"
  ylab           = "Number of elements"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ####################################
  # CELL'S METABOLIC NODES AND EDGES #
  ####################################
  variables_list = ["t", "metabolic_nb_nodes", "metabolic_nb_edges"]
  js_filename    = path+"viewer/src/js/mean_metabolic_nodes_edges.js"
  div_id         = "mean_metabolic_nodes_edges_div"
  title          = "Metabolic network"
  xlab           = "Simulation time"
  ylab           = "Number of elements"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  #########################
  # REGULATION REDUNDANCY #
  #########################
  variables_list = ["t", "regulation_redundancy"]
  js_filename    = path+"viewer/src/js/mean_regulation_redundancy.js"
  div_id         = "mean_regulation_redundancy_div"
  title          = "Regulation redundancy"
  xlab           = "Simulation time"
  ylab           = "Number of redundant units"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ########################
  # METABOLIC REDUNDANCY #
  ########################
  variables_list = ["t", "metabolic_redundancy"]
  js_filename    = path+"viewer/src/js/mean_metabolic_redundancy.js"
  div_id         = "mean_metabolic_redundancy_div"
  title          = "Metabolic redundancy"
  xlab           = "Simulation time"
  ylab           = "Number of redundant units"
  dygraph_name   = "mean_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1

### WRITE BEST LINEAGE PHENOTYPE FIGURES ###
def write_best_lineage_phenotype_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/phenotype_lineage_best.txt"
  REVERTED = get_file_data(filename, vector_length, "REVERTED")
  counter  = 1
  ###################
  # PROTEIN AMOUNTS #
  ###################
  variables_list = ["generation", "inherited_TF_amount", "inherited_E_amount", "TF_amount", "E_amount"]
  js_filename    = path+"viewer/src/js/best_lineage_protein_amounts.js"
  div_id         = "best_lineage_protein_amounts_div"
  title          = "Protein amounts"
  xlab           = "Generations"
  ylab           = "Concentration"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  #####################
  # METABOLIC AMOUNTS #
  #####################
  variables_list = ["generation", "inherited_metabolic_amount", "metabolic_amount"]
  js_filename    = path+"viewer/src/js/best_lineage_metabolic_amounts.js"
  div_id         = "best_lineage_metabolic_amounts_div"
  title          = "Metabolic concentrations"
  xlab           = "Generations"
  ylab           = "Concentration"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  #########
  # SCORE #
  #########
  variables_list = ["generation", "min_score", "mean_score", "max_score"]
  js_filename    = path+"viewer/src/js/best_lineage_score.js"
  div_id         = "best_lineage_score_div"
  title          = "Score"
  xlab           = "Generations"
  ylab           = "Score"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ##########
  # ENERGY #
  ##########
  variables_list = ["generation", "min_energy", "mean_energy", "max_energy"]
  js_filename    = path+"viewer/src/js/best_lineage_energy.js"
  div_id         = "best_lineage_energy_div"
  title          = "Energy"
  xlab           = "Generations"
  ylab           = "Energy"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ############
  # LIFESPAN #
  ############
  variables_list = ["generation", "lifespan"]
  js_filename    = path+"viewer/src/js/best_lineage_lifespan.js"
  div_id         = "best_lineage_lifespan_div"
  title          = "Lifespan"
  xlab           = "Generations"
  ylab           = "Lifespan (simulation steps)"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  #######################
  # NUMBER OF DIVISIONS #
  #######################
  variables_list = ["generation", "number_of_divisions"]
  js_filename    = path+"viewer/src/js/best_lineage_nb_divisions.js"
  div_id         = "best_lineage_nb_divisions_div"
  title          = "Number of divisions"
  xlab           = "Generations"
  ylab           = "Number of divisions"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  #########################
  # TOXICITY ACCUMULATION #
  #########################
  variables_list = ["generation", "toxicity"]
  js_filename    = path+"viewer/src/js/best_lineage_toxicity.js"
  div_id         = "best_lineage_toxicity_div"
  title          = "Toxicity accumulation"
  xlab           = "Generations"
  ylab           = "Toxicity accumulation"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ####################
  # CELL'S EXCHANGES #
  ####################
  variables_list = ["generation", "metabolic_uptake", "metabolic_release"]
  js_filename    = path+"viewer/src/js/best_lineage_exchanges.js"
  div_id         = "best_lineage_exchanges_div"
  title          = "Exchanges with environment"
  xlab           = "Generations"
  ylab           = "Concentrations"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  ################################
  # CELL'S METABOLIC GROWTH RATE #
  ################################
  variables_list = ["generation", "metabolic_growth_rate"]
  js_filename    = path+"viewer/src/js/best_lineage_metabolic_growth_rate.js"
  div_id         = "best_lineage_metabolic_growth_rate_div"
  title          = "Metabolic growth rate"
  xlab           = "Generations"
  ylab           = "Growth rate"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ####################################
  # CELL'S METABOLIC NODES AND EDGES #
  ####################################
  variables_list = ["generation", "grn_nb_nodes", "grn_nb_edges"]
  js_filename    = path+"viewer/src/js/best_lineage_grn_nodes_edges.js"
  div_id         = "best_lineage_grn_nodes_edges_div"
  title          = "Genetic regulation network structure"
  xlab           = "Generations"
  ylab           = "Number of elements"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ####################################
  # CELL'S METABOLIC NODES AND EDGES #
  ####################################
  variables_list = ["generation", "metabolic_nb_nodes", "metabolic_nb_edges"]
  js_filename    = path+"viewer/src/js/best_lineage_metabolic_nodes_edges.js"
  div_id         = "best_lineage_metabolic_nodes_edges_div"
  title          = "Metabolic network structure"
  xlab           = "Generations"
  ylab           = "Number of elements"
  dygraph_name   = "best_lineage_phenotype_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1

### WRITE MEAN GENOME STRUCTURE FIGURES ###
def write_mean_genome_structure_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/genome_structure_mean.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter  = 1
  ###############
  # GENOME SIZE #
  ###############
  variables_list = ["t", "genome_size", "functional_size"]
  js_filename    = path+"viewer/src/js/mean_genome_size.js"
  div_id         = "mean_genome_size_div"
  title          = "Genome size"
  xlab           = "Simulation time"
  ylab           = "Number of genetic units"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ###########################
  # NUMBER OF GENETIC UNITS #
  ###########################
  variables_list = ["t", "nb_NC", "nb_E", "nb_TF", "nb_BS", "nb_P"]
  js_filename    = path+"viewer/src/js/mean_nb_genetic_units.js"
  div_id         = "mean_nb_genetic_units_div"
  title          = "Genome composition"
  xlab           = "Simulation time"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ####################################
  # ENZYME GENETIC UNITS COMPOSITION #
  ####################################
  variables_list = ["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"]
  js_filename    = path+"viewer/src/js/mean_E_composition.js"
  div_id         = "mean_E_composition_div"
  title          = "Enzyme unit classes"
  xlab           = "Simulation time"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ####################################
  # NUMBER OF ENHANCERS OR OPERATORS #
  ####################################
  variables_list = ["t", "nb_functional_regions", "nb_enhancers", "nb_operators"]
  js_filename    = path+"viewer/src/js/mean_nb_enhancers_operators.js"
  div_id         = "mean_nb_enhancers_operators_div"
  title          = "Number of enhancers and operators"
  xlab           = "Simulation time"
  ylab           = "Number of sites"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ##################################
  # SIZE OF ENHANCERS OR OPERATORS #
  ##################################
  variables_list = ["t", "enhancer_size", "operator_size"]
  js_filename    = path+"viewer/src/js/mean_size_enhancers_operators.js"
  div_id         = "mean_size_enhancers_operators_div"
  title          = "Size of enhancers and operators"
  xlab           = "Simulation time"
  ylab           = "Size"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ####################
  # TYPES OF REGIONS #
  ####################
  variables_list = ["t", "nb_E_regions", "nb_TF_regions", "nb_mixed_regions"]
  js_filename    = path+"viewer/src/js/mean_region_types.js"
  div_id         = "mean_region_types_div"
  title          = "Types of functional regions"
  xlab           = "Simulation time"
  ylab           = "Number by type"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ###########################
  # SIZE OF REGIONS BY TYPE #
  ###########################
  variables_list = ["t", "E_region_size", "TF_region_size", "mixed_region_size"]
  js_filename    = path+"viewer/src/js/mean_region_size.js"
  div_id         = "mean_region_size_div"
  title          = "Size of functional regions"
  xlab           = "Simulation time"
  ylab           = "Mean size by type"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1
  ###########################
  # SIZE OF OPERONS BY TYPE #
  ###########################
  variables_list = ["t", "E_operon_size", "TF_operon_size", "mixed_operon_size"]
  js_filename    = path+"viewer/src/js/mean_operon_size.js"
  div_id         = "mean_operon_size_div"
  title          = "Size of operons"
  xlab           = "Simulation time"
  ylab           = "Mean size by type"
  dygraph_name   = "mean_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1

### WRITE BEST LINEAGE GENOME STRUCTURE FIGURES ###
def write_best_lineage_genome_structure_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/genome_structure_lineage_best.txt"
  REVERTED = get_file_data(filename, vector_length, "REVERTED")
  counter  = 1
  ###############
  # GENOME SIZE #
  ###############
  variables_list = ["generation", "genome_size", "genome_functional_size"]
  js_filename    = path+"viewer/src/js/best_lineage_genome_size.js"
  div_id         = "best_lineage_genome_size_div"
  title          = "Genome size"
  xlab           = "Generations"
  ylab           = "Number of genetic units"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ###########################
  # NUMBER OF GENETIC UNITS #
  ###########################
  variables_list = ["generation", "genome_nb_NC", "genome_nb_E", "genome_nb_TF", "genome_nb_BS", "genome_nb_P"]
  js_filename    = path+"viewer/src/js/best_lineage_nb_genetic_units.js"
  div_id         = "best_lineage_nb_genetic_units_div"
  title          = "Genome composition"
  xlab           = "Generations"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  genome_option = ["stackedGraph: true"]
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  ####################################
  # ENZYME GENETIC UNITS COMPOSITION #
  ####################################
  variables_list = ["generation", "genome_nb_inner_enzymes", "genome_nb_inflow_pumps", "genome_nb_outflow_pumps"]
  js_filename    = path+"viewer/src/js/best_lineage_E_composition.js"
  div_id         = "best_lineage_E_composition_div"
  title          = "Enzyme unit classes"
  xlab           = "Generations"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  #########################
  # REGULATION REDUNDANCY #
  #########################
  variables_list = ["generation", "mean_regulation_redundancy"]
  js_filename    = path+"viewer/src/js/best_lineage_regulation_redundancy.js"
  div_id         = "best_lineage_regulation_redundancy_div"
  title          = "Regulation redundancy"
  xlab           = "Generations"
  ylab           = "Number of redundant genetic units"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ########################
  # METABOLIC REDUNDANCY #
  ########################
  variables_list = ["generation", "mean_metabolic_redundancy"]
  js_filename    = path+"viewer/src/js/best_lineage_metabolic_redundancy.js"
  div_id         = "best_lineage_metabolic_redundancy_div"
  title          = "Metabolic redundancy"
  xlab           = "Generations"
  ylab           = "Number of redundant genetic units"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ####################################
  # NUMBER OF ENHANCERS OR OPERATORS #
  ####################################
  variables_list = ["generation", "nb_functional_regions", "nb_enhancers", "nb_operators"]
  js_filename    = path+"viewer/src/js/best_lineage_nb_enhancers_operators.js"
  div_id         = "best_lineage_nb_enhancers_operators_div"
  title          = "Number of enhancers and operators"
  xlab           = "Generations"
  ylab           = "Number of sites"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ##################################
  # SIZE OF ENHANCERS OR OPERATORS #
  ##################################
  variables_list = ["generation", "mean_enhancer_size", "mean_operator_size"]
  js_filename    = path+"viewer/src/js/best_lineage_size_enhancers_operators.js"
  div_id         = "best_lineage_size_enhancers_operators_div"
  title          = "Size of enhancers and operators"
  xlab           = "Generations"
  ylab           = "Size"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ####################
  # TYPES OF REGIONS #
  ####################
  variables_list = ["generation", "nb_E_regions", "nb_TF_regions", "nb_mixed_regions"]
  js_filename    = path+"viewer/src/js/best_lineage_region_types.js"
  div_id         = "best_lineage_region_types_div"
  title          = "Types of functional regions"
  xlab           = "Generations"
  ylab           = "Number by type"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  ###########################
  # SIZE OF REGIONS BY TYPE #
  ###########################
  variables_list = ["generation", "mean_E_region_size", "mean_TF_region_size", "mean_mixed_region_size"]
  js_filename    = path+"viewer/src/js/best_lineage_region_size.js"
  div_id         = "best_lineage_region_size_div"
  title          = "Size of functional regions"
  xlab           = "Generations"
  ylab           = "Mean size by type"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1
  ###########################
  # SIZE OF OPERONS BY TYPE #
  ###########################
  variables_list = ["generation", "mean_E_operon_size", "mean_TF_operon_size", "mean_mixed_operon_size"]
  js_filename    = path+"viewer/src/js/best_lineage_operon_size.js"
  div_id         = "best_lineage_operon_size_div"
  title          = "Size of operons"
  xlab           = "Generations"
  ylab           = "Mean size by type"
  dygraph_name   = "best_lineage_genome_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, REVERTED)
  counter += 1

### WRITE MEAN INHERITED PROTEINS FIGURES ###
def write_mean_inherited_proteins_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/inherited_proteins_mean.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter  = 1
  ################################
  # NUMBER OF INHERITED PROTEINS #
  ################################
  variables_list = ["t", "nb_E", "nb_TF"]
  js_filename    = path+"viewer/src/js/mean_number_inherited_proteins.js"
  div_id         = "mean_number_inherited_proteins_div"
  title          = "Number of inherited proteins"
  xlab           = "Simulation time"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "mean_inherited_proteins_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ########################
  # ENZYME TUPLE CLASSES #
  ########################
  variables_list = ["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"]
  js_filename    = path+"viewer/src/js/mean_inherited_E_composition.js"
  div_id         = "mean_inherited_E_composition_div"
  title          = "Inherited enzyme unit classes"
  xlab           = "Simulation time"
  ylab           = "Number of genetic units by class"
  dygraph_name   = "mean_inherited_proteins_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1

### WRITE BEST LINEAGE INHERITED PROTEINS FIGURES ###
def write_best_lineage_inherited_proteins_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/inherited_proteins_lineage_best.txt"
  REVERTED = get_file_data(filename, vector_length, "REVERTED")
  counter  = 1
  ################################
  # NUMBER OF INHERITED PROTEINS #
  ################################
  variables_list = ["generation", "inherited_nb_E", "inherited_nb_TF"]
  js_filename    = path+"viewer/src/js/best_lineage_number_inherited_proteins.js"
  div_id         = "best_lineage_number_inherited_proteins_div"
  title          = "Number of inherited proteins"
  xlab           = "Generations"
  ylab           = "Number of genetic units by type"
  dygraph_name   = "best_lineage_inherited_proteins_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1
  ########################
  # ENZYME TUPLE CLASSES #
  ########################
  variables_list = ["generation", "inherited_nb_inner_enzymes", "inherited_nb_inflow_pumps", "inherited_nb_outflow_pumps"]
  js_filename    = path+"viewer/src/js/best_lineage_inherited_E_composition.js"
  div_id         = "best_lineage_inherited_E_composition_div"
  title          = "Inherited enzyme unit classes"
  xlab           = "Generations"
  ylab           = "Number of genetic units by class"
  dygraph_name   = "best_lineage_inherited_proteins_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, REVERTED)
  counter += 1

### WRITE BEST LINEAGE FIXED MUTATIONS FIGURES ###
def write_best_lineage_fixed_mutations_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/fixed_mutations_lineage_best.txt"
  CUMULATED = get_file_data(filename, vector_length, "CUMULATED")
  counter  = 1
  ###################
  # POINT MUTATIONS #
  ###################
  variables_list = ["generation", "nb_NC_point_mutations", "nb_E_point_mutations", "nb_TF_point_mutations", "nb_BS_point_mutations", "nb_P_point_mutations"]
  js_filename    = path+"viewer/src/js/best_lineage_point_mutations.js"
  div_id         = "best_lineage_point_mutations_div"
  title          = "Point mutations"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #######################
  # NC UNIT TRANSITIONS #
  #######################
  variables_list = ["generation", "nb_NC_to_E_transitions", "nb_NC_to_TF_transitions", "nb_NC_to_BS_transitions", "nb_NC_to_P_transitions"]
  js_filename    = path+"viewer/src/js/best_lineage_NC_transitions.js"
  div_id         = "best_lineage_NC_transitions_div"
  title          = "NC units transitions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  ######################
  # E UNIT TRANSITIONS #
  ######################
  variables_list = ["generation", "nb_E_to_NC_transitions", "nb_E_to_TF_transitions", "nb_E_to_BS_transitions", "nb_E_to_P_transitions"]
  js_filename    = path+"viewer/src/js/best_lineage_E_transitions.js"
  div_id         = "best_lineage_E_transitions_div"
  title          = "E units transitions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #######################
  # TF UNIT TRANSITIONS #
  #######################
  variables_list = ["generation", "nb_TF_to_NC_transitions", "nb_TF_to_E_transitions", "nb_TF_to_BS_transitions", "nb_TF_to_P_transitions"]
  js_filename    = path+"viewer/src/js/best_lineage_TF_transitions.js"
  div_id         = "best_lineage_TF_transitions_div"
  title          = "TF units transitions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, CUMULATED)
  counter += 1
  #######################
  # BS UNIT TRANSITIONS #
  #######################
  variables_list = ["generation", "nb_BS_to_NC_transitions", "nb_BS_to_E_transitions", "nb_BS_to_TF_transitions", "nb_BS_to_P_transitions"]
  js_filename    = path+"viewer/src/js/best_lineage_BS_transitions.js"
  div_id         = "best_lineage_BS_transitions_div"
  title          = "BS units transitions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  ######################
  # P UNIT TRANSITIONS #
  ######################
  variables_list = ["generation", "nb_P_to_NC_transitions", "nb_P_to_E_transitions", "nb_P_to_TF_transitions", "nb_P_to_BS_transitions"]
  js_filename    = path+"viewer/src/js/best_lineage_P_transitions.js"
  div_id         = "best_lineage_P_transitions_div"
  title          = "P units transitions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #################################
  # METABOLIC SPACE MUTATION SIZE #
  #################################
  variables_list = ["generation", "s_mutation_size", "p_mutation_size"]
  js_filename    = path+"viewer/src/js/best_lineage_metabolic_space_mutation_size.js"
  div_id         = "best_lineage_metabolic_space_mutation_size_div"
  title          = "Metabolic space mutation size"
  xlab           = "Generations"
  ylab           = "Mutation size"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  ###############################
  # RATE CONSTANT MUTATION SIZE #
  ###############################
  variables_list = ["generation", "kcat_km_ratio_mutation_size", "kcat_mutation_size", "BS_tag_mutation_size", "coE_tag_mutation_size", "beta_mutation_size"]
  js_filename    = path+"viewer/src/js/best_lineage_rate_constant_mutation_size.js"
  div_id         = "best_lineage_rate_constant_mutation_size_div"
  title          = "Rate constants mutation size"
  xlab           = "Generations"
  ylab           = "Mutation size"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  ################
  # DUPLICATIONS #
  ################
  variables_list = ["generation", "nb_duplicated_NC", "nb_duplicated_E", "nb_duplicated_TF", "nb_duplicated_BS", "nb_duplicated_P"]
  js_filename    = path+"viewer/src/js/best_lineage_duplications.js"
  div_id         = "best_lineage_duplications_div"
  title          = "Duplications"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #############
  # DELETIONS #
  #############
  variables_list = ["generation", "nb_deleted_NC", "nb_deleted_E", "nb_deleted_TF", "nb_deleted_BS", "nb_deleted_P"]
  js_filename    = path+"viewer/src/js/best_lineage_deletions.js"
  div_id         = "best_lineage_deletions_div"
  title          = "Deletions"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  ############################
  # NUMBER OF REARRANGEMENTS #
  ############################
  variables_list = ["generation", "nb_duplications", "nb_deletions", "nb_translocations", "nb_inversions"]
  js_filename    = path+"viewer/src/js/best_lineage_nb_rearrangements.js"
  div_id         = "best_lineage_nb_rearrangements_div"
  title          = "Number of rearrangements"
  xlab           = "Generations"
  ylab           = "Number of events"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #######################
  # REARRANGEMENTS SIZE #
  #######################
  variables_list = ["generation", "duplication_size", "deletion_size", "translocation_size", "inversion_size"]
  js_filename    = path+"viewer/src/js/best_lineage_rearrangement_size.js"
  div_id         = "best_lineage_rearrangement_size_div"
  title          = "Rearrangements size"
  xlab           = "Generations"
  ylab           = "Rearrangements size"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1
  #######################################
  # NUMBER OF HORIZONTAL GENE TRANSFERS #
  #######################################
  variables_list = ["generation", "nb_NC_HGT", "nb_E_HGT", "nb_TF_HGT", "nb_BS_HGT", "nb_P_HGT"]
  js_filename    = path+"viewer/src/js/best_lineage_nb_hgt.js"
  div_id         = "best_lineage_nb_hgt_div"
  title          = "Number of HGT"
  xlab           = "Generations"
  ylab           = "Number of HGT"
  dygraph_name   = "best_lineage_fixed_mutations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, CUMULATED)
  counter += 1

### WRITE TROPHIC NETWORK FIGURES ###
def write_trophic_network_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/trophic_network_profile.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter  = 1
  #############################################
  # NUMBER OF GROUPS BY TROPHIC NETWORK LEVEL #
  #############################################
  variables_list = ["t", "level0_nb_groups", "level1_nb_groups", "level2_nb_groups", "nolevel_nb_groups"]
  js_filename    = path+"viewer/src/js/trophic_network_nb_groups.js"
  div_id         = "trophic_network_nb_groups_div"
  title          = "Number of groups per trophic level"
  xlab           = "Simulation time"
  ylab           = "Number of groups"
  dygraph_name   = "trophic_network_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ############################################
  # NUMBER OF CELLS BY TROPHIC NETWORK LEVEL #
  ############################################
  variables_list = ["t", "level0_nb_cells", "level1_nb_cells", "level2_nb_cells", "nolevel_nb_cells"]
  js_filename    = path+"viewer/src/js/trophic_network_nb_cells.js"
  div_id         = "trophic_network_nb_cells_div"
  title          = "Number of cells per trophic level"
  xlab           = "Simulation time"
  ylab           = "Number of cells"
  dygraph_name   = "trophic_network_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ###############################################
  # NUMBER OF GROUP APPEARANCES AND EXTINCTIONS #
  ###############################################
  variables_list = ["t", "nb_group_appearances", "nb_group_extinctions"]
  js_filename    = path+"viewer/src/js/trophic_network_group_events.js"
  div_id         = "trophic_network_group_events_div"
  title          = "Number of trophic group new events"
  xlab           = "Simulation time"
  ylab           = "Number of groups"
  dygraph_name   = "trophic_network_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  #######################
  # MEAN GROUP LIFESPAN #
  #######################
  variables_list = ["t", "mean_group_lifespan"]
  js_filename    = path+"viewer/src/js/trophic_network_group_lifespan.js"
  div_id         = "trophic_network_group_lifespan_div"
  title          = "Mean trophic group lifespan"
  xlab           = "Simulation time"
  ylab           = "Lifespan"
  dygraph_name   = "trophic_network_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1

### WRITE TREE STRUCTURE FIGURES ###
def write_tree_structure_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/tree_structure.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter = 1
  #######################################
  # NUMBER OF NODES IN THE LINEAGE TREE #
  #######################################
  variables_list = ["t", "lineage_nb_nodes"]
  js_filename    = path+"viewer/src/js/lineage_nb_nodes.js"
  div_id         = "lineage_nb_nodes_div"
  title          = "Number of nodes in the lineage tree"
  xlab           = "Simulation time"
  ylab           = "Number of nodes"
  dygraph_name   = "tree_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ############################################
  # NUMBER OF NODES IN THE PHYLOGENETIC TREE #
  ############################################
  variables_list = ["t", "phylogeny_nb_nodes"]
  js_filename    = path+"viewer/src/js/phylogeny_nb_nodes.js"
  div_id         = "phylogeny_nb_nodes_div"
  title          = "Number of nodes in the phylogenetic tree"
  xlab           = "Simulation time"
  ylab           = "Number of nodes"
  dygraph_name   = "tree_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  #######################
  # COMMON ANCESTOR AGE #
  #######################
  variables_list = ["t", "phylogeny_ca_age"]
  js_filename    = path+"viewer/src/js/common_ancestor_age.js"
  div_id         = "common_ancestor_age_div"
  title          = "Common ancestor age"
  xlab           = "Simulation time"
  ylab           = "Age"
  dygraph_name   = "tree_structure_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1

### WRITE GLOBAL CONCENTRATIONS FIGURES ###
def write_global_concentrations_figures( path, stacked_option, filled_option ):
  filename = path+"statistics/global_concentrations.txt"
  DATA = get_file_data(filename, vector_length, "DATA")
  counter = 1
  ###########################################
  # REPARTITION OF THE MATTER IN THE SYSTEM #
  ###########################################
  variables_list = ["t", "matter_pop", "matter_env"]
  js_filename    = path+"viewer/src/js/matter_repartition.js"
  div_id         = "matter_repartition_div"
  title          = "Repartition of matter in the world"
  xlab           = "Simulation time"
  ylab           = "Amounts"
  dygraph_name   = "global_concentrations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, stacked_option, dygraph_name, DATA)
  counter += 1
  ################
  # FLUXES BILAN #
  ################
  variables_list = ["t", "matter_inflow", "matter_outflow"]
  js_filename    = path+"viewer/src/js/fluxes_bilan.js"
  div_id         = "fluxes_bilan_div"
  title          = "Bilan of matter fluxes"
  xlab           = "Simulation time"
  ylab           = "Amounts"
  dygraph_name   = "global_concentrations_"+str(counter)
  write_dygraph_js(variables_list, js_filename, div_id, title, xlab, ylab, filled_option, dygraph_name, DATA)
  counter += 1

############
#   MAIN   #
############
if __name__ == '__main__':
  
  PATH = read_args(sys.argv)

  #-------------------------------------------------------#
  # 1) DEFINE OPTIONS                                     #
  #-------------------------------------------------------#
  stacked_option = ["stackedGraph: true"]
  filled_option  = ["fillGraph: true"]
  vector_length  = 500

  # POPULATION FIGURES ---------------------------------------------------------------------#
  write_population_figures(PATH, stacked_option, filled_option)

  # PHENOTYPE FIGURES ----------------------------------------------------------------------#
  write_mean_phenotype_figures(PATH, stacked_option, filled_option)
  write_best_lineage_phenotype_figures(PATH, stacked_option, filled_option)

  # GENOME STRUCTURE FIGURES ---------------------------------------------------------------#
  write_mean_genome_structure_figures(PATH, stacked_option, filled_option)
  write_best_lineage_genome_structure_figures(PATH, stacked_option, filled_option)

  # INHERITED PROTEINS FIGURES -------------------------------------------------------------#
  write_mean_inherited_proteins_figures(PATH, stacked_option, filled_option)
  write_best_lineage_inherited_proteins_figures(PATH, stacked_option, filled_option)
  
  # FIXED MUTATIONS FIGURES ----------------------------------------------------------------#
  write_best_lineage_fixed_mutations_figures(PATH, stacked_option, filled_option)

  # TROPHIC NETWORK FIGURES ----------------------------------------------------------------#
  write_trophic_network_figures(PATH, stacked_option, filled_option)

  # TREE STRUCTURE FIGURES -----------------------------------------------------------------#
  write_tree_structure_figures(PATH, stacked_option, filled_option)

  # GLOBAL CONCENTRATIONS FIGURES ----------------------------------------------------------#
  write_global_concentrations_figures(PATH, stacked_option, filled_option)

