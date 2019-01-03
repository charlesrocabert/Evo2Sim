
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

### get index of the value in the vector ###
get_index <- function(vector, value)
{
  for (i in seq(1, length(vector)))
  {
    if (vector[i] == value)
    {
      return(i)
    }
  }
  return(-1)
}

### get index of the edge ending by value ###
get_edge_index <- function(edges, value)
{
  for (i in seq(1, length(edges[,1])))
  {
    if (edges[i,2] == value)
    {
      return(i)
    }
  }
  return(-1)
}

### get a specific color from the trophic level ###
get_color_from_trophic_level <- function(x)
{
  if (x == 0)
  {
    return("purple")
  }
  else if (x == 1)
  {
    return("cornflowerblue")
  }
  else if (x == 2)
  {
    return("darkolivegreen3")
  }
  else if (x == 4)
  {
    return("lightgrey")
  }
  else
  {
    return("black")
  }
}

#-------------------------------------#
# 1) Read command line arguments      #
#-------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0)
{
  stop("You must give the experiment path")
}
experiment_path = args[1]
setwd(experiment_path)

silence_warnings = F
if (length(args) == 2)
{
  if (args[2] == "-no-warnings")
  {
    silence_warnings = T
  }
}

#-------------------------------------#
# 2) Hide error messages and warnings #
#-------------------------------------#
if (silence_warnings == T)
{
  options(show.error.messages = FALSE, warn=-1)
}

#-------------------------------------#
# 3) Load libraries                   #
#-------------------------------------#
library("ape")
library("RColorBrewer")

#-------------------------------------#
# 4) Open phylogenetic tree and data  #
#-------------------------------------#
input_tree = read.tree("./statistics/phylogenetic_tree.phb")
phyData    = read.table("./statistics/phylogenetic_tree_trophic_data.txt", h=T, sep=" ")

### Load variables ###
variables = c("trophic_level")

#####################################################
# Plot colored phylogenetic tree for each character #
#####################################################
for (variable in variables)
{
  ###############################
  # If there are multiple roots #
  ###############################
  if (length(input_tree) != 6)
  {
    print("Evaluating multiple trees ...")
    for (tree_number in seq(1, length(input_tree)))
    {
      tree = input_tree[[tree_number]]
      
      ### Get the number of edges and create color palette ###
      Nedges = length(tree$edge[,1])
      
      ### Get data ###
      data = data.frame(phyData[,1], phyData[,variable])
      Nnodes = length(data[,1])
      
      ### Build colors vector to color tree edges ###
      colors = rep("black", Nedges)
      for (i in seq(1, Nnodes))
      {
        tip_index = get_index(tree$tip, data[i,1])
        node_index = get_index(tree$node, data[i,1])
        if (tip_index != -1 && node_index == -1)
        {
          edge_index = get_edge_index(tree$edge, tip_index)
          if (edge_index != -1)
          {
            colors[edge_index] = get_color_from_trophic_level(data[i,2])
          }
        }
        else if (tip_index == -1 && node_index != -1)
        {
          edge_index = get_edge_index(tree$edge, node_index+length(tree$tip))
          if (edge_index != -1)
          {
            colors[edge_index] = get_color_from_trophic_level(data[i,2])
          }
        }
      }
      
      ### Make figure ###
      
      #########################################################
      #     PHYLOGENETIC TREE IN PHYLOGRAM REPRESENTATION     #
      #########################################################
      png("./figures/phylogenetic_tree_phylogram_trophic_level.png", width=1200, height=1200, pointsize=20)
      plot(tree, type="phylogram", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", edge.color=colors, main=paste("Phylogenetic tree with ", variable, sep=""), edge.width=3)
      add.scale.bar(cex=1, font=2, length=1000)
      dev.off()
      
      ###################################################
      #     PHYLOGENETIC TREE IN FAN REPRESENTATION     #
      ###################################################
      png("./figures/phylogenetic_tree_fan_trophic_level.png", width=1200, height=1200, pointsize=20)
      plot(tree, type="fan", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", edge.color=colors, main=paste("Phylogenetic tree with ", variable, sep=""), edge.width=3)
      add.scale.bar(cex=1, font=2, length=1000)
      dev.off()
    }
  } else {
    print("Evaluating one tree for :")
    print(variable)
    #############################
    # Else if there is one root #
    #############################
    tree = input_tree
    ### Get the number of edges and create color palette ###
    Nedges = length(tree$edge[,1])
    
    ### Get data ###
    data = data.frame(phyData[,1], phyData[,variable])
    Nnodes = length(data[,1])
    
    ### Build colors vector to color tree edges ###
    colors = rep("black", Nedges)
    for (i in seq(1, Nnodes))
    {
      tip_index = get_index(tree$tip, data[i,1])
      node_index = get_index(tree$node, data[i,1])
      if (tip_index != -1 && node_index == -1)
      {
        edge_index = get_edge_index(tree$edge, tip_index)
        if (edge_index != -1)
        {
          colors[edge_index] = get_color_from_trophic_level(data[i,2])
        }
      }
      else if (tip_index == -1 && node_index != -1)
      {
        edge_index = get_edge_index(tree$edge, node_index+length(tree$tip))
        if (edge_index != -1)
        {
          colors[edge_index] = get_color_from_trophic_level(data[i,2])
        }
      }
    }
    
    ### Make figure ###
    
    #########################################################
    #     PHYLOGENETIC TREE IN PHYLOGRAM REPRESENTATION     #
    #########################################################
    png("./figures/phylogenetic_tree_phylogram_trophic_level.png", width=1200, height=1200, pointsize=20)
    plot(tree, type="phylogram", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", edge.color=colors, main=paste("Phylogenetic tree with ", variable, sep=""), edge.width=3)
    add.scale.bar(cex=1, font=2, length=1000)
    dev.off()
    
    ###################################################
    #     PHYLOGENETIC TREE IN FAN REPRESENTATION     #
    ###################################################
    png("./figures/phylogenetic_tree_fan_trophic_level.png", width=1200, height=1200, pointsize=20)
    plot(tree, type="fan", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", edge.color=colors, main=paste("Phylogenetic tree with ", variable, sep=""), edge.width=3)
    add.scale.bar(cex=1, font=2, length=1000)
    dev.off()
  }
}
