
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

#-------------------------------------#
# 4) Open data                        #
#-------------------------------------#
tree <- read.tree("./statistics/phylogenetic_tree.phb")

#########################################################
#     PHYLOGENETIC TREE IN PHYLOGRAM REPRESENTATION     #
#########################################################
png("./figures/phylogenetic_tree_phylogram.png", width=1200, height=1200, pointsize=20)
plot(tree, type="phylogram", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", main="Phylogenetic tree", edge.width=3)
add.scale.bar(cex=1, font=2, length=1000)
dev.off()

###################################################
#     PHYLOGENETIC TREE IN FAN REPRESENTATION     #
###################################################
png("./figures/phylogenetic_tree_fan.png", width=1200, height=1200, pointsize=20)
plot(tree, type="fan", show.tip.label=F, cex=0.5, label.offset=100, root.edge=T, direction="rightward", main="Phylogenetic tree", edge.width=3)
add.scale.bar(cex=1, font=2, length=1000)
dev.off()


