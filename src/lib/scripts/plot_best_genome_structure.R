
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
# 3) Open data                        #
#-------------------------------------#
data = read.table("./statistics/best_genome_structure.txt", h=T, sep=" ")

#######################################
#     CODING/NON CODING PIE CHART     #
#######################################
png("./figures/coding_piechart.png", width=800, height=800, pointsize=20)
par(mfrow=c(1,1), mar=c(5,5,5,5))

L = 1

# Pie Chart from data frame with Appended Sample Sizes
piedata = c(data$nb_NC[L], data$nb_E[L], data$nb_TF[L], data$nb_BS[L], data$nb_P[L])
pielabs = paste(c("Non coding", "Enzymes", "Transcription factors", "Binding sites", "Promoters"), "\n", piedata, sep="")
pie(piedata, labels=pielabs, main="Genome composition", col=c("lightgrey", "tomato", "purple", "darkolivegreen3", "cornflowerblue"))

dev.off()

##############################################
#     ENZYME UNITS COMPOSITION PIE CHART     #
##############################################
png("./figures/enzyme_units_composition_piechart.png", width=800, height=800, pointsize=20)
par(mfrow=c(1,1), mar=c(5,5,5,5))

# Pie Chart from data frame with Appended Sample Sizes
piedata = c(data$nb_inner_enzymes[L], data$nb_inflow_pumps[L], data$nb_outflow_pumps[L])
pielabs = paste(c("Inner enzymes", "Inflowing pumps", "Outflowing pumps"), "\n", piedata, sep="")
pie(piedata, labels=pielabs, main="Types of enzymes", col=c("tomato", "darkolivegreen3", "cornflowerblue"))

dev.off()
