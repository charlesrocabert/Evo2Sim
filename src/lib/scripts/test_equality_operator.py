
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

### Detect errors of equality operator usage ('=' instead of '==') ###
def detect_equality_operator_errors( file_to_explore ):
	f = open(file_to_explore, "r")
	l = f.readline()
	count = 1
	while l:
		l = l.strip("\n")
		l = l.strip("\t")
		l = l.strip(" ")
		if l.startswith("if") or l.startswith("else if") or l.startswith("else") or l.startswith("while") or l.startswith("assert"):
			if l.find(" = ") != -1:
				print "Equality operator error find in file "+file_to_explore+" at line "+str(count)
				print "code :'"+l+"'"
				sys.exit()
		l = f.readline()
	f.close()


############
#   MAIN   #
############
if __name__ == '__main__':
	files_to_explore = [
	"./src/create.cpp",
	"./src/bootstrap.cpp",
	"./src/run.cpp",
  "./src/generate_figures.cpp",
	"./src/recover_parameters.cpp",
	"./src/unitary_tests.cpp",
  "./src/integrated_tests.cpp",
  "./src/AB_post_treatments/AB_get_deepest_tree.cpp",
  "./src/AB_post_treatments/AB_trophic_groups_history.cpp",
  "./src/AB_post_treatments/AB_black_queen.cpp",
  "./src/AB_post_treatments/AB_create_competition_experiment.cpp",
  "./src/AB_post_treatments/AB_cross_feeding.cpp",
  "./src/AB_post_treatments/AB_frequency_dependence.cpp",
  "./src/AB_post_treatments/AB_heatmap.cpp",
  "./src/AB_post_treatments/AB_phylogenetic_diversity.cpp",
  "./src/AB_post_treatments/AB_recover_genome_statistics.cpp",
  "./src/AB_post_treatments/AB_recover_mutational_history.cpp",
  "./src/AB_post_treatments/AB_recover_redundancy.cpp",
  "./src/AB_post_treatments/AB_recover_statistics.cpp",
  "./src/AB_post_treatments/AB_run_competition_experiment.cpp",
  "src/lib/Cell.cpp",
  "src/lib/Cell.h",
  "src/lib/Enums.h",
  "src/lib/Environment.cpp",
  "src/lib/Environment.h",
  "src/lib/Genome.cpp",
  "src/lib/Genome.h",
  "src/lib/GraphicDisplay.cpp",
  "src/lib/GraphicDisplay.h",
  "src/lib/InheritedProteins.cpp",
  "src/lib/InheritedProteins.h",
  "src/lib/IntegratedTests.cpp",
  "src/lib/IntegratedTests.h",
  "src/lib/Macros.h",
  "src/lib/MutationEvent.cpp",
  "src/lib/MutationEvent.h",
  "src/lib/MutationVector.cpp",
  "src/lib/MutationVector.h",
  "src/lib/Node.cpp",
  "src/lib/Node.h",
  "src/lib/ODE.cpp",
  "src/lib/ODE.h",
  "src/lib/Parameters.cpp",
  "src/lib/Parameters.h",
  "src/lib/Population.cpp",
  "src/lib/Population.h",
  "src/lib/Prng.cpp",
  "src/lib/Prng.h",
  "src/lib/ReplicationReport.cpp",
  "src/lib/ReplicationReport.h",
  "src/lib/Simulation.cpp",
  "src/lib/Simulation.h",
  "src/lib/SpeciesList.cpp",
  "src/lib/SpeciesList.h",
  "src/lib/Statistics.cpp",
  "src/lib/Statistics.h",
  "src/lib/Structs.h",
  "src/lib/Tree.cpp",
  "src/lib/Tree.h",
  "src/lib/TrophicGroup.cpp",
  "src/lib/TrophicGroup.h",
  "src/lib/TrophicNetwork.cpp",
  "src/lib/TrophicNetwork.h",
  "src/lib/UnitaryTests.cpp",
  "src/lib/UnitaryTests.h"
	]

	for file_to_explore in files_to_explore:
		detect_equality_operator_errors(file_to_explore)

