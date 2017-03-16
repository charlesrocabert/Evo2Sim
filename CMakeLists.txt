
#****************************************************************************
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
#****************************************************************************

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the minimum required CMake version                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
cmake_minimum_required(VERSION 2.6)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the project name                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
project(EVO2SIM)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the package name                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(PACKAGE "\"Evo2Sim\"")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the version number                                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 2)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set other options                                                                #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(AB_POST_TREATMENTS 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Add definitions to CXX flags                                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ADD_DEFINITIONS(-std=c++11)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set DEBUG and RELEASE flags                                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -g -pg -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -O2 -Wall -Wextra -pedantic")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the modules path                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the executables output path                                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(EXECUTABLE_OUTPUT_PATH ../build/bin)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Define the executables                                                           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
set(CREATE_EXECUTABLE evo2sim_create)
add_executable(${CREATE_EXECUTABLE} src/evo2sim_create.cpp)
set(BOOTSTRAP_EXECUTABLE evo2sim_bootstrap)
add_executable(${BOOTSTRAP_EXECUTABLE} src/evo2sim_bootstrap.cpp)
set(RUN_EXECUTABLE evo2sim_run)
add_executable(${RUN_EXECUTABLE} src/evo2sim_run.cpp)
set(GENERATE_FIGURES_EXECUTABLE evo2sim_generate_figures)
add_executable(${GENERATE_FIGURES_EXECUTABLE} src/evo2sim_generate_figures.cpp)
set(PARAMETERS_EXECUTABLE evo2sim_recover_parameters)
add_executable(${PARAMETERS_EXECUTABLE} src/evo2sim_recover_parameters.cpp)
set(UNITARY_TESTS_EXECUTABLE evo2sim_unitary_tests)
add_executable(${UNITARY_TESTS_EXECUTABLE} src/evo2sim_unitary_tests.cpp)
set(INTEGRATED_TESTS_EXECUTABLE evo2sim_integrated_tests)
add_executable(${INTEGRATED_TESTS_EXECUTABLE} src/evo2sim_integrated_tests.cpp)

### AB post treatments ###

if(AB_POST_TREATMENTS)
  set(AB_STATISTICS_EXECUTABLE evo2sim_AB_recover_statistics)
  add_executable(${AB_STATISTICS_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_recover_statistics.cpp)
  set(AB_DIVERSITY_EXECUTABLE evo2sim_AB_phylogenetic_diversity)
  add_executable(${AB_DIVERSITY_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_phylogenetic_diversity.cpp)
  set(AB_GROUPS_HISTORY_EXECUTABLE evo2sim_AB_trophic_groups_history)
  add_executable(${AB_GROUPS_HISTORY_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_trophic_groups_history.cpp)
  set(AB_BLACK_QUEEN_EXECUTABLE evo2sim_AB_black_queen)
  add_executable(${AB_BLACK_QUEEN_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_black_queen.cpp)
  set(AB_FREQ_DEP_EXECUTABLE evo2sim_AB_frequency_dependence)
  add_executable(${AB_FREQ_DEP_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_frequency_dependence.cpp)
  set(AB_COMPETITION_EXECUTABLE evo2sim_AB_run_competition_experiment)
  add_executable(${AB_COMPETITION_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_run_competition_experiment.cpp)
  set(AB_LAST_INDIVIDUALS_EXECUTABLE evo2sim_AB_last_individuals)
  add_executable(${AB_LAST_INDIVIDUALS_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_last_individuals.cpp)
  set(AB_TROPHIC_STATISTICS_EXECUTABLE evo2sim_AB_recover_trophic_statistics)
  add_executable(${AB_TROPHIC_STATISTICS_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_recover_trophic_statistics.cpp)
  set(AB_HEATMAPS_EXECUTABLE evo2sim_AB_recover_metabolic_heatmaps)
  add_executable(${AB_HEATMAPS_EXECUTABLE} src/AB_post_treatments/evo2sim_AB_recover_metabolic_heatmaps.cpp)
endif(AB_POST_TREATMENTS)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Find and link external libraries                                                 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${GSL_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${GSL_LIBRARIES})

  if(AB_POST_TREATMENTS)
    target_link_libraries(${AB_STATISTICS_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_DIVERSITY_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_GROUPS_HISTORY_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_BLACK_QUEEN_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_FREQ_DEP_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_COMPETITION_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_LAST_INDIVIDUALS_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_TROPHIC_STATISTICS_EXECUTABLE} ${GSL_LIBRARIES})
    target_link_libraries(${AB_HEATMAPS_EXECUTABLE} ${GSL_LIBRARIES})
  endif(AB_POST_TREATMENTS)
endif(GSL_FOUND)

find_package(TBB REQUIRED)
if(TBB_FOUND)
  include_directories(${TBB_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${TBB_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${TBB_LIBRARIES})

  if(AB_POST_TREATMENTS)
    target_link_libraries(${AB_STATISTICS_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_DIVERSITY_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_GROUPS_HISTORY_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_BLACK_QUEEN_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_FREQ_DEP_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_COMPETITION_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_LAST_INDIVIDUALS_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_TROPHIC_STATISTICS_EXECUTABLE} ${TBB_LIBRARIES})
    target_link_libraries(${AB_HEATMAPS_EXECUTABLE} ${TBB_LIBRARIES})
  endif(AB_POST_TREATMENTS)
endif(TBB_FOUND)

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(${CREATE_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${RUN_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${PARAMETERS_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${ZLIB_LIBRARIES})
  target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${ZLIB_LIBRARIES})

  if(AB_POST_TREATMENTS)
    target_link_libraries(${AB_STATISTICS_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_DIVERSITY_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_GROUPS_HISTORY_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_BLACK_QUEEN_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_FREQ_DEP_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_COMPETITION_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_LAST_INDIVIDUALS_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_TROPHIC_STATISTICS_EXECUTABLE} ${ZLIB_LIBRARIES})
    target_link_libraries(${AB_HEATMAPS_EXECUTABLE} ${ZLIB_LIBRARIES})
  endif(AB_POST_TREATMENTS)
endif(ZLIB_FOUND)

if(WITH_GRAPHICS_CONTEXT)
  find_package(SFML 2 REQUIRED graphics window system)
  if(SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries(${CREATE_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${BOOTSTRAP_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${RUN_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${PARAMETERS_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${UNITARY_TESTS_EXECUTABLE} ${SFML_LIBRARIES})
    target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} ${SFML_LIBRARIES})

    if(AB_POST_TREATMENTS)
      target_link_libraries(${AB_STATISTICS_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_DIVERSITY_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_GROUPS_HISTORY_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_BLACK_QUEEN_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_FREQ_DEP_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_COMPETITION_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_LAST_INDIVIDUALS_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_TROPHIC_STATISTICS_EXECUTABLE} ${SFML_LIBRARIES})
      target_link_libraries(${AB_HEATMAPS_EXECUTABLE} ${SFML_LIBRARIES})
    endif(AB_POST_TREATMENTS)
  endif(SFML_FOUND)
endif(WITH_GRAPHICS_CONTEXT)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Create and link Evo2Sim library                                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if(WITH_GRAPHICS_CONTEXT)
  add_library(
    Evo2Sim
    src/lib/Cell.cpp
    src/lib/Cell.h
    src/lib/Enums.h
    src/lib/Environment.cpp
    src/lib/Environment.h
    src/lib/Genome.cpp
    src/lib/Genome.h
    src/lib/GraphicDisplay.cpp
    src/lib/GraphicDisplay.h
    src/lib/InheritedProteins.cpp
    src/lib/InheritedProteins.h
    src/lib/IntegratedTests.cpp
    src/lib/IntegratedTests.h
    src/lib/Macros.h
    src/lib/MutationEvent.cpp
    src/lib/MutationEvent.h
    src/lib/MutationVector.cpp
    src/lib/MutationVector.h
    src/lib/Node.cpp
    src/lib/Node.h
    src/lib/ODE.cpp
    src/lib/ODE.h
    src/lib/Parameters.cpp
    src/lib/Parameters.h
    src/lib/Population.cpp
    src/lib/Population.h
    src/lib/Prng.cpp
    src/lib/Prng.h
    src/lib/ReplicationReport.cpp
    src/lib/ReplicationReport.h
    src/lib/Simulation.cpp
    src/lib/Simulation.h
    src/lib/SpeciesList.cpp
    src/lib/SpeciesList.h
    src/lib/Statistics.cpp
    src/lib/Statistics.h
    src/lib/Structs.h
    src/lib/Tree.cpp
    src/lib/Tree.h
    src/lib/TrophicGroup.cpp
    src/lib/TrophicGroup.h
    src/lib/TrophicNetwork.cpp
    src/lib/TrophicNetwork.h
    src/lib/UnitaryTests.cpp
    src/lib/UnitaryTests.h
  )
elseif(NOT WITH_GRAPHICS_CONTEXT)
  add_library(
    Evo2Sim
    src/lib/Cell.cpp
    src/lib/Cell.h
    src/lib/Enums.h
    src/lib/Environment.cpp
    src/lib/Environment.h
    src/lib/Genome.cpp
    src/lib/Genome.h
    src/lib/InheritedProteins.cpp
    src/lib/InheritedProteins.h
    src/lib/IntegratedTests.cpp
    src/lib/IntegratedTests.h
    src/lib/Macros.h
    src/lib/MutationEvent.cpp
    src/lib/MutationEvent.h
    src/lib/MutationVector.cpp
    src/lib/MutationVector.h
    src/lib/Node.cpp
    src/lib/Node.h
    src/lib/ODE.cpp
    src/lib/ODE.h
    src/lib/Parameters.cpp
    src/lib/Parameters.h
    src/lib/Population.cpp
    src/lib/Population.h
    src/lib/Prng.cpp
    src/lib/Prng.h
    src/lib/ReplicationReport.cpp
    src/lib/ReplicationReport.h
    src/lib/Simulation.cpp
    src/lib/Simulation.h
    src/lib/SpeciesList.cpp
    src/lib/SpeciesList.h
    src/lib/Statistics.cpp
    src/lib/Statistics.h
    src/lib/Structs.h
    src/lib/Tree.cpp
    src/lib/Tree.h
    src/lib/TrophicGroup.cpp
    src/lib/TrophicGroup.h
    src/lib/TrophicNetwork.cpp
    src/lib/TrophicNetwork.h
    src/lib/UnitaryTests.cpp
    src/lib/UnitaryTests.h
  )
endif(WITH_GRAPHICS_CONTEXT)

target_link_libraries(Evo2Sim z tbb gsl gslcblas)

target_link_libraries(${CREATE_EXECUTABLE} Evo2Sim)
target_link_libraries(${BOOTSTRAP_EXECUTABLE} Evo2Sim)
target_link_libraries(${RUN_EXECUTABLE} Evo2Sim)
target_link_libraries(${GENERATE_FIGURES_EXECUTABLE} Evo2Sim)
target_link_libraries(${PARAMETERS_EXECUTABLE} Evo2Sim)
target_link_libraries(${UNITARY_TESTS_EXECUTABLE} Evo2Sim)
target_link_libraries(${INTEGRATED_TESTS_EXECUTABLE} Evo2Sim)

if(AB_POST_TREATMENTS)
  target_link_libraries(${AB_STATISTICS_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_DIVERSITY_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_GROUPS_HISTORY_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_BLACK_QUEEN_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_FREQ_DEP_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_COMPETITION_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_LAST_INDIVIDUALS_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_TROPHIC_STATISTICS_EXECUTABLE} Evo2Sim)
  target_link_libraries(${AB_HEATMAPS_EXECUTABLE} Evo2Sim)
endif(AB_POST_TREATMENTS)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Build the Config.h file                                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/modules/Config.h.in"
  "${PROJECT_BINARY_DIR}/Config.h"
)
include_directories("${PROJECT_BINARY_DIR}")

