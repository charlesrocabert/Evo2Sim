
/**
 * \file      evo2sim_AB_last_individuals.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      13-11-2016
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Save last best individual data of A and B populations
 */

/****************************************************************************
 * Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * Web: https://github.com/charlesrocabert/Evo2Sim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "../../cmake/Config.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <assert.h>

#include "../lib/Macros.h"
#include "../lib/Enums.h"
#include "../lib/Structs.h"
#include "../lib/Parameters.h"
#include "../lib/Population.h"
#include "../lib/Environment.h"
#include "../lib/Tree.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME          = "build/bin/evo2sim_AB_last_individuals";
const std::string DEFAULT_FILENAME         = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH  = "population_to_load";
const std::string DEFAULT_ENVIRONMENT_PATH = "environment_to_load";

void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path, std::string& optional_environment_path );
void printUsage( void );
void printHeader( void );


/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  printHeader();
  
  /***********************************/
  /* 1) Read command line arguments  */
  /***********************************/
  size_t      rep                       = 0;
  std::string optional_filename         = "";
  std::string optional_population_path  = "";
  std::string optional_environment_path = "";
  readArgs(argc, argv, rep, optional_filename, optional_population_path, optional_environment_path);
  
  /***********************************/
  /* 2) Load parameters from file    */
  /***********************************/
  Parameters* parameters = new Parameters();
  if (strcmp(optional_filename.c_str(), "") != 0)
  {
    parameters->load_parameters_from_file(optional_filename);
  }
  else
  {
    parameters->load_parameters_from_file(DEFAULT_FILENAME);
  }
  
  /***********************************/
  /* 3) Load the evolved population  */
  /***********************************/
  Population* evolved_population = NULL;
  if (strcmp(optional_population_path.c_str(), "") != 0)
  {
    gzFile pop_file    = gzopen(optional_population_path.c_str(), "r");
    evolved_population = new Population(parameters, pop_file);
    gzclose(pop_file);
  }
  else
  {
    gzFile pop_file    = gzopen(DEFAULT_POPULATION_PATH.c_str(), "r");
    evolved_population = new Population(parameters, pop_file);
    gzclose(pop_file);
  }
  
  /***********************************/
  /* 4) Load the evolved environment */
  /***********************************/
  Environment* evolved_environment = NULL;
  if (strcmp(optional_environment_path.c_str(), "") != 0)
  {
    gzFile env_file     = gzopen(optional_environment_path.c_str(), "r");
    evolved_environment = new Environment(parameters, env_file);
    gzclose(env_file);
  }
  else
  {
    gzFile env_file     = gzopen(DEFAULT_ENVIRONMENT_PATH.c_str(), "r");
    evolved_environment = new Environment(parameters, env_file);
    gzclose(env_file);
  }
  
  /***********************************/
  /* 5) Run the post-treatment       */
  /***********************************/
  double best_A_score = 0.0;
  double best_B_score = 0.0;
  size_t best_A_index = 0;
  size_t best_B_index = 0;
  for (size_t i = 0; i < evolved_population->get_width()*evolved_population->get_height(); i++)
  {
    Cell* cell = evolved_population->get_cell(i);
    if (cell->isActive() && cell->isAlive())
    {
      if ((cell->get_trophic_level() == 0 || cell->get_trophic_level() == 1) && best_A_score < cell->get_score())
      {
        best_A_score = cell->get_score();
        best_A_index = i;
      }
      else if ((cell->get_trophic_level() == 2 || cell->get_trophic_level() == 3) && best_B_score < cell->get_score())
      {
        best_B_score = cell->get_score();
        best_B_index = i;
      }
    }
  }
  Cell* A_cell = evolved_population->get_cell(best_A_index);
  A_cell->load_genome_in_ODE_system(evolved_environment, true, false);
  Cell* B_cell = evolved_population->get_cell(best_B_index);
  B_cell->load_genome_in_ODE_system(evolved_environment, true, false);
  
  /* Save best A individual */
  std::ofstream genome_file("best_A_genome.txt", std::ios::out | std::ios::trunc);
  std::ofstream nodes_file("best_A_nodes.txt", std::ios::out | std::ios::trunc);
  std::ofstream edges_file("best_A_edges.txt", std::ios::out | std::ios::trunc);
  std::ofstream state_file("best_A_state.txt", std::ios::out | std::ios::trunc);
  A_cell->write_genome(genome_file);
  A_cell->write_metabolic_network(nodes_file, edges_file);
  A_cell->write_metabolic_amounts(state_file);
  genome_file.close();
  nodes_file.close();
  edges_file.close();
  state_file.close();
  
  /* Save best B individual */
  genome_file.open("best_B_genome.txt", std::ios::out | std::ios::trunc);
  nodes_file.open("best_B_nodes.txt", std::ios::out | std::ios::trunc);
  edges_file.open("best_B_edges.txt", std::ios::out | std::ios::trunc);
  state_file.open("best_B_state.txt", std::ios::out | std::ios::trunc);
  B_cell->write_genome(genome_file);
  B_cell->write_metabolic_network(nodes_file, edges_file);
  B_cell->write_metabolic_amounts(state_file);
  genome_file.close();
  nodes_file.close();
  edges_file.close();
  state_file.close();
  
  /***********************************/
  /* 6) Free the memory              */
  /***********************************/
  delete evolved_population;
  evolved_population = NULL;
  delete parameters;
  parameters = NULL;
  
  return EXIT_SUCCESS;
}


/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& rep
 * \param    std::string& optional_filename
 * \param    string& optional_population_path
 * \param    string& optional_environment_path
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path, std::string& optional_environment_path )
{
  for (int i = 0; i < argc; i++)
  {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
    {
      printUsage();
      exit(EXIT_SUCCESS);
    }
    if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0)
    {
      std::cout << PACKAGE << " (" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << ")\n";
      exit(EXIT_SUCCESS);
    }
    if (strcmp(argv[i], "-rep") == 0 || strcmp(argv[i], "--rep") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        rep = (size_t)atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--file") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_filename = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-pop") == 0 || strcmp(argv[i], "--population-path") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_population_path = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-env") == 0 || strcmp(argv[i], "--environment-path") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_environment_path = argv[i+1];
      }
    }
  }
}

/**
 * \brief    Print usage
 * \details  --
 * \param    void
 * \return   \e void
 */
void printUsage( void )
{
  std::cout << "\n";
  std::cout << "***************************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << "                                                                           \n";
  std::cout << " Multi-scale and individual-based computational model dedicated            \n";
  std::cout << " to in silico experimental evolution.                                      \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
  std::cout << "Save last best individual data of A and B populations:\n";
  std::cout << "------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_last_individuals -h or --help\n";
  std::cout << "   or: evo2sim_AB_last_individuals [-rep <nb-repetitions>] [-f <param-file>] [-pop <population-file>]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -rep, --rep\n";
  std::cout << "        specify the number of repetitions (mandatory)\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify parameters file (default: parameters.txt)\n";
  std::cout << "  -pop, --population-path\n";
  std::cout << "        specify the path of the population backup to load (default: population_to_load)\n";
  std::cout << "  -env, --environment-path\n";
  std::cout << "        specify the path of the environment backup to load (default: environment_to_load)\n";
  std::cout << "\n";
}

/**
 * \brief    Print header
 * \details  --
 * \param    void
 * \return   \e void
 */
void printHeader( void )
{
  std::cout << "\n";
  std::cout << "***************************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << "                                                                           \n";
  std::cout << " Multi-scale and individual-based computational model dedicated            \n";
  std::cout << " to in silico experimental evolution.                                      \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Save last best individual data of A and B populations.\n";
  std::cout << "\n";
}

