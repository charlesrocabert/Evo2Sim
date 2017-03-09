
/**
 * \file      evo2sim_AB_black_queen.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      27-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Measure the black queen effect on A/B ecotypes
 */

/****************************************************************************
 * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * E-mail: charles.rocabert@gmail.com
 * Web: http://www.evoevo.eu/
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
#include "../lib/Tree.h"
#include "../lib/Node.h"

const std::string EXECUTABLE_NAME  = "build/bin/evo2sim_AB_black_queen";
const std::string DEFAULT_FILENAME = "parameters.txt";

void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& backup );
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
  
  /**********************************/
  /* 1) Read command line arguments */
  /**********************************/
  
  std::string optional_filename = "";
  size_t      backup            = 0;
  readArgs(argc, argv, optional_filename, backup);
  
  /**********************************/
  /* 2) Load parameters from file   */
  /**********************************/
  
  Parameters* parameters = new Parameters();
  if (strcmp(optional_filename.c_str(), "") != 0)
  {
    parameters->load_parameters_from_file(optional_filename);
  }
  else
  {
    parameters->load_parameters_from_file(DEFAULT_FILENAME);
  }
  
  /**********************************/
  /* 3) Load backups                */
  /**********************************/
  
  /*---------------------------------*/
  /* 3.1) Load the population        */
  /*---------------------------------*/
  std::stringstream population_filename;
  population_filename << "population/population_" << backup;
  gzFile pop_file = gzopen(population_filename.str().c_str(), "r");
  Population* evolved_population = new Population(parameters, pop_file);
  gzclose(pop_file);
  
  /*---------------------------------*/
  /* 3.2) Load the lineage tree      */
  /*---------------------------------*/
  std::stringstream lineage_tree_filename;
  lineage_tree_filename << "tree/lineage_tree_" << backup;
  gzFile lineage_tree_file = gzopen(lineage_tree_filename.str().c_str(), "r");
  Tree* evolved_lineage = new Tree(parameters, evolved_population, lineage_tree_file);
  gzclose(lineage_tree_file);
  
  /*---------------------------------*/
  /* 3.3) Load the phylogenetic tree */
  /*---------------------------------*/
  std::stringstream phylogenetic_tree_filename;
  phylogenetic_tree_filename << "tree/phylogenetic_tree_" << backup;
  gzFile phylogenetic_tree_file = gzopen(phylogenetic_tree_filename.str().c_str(), "r");
  Tree* evolved_phylogeny = new Tree(parameters, evolved_population, phylogenetic_tree_file);
  gzclose(phylogenetic_tree_file);
  
  
  /**********************************/
  /* 4) Run the post-treatment      */
  /**********************************/
  
  /*-------------------------------*/
  /* 4.1) Get the master root      */
  /*-------------------------------*/
  Node* master_root = evolved_lineage->get_node(0);
  
  /*-------------------------------*/
  /* 4.2) Get the common ancestor  */
  /*-------------------------------*/
  std::ofstream output("rootedge.txt", std::ios::out | std::ios::trunc);
  Node* CA = master_root;
  while (CA->get_number_of_children() == 1)
  {
    if (CA->get_id() != 0)
    {
      output << CA->get_replication_report()->get_birth_time() << " ";
      output << CA->get_replication_report()->get_generation() << " ";
      output << CA->get_replication_report()->get_trophic_level() << "\n";
    }
    CA = CA->get_child(0);
  }
  if (CA->get_id() != 0)
  {
    output << CA->get_replication_report()->get_birth_time() << " ";
    output << CA->get_replication_report()->get_generation() << " ";
    output << CA->get_replication_report()->get_trophic_level() << "\n";
  }
  output.close();
  
  /*-------------------------------*/
  /* 4.3) Get nascent A/B subtrees */
  /*-------------------------------*/
  Node* subA1 = CA->get_child(0);
  Node* subA2 = CA->get_child(1);
  
  /*-------------------------------*/
  /* 4.4) Climb subbranches        */
  /*-------------------------------*/
  output.open("subbranch1.txt", std::ios::out | std::ios::trunc);
  while (subA1->get_number_of_children() == 1)
  {
    output << subA1->get_replication_report()->get_birth_time() << " ";
    output << subA1->get_replication_report()->get_generation() << " ";
    output << subA1->get_replication_report()->get_trophic_level() << "\n";
    subA1 = subA1->get_child(0);
  }
  output.close();
  
  output.open("subbranch2.txt", std::ios::out | std::ios::trunc);
  while (subA2->get_number_of_children() == 1)
  {
    output << subA2->get_replication_report()->get_birth_time() << " ";
    output << subA2->get_replication_report()->get_generation() << " ";
    output << subA2->get_replication_report()->get_trophic_level() << "\n";
    subA2 = subA2->get_child(0);
  }
  output.close();
  
  /*****************************************/
  /* 5) Free the memory                    */
  /*****************************************/
  delete evolved_phylogeny;
  evolved_phylogeny = NULL;
  delete evolved_lineage;
  evolved_lineage = NULL;
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
 * \param    std::string& optional_filename
 * \param    size_t& backup
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& backup )
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
    if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--backup-date") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup = (size_t)atoi(argv[i+1]);
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
  std::cout << " This software is dedicated to EvoEvo WP2 models development               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " E-mail: charles.rocabert@gmail.com                                        \n";
  std::cout << " Web: http://www.evoevo.eu/                                                \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
  std::cout << "Recover the evolution of ecotypes A/B along the lineage tree:\n";
  std::cout << "-------------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_black_queen -h or --help\n";
  std::cout << "   or: evo2sim_AB_black_queen [-f <param-file>] [-b <backup-date>]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify the parameters file (default: parameters.txt)\n";
  std::cout << "  -b, --backup\n";
  std::cout << "        specify the backup date (mandatory)\n";
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
  std::cout << " This software is dedicated to EvoEvo WP2 models development               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " E-mail: charles.rocabert@gmail.com                                        \n";
  std::cout << " Web: http://www.evoevo.eu/                                                \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Recover the evolution of ecotypes A/B along the lineage tree.\n";
  std::cout << "\n";
}
