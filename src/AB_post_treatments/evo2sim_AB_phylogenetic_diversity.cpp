
/**
 * \file      evo2sim_AB_phylogenetic_diversity.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      06-01-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Compute the AB strains phylogenetic diversity
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
#include <unordered_map>

#include "../lib/Macros.h"
#include "../lib/Enums.h"
#include "../lib/Structs.h"
#include "../lib/Parameters.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME = "build/bin/evo2sim_AB_phylogenetic_diversity";

void readArgs( int argc, char const** argv, size_t& backup_start, size_t& backup_end, size_t& backup_step );
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
  
  /**************************************************/
  /* 1) Read command line arguments                 */
  /**************************************************/
  size_t backup_start = 0;
  size_t backup_end   = 0;
  size_t backup_step  = 0;
  readArgs(argc, argv, backup_start, backup_end, backup_step);
  
  /**************************************************/
  /* 2) Recover the simulation state at each backup */
  /**************************************************/
  std::ofstream file("./recovered_phylogenetic_diversity.txt", std::ios::out | std::ios::trunc);
  file << "backup Sp1 Sp2 depth\n";
  size_t current_step = backup_start;
  while (current_step <= backup_end)
  {
    std::cout << "> Working with backup " << current_step << " ...\n";
    
    /*----------------------------------------*/
    /* 2.1) Open the backups                  */
    /*----------------------------------------*/
    Parameters* parameters = new Parameters(current_step);
    Simulation* simulation = new Simulation(parameters, current_step, false);
    simulation->initialize(FROM_BACKUP);
    
    /*----------------------------------------*/
    /* 2.2) Compute SL diversity and write it */
    /*----------------------------------------*/
    double Sp_1  = 0.0;
    double Sp_2  = 0.0;
    double depth = (double)current_step-simulation->get_phylogenetic_tree()->get_common_ancestor_age();
    simulation->get_phylogenetic_tree()->compute_common_ancestor_SL_repartition(Sp_1, Sp_2);
    if (Sp_1 != -1.0 && Sp_2 != -1.0)
    {
      file << current_step << " " << Sp_1 << " " << Sp_2 << " " << depth << "\n";
    }
    
    /*----------------------------------------*/
    /* 2.3) Free memory                       */
    /*----------------------------------------*/
    delete simulation;
    simulation = NULL;
    delete parameters;
    parameters = NULL;
    current_step += backup_step;
  }
  file.close();
  return EXIT_SUCCESS;
}


/**
 * \brief    Read arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& backup_start
 * \param    size_t& backup_end
 * \param    size_t& backup_step
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& backup_start, size_t& backup_end, size_t& backup_step )
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
    if (strcmp(argv[i], "-start") == 0 || strcmp(argv[i], "--backup-start") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup_start = (size_t)atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-end") == 0 || strcmp(argv[i], "--backup-end") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup_end = (size_t)atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-step") == 0 || strcmp(argv[i], "--backup-step") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup_step = (size_t)atoi(argv[i+1]);
      }
    }
  }
  if (backup_end < backup_start || backup_step <= 0)
  {
    std::cout << "Error: incorrect command line structure (missing arguments).\n";
    exit(EXIT_FAILURE);
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
  std::cout << "Compute the AB strains phylogenetic diversity:\n";
  std::cout << "----------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_phylogenetic_diversity -h or --help\n";
  std::cout << "   or: evo2sim_AB_phylogenetic_diversity [-start <first-backup>] [-end <last-backup>] [-step <backup-step>]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -start, --backup-start\n";
  std::cout << "        specify the starting date of the backup to load (mandatory)\n";
  std::cout << "  -end, --backup-end\n";
  std::cout << "        specify the ending date of the backup to load (mandatory)\n";
  std::cout << "  -step, --backup-step\n";
  std::cout << "        specify the step size of backup saving (mandatory)\n";
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
  std::cout << "Compute the AB strains phylogenetic diversity.\n";
  std::cout << "\n";
}
