
/**
 * \file      evo2sim_AB_recover_trophic_statistics.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      06-05-2016
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Recover some trophic network statistics from AB populations
 */

/****************************************************************************
 * Evo2Sim (Evolution of Evolution Simulator)
 * -------------------------------------------
 * Digital evolution model dedicated to
 * bacterial in silico experimental evolution.
 *
 * Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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
#include <unordered_map>

#include "../lib/Macros.h"
#include "../lib/Enums.h"
#include "../lib/Structs.h"
#include "../lib/Parameters.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME = "build/bin/evo2sim_AB_recover_trophic_statistics";

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
  
  /**********************************/
  /* 1) Read command line arguments */
  /**********************************/
  size_t backup_start = 0;
  size_t backup_end   = 0;
  size_t backup_step  = 0;
  readArgs(argc, argv, backup_start, backup_end, backup_step);
  
  /**********************************/
  /* 2) Recover statistics          */
  /**********************************/
  std::ofstream history("recovered_trophic_stats.txt", std::ios::out | std::ios::trunc);
  
  history << "time" << " ";
  history << "A_count" << " ";
  history << "A_genome_size" << " ";
  history << "A_non_coding" << " ";
  history << "A_redundancy" << " ";
  history << "A_pump" << " ";
  history << "A_prod" << " ";
  history << "A_pm_pump" << " ";
  history << "A_pm_prod" << " ";
  history << "B_count" << " ";
  history << "B_genome_size" << " ";
  history << "B_non_coding" << " ";
  history << "B_redundancy" << " ";
  history << "B_pump" << " ";
  history << "B_prod" << " ";
  history << "B_pm_pump" << " ";
  history << "B_pm_prod" << "\n";
  
  size_t current_step = backup_start;
  while (current_step <= backup_end)
  {
    std::cout << "> working with backup " << current_step << " ...\n";
    Parameters* parameters = new Parameters(current_step);
    Simulation* simulation = new Simulation(parameters, current_step, false);
    simulation->initialize(FROM_BACKUP);
    
    double A_count       = 0.0;
    double A_genome_size = 0.0;
    double A_non_coding  = 0.0;
    double A_redundancy  = 0.0;
    double A_pump        = 0.0;
    double A_prod        = 0.0;
    double A_pm_pump     = 0.0;
    double A_pm_prod     = 0.0;
    
    double B_count       = 0.0;
    double B_genome_size = 0.0;
    double B_non_coding  = 0.0;
    double B_redundancy  = 0.0;
    double B_pump        = 0.0;
    double B_prod        = 0.0;
    double B_pm_pump     = 0.0;
    double B_pm_prod     = 0.0;
    
    for (size_t x = 0; x < parameters->get_width(); x++)
    {
      for (size_t y = 0; y < parameters->get_height(); y++)
      {
        Cell* cell = simulation->get_population()->get_cell(x, y);
        if (cell->isActive() && cell->isAlive())
        {
          if (cell->get_trophic_level() == LEVEL_0 || cell->get_trophic_level() == LEVEL_1)
          {
            A_count       += 1.0;
            A_genome_size += cell->get_genome()->get_size();
            A_non_coding  += cell->get_genome()->get_size()-cell->get_genome()->get_coding_size();
            A_redundancy  += cell->get_replication_report()->get_mean_metabolic_redundancy();
            std::vector<int> pumped, produced, pm_pumped, pm_produced;
            cell->get_pumped_and_produced_metabolites(pumped, produced, pm_pumped, pm_produced);
            A_pump        += (double)pumped.size();
            A_prod        += (double)produced.size();
            A_pm_pump     += (double)pm_pumped.size();
            A_pm_prod     += (double)pm_produced.size();
          }
          else if (cell->get_trophic_level() == LEVEL_2 || cell->get_trophic_level() == NO_LEVEL)
          {
            B_count       += 1.0;
            B_genome_size += cell->get_genome()->get_size();
            B_non_coding  += cell->get_genome()->get_size()-cell->get_genome()->get_coding_size();
            B_redundancy  += cell->get_replication_report()->get_mean_metabolic_redundancy();
            std::vector<int> pumped, produced, pm_pumped, pm_produced;
            cell->get_pumped_and_produced_metabolites(pumped, produced, pm_pumped, pm_produced);
            B_pump        += (double)pumped.size();
            B_prod        += (double)produced.size();
            B_pm_pump     += (double)pm_pumped.size();
            B_pm_prod     += (double)pm_produced.size();
          }
        }
      }
    }
    if (A_count > 0.0)
    {
      A_genome_size /= A_count;
      A_non_coding  /= A_count;
      A_redundancy  /= A_count;
      A_pump        /= A_count;
      A_prod        /= A_count;
      A_pm_pump     /= A_count;
      A_pm_prod     /= A_count;
    }
    if (B_count > 0.0)
    {
      B_genome_size /= B_count;
      B_non_coding  /= B_count;
      B_redundancy  /= B_count;
      B_pump        /= B_count;
      B_prod        /= B_count;
      B_pm_pump     /= B_count;
      B_pm_prod     /= B_count;
    }
    
    history << current_step << " ";
    history << A_count << " ";
    history << A_genome_size << " ";
    history << A_non_coding << " ";
    history << A_redundancy << " ";
    history << A_pump << " ";
    history << A_prod << " ";
    history << A_pm_pump << " ";
    history << A_pm_prod << " ";
    history << B_count << " ";
    history << B_genome_size << " ";
    history << B_non_coding << " ";
    history << B_redundancy << " ";
    history << B_pump << " ";
    history << B_prod << " ";
    history << B_pm_pump << " ";
    history << B_pm_prod << "\n";
    
    delete simulation;
    simulation = NULL;
    delete parameters;
    parameters = NULL;
    current_step += backup_step;
  }
  
  history.close();
  
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
  std::cout << "                                                                           \n";
  std::cout << " Digital evolution model dedicated to                                      \n";
  std::cout << " bacterial in silico experimental evolution.                               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
  std::cout << "Recover some trophic network statistics from AB populations (at backups resolution):\n";
  std::cout << "------------------------------------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_recover_trophic_statistics -h or --help\n";
  std::cout << "   or: evo2sim_AB_recover_trophic_statistics [-start <first-backup>] [-end <last-backup>] [-step <backup-step>]\n";
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
  std::cout << "                                                                           \n";
  std::cout << " Digital evolution model dedicated to                                      \n";
  std::cout << " bacterial in silico experimental evolution.                               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Recover some trophic network statistics from AB populations (at backups resolution).\n";
  std::cout << "\n";
}
