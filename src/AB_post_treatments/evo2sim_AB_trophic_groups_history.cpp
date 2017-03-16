
/**
 * \file      evo2sim_AB_trophic_groups_history.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      06-01-2016
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Recover the history of AB trophic groups
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
#include <unordered_map>

#include "../lib/Macros.h"
#include "../lib/Enums.h"
#include "../lib/Structs.h"
#include "../lib/Parameters.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME = "build/bin/evo2sim_AB_trophic_groups_history";

void   readArgs( int argc, char const** argv, size_t& backup_start, size_t& backup_end, size_t& backup_step );
void   printUsage( void );
void   printHeader( void );
size_t get_uptake_profile_max_length( size_t backup_start, size_t backup_end, size_t backup_step );
void   get_profiles_map( size_t L, size_t backup_start, size_t backup_end, size_t backup_step, std::unordered_map<std::string, size_t>* profiles, std::unordered_map<size_t, int>* groups );


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
  
  /*****************************************/
  /* 1) Read command line arguments        */
  /*****************************************/
  size_t backup_start = 0;
  size_t backup_end   = 0;
  size_t backup_step  = 0;
  readArgs(argc, argv, backup_start, backup_end, backup_step);
  
  /*****************************************/
  /* 2) Recover the trophic groups history */
  /*****************************************/
  
  /*------------------------*/
  /* 2.1) Recover meta data */
  /*------------------------*/
  size_t L = get_uptake_profile_max_length(backup_start, backup_end, backup_step);
  std::unordered_map<std::string, size_t> profiles;
  std::unordered_map<size_t, int> groups;
  std::unordered_map<std::string, size_t>::iterator it;
  get_profiles_map(L, backup_start, backup_end, backup_step, &profiles, &groups);
  
  /*------------------------*/
  /* 2.2) Write headers     */
  /*------------------------*/
  std::cout << "Writing the data:\n";
  std::ofstream history("trophic_groups_history.txt", std::ios::out | std::ios::trunc);
  for (it = profiles.begin(); it != profiles.end(); ++it)
  {
    if (it == profiles.begin())
    {
      history << it->second;
    }
    else
    {
      history << " " << it->second;
    }
  }
  history << "\n";
  
  std::ofstream levels("trophic_groups_levels.txt", std::ios::out | std::ios::trunc);
  for (it = profiles.begin(); it != profiles.end(); ++it)
  {
    if (it == profiles.begin())
    {
      levels << groups[it->second];
    }
    else
    {
      levels << " " << groups[it->second];
    }
  }
  levels << "\n";
  levels.close();
  
  /*------------------------*/
  /* 2.3) Write data        */
  /*------------------------*/
  size_t current_step = backup_start;
  while (current_step <= backup_end)
  {
    std::cout << "> working with backup " << current_step << " ...\n";
    Parameters* parameters = new Parameters(current_step);
    Simulation* simulation = new Simulation(parameters, current_step, false);
    simulation->initialize(FROM_BACKUP);
    
    for (it = profiles.begin(); it != profiles.end(); ++it)
    {
      size_t counter = 0;
      TrophicGroup* group = simulation->get_trophic_network()->get_first_group();
      while (group != NULL)
      {
        std::string profile = group->_uptake_profile;
        if (profile.size() < L)
        {
          std::string seq;
          seq.assign(L-profile.size(), '0');
          profile += seq;
        }
        if (it->first.compare(profile) == 0)
        {
          counter += group->get_number_of_cells();
        }
        group = simulation->get_trophic_network()->get_next_group();
      }
      if (it == profiles.begin())
      {
        history << counter;
      }
      else
      {
        history << " " << counter;
      }
    }
    history << "\n";
    
    delete simulation;
    simulation = NULL;
    delete parameters;
    parameters = NULL;
    current_step += backup_step;
  }
  
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
  std::cout << "Recover the AB trophic groups history (at backups resolution):\n";
  std::cout << "--------------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_trophic_groups_history -h or --help\n";
  std::cout << "   or: evo2sim_AB_trophic_groups_history [-start <first-backup>] [-end <last-backup>] [-step <backup-step>]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -start, --backup-start\n";
  std::cout << "        specify the starting date of the backup to load\n";
  std::cout << "  -end, --backup-end\n";
  std::cout << "        specify the ending date of the backup to load\n";
  std::cout << "  -step, --backup-step\n";
  std::cout << "        specify the step size of backup saving\n";
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
  std::cout << "Recover the AB trophic groups history (at backups resolution).\n";
  std::cout << "\n";
}

/**
 * \brief    Get the maximum length of uptake profiles in the backups
 * \details  --
 * \param    size_t backup_start
 * \param    size_t backup_end
 * \param    size_t backup_step
 * \return   \e size_t
 */
size_t get_uptake_profile_max_length( size_t backup_start, size_t backup_end, size_t backup_step )
{
  std::cout << "Recovering the uptake profiles maximum length:\n";
  size_t max_length   = 0;
  size_t current_step = backup_start;
  while (current_step <= backup_end)
  {
    std::cout << "> working with backup " << current_step << " ...\n";
    Parameters* parameters = new Parameters(current_step);
    Simulation* simulation = new Simulation(parameters, current_step, false);
    simulation->initialize(FROM_BACKUP);
    size_t L = simulation->get_trophic_network()->get_first_group()->_uptake_profile.size();
    if (max_length < L)
    {
      max_length = L;
    }
    delete simulation;
    simulation = NULL;
    delete parameters;
    parameters = NULL;
    current_step += backup_step;
  }
  return max_length;
}

/**
 * \brief    Get the profiles map
 * \details  --
 * \param    size_t L
 * \param    size_t backup_start
 * \param    size_t backup_end
 * \param    size_t backup_step
 * \param    std::unordered_map<std::string, size_t>* profiles
 * \param    std::unordered_map<size_t, int>* groups
 * \return   \e void
 */
void get_profiles_map( size_t L, size_t backup_start, size_t backup_end, size_t backup_step, std::unordered_map<std::string, size_t>* profiles, std::unordered_map<size_t, int>* groups )
{
  std::cout << "Recovering the profiles map:\n";
  profiles->clear();
  groups->clear();
  size_t counter      = 0;
  size_t current_step = backup_start;
  while (current_step <= backup_end)
  {
    std::cout << "> working with backup " << current_step << " ...\n";
    Parameters* parameters = new Parameters(current_step);
    Simulation* simulation = new Simulation(parameters, current_step, false);
    simulation->initialize(FROM_BACKUP);
    TrophicGroup* group = simulation->get_trophic_network()->get_first_group();
    while (group != NULL)
    {
      if (group->get_identifier() != 0)
      {
        std::string profile = group->_uptake_profile;
        if (profile.size() < L)
        {
          std::string seq;
          seq.assign(L-profile.size(), '0');
          profile += seq;
        }
        if (profiles->find(profile) == profiles->end())
        {
          profiles->insert({profile, counter});
          groups->insert({counter, group->get_trophic_level()});
          counter += 1;
        }
      }
      group = simulation->get_trophic_network()->get_next_group();
    }
    delete simulation;
    simulation = NULL;
    delete parameters;
    parameters = NULL;
    current_step += backup_step;
  }
}
