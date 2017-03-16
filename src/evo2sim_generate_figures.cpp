
/**
 * \file      evo2sim_generate_figures.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      05-04-2015
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Generate figures from backup files
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

#include "../cmake/Config.h"

#include <iostream>
#include <cstring>
#include <assert.h>

#include "./lib/Macros.h"
#include "./lib/Parameters.h"
#include "./lib/Simulation.h"

const std::string EXECUTABLE_NAME = "build/bin/evo2sim_generate_figures";

void readArgs( int argc, char const** argv, size_t& backup_time );
void printUsage( void );
void printHeader( void );
void write_simulation_data( Simulation* simulation );
void generate_figures( char const** argv, Simulation* simulation );


/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  /**********************************/
  /* 1) Read command line arguments */
  /**********************************/
  size_t  backup_time = 0;
  readArgs(argc, argv, backup_time);
  
  /**********************************/
  /* 2) Load parameters from backup */
  /**********************************/
  printHeader();
  Parameters* parameters = new Parameters(backup_time);
  
  /**********************************/
  /* 3) Load simulation from backup */
  /**********************************/
  Simulation* simulation = new Simulation(parameters, backup_time, false);
  
  /**********************************/
  /* 4) Initialize simulation       */
  /**********************************/
  simulation->initialize( FROM_BACKUP );
  
  /**********************************/
  /* 5) Generate figures            */
  /**********************************/
  printf("> Generate figures ...\n");
  write_simulation_data(simulation);
  generate_figures(argv, simulation);
  
  /**********************************/
  /* 6) Close statistics files and  */
  /*    free memory                 */
  /**********************************/
  printf("> Figures generated.\n\nTo see the figures, open viewer/viewer.html in an internet browser.\n");
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
  delete parameters;
  parameters = NULL;
  return EXIT_SUCCESS;
}


/**
 * \brief    Read arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& backup_time
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& backup_time )
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
    if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--backup-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: -b command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup_time = atoi(argv[i+1]);
      }
    }
  }
  if (backup_time <= 0)
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
  std::cout << "Generate figures from backup files:\n";
  std::cout << "-----------------------------------\n";
  std::cout << "Usage: generate_figures -h or --help\n";
  std::cout << "   or: generate_figures [-b <backup-time>] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit (optional)\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit (optional)\n";
  std::cout << "  -b, --backup-time\n";
  std::cout << "        set the date of the backup to load (mandatory)\n";
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
  std::cout << "Generate figures from backup files.\n";
  std::cout << "\n";
}

/**
 * \brief    Write simulation data
 * \details  --
 * \param    Simulation* simulation
 * \return   \e void
 */
void write_simulation_data( Simulation* simulation )
{
  if (simulation->get_population()->get_population_size() == 0)
  {
    return;
  }
  if (!simulation->get_parameters()->get_parallel_computing())
  {
    simulation->get_statistics()->write_best_genome_structure_file();
    simulation->get_statistics()->write_best_genome_file();
    simulation->get_statistics()->write_best_inherited_proteins_file();
    simulation->get_statistics()->write_best_genetic_regulation_network_file();
    simulation->get_statistics()->write_best_metabolic_network_file();
    simulation->get_statistics()->write_best_metabolic_amounts_file();
    simulation->get_statistics()->write_last_environment_metabolic_amounts_file();
    simulation->get_statistics()->write_last_trophic_network_file();
    if (simulation->get_parameters()->get_simulation_backup_step() > 0 && simulation->get_lineage_tree()->get_number_of_nodes() > 0)
    {
      simulation->get_statistics()->write_last_lineage_tree_statistics();
    }
    if (simulation->get_parameters()->get_simulation_backup_step() > 0 && simulation->get_phylogenetic_tree()->get_number_of_nodes() > 0)
    {
      simulation->get_statistics()->write_last_phylogenetic_tree_statistics();
    }
  }
  else if (simulation->get_parameters()->get_parallel_computing())
  {
    tbb::task_group tasks;
    tasks.run([=]{simulation->get_statistics()->write_best_genome_structure_file();});
    tasks.run([=]{simulation->get_statistics()->write_best_genome_file();});
    tasks.run([=]{simulation->get_statistics()->write_best_inherited_proteins_file();});
    tasks.run([=]{simulation->get_statistics()->write_best_genetic_regulation_network_file();});
    tasks.run([=]{simulation->get_statistics()->write_best_metabolic_network_file();});
    tasks.run([=]{simulation->get_statistics()->write_best_metabolic_amounts_file();});
    tasks.run([=]{simulation->get_statistics()->write_last_environment_metabolic_amounts_file();});
    tasks.run([=]{simulation->get_statistics()->write_last_trophic_network_file();});
    if (simulation->get_parameters()->get_simulation_backup_step() > 0 && simulation->get_lineage_tree()->get_number_of_nodes() > 0)
    {
      tasks.run([=]{simulation->get_statistics()->write_last_lineage_tree_statistics();});
    }
    if (simulation->get_parameters()->get_simulation_backup_step() > 0 && simulation->get_phylogenetic_tree()->get_number_of_nodes() > 0)
    {
      tasks.run([=]{simulation->get_statistics()->write_last_phylogenetic_tree_statistics();});
    }
    tasks.wait();
  }
}

/**
 * \brief    Generate figures
 * \details  --
 * \param    char const** argv
 * \param    Simulation* simulation
 * \return   \e void
 */
void generate_figures( char const** argv, Simulation* simulation )
{
  if (simulation->get_population()->get_population_size() == 0)
  {
    return;
  }
  simulation->get_statistics()->plot_population_figures(std::string(argv[0]), EXECUTABLE_NAME);
  simulation->get_statistics()->plot_trophic_network_figures(std::string(argv[0]), EXECUTABLE_NAME);
  if (simulation->get_parameters()->get_simulation_backup_step() > 0 && simulation->get_phylogenetic_tree()->get_number_of_nodes() > 0)
  {
    simulation->get_statistics()->plot_phylogenetic_tree_figures(std::string(argv[0]), EXECUTABLE_NAME);
  }
}
