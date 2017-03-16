
/**
 * \file      evo2sim_run.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Run a simulation from backup files
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
#include <tbb/tbb.h>

#include "./lib/Macros.h"
#include "./lib/Parameters.h"
#include "./lib/Simulation.h"

#if WITH_GRAPHICS_CONTEXT
  #include <SFML/Graphics.hpp>
  #include "./lib/GraphicDisplay.h"
#endif

const std::string EXECUTABLE_NAME         = "build/bin/evo2sim_run";
const size_t      DEFAULT_BACKUP_TIME     = 0;
const size_t      DEFAULT_SIMULATION_TIME = 10000;

void readArgs( int argc, char const** argv, size_t& backup_time, size_t& simulation_time, bool& activate_graphic_display );
void printUsage( void );
void printHeader( void );
void write_last_backup_step( size_t last_backup, bool simulation_end );
void write_simulation_data( Simulation* simulation );
void generate_figures( char const** argv, Simulation* simulation );
#if WITH_GRAPHICS_CONTEXT
  void create_graphic_context( Parameters* parameters, sf::RenderWindow* pop_window, sf::RenderWindow* env_window );
#endif


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
  size_t backup_time              = 0;
  size_t simulation_time          = 0;
  bool   activate_graphic_display = false;
  readArgs(argc, argv, backup_time, simulation_time, activate_graphic_display);
  if (backup_time == 0)
  {
    backup_time = DEFAULT_BACKUP_TIME;
  }
  if (simulation_time == 0)
  {
    simulation_time = DEFAULT_SIMULATION_TIME;
  }
  
  /**********************************/
  /* 2) Load parameters from backup */
  /**********************************/
  printHeader();
  Parameters* parameters = new Parameters(backup_time);
  parameters->write("parameters.out");
  
  /**********************************/
  /* 3) Load simulation from backup */
  /*    and graphic display         */
  /**********************************/
  Simulation* simulation = new Simulation(parameters, backup_time, true);
#if WITH_GRAPHICS_CONTEXT
  GraphicDisplay* graphic_display = NULL;
  sf::RenderWindow pop_window;
  sf::RenderWindow env_window;
  std::string font_path = std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/fonts/andale_mono.ttf";
#endif
  
  /**********************************/
  /* 4) Open render windows         */
  /**********************************/
#if WITH_GRAPHICS_CONTEXT
  if (activate_graphic_display)
  {
    create_graphic_context(parameters, &pop_window, &env_window);
    graphic_display = new GraphicDisplay(parameters, simulation, &pop_window, &env_window, font_path);
  }
#endif
  
  /**********************************/
  /* 5) Initialize simulation       */
  /**********************************/
  simulation->initialize( FROM_BACKUP );
  write_simulation_data(simulation);
  if (parameters->get_figures_generation_step() > 0)
  {
    generate_figures(argv, simulation);
  }
  
  /**********************************/
  /* 6) Run simulation              */
  /**********************************/
  double TIME        = backup_time + simulation_time;
  bool   extinction  = true;
  bool   graphics    = true;
  size_t last_backup = backup_time;
  write_last_backup_step(last_backup, false);
  
  while (simulation->get_population()->get_time() < TIME && extinction && graphics)
  {
    /*---------------------------------------------------*/
    /* 6.1) Save best cell state in track.txt            */
    /*---------------------------------------------------*/
    simulation->get_statistics()->write_best_cell_state();
    
    /*---------------------------------------------------*/
    /* 6.2) Update simulation                            */
    /*---------------------------------------------------*/
    simulation->update();
    simulation->get_statistics()->write_stats();
    if (simulation->get_population()->get_population_size() == 0)
    {
      extinction = false;
    }
    
    /*---------------------------------------------------*/
    /* 6.3) Display graphics                             */
    /*---------------------------------------------------*/
#if WITH_GRAPHICS_CONTEXT
    if (activate_graphic_display)
    {
      graphics = graphic_display->display();
    }
#endif
    
    /*---------------------------------------------------*/
    /* 6.4) Save simulation                              */
    /*---------------------------------------------------*/
    if (parameters->get_simulation_backup_step() > 0 && simulation->get_population()->get_time() % parameters->get_simulation_backup_step() == 0)
    {
      simulation->save_simulation();
      last_backup = simulation->get_population()->get_time();
    }
    
    /*---------------------------------------------------*/
    /* 6.5) Print simulation report and generate figures */
    /*---------------------------------------------------*/
    if (simulation->get_population()->get_time()%500 == 0)
    {
      std::cout << "-------------- Simulation report ----------------\n";
      std::cout << " Elapsed time        : " << simulation->get_population()->get_time() << "\n";
      std::cout << " Elapsed generations : " << simulation->get_statistics()->_mean_generations << "\n";
      std::cout << " Population size     : " << simulation->get_population()->get_population_size() << "\n";
      std::cout << " Best score          : " << simulation->get_max_score() << "\n";
      std::cout << "-------------------------------------------------\n";
    }
    if (parameters->get_figures_generation_step() > 0 && simulation->get_population()->get_time()%parameters->get_figures_generation_step() == 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
      generate_figures(argv, simulation);
    }
    write_last_backup_step(last_backup, false);
  }
  
  /**********************************/
  /* 7) final output                */
  /**********************************/
  std::cout << "-------------- End of simulation ----------------\n";
  std::cout << "                                                 \n";
  if (simulation->get_population()->get_population_size() > 0)
  {
    std::cout << " End of simulation at " << simulation->get_population()->get_time() << " timesteps.\n";
  }
  else
  {
    std::cout << " Population extinction at " << simulation->get_population()->get_time() << " timesteps.\n";
  }
  std::cout << "                                                 \n";
  std::cout << "-------------------------------------------------\n\n";
  if (parameters->get_figures_generation_step() > 0)
  {
    simulation->get_statistics()->flush_files();
    write_simulation_data(simulation);
    generate_figures(argv, simulation);
  }
  write_last_backup_step(last_backup, true);
  
  /**********************************/
  /* 8) close statistics files and  */
  /*    free memory                 */
  /**********************************/
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
#if WITH_GRAPHICS_CONTEXT
  if (activate_graphic_display)
  {
    delete graphic_display;
    graphic_display = NULL;
  }
#endif
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
 * \param    size_t& simulation_time
 * \param    bool& activate_graphic_display
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& backup_time, size_t& simulation_time, bool& activate_graphic_display )
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
    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--simulation-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: -t command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        simulation_time = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graphics") == 0)
    {
      activate_graphic_display = true;
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
  std::cout << " Evo2Sim is a multi-scale, individual-based computational model dedicated  \n";
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
  std::cout << "Run a simulation from backup files:\n";
  std::cout << "-----------------------------------\n";
  std::cout << "Usage: run -h or --help\n";
  std::cout << "   or: run [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit (optional)\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit (optional)\n";
  std::cout << "  -b, --backup-time\n";
  std::cout << "        set the date of the backup to load (default: 0)\n";
  std::cout << "  -t, --simulation-time\n";
  std::cout << "        set the duration of the simulation (default: 10000)\n";
  std::cout << "  -g, --graphics\n";
  std::cout << "        activate graphic display (optional)\n\n";
  std::cout << "Statistic files content is automatically managed when a simulation is reloaded from backup to avoid data loss.\n";
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
  std::cout << " Evo2Sim is a multi-scale, individual-based computational model dedicated  \n";
  std::cout << " to in silico experimental evolution.                                      \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Run a simulation from backup files.\n";
  std::cout << "\n";
}

/**
 * \brief    Write the last backup step and the simulation state
 * \details  --
 * \param    size_t last_backup
 * \param    bool simulation_end
 * \return   \e void
 */
void write_last_backup_step( size_t last_backup, bool simulation_end )
{
  std::ofstream file("last_backup.txt", std::ios::out | std::ios::trunc);
  if (!simulation_end)
  {
    file << last_backup << " wait\n";
  }
  else
  {
    file << last_backup << " end\n";
  }
  file.close();
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

/**
 * \brief    Create the graphic context
 * \details  Basically, open graphic windows
 * \param    Parameters* parameters
 * \param    sf::RenderWindow* pop_window
 * \param    sf::RenderWindow* env_window
 * \return   \e void
 */
#if WITH_GRAPHICS_CONTEXT
void create_graphic_context( Parameters* parameters, sf::RenderWindow* pop_window, sf::RenderWindow* env_window )
{
  int x_inter = ((int)parameters->get_width()-1)*CELL_SPACE;
  int y_inter = ((int)parameters->get_height()-1)*CELL_SPACE;
  int width   = (int)parameters->get_width()*CELL_SCALE + x_inter + 4*SPAN + GRADIENT_SCALE + TEXT_SCALE;
  int height  = (int)parameters->get_height()*CELL_SCALE + y_inter + 2*SPAN;
  
  pop_window->create(sf::VideoMode(width, height), "Population");
  env_window->create(sf::VideoMode(width, height), "Environment");
  
  pop_window->setPosition(sf::Vector2i(100, 100));
  env_window->setPosition(sf::Vector2i(100+width+50, 100));
  
  if (FRAMERATE > 0)
  {
    pop_window->setFramerateLimit(FRAMERATE);
    env_window->setFramerateLimit(FRAMERATE);
  }
}
#endif
