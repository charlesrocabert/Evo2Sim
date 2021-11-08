
/**
 * \file      evo2sim_AB_run_competition_experiment.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      25-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Run a long-term competition experiments with A/B ecotypes
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
/* #include <tbb/tbb.h> */

#include "../lib/Macros.h"
#include "../lib/Enums.h"
#include "../lib/Structs.h"
#include "../lib/Parameters.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME         = "build/bin/evo2sim_AB_run_competition_experiment";
const size_t      DEFAULT_BACKUP_TIME     = 0;
const size_t      DEFAULT_SIMULATION_TIME = 10000;

void readArgs( int argc, char const** argv, size_t& backup, size_t& simulation_time, unsigned long int& seed, bool& mutations, bool& periodic, bool& continuous, bool& activate_graphic_display, bool& parallelism );
void printUsage( void );
void printHeader( void );
void write_last_backup_step( size_t last_backup, bool simulation_end );
void write_simulation_data( Simulation* simulation );
void write_SL_proportion( Simulation* simulation, std::ofstream& flux );
void generate_periodic_environment_properties( environment_properties* properties );
void generate_continuous_environment_properties( environment_properties* properties );


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
  
  /************************************/
  /* 1) Read command line arguments   */
  /************************************/
  size_t            backup                   = 0;
  size_t            simulation_time          = 0;
  unsigned long int seed                     = 0;
  bool              mutations                = false;
  bool              periodic                 = false;
  bool              continuous               = false;
  bool              activate_graphic_display = false;
  bool              parallelism              = false;
  readArgs(argc, argv, backup, simulation_time, seed, mutations, periodic, continuous, activate_graphic_display, parallelism);
  if (backup == 0)
  {
    backup = DEFAULT_BACKUP_TIME;
  }
  if (simulation_time == 0)
  {
    simulation_time = DEFAULT_SIMULATION_TIME;
  }
  if (seed == 0)
  {
    printf("You must give a seed. Exit.\n");
    exit(EXIT_FAILURE);
  }
  if ((periodic && continuous) || (!periodic && !continuous))
  {
    printf("You must choose between periodic environment OR continuous environment. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  /************************************/
  /* 2) Load parameters from backup   */
  /************************************/
  Parameters* parameters = new Parameters(backup);
  parameters->write("parameters.out");
  parameters->set_seed(seed);
  
  /************************************/
  /* 3) Load simulation from backup   */
  /************************************/
  Simulation* simulation = new Simulation(parameters, backup, true);
  
  /************************************/
  /* 4) Modify parameters             */
  /************************************/
  
  /*-------------------------*/
  /* 4.1) Set mutation rates */
  /*-------------------------*/
  if (!mutations)
  {
    parameters->set_point_mutation_rate(0.0);
    parameters->set_duplication_rate(0.0);
    parameters->set_deletion_rate(0.0);
    parameters->set_translocation_rate(0.0);
    parameters->set_inversion_rate(0.0);
    parameters->set_transition_rate(0.0);
    parameters->set_breakpoint_rate(0.0);
    for (size_t x = 0; x < parameters->get_width(); x++)
    {
      for (size_t y = 0; y < parameters->get_height(); y++)
      {
        Cell* cell = simulation->get_population()->get_cell(x, y);
        cell->set_point_mutation_rate(0.0);
        cell->set_duplication_rate(0.0);
        cell->set_deletion_rate(0.0);
        cell->set_translocation_rate(0.0);
        cell->set_inversion_rate(0.0);
        cell->set_transition_rate(0.0);
        cell->set_breakpoint_rate(0.0);
      }
    }
  }
  
  /*-------------------------*/
  /* 4.2) Set environment    */
  /*-------------------------*/
  if (periodic && !continuous)
  {
    environment_properties properties;
    generate_periodic_environment_properties(&properties);
    parameters->set_environment_properties(&properties);
  }
  else if (!periodic && continuous)
  {
    environment_properties properties;
    generate_continuous_environment_properties(&properties);
    parameters->set_environment_properties(&properties);
  }
  else
  {
    printf("Error cith periodic and continuous parameters. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  /*-------------------------*/
  /* 4.3) Set parallelism    */
  /*-------------------------*/
  parameters->set_parallel_computing(parallelism);
  
  /************************************/
  /* 5) Initialize simulation         */
  /************************************/
  simulation->initialize(FROM_BACKUP);
  simulation->get_phylogenetic_tree()->prune();
  simulation->get_phylogenetic_tree()->shorten();
  simulation->get_lineage_tree()->prune();
  simulation->get_lineage_tree()->shorten();
  
  std::ofstream AB_file("AB_proportion.txt", std::ios::out | std::ios::trunc);
  std::ofstream TN_file("AB_trophic_network.txt", std::ios::out | std::ios::trunc);
  write_SL_proportion(simulation, AB_file);
  simulation->get_trophic_network()->write_trophic_network(simulation->get_population()->get_time(), TN_file);
  
  write_simulation_data(simulation);
  
  /************************************/
  /* 6) Run simulation                */
  /************************************/
  double TIME = backup + simulation_time;
  bool   run  = true;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*---------------------------------------------------*/
    /* 6.1) Update simulation                            */
    /*---------------------------------------------------*/
    simulation->update();
    write_SL_proportion(simulation, AB_file);
    simulation->get_trophic_network()->write_trophic_network(simulation->get_population()->get_time(), TN_file);
    
    /*---------------------------------------------------*/
    /* 6.2) Save simulation                              */
    /*---------------------------------------------------*/
    if (parameters->get_simulation_backup_step() > 0 && simulation->get_population()->get_time() % parameters->get_simulation_backup_step() == 0)
    {
      simulation->save_simulation();
    }
    
    /*---------------------------------------------------*/
    /* 6.3) Stop if popsize = 0, else write stats        */
    /*---------------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      run = false;
    }
    else
    {
      simulation->get_statistics()->write_stats();
    }
    
    /*---------------------------------------------------*/
    /* 6.4) Print simulation report and generate figures */
    /*---------------------------------------------------*/
    size_t step = 500;
    if (parameters->get_figures_generation_step() > 0 && parameters->get_figures_generation_step() < 500)
    {
      step = parameters->get_figures_generation_step();
    }
    if (simulation->get_population()->get_time()%step == 0 && step == parameters->get_figures_generation_step())
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
    }
    
    /*---------------------------------------------------*/
    /* 6.5) Write signal file                            */
    /*---------------------------------------------------*/
    write_last_backup_step(simulation->get_population()->get_time(), false);
  }
  
  /************************************/
  /* 7) final output                  */
  /************************************/
  if (simulation->get_population()->get_population_size() > 0)
  {
    if (parameters->get_figures_generation_step() > 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
    }
  }
  else
  {
    if (parameters->get_figures_generation_step() > 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
    }
  }
  
  /************************************/
  /* 8) Write signal file             */
  /************************************/
  write_last_backup_step(simulation->get_population()->get_time(), true);
  
  /************************************/
  /* 9) close statistics files and    */
  /*     free memory                  */
  /************************************/
  AB_file.close();
  TN_file.close();
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
  delete parameters;
  parameters = NULL;
  
  int success = 0;
  success = system("zip statistics/AB_trophic_network.zip statistics/AB_trophic_network.txt");
  success = system("rm statistics/AB_trophic_network.txt");
  (void)success;
  return EXIT_SUCCESS;
}


/**
 * \brief    Read arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& backup
 * \param    size_t& simulation_time
 * \param    unsigned long int& seed
 * \param    bool& mutations
 * \param    bool& periodic
 * \param    bool& continuous
 * \param    bool& activate_graphic_display
 * \param    bool& parallelism
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& backup, size_t& simulation_time, unsigned long int& seed, bool& mutations, bool& periodic, bool& continuous, bool& activate_graphic_display, bool& parallelism )
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
    if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--backup") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--simulation-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        simulation_time = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-seed") == 0 || strcmp(argv[i], "--seed") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        seed = (unsigned long int)atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-mutations") == 0 || strcmp(argv[i], "--mutations") == 0)
    {
      mutations = true;
    }
    if (strcmp(argv[i], "-periodic") == 0 || strcmp(argv[i], "--periodic") == 0)
    {
      periodic = true;
    }
    if (strcmp(argv[i], "-continuous") == 0 || strcmp(argv[i], "--continuous") == 0)
    {
      continuous = true;
    }
    if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graphics") == 0)
    {
      activate_graphic_display = true;
    }
    if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--parallelism") == 0)
    {
      parallelism = true;
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
  std::cout << "Run a competition experiment between ecotypes A and B:\n";
  std::cout << "------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_run_competition_experiment -h or --help\n";
  std::cout << "   or: evo2sim_AB_run_competition_experiment [-b <backup-time>] [-t <simulation-time>] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -b, --backup\n";
  std::cout << "        specify the backup date (default: 0)\n";
  std::cout << "  -t, --simulation-time\n";
  std::cout << "        specify the duration of the simulation (default: 10000)\n";
  std::cout << "  -seed, --seed\n";
  std::cout << "        specify the prng seed (mandatory)\n";
  std::cout << "  -mutations, --mutations\n";
  std::cout << "        activate mutation rates\n";
  std::cout << "  -periodic, --periodic\n";
  std::cout << "        activate the periodic environment\n";
  std::cout << "  -continuous, --continuous\n";
  std::cout << "        activate the continuous environment\n";
  std::cout << "  -g, --graphics\n";
  std::cout << "        activate graphic display\n";
  std::cout << "  -p, --parallelism\n";
  std::cout << "        activate parallel computation\n";
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
  std::cout << "Run a competition experiment between ecotypes A and B.\n";
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
    /*
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
     */
  }
}

/**
 * \brief    Write the proportion of A/B ecotypes
 * \details  --
 * \param    Simulation* simulation
 * \param    std::ofstream& flux
 * \return   \e void
 */
void write_SL_proportion( Simulation* simulation, std::ofstream& flux )
{
  size_t nbA = simulation->get_trophic_network()->get_nb_level_0_cells()+simulation->get_trophic_network()->get_nb_level_1_cells();;
  size_t nbB = simulation->get_trophic_network()->get_nb_level_2_cells()+simulation->get_trophic_network()->get_nb_no_level_cells();
  flux << nbA << " " << nbB << "\n";
}

/**
 * \brief    Generate periodic environment properties
 * \details  --
 * \param    environment_properties* properties
 * \return   \e void
 */
void generate_periodic_environment_properties( environment_properties* properties )
{
  properties->number_of_init_cycles = 0;

  properties->species_tag_range.min = 10.0;
  properties->species_tag_range.max = 10.0;
  
  properties->concentration_range.min = 10.0;
  properties->concentration_range.max = 10.0;
  
  properties->number_of_species_range.min = 1.0;
  properties->number_of_species_range.max = 1.0;
    
  properties->interaction_scheme  = INTERACTION;
  properties->renewal_scheme      = CLEAR_MATTER;
  properties->variation_scheme    = PERIODIC_SCHEME;
  properties->localization_scheme = GLOBAL_LOCALIZATION;
  properties->metabolic_scheme    = UNIQUE_METABOLITE;
    
  properties->introduction_rate     = 0.003;
  properties->diffusion_coefficient = 0.1;
  properties->degradation_rate      = 0.0;
}

/**
 * \brief    Generate continuous environment properties
 * \details  --
 * \param    environment_properties* properties
 * \return   \e void
 */
void generate_continuous_environment_properties( environment_properties* properties )
{
  properties->number_of_init_cycles = 0;
  
  properties->species_tag_range.min = 10.0;
  properties->species_tag_range.max = 10.0;
  
  properties->concentration_range.min = 0.03;
  properties->concentration_range.max = 0.03;
  
  properties->number_of_species_range.min = 1.0;
  properties->number_of_species_range.max = 1.0;
  
  properties->interaction_scheme  = INTERACTION;
  properties->renewal_scheme      = KEEP_MATTER;
  properties->variation_scheme    = PERIODIC_SCHEME;
  properties->localization_scheme = GLOBAL_LOCALIZATION;
  properties->metabolic_scheme    = UNIQUE_METABOLITE;
  
  properties->introduction_rate     = 1.0;
  properties->diffusion_coefficient = 0.1;
  properties->degradation_rate      = 0.003;
}
