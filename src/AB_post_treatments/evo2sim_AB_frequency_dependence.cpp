
/**
 * \file      evo2sim_AB_frequency_dependence.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      18-03-2016
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Measure the AB frequency dependent fitness by competition experiments
 */

/****************************************************************************
 * Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon
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

const std::string EXECUTABLE_NAME         = "build/bin/evo2sim_AB_frequency_dependence";
const std::string DEFAULT_FILENAME        = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH = "population_to_load";

void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path, size_t& simulation_time );
void printUsage( void );
void printHeader( void );
void create_folders( void );
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, bool identical, double B_prop, unsigned long int seed );
Simulation* create_simulation( Parameters* parameters, Population* evolved_population, unsigned long int seed, bool identical, double B_prop );
bool run_simulation( Simulation* simulation, size_t simulation_time );
void run_competition_experiment( Prng* prng, Parameters* parameters, Population* evolved_population, size_t reps, double B_prop, size_t simulation_time );
void measure_frequency_dependent_fitness( Parameters* parameters, Population* evolved_population, size_t reps, size_t simulation_time );


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
  size_t      rep                      = 0;
  std::string optional_filename        = "";
  std::string optional_population_path = "";
  size_t      simulation_time          = 0;
  readArgs(argc, argv, rep, optional_filename, optional_population_path, simulation_time);
  if (simulation_time == 0)
  {
    printf("--simulation-time (or -time) parameter is mandatory. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
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
  /* 3) Load the evolved population */
  /**********************************/
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
  
  /**********************************/
  /* 4) Run the post-treatment      */
  /**********************************/
  measure_frequency_dependent_fitness(parameters, evolved_population, rep, simulation_time);
  
  /**********************************/
  /* 5) Free the memory             */
  /**********************************/
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
 * \param    size_t& simulation_time
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path, size_t& simulation_time )
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
    if (strcmp(argv[i], "-time") == 0 || strcmp(argv[i], "--simulation-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        simulation_time = (size_t)atoi(argv[i+1]);
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
  std::cout << " Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
  std::cout << "Analyze the frequency-dependent fitness between ecotypes A and B:\n";
  std::cout << "-----------------------------------------------------------------\n";
  std::cout << "Usage: evo2sim_AB_frequency_dependence -h or --help\n";
  std::cout << "   or: evo2sim_AB_frequency_dependence [-rep <nb-repetitions>] [-f <param-file>] [-pop <population-file>] [-time <simulation-time>]\n";
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
  std::cout << "  -time, --simulation-time\n";
  std::cout << "        specify the simulation time(mandatory)\n";
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
  std::cout << " Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Analyze the frequency-dependent fitness between ecotypes A and B.\n";
  std::cout << "\n";
}

/**
 * \brief    Create folders
 * \details  --
 * \param    void
 * \return   \e void
 */
void create_folders( void )
{
  int success = 0;
  success = system("rm -rf statistics");
  mkdir("statistics", 0777);
  success = system("rm -rf environment");
  mkdir("environment", 0777);
  success = system("rm -rf population");
  mkdir("population", 0777);
  success = system("rm -rf trophic_network");
  mkdir("trophic_network", 0777);
  success = system("rm -rf tree");
  mkdir("tree", 0777);
  success = system("rm -rf parameters");
  mkdir("parameters", 0777);
  success = system("rm -rf prng");
  mkdir("prng", 0777);
  (void)success;
}

/**
 * \brief    Replace the individuals of the simulation by evolved individuals in the right proportions
 * \details  --
 * \param    Simulation* simulation
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    bool identical
 * \param    double B_prop
 * \param    unsigned long int seed
 * \return   \e void
 */
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, bool identical, double B_prop, unsigned long int seed )
{
  assert(B_prop >= 0.0);
  assert(B_prop <= 1.0);
  if (parameters->get_width() != evolved_population->get_width() || parameters->get_height() != evolved_population->get_height())
  {
    printf("Error: evolved population grid size is different from parameters file grid size. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  /*--------------------------------------------------------------*/
  /* 1) If identical boolean is true, just replace every cells    */
  /*--------------------------------------------------------------*/
  if (identical)
  {
    for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
    {
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(i));
    }
    return;
  }
  
  /*--------------------------------------------------------------*/
  /* 2) Isolate evolved identifiers thanks to their trophic level */
  /*--------------------------------------------------------------*/
  std::vector<size_t> A_indexes;
  std::vector<size_t> B_indexes;
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    if (evolved_population->get_cell(i)->isAlive() && (evolved_population->get_cell(i)->get_trophic_level() == LEVEL_0 || evolved_population->get_cell(i)->get_trophic_level() == LEVEL_1))
    {
      A_indexes.push_back(i);
    }
    else if (evolved_population->get_cell(i)->isAlive() && (evolved_population->get_cell(i)->get_trophic_level() == LEVEL_2 || evolved_population->get_cell(i)->get_trophic_level() == NO_LEVEL))
    {
      B_indexes.push_back(i);
    }
  }
  
  /*--------------------------------------------------------------*/
  /* 3) Replace cells in the right proportion                     */
  /*--------------------------------------------------------------*/
  Prng* prng = new Prng();
  prng->set_seed(seed);
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    if (prng->uniform() < B_prop)
    {
      size_t draw = (size_t)prng->uniform(0, (int)B_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(B_indexes[draw]));
    }
    else
    {
      size_t draw = (size_t)prng->uniform(0, (int)A_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(A_indexes[draw]));
    }
  }
  A_indexes.clear();
  B_indexes.clear();
  delete prng;
  prng = NULL;
}

/**
 * \brief    Create a simulation
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    unsigned long int seed
 * \param    bool identical
 * \param    double B_prop
 * \return   \e Simulation*
 */
Simulation* create_simulation( Parameters* parameters, Population* evolved_population, unsigned long int seed, bool identical, double B_prop )
{
  /*----------------------------------*/
  /* 1) Create simulation folders     */
  /*----------------------------------*/
  int success = 0;
  success = system("rm -rf environment");
  mkdir("environment", 0777);
  success = system("rm -rf population");
  mkdir("population", 0777);
  success = system("rm -rf trophic_network");
  mkdir("trophic_network", 0777);
  success = system("rm -rf tree");
  mkdir("tree", 0777);
  success = system("rm -rf parameters");
  mkdir("parameters", 0777);
  success = system("rm -rf prng");
  mkdir("prng", 0777);
  success = system("rm -rf statistics");
  mkdir("statistics", 0777);
  (void)success;
  
  /*----------------------------------*/
  /* 2) Initialize PRNG seed          */
  /*----------------------------------*/
  parameters->set_seed(seed);
  
  /*----------------------------------*/
  /* 3) Create a new simulation       */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters);
  
  /*----------------------------------*/
  /* 4) Load the evolved population   */
  /*----------------------------------*/
  replace_population(simulation, parameters, evolved_population, identical, B_prop, seed);
  
  /*----------------------------------*/
  /* 5) Set mutation rates            */
  /*----------------------------------*/
  for (size_t x = 0; x < parameters->get_width(); x++)
  {
    for (size_t y = 0; y < parameters->get_height(); y++)
    {
      Cell* cell = simulation->get_population()->get_cell(x, y);
      if (cell->isAlive())
      {
        cell->set_point_mutation_rate(parameters->get_point_mutation_rate());
        cell->set_duplication_rate(parameters->get_duplication_rate());
        cell->set_deletion_rate(parameters->get_deletion_rate());
        cell->set_translocation_rate(parameters->get_translocation_rate());
        cell->set_inversion_rate(parameters->get_inversion_rate());
        cell->set_transition_rate(parameters->get_transition_rate());
      }
    }
  }
  
  /*----------------------------------*/
  /* 6) initialize the population     */
  /*----------------------------------*/
  simulation->initialize_from_evolved_population();
  
  return simulation;
}

/**
 * \brief    Run a simulation
 * \details  --
 * \param    Simulation* simulation
 * \param    size_t simulation_time
 * \param    std::string filename
 * \return   \e bool
 */
bool run_simulation( Simulation* simulation, size_t simulation_time, std::string filename )
{
  /*---------------------------*/
  /* 1) Run simulation         */
  /*---------------------------*/
  std::ofstream AB_frequency(filename.c_str(), std::ios::out | std::ios::trunc);
  AB_frequency << "t A B\n";
  AB_frequency << simulation->get_trophic_network()->get_nb_level_0_cells()+simulation->get_trophic_network()->get_nb_level_1_cells() << " ";
  AB_frequency << simulation->get_trophic_network()->get_nb_level_2_cells()+simulation->get_trophic_network()->get_nb_no_level_cells() << "\n";
  double TIME    = simulation->get_population()->get_time()+simulation_time;
  bool   run     = true;
  bool   extinct = false;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    simulation->update();
    AB_frequency << simulation->get_trophic_network()->get_nb_level_0_cells()+simulation->get_trophic_network()->get_nb_level_1_cells() << " ";
    AB_frequency << simulation->get_trophic_network()->get_nb_level_2_cells()+simulation->get_trophic_network()->get_nb_no_level_cells() << "\n";
    if (simulation->get_population()->get_population_size() == 0)
    {
      extinct = true;
      run     = false;
    }
  }
  
  /*---------------------------*/
  /* 2) Evaluate final state   */
  /*---------------------------*/
  if (simulation->get_trophic_network()->get_nb_level_0_cells() == 0 && simulation->get_trophic_network()->get_nb_level_1_cells() == 0)
  {
    extinct = true;
  }
  if (simulation->get_trophic_network()->get_nb_level_2_cells() == 0 && simulation->get_trophic_network()->get_nb_no_level_cells() == 0)
  {
    extinct = true;
  }
  
  /*---------------------------*/
  /* 3) close statistics files */
  /*---------------------------*/
  simulation->get_statistics()->close_files();
  AB_frequency.close();
  
  return extinct;
}

/**
 * \brief    Run a competition experiment with repetitions
 * \details  --
 * \param    Prng* prng
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    size_t reps
 * \param    double B_prop
 * \param    size_t simulation_time
 * \return   \e void
 */
void run_competition_experiment( Prng* prng, Parameters* parameters, Population* evolved_population, size_t reps, double B_prop, size_t simulation_time )
{
  std::cout << "Running competition experiment with initial B frequency of " << B_prop << " ...\n";
  
  size_t rep  = 0;
  while (rep < reps)
  {
    std::cout << "> Rep " << rep << " ...\n";
    std::stringstream filename;
    filename << "AB_frequency_" << B_prop << "_" << rep << ".txt";
    Simulation* simulation = create_simulation(parameters, evolved_population, (unsigned long int)prng->uniform(1, 10000000), false, B_prop);
    run_simulation(simulation, simulation_time, filename.str());
    assert(simulation->get_population()->get_time() == simulation_time);
    delete simulation;
    simulation = NULL;
    rep++;
  }
}

/**
 * \brief    Measure the frequency-dependent fitness of S/L strains
 * \details  --
 * \param    Parameters* parameter
 * \param    Population* evolved_population
 * \param    size_t reps
 * \param    size_t simulation_time
 * \return   \e void
 */
void measure_frequency_dependent_fitness( Parameters* parameters, Population* evolved_population, size_t reps, size_t simulation_time )
{
  Prng* prng = new Prng();
  prng->set_seed(parameters->get_seed());
  
  /*-------------------------------------------------*/
  /* 1) evaluate B relative fitness at frequency 0.1 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.1, simulation_time);
  
  /*-------------------------------------------------*/
  /* 2) evaluate B relative fitness at frequency 0.2 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.2, simulation_time);
  
  /*-------------------------------------------------*/
  /* 3) evaluate B relative fitness at frequency 0.3 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.3, simulation_time);
  
  /*-------------------------------------------------*/
  /* 4) evaluate B relative fitness at frequency 0.4 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.4, simulation_time);
  
  /*-------------------------------------------------*/
  /* 5) evaluate B relative fitness at frequency 0.5 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.5, simulation_time);
  
  /*-------------------------------------------------*/
  /* 6) evaluate B relative fitness at frequency 0.6 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.6, simulation_time);
  
  /*-------------------------------------------------*/
  /* 7) evaluate B relative fitness at frequency 0.7 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.7, simulation_time);
  
  /*-------------------------------------------------*/
  /* 8) evaluate B relative fitness at frequency 0.8 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.8, simulation_time);
  
  /*-------------------------------------------------*/
  /* 9) evaluate B relative fitness at frequency 0.9 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, reps, 0.9, simulation_time);
  
  delete prng;
  prng = NULL;
}
