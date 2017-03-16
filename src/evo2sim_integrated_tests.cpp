
/**
 * \file      evo2sim_integrated_tests.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Run integrated tests
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
#include "./lib/IntegratedTests.h"

const std::string EXECUTABLE_NAME  = "build/bin/evo2sim_integrated_tests";
const std::string DEFAULT_FILENAME = "parameters.txt";

void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& number_of_tests, size_t& number_of_steps, bool& random_seed, bool& random_parameters );
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
  /************************************/
  /* 1) Read command line arguments   */
  /************************************/
  std::string optional_filename = "";
  size_t      number_of_tests   = 1;
  size_t      number_of_steps   = 1;
  bool        random_seed       = false;
  bool        random_parameters = false;
  readArgs(argc, argv, optional_filename, number_of_tests, number_of_steps, random_seed, random_parameters);
  
  /************************************/
  /* 2) Load parameters from file     */
  /************************************/
  printHeader();
  Parameters* parameters1 = new Parameters();
  Parameters* parameters2 = new Parameters();
  if (strcmp(optional_filename.c_str(), "") != 0)
  {
    parameters1->load_parameters_from_file(optional_filename);
    parameters2->load_parameters_from_file(optional_filename);
  }
  else
  {
    parameters1->load_parameters_from_file(DEFAULT_FILENAME);
    parameters2->load_parameters_from_file(DEFAULT_FILENAME);
  }
  
  /************************************/
  /* 3) Load and run integrated tests */
  /************************************/
  IntegratedTests* integrated_tests = new IntegratedTests(parameters1, parameters2);
  integrated_tests->run_integrated_tests(number_of_tests, number_of_steps, random_seed, random_parameters);
  
  /************************************/
  /* 4) Free memory                   */
  /************************************/
  delete integrated_tests;
  integrated_tests = NULL;
  delete parameters1;
  parameters1 = NULL;
  delete parameters2;
  parameters2 = NULL;
  return EXIT_SUCCESS;
}


/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    string& optional_file
 * \param    size_t& number_of_tests
 * \param    size_t& number_of_steps
 * \param    bool& random_seed
 * \param    bool& random_parameters
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& number_of_tests, size_t& number_of_steps, bool& random_seed, bool& random_parameters )
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
        std::cout << "Error: -f command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_filename = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-tests") == 0 || strcmp(argv[i], "--number-of-tests") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: -tests command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        number_of_tests = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-steps") == 0 || strcmp(argv[i], "--number-of-steps") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: -steps command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        number_of_steps = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-rs") == 0 || strcmp(argv[i], "--random-seed") == 0)
    {
      random_seed = true;
    }
    if (strcmp(argv[i], "-rp") == 0 || strcmp(argv[i], "--random-parameters") == 0)
    {
      random_parameters = true;
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
  std::cout << "Run integrated tests:\n";
  std::cout << "---------------------\n";
  std::cout << "Usage: integrated_tests -h or --help\n";
  std::cout << "   or: integrated_tests [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit (optional)\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit (optional)\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify the parameters file (default: parameters.txt)\n";
  std::cout << "  -tests, --number-of-tests\n";
  std::cout << "        specify the number of tests with different seeds (default: 1)\n";
  std::cout << "  -steps, --number-of-steps\n";
  std::cout << "        specify the number of steps by test (default: 1)\n";
  std::cout << "  -rs, --random-seed\n";
  std::cout << "        the prng seed is drawn at random for each test (optional)\n";
  std::cout << "  -rp, --random-parameters\n";
  std::cout << "        the parameters are drawn at random for each test (optional)\n\n";
  std::cout << "To use the integrated tests, the software must be compiled in DEBUG mode (see INSTALL).\n";
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
  std::cout << "Run integrated tests.\n";
  std::cout << "\n";
}
