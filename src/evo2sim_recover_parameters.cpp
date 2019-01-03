
/**
 * \file      evo2sim_recover_parameters.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      05-04-2015
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Recover the parameters file from backup files
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

#include "../cmake/Config.h"

#include <iostream>
#include <cstring>
#include <assert.h>

#include "./lib/Macros.h"
#include "./lib/Parameters.h"
#include "./lib/Simulation.h"

const std::string EXECUTABLE_NAME = "build/bin/evo2sim_recover_parameters";

void readArgs( int argc, char const** argv, std::string filename );
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
  /**********************************/
  /* 1) Read command line arguments */
  /**********************************/
  std::string filename = "";
  readArgs(argc, argv, filename);
  
  /**********************************/
  /* 2) Recover the parameters file */
  /**********************************/
  printHeader();
  Parameters* parameters = new Parameters(0);
  parameters->write(filename);
  delete parameters;
  parameters = NULL;
  return EXIT_SUCCESS;
}


/**
 * \brief    Read arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    std::string filename
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string filename )
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
        filename = argv[i+1];
      }
    }
  }
  if (filename.compare("") == 0)
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
  std::cout << " Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " Web: https://github.com/charlesrocabert/Evo2Sim                           \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
  std::cout << "Recover the parameters file from backup files:\n";
  std::cout << "----------------------------------------------\n";
  std::cout << "Usage: recover_parameters -h or --help\n";
  std::cout << "   or: recover_parameters [-f <filename>] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit (optional)\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit (optional)\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify the name of the parameters file to save (mandatory)\n";
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
  std::cout << "Recover the parameters file from backup files.\n";
  std::cout << "\n";
}
