
/**
 * \file      SpeciesList.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     SpeciesList class definition
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

#include "SpeciesList.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
SpeciesList::SpeciesList( void )
{
  _X = new double[SPECIES_LIST_BUFFER];
  for (size_t i = 0; i < SPECIES_LIST_BUFFER; i++)
  {
    _X[i] = 0.0;
  }
  _size   = SPECIES_LIST_BUFFER;
  _amount = 0.0;
}

/**
 * \brief    Constructor from backup file
 * \details  Load SpeciesList class from backup file
 * \param    gzFile backup_file
 * \return   \e void
 */
SpeciesList::SpeciesList( gzFile backup_file )
{
  gzread( backup_file, &_size, sizeof(_size) );
  _X = new double[_size];
  for (size_t i = 0; i < _size; i++)
  {
    gzread( backup_file, &_X[i], sizeof(_X[i]) );
  }
  gzread( backup_file, &_amount, sizeof(_amount) );
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const SpeciesList& list
 * \return   \e void
 */
SpeciesList::SpeciesList( SpeciesList const& species_list )
{
  _X = new double[species_list._size];
  memcpy(_X, species_list._X, sizeof(double)*species_list._size);
  _size   = species_list._size;
  _amount = species_list._amount;
}

/*----------------------------
 * DESTRUCTORS
 *----------------------------*/

/**
 * \brief    Destructor
 * \details  --
 * \param    void
 * \return   \e void
 */
SpeciesList::~SpeciesList( void )
{
  delete[] _X;
  _X = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void SpeciesList::save( gzFile backup_file )
{
  gzwrite( backup_file, &_size, sizeof(_size) );
  for (size_t i = 0; i < _size; i++)
  {
    gzwrite( backup_file, &_X[i], sizeof(_X[i]) );
  }
  gzwrite( backup_file, &_amount, sizeof(_amount) );
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
