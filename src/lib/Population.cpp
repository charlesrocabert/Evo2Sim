
/**
 * \file      Population.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Population class definition
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

#include "Population.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
Population::Population( Parameters* parameters )
{
  _parameters = parameters;
  _width      = parameters->get_width();
  _height     = parameters->get_height();
  _grid       = new Cell*[_width*_height];
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _grid[pos] = new Cell(_parameters, _parameters->get_population_prng(pos));
  }
  _current_id      = 0;
  _best_id         = 0;
  _best_pos        = 0;
  _time            = 0;
  _population_size = 0;
  _previous_size   = 0;
  _growth_rate     = 0.0;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Population class from backup file
 * \param    Parameters* parameters
 * \param    gzFile backup_file
 * \return   \e void
 */
Population::Population( Parameters* parameters, gzFile backup_file )
{
  _parameters = parameters;
  gzread( backup_file, &_width,  sizeof(_width) );
  gzread( backup_file, &_height, sizeof(_height) );
  _grid   = new Cell*[_width*_height];
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    _grid[pos] = new Cell(_parameters, _parameters->get_population_prng(pos), backup_file);
  }
  gzread( backup_file, &_current_id,      sizeof(_current_id) );
  gzread( backup_file, &_best_id,         sizeof(_best_id) );
  gzread( backup_file, &_best_pos,        sizeof(_best_pos) );
  gzread( backup_file, &_time,            sizeof(_time) );
  gzread( backup_file, &_population_size, sizeof(_population_size) );
  gzread( backup_file, &_previous_size,   sizeof(_previous_size) );
  gzread( backup_file, &_growth_rate,     sizeof(_growth_rate) );
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
Population::~Population( void )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    delete _grid[i];
    _grid[i] = NULL;
  }
  delete[] _grid;
  _grid = NULL;
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
void Population::save( gzFile backup_file )
{
  gzwrite( backup_file, &_width,  sizeof(_width) );
  gzwrite( backup_file, &_height, sizeof(_height) );
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i]->save(backup_file);
  }
  gzwrite( backup_file, &_current_id,      sizeof(_current_id) );
  gzwrite( backup_file, &_best_id,         sizeof(_best_id) );
  gzwrite( backup_file, &_best_pos,        sizeof(_best_pos) );
  gzwrite( backup_file, &_time,            sizeof(_time) );
  gzwrite( backup_file, &_population_size, sizeof(_population_size) );
  gzwrite( backup_file, &_previous_size,   sizeof(_previous_size) );
  gzwrite( backup_file, &_growth_rate,     sizeof(_growth_rate) );
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
