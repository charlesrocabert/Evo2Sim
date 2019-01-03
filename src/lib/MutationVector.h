
/**
 * \file      MutationVector.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MutationVector class declaration
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

#ifndef __EVOEVO__MutationVector__
#define __EVOEVO__MutationVector__

#include <iostream>
#include <zlib.h>
#include <cmath>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Prng.h"


class MutationVector
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  MutationVector( void );
  MutationVector( gzFile backup_file );
  MutationVector( const MutationVector& vector );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~MutationVector( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline genetic_unit* get_dX( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  bool draw( Prng* prng, const double* mutation_rates );
  bool breakpoint_draw( Prng* prng, const double* mutation_rates );
  void save( gzFile backup_file );
  void clear( void );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void load_genetic_unit( gzFile backup_file, genetic_unit* unit );
  void save_genetic_unit( gzFile backup_file, genetic_unit* unit );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  genetic_unit* _dX; /*!< Mutation vector */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get dX vector
 * \details  --
 * \param    void
 * \return   \e genetic_unit*
 */
inline genetic_unit* MutationVector::get_dX( void ) const
{
  return _dX;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__MutationVector__) */
