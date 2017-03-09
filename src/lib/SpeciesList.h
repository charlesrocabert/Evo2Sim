
/**
 * \file      SpeciesList.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     SpeciesList class declaration
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

#ifndef __EVOEVO__SpeciesList__
#define __EVOEVO__SpeciesList__

#include <iostream>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <cstring>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"


class SpeciesList
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  SpeciesList( void );
  SpeciesList( gzFile backup_file );
  SpeciesList( const SpeciesList& species_list );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~SpeciesList( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline double* get_X( void );
  inline double  get( int tag );
  inline size_t  get_size( void );
  inline double  get_amount( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  inline void set( int tag, double x );
  inline void add( int tag, double x );
  inline void remove( int tag, double x );
  inline void clear( void );
  inline void reset( bool metabolic_inheritance );
  
  inline void increase_size( size_t size );
  inline void decrease_size( size_t size );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  double* _X;      /*!< Vector of concentrations */
  size_t  _size;   /*!< Size                     */
  double  _amount; /*!< Total amount             */
  
};

/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get state vector X
 * \details  --
 * \param    void
 * \return   \e double*
 */
inline double* SpeciesList::get_X( void )
{
  return _X;
}

/**
 * \brief    Get
 * \details  Check if tag in [1,size] and return concentration
 * \param    int tag
 * \return   \e double
 */
inline double SpeciesList::get( int tag )
{
  assert(tag > 0);
  if (tag <= (int)_size)
  {
    return _X[tag-1];
  }
  return 0.0;
}

/**
 * \brief    Get species list size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t SpeciesList::get_size( void )
{
  return _size;
}

/**
 * \brief    Get amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double SpeciesList::get_amount( void )
{
  return _amount;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set
 * \details  Set species tag concentration. If tag is higher than species list size, increase it
 * \param    int tag
 * \param    double x
 * \return   \e void
 */
inline void SpeciesList::set( int tag, double x )
{
  assert(tag > 0);
  if (_size < (size_t)tag)
  {
    increase_size((size_t)tag);
  }
  _amount   -= _X[tag-1];
  _X[tag-1] = x;
  _amount   += x;
}

/**
 * \brief    Add
 * \details  Add concentration to species. If tag is higher than species list size, increase it
 * \param    int tag
 * \param    double x
 * \return   \e void
 */
inline void SpeciesList::add( int tag, double x )
{
  assert(tag > 0);
  if (_size < (size_t)tag)
  {
    increase_size((size_t)tag);
  }
  _amount   -= _X[tag-1];
  _X[tag-1] += x;
  _amount   += _X[tag-1];
}

/**
 * \brief    Remove
 * \details  Remove concentration from species. If tag is higher than species list size, increase it
 * \param    int tag
 * \param    double x
 * \return   \e void
 */
inline void SpeciesList::remove( int tag, double x )
{
  assert(tag > 0);
  if (_size < (size_t)tag)
  {
    increase_size((size_t)tag);
  }
  _amount   -= _X[tag-1];
  _X[tag-1] -= x;
  _amount   += _X[tag-1];
}

/**
 * \brief    Clear species list
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void SpeciesList::clear( void )
{
  delete[] _X;
  _X = new double[SPECIES_LIST_BUFFER];
  for (size_t i = 0; i < SPECIES_LIST_BUFFER; i++)
  {
    _X[i] = 0.0;
  }
  _size    = SPECIES_LIST_BUFFER;
  _amount  = 0.0;
}

/**
 * \brief    Reset the species list to zero
 * \details  --
 * \param    bool metabolic_inheritance
 * \return   \e void
 */
inline void SpeciesList::reset( bool metabolic_inheritance )
{
  /*---------------------------------------------------------*/
  /* A) If daughter cell inherits mother cytoplasmic content */
  /*---------------------------------------------------------*/
  if (metabolic_inheritance)
  {
    for (size_t i = 0; i < _size; i++)
    {
      _amount  -= _X[i];
      _X[i]    *= 0.5;
      _amount  += _X[i];
    }
  }
  /*---------------------------------------------------------*/
  /* B) Else if cytoplasm is reset to 0 at division          */
  /*---------------------------------------------------------*/
  else
  {
    for (size_t i = 0; i < _size; i++)
    {
      _X[i] = 0.0;
    }
    _amount  = 0.0;
  }
}

/**
 * \brief    Increase the size of the species list
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void SpeciesList::increase_size( size_t size )
{
  assert(size > _size);
  double* new_x = new double[size];
  memcpy(new_x, _X, sizeof(double)*_size);
  for (size_t i = _size; i < size; i++)
  {
    new_x[i] = 0.0;
  }
  delete[] _X;
  _X = new_x;
  _size = size;
}

/**
 * \brief    Decrease the size of the species list
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void SpeciesList::decrease_size( size_t size )
{
  assert(size < _size);
  double* new_x = new double[size];
  memcpy(new_x, _X, sizeof(double)*size);
  delete[] _X;
  _X = new_x;
  _size = size;
}


#endif /* defined(__EVOEVO__SpeciesList__) */
