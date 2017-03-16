
/**
 * \file      Environment.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Environment class declaration
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

#ifndef __EVOEVO__Environment__
#define __EVOEVO__Environment__

#include <iostream>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "SpeciesList.h"
#include "Prng.h"


class Environment
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Environment( void ) = delete;
  Environment( Parameters* parameters );
  Environment( Parameters* parameters, gzFile backup_file );
  Environment( const Environment& environment );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Environment( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline double       get( size_t x, size_t y, int tag );
  inline SpeciesList* get_species_list( size_t x, size_t y );
  inline SpeciesList* get_X( void );
  inline size_t       get_width( void ) const;
  inline size_t       get_height( void ) const;
  inline size_t       get_size( void ) const;
  inline size_t       get_size( size_t x, size_t y ) const;
  inline double       get_total_amount( void ) const;
  inline double       get_min_amount( void ) const;
  inline double       get_max_amount( void ) const;
  inline double       get_inflowing_amount( void ) const;
  inline double       get_outflowing_amount( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  inline void set( size_t x, size_t y, int tag, double concentration );
  inline void add( size_t x, size_t y, int tag, double concentration );
  inline void remove( size_t x, size_t y, int tag, double concentration );
  inline void increase_size( size_t x, size_t y, int tag );
  inline void set_min_amount( double min_amount );
  inline void set_max_amount( double max_amount );
  inline void set_inflowing_amount( double inflowing_amount );
  inline void set_outflowing_amount( double outflowing_amount );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void actualize_environment_state( void );
  void compute_diffusion_and_degradation( void );
  void save( gzFile backup_file );
  void write_metabolic_state_vector( std::string filename );
  void write_local_metabolic_state_vector( std::string filename, size_t pos );
  void update( bool initialize, size_t time );
  void clear( void );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void compute_classic_diffusion_and_degradation( void );
  void compute_perfect_diffusion_and_degradation( void );
  
  void update_global_unique( bool initialize, double factor );
  void update_global_multiple( bool initialize, double factor );
  void update_random_unique( bool initialize, double factor );
  void update_random_multiple( bool initialize, double factor );
  void update_spot_unique( bool initialize, double factor );
  void update_spot_multiple( bool initialize, double factor );
  void update_center_unique( bool initialize, double factor );
  void update_center_multiple( bool initialize, double factor );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ parameters */
  
  Parameters*             _parameters; /*!< Simulation parameters  */
  Prng*                   _prng;       /*!< Prng                   */
  size_t                  _width;      /*!< Grid's width           */
  size_t                  _height;     /*!< Grid's height          */
  environment_properties* _properties; /*!< Environment properties */
  
  /*------------------------------------------------------------------ environment structure */
  
  SpeciesList** _grid;              /*!< Grid                        */
  SpeciesList*  _X;                 /*!< Metabolic state vector      */
  size_t        _size;              /*!< Species lists common size   */
  double        _total_amount;      /*!< Total amount                */
  double        _min_amount;        /*!< Minimum amount              */
  double        _max_amount;        /*!< Maximum amount              */
  double        _inflowing_amount;  /*!< Amount of inflowing matter  */
  double        _outflowing_amount; /*!< Amount of outflowing matter */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get specified metabolite tag concentration
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \param    int tag
 * \return   \e double
 */
inline double Environment::get( size_t x, size_t y, int tag )
{
  assert(x < _width);
  assert(y < _height);
  assert(tag > 0);
  if (tag > (int)_grid[x*_height+y]->get_size())
  {
    return 0.0;
  }
  return _grid[x*_height+y]->get(tag);
}

/**
 * \brief    Get species list at position (x, y)
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \return   \e SpeciesList*
 */
inline SpeciesList* Environment::get_species_list( size_t x, size_t y )
{
  assert(x < _width);
  assert(y < _height);
  return _grid[x*_height+y];
}

/**
 * \brief    Get species lists common size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Environment::get_size( void ) const
{
  return _size;
}

/**
 * \brief    Get species lists size
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \return   \e size_t
 */
inline size_t Environment::get_size( size_t x, size_t y ) const
{
  return _grid[x*_height+y]->get_size();
}

/**
 * \brief    Get metabolic state vector X
 * \details  --
 * \param    void
 * \return   \e SpeciesList*
 */
inline SpeciesList* Environment::get_X( void )
{
  return _X;
}

/**
 * \brief    Get grid width
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Environment::get_width( void ) const
{
  return _width;
}

/**
 * \brief    Get grid height
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Environment::get_height( void ) const
{
  return _height;
}

/**
 * \brief    Get total amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Environment::get_total_amount( void ) const
{
  return _total_amount;
}

/**
 * \brief    Get minimum amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Environment::get_min_amount( void ) const
{
  return _min_amount;
}

/**
 * \brief    Get maximum amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Environment::get_max_amount( void ) const
{
  return _max_amount;
}

/**
 * \brief    Get inflowing amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Environment::get_inflowing_amount( void ) const
{
  return _inflowing_amount;
}

/**
 * \brief    Get outflowing amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Environment::get_outflowing_amount( void ) const
{
  return _outflowing_amount;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \param    int tag
 * \param    double concentration
 * \return   \e void
 */
inline void Environment::set( size_t x, size_t y, int tag, double concentration )
{
  assert(x < _width);
  assert(y < _height);
  assert(tag > 0);
  if (_grid[x*_height+y]->get_size() < (size_t)tag)
  {
    _grid[x*_height+y]->increase_size((size_t)tag);
  }
  _grid[x*_height+y]->set(tag, concentration);
}

/**
 * \brief    Add
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \param    int tag
 * \param    double concentration
 * \return   \e void
 */
inline void Environment::add( size_t x, size_t y, int tag, double concentration )
{
  assert(x < _width);
  assert(y < _height);
  assert(tag > 0);
  if (_grid[x*_height+y]->get_size() < (size_t)tag)
  {
    _grid[x*_height+y]->increase_size((size_t)tag);
  }
  _grid[x*_height+y]->add(tag, concentration);
}

/**
 * \brief    Remove
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \param    int tag
 * \param    double concentration
 * \return   \e void
 */
inline void Environment::remove( size_t x, size_t y, int tag, double concentration )
{
  assert(x < _width);
  assert(y < _height);
  assert(tag > 0);
  if (_grid[x*_height+y]->get_size() < (size_t)tag)
  {
    _grid[x*_height+y]->increase_size((size_t)tag);
  }
  _grid[x*_height+y]->remove(tag, concentration);
}

/**
 * \brief    Increase species list size
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \param    int tag
 * \return   \e void
 */
inline void Environment::increase_size( size_t x, size_t y, int tag )
{
  assert(x < _width);
  assert(y < _height);
  assert(tag > 0);
  if (_grid[x*_height+y]->get_size() < (size_t)tag)
  {
    _grid[x*_height+y]->increase_size((size_t)tag);
  }
}

/**
 * \brief    Set the miminum amount
 * \details  --
 * \param    double min_amount
 * \return   \e void
 */
inline void Environment::set_min_amount( double min_amount )
{
  assert(min_amount >= 0.0);
  _min_amount = min_amount;
}

/**
 * \brief    Set the maximum amount
 * \details  --
 * \param    double max_amount
 * \return   \e void
 */
inline void Environment::set_max_amount( double max_amount )
{
  assert(max_amount >= 0.0);
  _max_amount = max_amount;
}

/**
 * \brief    Set the inflowing amount
 * \details  --
 * \param    double inflowing_amount
 * \return   \e void
 */
inline void Environment::set_inflowing_amount( double inflowing_amount )
{
  assert(inflowing_amount >= 0.0);
  _inflowing_amount = inflowing_amount;
}

/**
 * \brief    Set the outflowing amount
 * \details  --
 * \param    double outflowing_amount
 * \return   \e void
 */
inline void Environment::set_outflowing_amount( double outflowing_amount )
{
  assert(outflowing_amount >= 0.0);
  _outflowing_amount = outflowing_amount;
}


#endif /* defined(__EVOEVO__Environment__) */
