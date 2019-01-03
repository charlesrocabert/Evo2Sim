
/**
 * \file      Population.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Population class declaration
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

#ifndef __EVOEVO__Population__
#define __EVOEVO__Population__

#include <iostream>
#include <zlib.h>
#include <vector>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Cell.h"


class Population
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Population( void ) = delete;
  Population( Parameters* parameters );
  Population( Parameters* parameters, gzFile backup_file );
  Population( const Population& population ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Population( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline size_t                 get_width( void ) const;
  inline size_t                 get_height( void ) const;
  inline Cell*                  get_cell( size_t i );
  inline Cell*                  get_cell( size_t x, size_t y );
  inline Cell*                  get_cell_by_id( unsigned long long int identifier );
  inline Cell*                  get_best_cell( void );
  inline size_t                 get_time( void ) const;
  inline size_t                 get_population_size( void ) const;
  inline double                 get_growth_rate( void ) const;
  inline unsigned long long int get_new_id( void );
  inline unsigned long long int get_best_id( void ) const;
  inline size_t                 get_best_position( void ) const;
  inline double                 get_total_amount( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Population& operator=(const Population&) = delete;
  
  inline void set_cell( size_t i, Cell* cell );
  inline void new_cell( size_t parent_position, size_t child_position );
  inline void set_population_size( size_t size );
  inline void set_previous_size( void );
  inline void compute_growth_rate( void );
  inline void increase_population_size( size_t size );
  inline void decrease_population_size( size_t size );
  inline void set_time( size_t time );
  inline void update_time( void );
  inline void set_best( unsigned long long int identifier, size_t i );

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
  Parameters*            _parameters;      /*!< Simulation parameters           */
  Cell**                 _grid;            /*!< population's grid               */
  size_t                 _width;           /*!< Grid's width                    */
  size_t                 _height;          /*!< Grid's height                   */
  unsigned long long int _current_id;      /*!< Current cell id                 */
  unsigned long long int _best_id;         /*!< Best identifier                 */
  size_t                 _best_pos;        /*!< Position of the best individual */
  size_t                 _time;            /*!< Evolution time                  */
  size_t                 _population_size; /*!< Population size                 */
  size_t                 _previous_size;   /*!< Previous population size        */
  double                 _growth_rate;     /*!< Population growth rate          */
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get grid width
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Population::get_width( void ) const
{
  return _width;
}

/**
 * \brief    Get grid height
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Population::get_height( void ) const
{
  return _height;
}

/**
 * \brief    Get cell on position i
 * \details  --
 * \param    size_t i
 * \return   \e Cell*
 */
inline Cell* Population::get_cell( size_t i )
{
  assert(i < _width*_height);
  return _grid[i];
}

/**
 * \brief    Get cell on coordinates (x, y)
 * \details  --
 * \param    size_t x
 * \param    size_t y
 * \return   \e Cell*
 */
inline Cell* Population::get_cell( size_t x, size_t y )
{
  assert(x < _width);
  assert(y < _height);
  return _grid[x*_height+y];
}

/**
 * \brief    Get cell by its identifier
 * \details  This method should only be used during tree loading from backup
 * \param    unsigned long long int identifier
 * \return   \e Cell*
 */
inline Cell* Population::get_cell_by_id( unsigned long long int identifier )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    if (_grid[i]->get_id() == identifier)
    {
      assert(_grid[i]->isAlive());
      return _grid[i];
    }
  }
  return NULL;
}

/**
 * \brief    Get best cell
 * \details  --
 * \param    void
 * \return   \e Cell*
 */
inline Cell* Population::get_best_cell( void )
{
  return _grid[_best_pos];
}

/**
 * \brief    Get time
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Population::get_time( void ) const
{
  return _time;
}

/**
 * \brief    Get population size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Population::get_population_size( void ) const
{
  return _population_size;
}

/**
 * \brief    Get population growth rate
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline double Population::get_growth_rate( void ) const
{
  return _growth_rate;
}

/**
 * \brief    Get new id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Population::get_new_id( void )
{
  _current_id++;
  return _current_id;
}

/**
 * \brief    Get best id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Population::get_best_id( void ) const
{
  return _best_id;
}

/**
 * \brief    Get best position
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Population::get_best_position( void ) const
{
  return _best_pos;
}

/**
 * \brief    Get the total amount of metabolites in the population
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Population::get_total_amount( void ) const
{
  double total_amount = 0.0;
  for (size_t i = 0; i < _width*_height; i++)
  {
    total_amount += _grid[i]->get_species_list()->get_amount();
  }
  return total_amount;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Set cell at position i
 * \details  --
 * \param    size_t i
 * \param    Cell* cell
 * \return   \e void
 */
inline void Population::set_cell( size_t i, Cell *cell )
{
  assert(i < _width*_height);
  _grid[i] = cell;
}

/**
 * \brief    Create a new cell from its parent on position child_position
 * \details  --
 * \param    size_t parent_position
 * \param    size_t child_position
 * \return   \e void
 */
inline void Population::new_cell( size_t parent_position, size_t child_position )
{
  assert(parent_position < _width*_height);
  assert(child_position < _width*_height);
  assert(parent_position != child_position);
  size_t y = child_position%_height;
  size_t x = (child_position-y)/_height;
  delete _grid[child_position];
  _grid[child_position] = NULL;
  _grid[child_position] = new Cell(*_grid[parent_position], _parameters->get_population_prng(child_position), get_new_id(), x, y, _time);
}

/**
 * \brief    Set population size
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void Population::set_population_size( size_t size )
{
  _population_size = size;
}

/**
 * \brief    Set previous population size
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Population::set_previous_size( void )
{
  _previous_size = _population_size;
}

/**
 * \brief    Compute growth rate
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Population::compute_growth_rate( void )
{
  _growth_rate = (double)_population_size-(double)_previous_size;
}

/**
 * \brief    Increase population size by "size"
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void Population::increase_population_size( size_t size )
{
  _population_size += size;
}

/**
 * \brief    Decrease population size by "size"
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void Population::decrease_population_size( size_t size )
{
  _population_size -= size;
}

/**
 * \brief    Set time
 * \details  --
 * \param    size_t time
 * \return   \e void
 */
inline void Population::set_time( size_t time )
{
  _time = time;
}

/**
 * \brief    Update time
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Population::update_time( void )
{
  _time++;
}

/**
 * \brief    Set best id
 * \details  --
 * \param    unsigned long long int identifier
 * \param    size_t i
 * \return   \e void
 */
inline void Population::set_best( unsigned long long int identifier, size_t i )
{
  _best_id  = identifier;
  _best_pos = i;
}


#endif /* defined(__EVOEVO__Population__) */
