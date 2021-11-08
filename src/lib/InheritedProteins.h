
/**
 * \file      InheritedProteins.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     InheritedProteins class declaration
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

#ifndef __Evo2Sim__InheritedProteins__
#define __Evo2Sim__InheritedProteins__

#include <iostream>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <cstring>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Parameters.h"
#include "Genome.h"


class InheritedProteins
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  InheritedProteins( void ) = delete;
  InheritedProteins( Parameters* parameters );
  InheritedProteins( Parameters* parameters, gzFile backup_file );
  InheritedProteins( const InheritedProteins& inherited_proteins );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~InheritedProteins( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline genetic_sequence* get_genetic_sequence( void );
  inline genetic_unit*     get_genetic_unit( size_t pos );
  inline size_t            get_size( void ) const;
  inline size_t            get_buffer_size( void ) const;
  inline double*           get_concentration_vector( void );
  inline size_t            get_nb_E( void ) const;
  inline size_t            get_nb_TF( void ) const;
  inline size_t            get_nb_inner_enzymes( void ) const;
  inline size_t            get_nb_inflow_pumps( void ) const;
  inline size_t            get_nb_outflow_pumps( void ) const;
  inline size_t*           get_Ei( void );
  inline size_t*           get_TFi( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  inline void add_genetic_unit( genetic_unit* unit );
  inline void clear( void );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void initialize_concentration_vector( void );
  void build_index_list( void );
  void save( gzFile backup_file );
  void replace_data( InheritedProteins* inherited_proteins );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void create_genetic_sequence( void );
  void copy_genetic_sequence( const genetic_sequence* model );
  void delete_genetic_sequence( void );
  void load_genetic_sequence( gzFile backup_file );
  void save_genetic_sequence( gzFile backup_file );
  void load_genetic_unit( gzFile backup_file, genetic_unit& unit );
  void save_genetic_unit( gzFile backup_file, genetic_unit& unit );
  
  void increase_buffer_size( size_t new_size );
  void decrease_buffer_size( void );
  void check_inherited_proteins_size( void );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ simulation parameters */
  
  Parameters* _parameters; /*!< Simulation parameters */
  
  /*------------------------------------------------------------------ genetic sequence */
  
  genetic_sequence* _genetic_sequence; /*!< Genetic sequence */
  
  /*------------------------------------------------------------------ concentration vector */
  
  double* _concentration_vector; /*!< Concentration vector */
  
  /*------------------------------------------------------------------ statistical data */
  
  size_t _nb_E;             /*!< Number of enzyme types (E)                */
  size_t _nb_TF;            /*!< Number of transcription factor types (TF) */
  size_t _nb_inner_enzymes; /*!< Number of inner metabolism enzymes        */
  size_t _nb_inflow_pumps;  /*!< Number of inflowing pumps                 */
  size_t _nb_outflow_pumps; /*!< Number of outflowing pumps                */
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  size_t* _Ei;  /*!< Vector of E positions   */
  size_t* _TFi; /*!< Vector of TF positions  */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ genetic sequence */

/*
 * \brief    Get the genetic sequence
 * \details  --
 * \param    void
 * \return   \e genetic_sequence*
 */
inline genetic_sequence* InheritedProteins::get_genetic_sequence( void )
{
  return _genetic_sequence;
}

/*
 * \brief    Get genetic unit at position 'pos'
 * \details  --
 * \param    size_t pos
 * \return   \e genetic_unit*
 */
inline genetic_unit* InheritedProteins::get_genetic_unit( size_t pos )
{
  assert(pos < _genetic_sequence->size);
  return &_genetic_sequence->x[pos];
}

/**
 * \brief    Get inherited proteins size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_size( void ) const
{
  return _genetic_sequence->size;
}

/**
 * \brief    Get buffer size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_buffer_size( void ) const
{
  return _genetic_sequence->buffer_size;
}

/*------------------------------------------------------------------ concentration vector */

/**
 * \brief    Get the concentration vector
 * \details  --
 * \param    void
 * \return   \e double*
 */
inline double* InheritedProteins::get_concentration_vector( void )
{
  return _concentration_vector;
}

/*------------------------------------------------------------------ statistical data */

/**
 * \brief    Get number of enzyme types (E)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_nb_E( void ) const
{
  return _nb_E;
}

/**
 * \brief    Get number of transcription factor types (TF)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_nb_TF( void ) const
{
  return _nb_TF;
}

/**
 * \brief    Get number of inner metabolism enzymes
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_nb_inner_enzymes( void ) const
{
  return _nb_inner_enzymes;
}

/**
 * \brief    Get number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_nb_inflow_pumps( void ) const
{
  return _nb_inflow_pumps;
}

/**
 * \brief    Get number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t InheritedProteins::get_nb_outflow_pumps( void ) const
{
  return _nb_outflow_pumps;
}

/*------------------------------------------------------------------ genetic unit positions */

/**
 * \brief    Get the list of E type positions
 * \details  --
 * \param    void
 * \return   \e size_t*
 */
inline size_t* InheritedProteins::get_Ei( void )
{
  return _Ei;
}

/**
 * \brief    Get the list of TF type positions
 * \details  --
 * \param    void
 * \return   \e size_t*
 */
inline size_t* InheritedProteins::get_TFi( void )
{
  return _TFi;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Add a genetic unit at the end of the sequene
 * \details  --
 * \param    genetic_unit* unit
 * \return   \e void
 */
inline void InheritedProteins::add_genetic_unit( genetic_unit* unit )
{
  if (_genetic_sequence->size+1 > _genetic_sequence->buffer_size)
  {
    increase_buffer_size(_genetic_sequence->size+1);
  }
  memcpy(&_genetic_sequence->x[_genetic_sequence->size], unit, sizeof(genetic_unit));
  _genetic_sequence->size++;
}

/**
 * \brief    Clear inherited proteins sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void InheritedProteins::clear( void )
{
  delete[] _genetic_sequence->x;
  _genetic_sequence->x           = NULL;
  _genetic_sequence->x           = new genetic_unit[INHERITED_PROTEINS_BUFFER];
  _genetic_sequence->size        = 0;
  _genetic_sequence->buffer_size = INHERITED_PROTEINS_BUFFER;
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  delete[] _Ei;
  _Ei  = NULL;
  delete[] _TFi;
  _TFi = NULL;
}


#endif /* defined(__Evo2Sim__InheritedProteins__) */
