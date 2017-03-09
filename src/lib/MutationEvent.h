
/**
 * \file      MutationEvent.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MutationEvent class declaration
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

#ifndef __EVOEVO__MutationEvent__
#define __EVOEVO__MutationEvent__

#include <iostream>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "MutationVector.h"


class MutationEvent
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  MutationEvent( void ) = delete;
  MutationEvent( mutation_type type, size_t point_mutation_location, MutationVector* vector );
  MutationEvent( mutation_type type, size_t hgt_insert, size_t size, size_t nb_NC, size_t nb_E, size_t nb_TF, size_t nb_BS, size_t nb_P );
  MutationEvent( mutation_type type, size_t src1, size_t src2, size_t tgt, size_t size, size_t nb_NC, size_t nb_E, size_t nb_TF, size_t nb_BS, size_t nb_P );
  MutationEvent( gzFile backup_file );
  MutationEvent( const MutationEvent& event );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~MutationEvent( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline mutation_type   get_mutation_type( void ) const;
  inline size_t          get_point_mutation_location( void ) const;
  inline size_t          get_hgt_insert( void ) const;
  inline size_t          get_nb_NC( void ) const;
  inline size_t          get_nb_E( void ) const;
  inline size_t          get_nb_TF( void ) const;
  inline size_t          get_nb_BS( void ) const;
  inline size_t          get_nb_P( void ) const;
  inline size_t          get_src_breakpoint1( void ) const;
  inline size_t          get_src_breakpoint2( void ) const;
  inline size_t          get_tgt_breakpoint( void ) const;
  inline size_t          get_size( void ) const;
  inline MutationVector* get_mutation_vector( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
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
  mutation_type   _mutation_type;           /*!< Mutation type                     */
  size_t          _point_mutation_location; /*!< Point mutation location on genome */
  size_t          _hgt_insert;              /*!< HGT insertion point               */
  size_t          _nb_NC;                   /*!< Number of NC types                */
  size_t          _nb_E;                    /*!< Number of E types                 */
  size_t          _nb_TF;                   /*!< Number of TF types                */
  size_t          _nb_BS;                   /*!< Number of BS types                */
  size_t          _nb_P;                    /*!< Number of P types                 */
  size_t          _src_breakpoint1;         /*!< First source breakpoint           */
  size_t          _src_breakpoint2;         /*!< Second source breakpoint          */
  size_t          _tgt_breakpoint;          /*!< Target breakpoint                 */
  size_t          _size;                    /*!< Rearrangement size                */
  MutationVector* _mutation_vector;         /*!< Point mutation vector             */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/**
 * \brief    Get mutation type
 * \details  Mutation type is POINT_MUTATION, DUPLICATION, DELETION, TRANSLOCATION or INVERSION
 * \param    void
 * \return   \e mutation_type
 */
inline mutation_type MutationEvent::get_mutation_type( void ) const
{
  return _mutation_type;
}

/**
 * \brief    Get location of point mutation
 * \details  Location is position of a genetic unit in the genome
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_point_mutation_location( void ) const
{
  return _point_mutation_location;
}

/**
 * \brief    Get insertion point of hgt
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_hgt_insert( void ) const
{
  return _hgt_insert;
}

/**
 * \brief    Get the number of NC types mutated, duplicated or deleted
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_nb_NC( void ) const
{
  return _nb_NC;
}

/**
 * \brief    Get the number of E types mutated, duplicated or deleted
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_nb_E( void ) const
{
  return _nb_E;
}

/**
 * \brief    Get the number of TF types mutated, duplicated or deleted
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_nb_TF( void ) const
{
  return _nb_TF;
}

/**
 * \brief    Get the number of BS types mutated, duplicated or deleted
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_nb_BS( void ) const
{
  return _nb_BS;
}

/**
 * \brief    Get the number of P types mutated, duplicated or deleted
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_nb_P( void ) const
{
  return _nb_P;
}

/**
 * \brief    Get first source breakpoint
 * \details  Location of the first source breakpoint of rearrangement
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_src_breakpoint1( void ) const
{
  return _src_breakpoint1;
}

/**
 * \brief    Get second source breakpoint
 * \details  Location of the second source breakpoint of rearrangement
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_src_breakpoint2( void ) const
{
  return _src_breakpoint2;
}

/**
 * \brief    Get target breakpoint
 * \details  Location of the target breakpoint of rearrangement
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_tgt_breakpoint( void ) const
{
  return _tgt_breakpoint;
}

/**
 * \brief    Get rearrangement size
 * \details  Get the length of the rearrangement
 * \param    void
 * \return   \e size_t
 */
inline size_t MutationEvent::get_size( void ) const
{
  return _size;
}

/**
 * \brief    Get point mutation vector
 * \details  Return the point mutation random vector
 * \param    void
 * \return   \e MutationVector*
 */
inline MutationVector* MutationEvent::get_mutation_vector( void )
{
  return _mutation_vector;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__MutationEvent__) */
