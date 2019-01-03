
/**
 * \file      Genome.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Genome class declaration
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

#ifndef __Evo2Sim__Genome__
#define __Evo2Sim__Genome__

#include <iostream>
#include <zlib.h>
#include <stdlib.h>
#include <assert.h>
#include <cstring>
#include <cmath>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Parameters.h"
#include "MutationVector.h"
#include "MutationEvent.h"
#include "ReplicationReport.h"
#include "Prng.h"


class Genome
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Genome( void ) = delete;
  Genome( const Genome& genome ) = delete;
  Genome( Parameters* parameters, Prng* prng, ReplicationReport* replication_report );
  Genome( Parameters* parameters, Prng* prng, ReplicationReport* replication_report, gzFile backup_file );
  Genome( const Genome& genome, Prng* prng, ReplicationReport* replication_report );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Genome( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  inline genetic_sequence* get_genetic_sequence( void );
  inline genetic_unit*     get_genetic_unit( size_t pos );
  inline size_t            get_size( void ) const;
  inline size_t            get_buffer_size( void ) const;
  inline size_t            get_coding_size( void ) const;
  inline size_t            get_non_coding_size( void ) const;
  inline double*           get_concentration_vector( void );
  inline size_t            get_nb_NC( void ) const;
  inline size_t            get_nb_E( void ) const;
  inline size_t            get_nb_TF( void ) const;
  inline size_t            get_nb_BS( void ) const;
  inline size_t            get_nb_P( void ) const;
  inline size_t            get_nb_inner_enzymes( void ) const;
  inline size_t            get_nb_inflow_pumps( void ) const;
  inline size_t            get_nb_outflow_pumps( void ) const;
  inline size_t*           get_TFi( void );
  inline size_t*           get_Pi( void );
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Genome& operator=(const Genome&) = delete;
  inline void add_genetic_unit( genetic_unit* unit );
  inline void clear( void );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void initialize_concentration_vector( void );
  void mutate( const double* mutation_rates );
  void shuffle( void );
  void save( gzFile backup_file );
  void replace_data( Genome* genome );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void do_point_mutations( const double* mutation_rates );
  void do_rearrangements( const double* mutation_rates );
  
  void do_duplication( const double* mutation_rates );
  void do_deletion( const double* mutation_rates );
  void do_translocation( const double* mutation_rates );
  void do_inversion( const double* mutation_rates );
  
  void do_hgt( const double* mutation_rates );
  void do_genetic_unit_mutation( const double* mutation_rates, size_t pos );
  void do_genetic_unit_mutation_at_breakpoints( const double* mutation_rates, size_t pos );
  void do_permutation( const double* mutation_rates, size_t pos1, size_t pos2 );
  
  void create_genetic_sequence( void );
  void copy_genetic_sequence( const genetic_sequence* model );
  void delete_genetic_sequence( void );
  void load_genetic_sequence( gzFile backup_file );
  void save_genetic_sequence( gzFile backup_file );
  void load_genetic_unit( gzFile backup_file, genetic_unit& unit );
  void save_genetic_unit( gzFile backup_file, genetic_unit& unit );
  void draw_random_genetic_unit( genetic_unit& unit );
  
  void increase_buffer_size( size_t new_size );
  void decrease_buffer_size( void );
  void check_genome_size( void );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  Prng* _prng; /*!< Pseudorandom numbers generator */
  
  /*------------------------------------------------------------------ simulation parameters */
  
  Parameters* _parameters; /*!< Simulation parameters */
  
  /*------------------------------------------------------------------ genetic sequence */
  
  genetic_sequence* _genetic_sequence; /*!< Genetic sequence */
  
  /*------------------------------------------------------------------ concentration vector */
  
  double* _concentration_vector; /*!< Concentration vector */
  
  /*------------------------------------------------------------------ statistical data */
  
  size_t _nb_NC;            /*!< Number of non coding types (NC)           */
  size_t _nb_E;             /*!< Number of enzyme types (E)                */
  size_t _nb_TF;            /*!< Number of transcription factor types (TF) */
  size_t _nb_BS;            /*!< Number of binding site types (BS)         */
  size_t _nb_P;             /*!< Number of promoter types (P)              */
  size_t _nb_inner_enzymes; /*!< Number of inner metabolism enzymes        */
  size_t _nb_inflow_pumps;  /*!< Number of inflowing pumps                 */
  size_t _nb_outflow_pumps; /*!< Number of outflowing pumps                */
  
  /*------------------------------------------------------------------ replication report */
  
  ReplicationReport* _replication_report; /*!< Replication report */
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  size_t* _TFi; /*!< Vector of TF positions */
  size_t* _Pi;  /*!< Vector of P positions  */
  
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
inline genetic_sequence* Genome::get_genetic_sequence( void )
{
  return _genetic_sequence;
}

/*
 * \brief    Get genetic unit at position 'pos'
 * \details  --
 * \param    size_t pos
 * \return   \e genetic_unit*
 */
inline genetic_unit* Genome::get_genetic_unit( size_t pos )
{
  assert(pos < _genetic_sequence->size);
  return &_genetic_sequence->x[pos];
}

/**
 * \brief    Get genome size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_size( void ) const
{
  return _genetic_sequence->size;
}

/**
 * \brief    Get buffer size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_buffer_size( void ) const
{
  return _genetic_sequence->buffer_size;
}

/**
 * \brief    Get coding size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_coding_size( void ) const
{
  return _genetic_sequence->size - _nb_NC;
}

/**
 * \brief    Get non coding size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_non_coding_size( void ) const
{
  return _nb_NC;
}

/*------------------------------------------------------------------ concentration vector */

/**
 * \brief    Get the concentration vector
 * \details  --
 * \param    void
 * \return   \e double*
 */
inline double* Genome::get_concentration_vector( void )
{
  return _concentration_vector;
}

/*------------------------------------------------------------------ statistical data */

/**
 * \brief    Get number of non coding types (NC)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_NC( void ) const
{
  return _nb_NC;
}

/**
 * \brief    Get number of enzyme types (E)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_E( void ) const
{
  return _nb_E;
}

/**
 * \brief    Get number of transcription factor types (TF)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_TF( void ) const
{
  return _nb_TF;
}

/**
 * \brief    Get number of binding site types (BS)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_BS( void ) const
{
  return _nb_BS;
}

/**
 * \brief    Get number of promoter types (P)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_P( void ) const
{
  return _nb_P;
}

/**
 * \brief    Get number of inner metabolism enzymes
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_inner_enzymes( void ) const
{
  return _nb_inner_enzymes;
}

/**
 * \brief    Get number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_inflow_pumps( void ) const
{
  return _nb_inflow_pumps;
}

/**
 * \brief    Get number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Genome::get_nb_outflow_pumps( void ) const
{
  return _nb_outflow_pumps;
}

/*------------------------------------------------------------------ genetic unit positions */

/**
 * \brief    Get the list of TF type positions
 * \details  --
 * \param    void
 * \return   \e size_t*
 */
inline size_t* Genome::get_TFi( void )
{
  return _TFi;
}

/**
 * \brief    Get the list of P type positions
 * \details  --
 * \param    void
 * \return   \e size_t*
 */
inline size_t* Genome::get_Pi( void )
{
  return _Pi;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/**
 * \brief    Add a genetic unit at the end of the sequence
 * \details  --
 * \param    genetic_unit* unit
 * \return   \e void
 */
inline void Genome::add_genetic_unit( genetic_unit* unit )
{
  if (_genetic_sequence->size+1 > _genetic_sequence->buffer_size)
  {
    increase_buffer_size(_genetic_sequence->size+1);
  }
  _genetic_sequence->x[_genetic_sequence->size] = *(unit);
  _genetic_sequence->size++;
}

/**
 * \brief    Clear genome sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Genome::clear( void )
{
  delete[] _genetic_sequence->x;
  _genetic_sequence->x           = NULL;
  _genetic_sequence->x           = new genetic_unit[GENOME_BUFFER];
  _genetic_sequence->size        = 0;
  _genetic_sequence->buffer_size = GENOME_BUFFER;
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  _nb_NC            = 0;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_BS            = 0;
  _nb_P             = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  delete[] _TFi;
  _TFi = NULL;
  delete[] _Pi;
  _Pi  = NULL;
}


#endif /* defined(__Evo2Sim__Genome__) */
