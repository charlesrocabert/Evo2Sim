
/**
 * \file      MutationEvent.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MutationEvent class definition
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

#include "MutationEvent.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor for point mutation event
 * \details  This constructor must be used to create a point mutation event
 * \param    mutation_type type
 * \param    size_t point_mutation_location
 * \param    MutationVector* vector
 * \return   \e void
 */
MutationEvent::MutationEvent( mutation_type type, size_t point_mutation_location, MutationVector* vector )
{
  _mutation_type           = type;
  _point_mutation_location = point_mutation_location;
  _hgt_insert              = 0;
  if (vector->get_dX()->type == NON_CODING ||
      vector->get_dX()->type == E_TO_NC_TRANSITION ||
      vector->get_dX()->type == TF_TO_NC_TRANSITION ||
      vector->get_dX()->type == BS_TO_NC_TRANSITION ||
      vector->get_dX()->type == P_TO_NC_TRANSITION)
  {
    _nb_NC = 1;
    _nb_E  = 0;
    _nb_TF = 0;
    _nb_BS = 0;
    _nb_P  = 0;
  }
  else if (vector->get_dX()->type == ENZYME ||
           vector->get_dX()->type == NC_TO_E_TRANSITION ||
           vector->get_dX()->type == TF_TO_E_TRANSITION ||
           vector->get_dX()->type == BS_TO_E_TRANSITION ||
           vector->get_dX()->type == P_TO_E_TRANSITION)
  {
    _nb_NC = 0;
    _nb_E  = 1;
    _nb_TF = 0;
    _nb_BS = 0;
    _nb_P  = 0;
  }
  else if (vector->get_dX()->type == TRANSCRIPTION_FACTOR ||
           vector->get_dX()->type == NC_TO_TF_TRANSITION ||
           vector->get_dX()->type == E_TO_TF_TRANSITION ||
           vector->get_dX()->type == BS_TO_TF_TRANSITION ||
           vector->get_dX()->type == P_TO_TF_TRANSITION)
  {
    _nb_NC = 0;
    _nb_E  = 0;
    _nb_TF = 1;
    _nb_BS = 0;
    _nb_P  = 0;
  }
  else if (vector->get_dX()->type == BINDING_SITE ||
           vector->get_dX()->type == NC_TO_BS_TRANSITION ||
           vector->get_dX()->type == E_TO_BS_TRANSITION ||
           vector->get_dX()->type == TF_TO_BS_TRANSITION ||
           vector->get_dX()->type == P_TO_BS_TRANSITION)
  {
    _nb_NC = 0;
    _nb_E  = 0;
    _nb_TF = 0;
    _nb_BS = 1;
    _nb_P  = 0;
  }
  else if (vector->get_dX()->type == PROMOTER ||
           vector->get_dX()->type == NC_TO_P_TRANSITION ||
           vector->get_dX()->type == E_TO_P_TRANSITION ||
           vector->get_dX()->type == TF_TO_P_TRANSITION ||
           vector->get_dX()->type == BS_TO_P_TRANSITION)
  {
    _nb_NC = 0;
    _nb_E  = 0;
    _nb_TF = 0;
    _nb_BS = 0;
    _nb_P  = 1;
  }
  _src_breakpoint1 = 0;
  _src_breakpoint2 = 0;
  _tgt_breakpoint  = 0;
  _size            = 1;
  _mutation_vector = vector;
}

/**
 * \brief    Constructor for HGT
 * \details  This constructor must be used to create a HGT event
 * \param    mutation_type type
 * \param    size_t hgt_insert
 * \param    size_t size
 * \param    size_t nb_NC
 * \param    size_t nb_E
 * \param    size_t nb_TF
 * \param    size_t nb_BS
 * \param    size_t nb_P
 * \return   \e void
 */
MutationEvent::MutationEvent( mutation_type type, size_t hgt_insert, size_t size, size_t nb_NC, size_t nb_E, size_t nb_TF, size_t nb_BS, size_t nb_P )
{
  assert(size == nb_NC+nb_E+nb_TF+nb_BS+nb_P);
  _mutation_type           = type;
  _point_mutation_location = 0;
  _hgt_insert              = hgt_insert;
  _nb_NC                   = nb_NC;
  _nb_E                    = nb_E;
  _nb_TF                   = nb_TF;
  _nb_BS                   = nb_BS;
  _nb_P                    = nb_P;
  _src_breakpoint1         = 0;
  _src_breakpoint2         = 0;
  _tgt_breakpoint          = 0;
  _size                    = size;
  _mutation_vector         = NULL;
}

/**
 * \brief    Constructor for large rearrangement
 * \details  This constructor must be used to create a rearrangement event
 * \param    mutation_type type
 * \param    size_t src1
 * \param    size_t src2
 * \param    size_t tgt
 * \param    size_t size
 * \param    size_t nb_NC
 * \param    size_t nb_E
 * \param    size_t nb_TF
 * \param    size_t nb_BS
 * \param    size_t nb_P
 * \return   \e void
 */
MutationEvent::MutationEvent( mutation_type type, size_t src1, size_t src2, size_t tgt, size_t size, size_t nb_NC, size_t nb_E, size_t nb_TF, size_t nb_BS, size_t nb_P )
{
  if (type == DUPLICATION || type == DELETION)
  {
    assert(size == nb_NC+nb_E+nb_TF+nb_BS+nb_P);
  }
  _mutation_type           = type;
  _point_mutation_location = 0;
  _hgt_insert              = 0;
  _nb_NC                   = nb_NC;
  _nb_E                    = nb_E;
  _nb_TF                   = nb_TF;
  _nb_BS                   = nb_BS;
  _nb_P                    = nb_P;
  _src_breakpoint1         = src1;
  _src_breakpoint2         = src2;
  _tgt_breakpoint          = tgt;
  _size                    = size;
  _mutation_vector         = NULL;
}

/**
 * \brief    Constructor from backup file
 * \details  Load MutationEvent class from backup file
 * \param    gzFile backup_file
 * \return   \e void
 */
MutationEvent::MutationEvent( gzFile backup_file )
{
  gzread( backup_file, &_mutation_type,           sizeof(_mutation_type) );
  gzread( backup_file, &_point_mutation_location, sizeof(_point_mutation_location) );
  gzread( backup_file, &_hgt_insert,              sizeof(_hgt_insert) );
  gzread( backup_file, &_nb_NC,                   sizeof(_nb_NC) );
  gzread( backup_file, &_nb_E,                    sizeof(_nb_E) );
  gzread( backup_file, &_nb_TF,                   sizeof(_nb_TF) );
  gzread( backup_file, &_nb_BS,                   sizeof(_nb_BS) );
  gzread( backup_file, &_nb_P,                    sizeof(_nb_P) );
  gzread( backup_file, &_src_breakpoint1,         sizeof(_src_breakpoint1) );
  gzread( backup_file, &_src_breakpoint2,         sizeof(_src_breakpoint2) );
  gzread( backup_file, &_tgt_breakpoint,          sizeof(_tgt_breakpoint) );
  gzread( backup_file, &_size,                    sizeof(_size) );
  _mutation_vector = NULL;
  if (_mutation_type == POINT_MUTATION)
  {
    _mutation_vector = new MutationVector(backup_file);
  }
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const MutationEvent& event
 * \return   \e void
 */
MutationEvent::MutationEvent( const MutationEvent& event )
{
  _mutation_type           = event._mutation_type;
  _point_mutation_location = event._point_mutation_location;
  _hgt_insert              = event._hgt_insert;
  _nb_NC                   = event._nb_NC;
  _nb_E                    = event._nb_E;
  _nb_TF                   = event._nb_TF;
  _nb_BS                   = event._nb_BS;
  _nb_P                    = event._nb_P;
  _src_breakpoint1         = event._src_breakpoint1;
  _src_breakpoint2         = event._src_breakpoint2;
  _tgt_breakpoint          = event._tgt_breakpoint;
  _size                    = event._size;
  _mutation_vector = NULL;
  if (_mutation_type == POINT_MUTATION)
  {
    _mutation_vector = new MutationVector(*event._mutation_vector);
  }
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
MutationEvent::~MutationEvent( void )
{
  delete _mutation_vector;
  _mutation_vector = NULL;
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
void MutationEvent::save( gzFile backup_file )
{
  gzwrite( backup_file, &_mutation_type,           sizeof(_mutation_type) );
  gzwrite( backup_file, &_point_mutation_location, sizeof(_point_mutation_location) );
  gzwrite( backup_file, &_hgt_insert,              sizeof(_hgt_insert) );
  gzwrite( backup_file, &_nb_NC,                   sizeof(_nb_NC) );
  gzwrite( backup_file, &_nb_E,                    sizeof(_nb_E) );
  gzwrite( backup_file, &_nb_TF,                   sizeof(_nb_TF) );
  gzwrite( backup_file, &_nb_BS,                   sizeof(_nb_BS) );
  gzwrite( backup_file, &_nb_P,                    sizeof(_nb_P) );
  gzwrite( backup_file, &_src_breakpoint1,         sizeof(_src_breakpoint1) );
  gzwrite( backup_file, &_src_breakpoint2,         sizeof(_src_breakpoint2) );
  gzwrite( backup_file, &_tgt_breakpoint,          sizeof(_tgt_breakpoint) );
  gzwrite( backup_file, &_size,                    sizeof(_size) );
  if (_mutation_type == POINT_MUTATION)
  {
    _mutation_vector->save(backup_file);
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
