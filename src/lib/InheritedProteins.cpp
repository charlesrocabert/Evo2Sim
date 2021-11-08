
/**
 * \file      InheritedProteins.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     InheritedProteins class definition
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

#include "InheritedProteins.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
InheritedProteins::InheritedProteins( Parameters* parameters )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  create_genetic_sequence();
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _Ei  = NULL;
  _TFi = NULL;
}

/**
 * \brief    Constructor from backup file
 * \details  Load InheritedProteins class from backup file
 * \param    Parameters* parameters
 * \param    gzFile backup_file
 * \return   \e void
 */
InheritedProteins::InheritedProteins( Parameters* parameters, gzFile backup_file )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  load_genetic_sequence(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_genetic_sequence->size > 0)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    for (size_t i = 0; i < _genetic_sequence->size; i++)
    {
      gzread( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzread( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzread( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzread( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzread( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzread( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    for (size_t i = 0; i < _nb_E; i++)
    {
      gzread( backup_file, &_Ei[i], sizeof(_Ei[i]) );
    }
  }
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzread( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const InheritedProteins& inherited_proteins
 * \return   \e void
 */
InheritedProteins::InheritedProteins( const InheritedProteins& inherited_proteins )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = inherited_proteins._parameters;
  
  /*------------------------------------------------------------------ genetic sequence */
  
  copy_genetic_sequence(inherited_proteins._genetic_sequence);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_genetic_sequence->size > 0)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    memcpy(_concentration_vector, inherited_proteins._concentration_vector, sizeof(double)*_genetic_sequence->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = inherited_proteins._nb_E;
  _nb_TF            = inherited_proteins._nb_TF;
  _nb_inner_enzymes = inherited_proteins._nb_inner_enzymes;
  _nb_inflow_pumps  = inherited_proteins._nb_inflow_pumps;
  _nb_outflow_pumps = inherited_proteins._nb_outflow_pumps;
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    memcpy(_Ei,  inherited_proteins._Ei,  sizeof(size_t)*_nb_E);
  }
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, inherited_proteins._TFi, sizeof(size_t)*_nb_TF);
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
InheritedProteins::~InheritedProteins( void )
{
  delete_genetic_sequence();
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  delete[] _Ei;
  _Ei  = NULL;
  delete[] _TFi;
  _TFi = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Initialize the concentration vector at zero (after mutation)
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::initialize_concentration_vector( void )
{
  delete[] _concentration_vector;
  _concentration_vector = new double[_genetic_sequence->size];
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    _concentration_vector[i] = 0.0;
  }
}

/**
 * \brief    Build the indexes list
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::build_index_list( void )
{
  /*---------------------------------------------------------*/
  /* 1) count genetic unit types                             */
  /*---------------------------------------------------------*/
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  genetic_unit* str = _genetic_sequence->x;
  for (size_t pos = 0; pos < _genetic_sequence->size; pos++)
  {
    if ( str[pos].type == ENZYME )
    {
      _nb_E++;
      if (str[pos].kcat > 0.0 && str[pos].s == str[pos].p)
      {
        _nb_inflow_pumps++;
      }
      else if (str[pos].kcat < 0.0 && str[pos].s == str[pos].p)
      {
        _nb_outflow_pumps++;
      }
      else if (str[pos].kcat != 0.0 && str[pos].s != str[pos].p)
      {
        _nb_inner_enzymes++;
      }
    }
    else if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _nb_TF++;
    }
  }
  
  /*---------------------------------------------------------*/
  /* 2) compute the indexes list for each genetic unit types */
  /*---------------------------------------------------------*/
  delete[] _Ei;
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
  }
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
  }
  size_t Ecount  = 0;
  size_t TFcount = 0;
  for (size_t pos = 0; pos < _genetic_sequence->size; pos++)
  {
    if (str[pos].type == ENZYME)
    {
      _Ei[Ecount] = pos;
      Ecount++;
    }
    else if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _TFi[TFcount] = pos;
      TFcount++;
    }
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ geneti sequence */
  
  save_genetic_sequence(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  if (_genetic_sequence->size > 0)
  {
    for (size_t i = 0; i < _genetic_sequence->size; i++)
    {
      gzwrite( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzwrite( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzwrite( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzwrite( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzwrite( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzwrite( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  if (_nb_E > 0)
  {
    for (size_t i = 0; i < _nb_E; i++)
    {
      gzwrite( backup_file, &_Ei[i], sizeof(_Ei[i]) );
    }
  }
  if (_nb_TF > 0)
  {
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzwrite( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
}

/**
 * \brief    Replace inherited proteins data
 * \details  --
 * \param    InheritedProteins* inherited_proteins
 * \return   \e void
 */
void InheritedProteins::replace_data( InheritedProteins* inherited_proteins )
{
  /*------------------------------------------------------------------ genetic sequence */
  
  delete_genetic_sequence();
  copy_genetic_sequence(inherited_proteins->get_genetic_sequence());
  
  /*------------------------------------------------------------------ concentration vector */
  
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  if (_genetic_sequence->size > 0)
  {
    _concentration_vector = new double[_genetic_sequence->size];
    memcpy(_concentration_vector, inherited_proteins->get_concentration_vector(), sizeof(double)*_genetic_sequence->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = inherited_proteins->get_nb_E();
  _nb_TF            = inherited_proteins->get_nb_TF();
  _nb_inner_enzymes = inherited_proteins->get_nb_inner_enzymes();
  _nb_inflow_pumps  = inherited_proteins->get_nb_inflow_pumps();
  _nb_outflow_pumps = inherited_proteins->get_nb_outflow_pumps();
  
  /*------------------------------------------------------------------ genetic unit positions */
  
  delete[] _Ei;
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    memcpy(_Ei,  inherited_proteins->get_Ei(),  sizeof(size_t)*_nb_E);
  }
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, inherited_proteins->get_TFi(), sizeof(size_t)*_nb_TF);
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Create a default genetic sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::create_genetic_sequence( void )
{
  _genetic_sequence              = new genetic_sequence;
  _genetic_sequence->x           = new genetic_unit[INHERITED_PROTEINS_BUFFER];
  _genetic_sequence->size        = 0;
  _genetic_sequence->buffer_size = INHERITED_PROTEINS_BUFFER;
}

/**
 * \brief    Copy a genetic sequence
 * \details  --
 * \param    const genetic_sequence* model
 * \return   \e void
 */
void InheritedProteins::copy_genetic_sequence( const genetic_sequence* model )
{
  _genetic_sequence              = new genetic_sequence;
  _genetic_sequence->x           = new genetic_unit[model->buffer_size];
  memcpy(_genetic_sequence->x, model->x, sizeof(genetic_unit)*model->buffer_size);
  _genetic_sequence->size        = model->size;
  _genetic_sequence->buffer_size = model->buffer_size;
}

/**
 * \brief    Delete the genetic sequence
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::delete_genetic_sequence( void )
{
  delete[] _genetic_sequence->x;
  _genetic_sequence->x = NULL;
  delete _genetic_sequence;
  _genetic_sequence = NULL;
}

/**
 * \brief    Load the genetic sequence from backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::load_genetic_sequence( gzFile backup_file )
{
  _genetic_sequence = new genetic_sequence;
  gzread( backup_file, &_genetic_sequence->size,        sizeof(_genetic_sequence->size) );
  gzread( backup_file, &_genetic_sequence->buffer_size, sizeof(_genetic_sequence->buffer_size) );
  _genetic_sequence->x = new genetic_unit[_genetic_sequence->buffer_size];
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    load_genetic_unit(backup_file, _genetic_sequence->x[i]);
  }
}

/**
 * \brief    Save the genetic sequence in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::save_genetic_sequence( gzFile backup_file )
{
  gzwrite( backup_file, &_genetic_sequence->size,        sizeof(_genetic_sequence->size) );
  gzwrite( backup_file, &_genetic_sequence->buffer_size, sizeof(_genetic_sequence->buffer_size) );
  for (size_t i = 0; i < _genetic_sequence->size; i++)
  {
    save_genetic_unit(backup_file, _genetic_sequence->x[i]);
  }
}

/**
 * \brief    Load a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit& unit
 * \return   \e void
 */
void InheritedProteins::load_genetic_unit( gzFile backup_file, genetic_unit& unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzread( backup_file, &unit.type,              sizeof(unit.type) );
  gzread( backup_file, &unit.identifier,        sizeof(unit.identifier) );
  gzread( backup_file, &unit.parent_identifier, sizeof(unit.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzread( backup_file, &unit.s,             sizeof(unit.s) );
  gzread( backup_file, &unit.p,             sizeof(unit.p) );
  gzread( backup_file, &unit.kcat,          sizeof(unit.kcat) );
  gzread( backup_file, &unit.kcat_km_ratio, sizeof(unit.kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzread( backup_file, &unit.BS_tag,         sizeof(unit.BS_tag) );
  gzread( backup_file, &unit.coE_tag,        sizeof(unit.coE_tag) );
  gzread( backup_file, &unit.free_activity,  sizeof(unit.free_activity) );
  gzread( backup_file, &unit.bound_activity, sizeof(unit.bound_activity) );
  gzread( backup_file, &unit.binding_window, sizeof(unit.binding_window));
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &unit.TF_tag, sizeof(unit.TF_tag));
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzread( backup_file, &unit.basal_expression_level, sizeof(unit.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzread( backup_file, &unit.functional, sizeof(unit.functional) );
}

/**
 * \brief    Save a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit& unit
 * \return   \e void
 */
void InheritedProteins::save_genetic_unit( gzFile backup_file, genetic_unit& unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzwrite( backup_file, &unit.type,              sizeof(unit.type) );
  gzwrite( backup_file, &unit.identifier,        sizeof(unit.identifier) );
  gzwrite( backup_file, &unit.parent_identifier, sizeof(unit.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzwrite( backup_file, &unit.s,             sizeof(unit.s) );
  gzwrite( backup_file, &unit.p,             sizeof(unit.p) );
  gzwrite( backup_file, &unit.kcat,          sizeof(unit.kcat) );
  gzwrite( backup_file, &unit.kcat_km_ratio, sizeof(unit.kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzwrite( backup_file, &unit.BS_tag,         sizeof(unit.BS_tag) );
  gzwrite( backup_file, &unit.coE_tag,        sizeof(unit.coE_tag) );
  gzwrite( backup_file, &unit.free_activity,  sizeof(unit.free_activity) );
  gzwrite( backup_file, &unit.bound_activity, sizeof(unit.bound_activity) );
  gzwrite( backup_file, &unit.binding_window, sizeof(unit.binding_window));
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &unit.TF_tag, sizeof(unit.TF_tag));
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &unit.basal_expression_level, sizeof(unit.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &unit.functional, sizeof(unit.functional) );
}

/**
 * \brief    Increase buffer size relatively to the new inherited proteins size
 * \details  --
 * \param    size_t new_size
 * \return   \e void
 */
void InheritedProteins::increase_buffer_size( size_t new_size )
{
  assert(new_size <= _genetic_sequence->size*2);
  _genetic_sequence->buffer_size = (new_size/INHERITED_PROTEINS_BUFFER+1)*INHERITED_PROTEINS_BUFFER;
  genetic_unit* new_x = new genetic_unit[_genetic_sequence->buffer_size];
  memcpy(new_x, _genetic_sequence->x, sizeof(genetic_unit)*_genetic_sequence->size);
  delete[] _genetic_sequence->x;
  _genetic_sequence->x = new_x;
}

/**
 * \brief    Decrease buffer size relatively to the inherited proteins size
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::decrease_buffer_size( void )
{
  if (_genetic_sequence->buffer_size > INHERITED_PROTEINS_BUFFER)
  {
    _genetic_sequence->buffer_size = (_genetic_sequence->size/INHERITED_PROTEINS_BUFFER+1)*INHERITED_PROTEINS_BUFFER;
    assert(_genetic_sequence->buffer_size >= _genetic_sequence->size);
    genetic_unit* new_x = new genetic_unit[_genetic_sequence->buffer_size];
    memcpy(new_x, _genetic_sequence->x, sizeof(genetic_unit)*_genetic_sequence->size);
    delete[] _genetic_sequence->x;
    _genetic_sequence->x = new_x;
  }
}
