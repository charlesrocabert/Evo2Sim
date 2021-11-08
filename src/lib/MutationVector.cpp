
/**
 * \file      MutationVector.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2021 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MutationVector class definition
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

#include "MutationVector.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
MutationVector::MutationVector( void )
{
  _dX = new genetic_unit;
  
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = NON_CODING;
  _dX->identifier        = 0;
  _dX->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s             = 0;
  _dX->p             = 0;
  _dX->kcat          = 0.0;
  _dX->kcat_km_ratio = 0.0;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = 0;
  _dX->coE_tag        = 0;
  _dX->free_activity  = false;
  _dX->bound_activity = false;
  _dX->binding_window = 0;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = false;
}

/**
 * \brief    Constructor from backup file
 * \details  Load MutationVector class from backup file
 * \param    gzFile backup_file
 * \return   \e void
 */
MutationVector::MutationVector( gzFile backup_file )
{
  _dX = new genetic_unit;
  load_genetic_unit(backup_file, _dX);
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const MutationVector& vector
 * \return   \e void
 */
MutationVector::MutationVector( const MutationVector& vector )
{
  _dX = new genetic_unit;
  
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = vector._dX->type;
  _dX->identifier        = vector._dX->identifier;
  _dX->parent_identifier = vector._dX->parent_identifier;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s             = vector._dX->s;
  _dX->p             = vector._dX->p;
  _dX->kcat          = vector._dX->kcat;
  _dX->kcat_km_ratio = vector._dX->kcat_km_ratio;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = vector._dX->BS_tag;
  _dX->coE_tag        = vector._dX->coE_tag;
  _dX->free_activity  = vector._dX->free_activity;
  _dX->bound_activity = vector._dX->bound_activity;
  _dX->binding_window = vector._dX->binding_window;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = vector._dX->TF_tag;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = vector._dX->basal_expression_level;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = vector._dX->functional;
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
MutationVector::~MutationVector( void )
{
  delete _dX;
  _dX = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/*
 * \brief    Draw vector
 * \details  --
 * \param    Prng* prng
 * \param    const double* mutation_rates
 * \return   \e bool
 */
bool MutationVector::draw( Prng* prng, const double* mutation_rates )
{
  bool mutate = false;
  
  /*-----------------------------------------------------*/
  /* 1) Mutate genetic unit type                         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[TRANSITION_RATE])
  {
    /* Draw uniformly between the five types of genetic units */
    double probas[5]           = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
    genetic_unit_type new_type = (genetic_unit_type)prng->roulette_wheel(probas, 1.0, 5);
    if (new_type != _dX->type)
    {
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* A) The genetic unit was NON_CODING type           */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (_dX->type == NON_CODING && new_type == ENZYME)
      {
        _dX->type = NC_TO_E_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = NC_TO_TF_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == BINDING_SITE)
      {
        _dX->type = NC_TO_BS_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == PROMOTER)
      {
        _dX->type = NC_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* B) The genetic unit was ENZYME type               */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == ENZYME && new_type == NON_CODING)
      {
        _dX->type = E_TO_NC_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = E_TO_TF_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == BINDING_SITE)
      {
        _dX->type = E_TO_BS_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == PROMOTER)
      {
        _dX->type = E_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* C) The genetic unit was TRANSCRIPTION_FACTOR type */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == NON_CODING)
      {
        _dX->type = TF_TO_NC_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == ENZYME)
      {
        _dX->type = TF_TO_E_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == BINDING_SITE)
      {
        _dX->type = TF_TO_BS_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == PROMOTER)
      {
        _dX->type = TF_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* D) The genetic unit was BINDING_SITE type         */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == BINDING_SITE && new_type == NON_CODING)
      {
        _dX->type = BS_TO_NC_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == ENZYME)
      {
        _dX->type = BS_TO_E_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = BS_TO_TF_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == PROMOTER)
      {
        _dX->type = BS_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* E) The genetic unit was PROMOTER type             */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == PROMOTER && new_type == NON_CODING)
      {
        _dX->type = P_TO_NC_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == ENZYME)
      {
        _dX->type = P_TO_E_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = P_TO_TF_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == BINDING_SITE)
      {
        _dX->type = P_TO_BS_TRANSITION;
      }
      mutate = true;
    }
  }
  
  /*-----------------------------------------------------*/
  /* 2) Mutate enzyme type (E) attributes                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->s = prng->uniform(-mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE], mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->p = prng->uniform(-mutation_rates[PRODUCT_TAG_MUTATION_SIZE], mutation_rates[PRODUCT_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->kcat = prng->gaussian(0.0, mutation_rates[KCAT_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->kcat_km_ratio = prng->gaussian(0.0, mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 3) Mutate transcription factor type (TF) attributes */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->BS_tag = prng->uniform(-mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE], mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->coE_tag = prng->uniform(-mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE], mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->free_activity = !_dX->free_activity;
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->bound_activity = !_dX->bound_activity;
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 4) Mutate binding site type (BS) attributes         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->TF_tag = prng->uniform(-mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 5) Mutate promoter type (P) attributes              */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->basal_expression_level = prng->gaussian(0.0, mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]);
    mutate = true;
  }
  return mutate;
}

/*
 * \brief    Draw vector at breakpoints
 * \details  --
 * \param    Prng* prng
 * \param    const double* mutation_rates
 * \return   \e bool
 */
bool MutationVector::breakpoint_draw( Prng* prng, const double* mutation_rates )
{
  bool mutate = false;
  
  /*-----------------------------------------------------*/
  /* 1) Mutate genetic unit type                         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    /* Draw uniformly between the five types of genetic units */
    double probas[5]           = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
    genetic_unit_type new_type = (genetic_unit_type)prng->roulette_wheel(probas, 1.0, 5);
    if (new_type != _dX->type)
    {
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* A) The genetic unit was NON_CODING type           */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (_dX->type == NON_CODING && new_type == ENZYME)
      {
        _dX->type = NC_TO_E_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = NC_TO_TF_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == BINDING_SITE)
      {
        _dX->type = NC_TO_BS_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == PROMOTER)
      {
        _dX->type = NC_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* B) The genetic unit was ENZYME type               */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == ENZYME && new_type == NON_CODING)
      {
        _dX->type = E_TO_NC_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = E_TO_TF_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == BINDING_SITE)
      {
        _dX->type = E_TO_BS_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == PROMOTER)
      {
        _dX->type = E_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* C) The genetic unit was TRANSCRIPTION_FACTOR type */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == NON_CODING)
      {
        _dX->type = TF_TO_NC_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == ENZYME)
      {
        _dX->type = TF_TO_E_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == BINDING_SITE)
      {
        _dX->type = TF_TO_BS_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == PROMOTER)
      {
        _dX->type = TF_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* D) The genetic unit was BINDING_SITE type         */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == BINDING_SITE && new_type == NON_CODING)
      {
        _dX->type = BS_TO_NC_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == ENZYME)
      {
        _dX->type = BS_TO_E_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = BS_TO_TF_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == PROMOTER)
      {
        _dX->type = BS_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* E) The genetic unit was PROMOTER type             */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == PROMOTER && new_type == NON_CODING)
      {
        _dX->type = P_TO_NC_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == ENZYME)
      {
        _dX->type = P_TO_E_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = P_TO_TF_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == BINDING_SITE)
      {
        _dX->type = P_TO_BS_TRANSITION;
      }
      mutate = true;
    }
  }
  
  /*-----------------------------------------------------*/
  /* 2) Mutate enzyme type (E) attributes                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->s = prng->uniform(-mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE], mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->p = prng->uniform(-mutation_rates[PRODUCT_TAG_MUTATION_SIZE], mutation_rates[PRODUCT_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->kcat = prng->gaussian(0.0, mutation_rates[KCAT_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->kcat_km_ratio = prng->gaussian(0.0, mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 3) Mutate transcription factor type (TF) attributes */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->BS_tag = prng->uniform(-mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE], mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->coE_tag = prng->uniform(-mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE], mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->free_activity = !_dX->free_activity;
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->bound_activity = !_dX->bound_activity;
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 4) Mutate binding site type (BS) attributes         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->TF_tag = prng->uniform(-mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 5) Mutate promoter type (P) attributes              */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[BREAKPOINT_RATE])
  {
    _dX->basal_expression_level = prng->gaussian(0.0, mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]);
    mutate = true;
  }
  return mutate;
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void MutationVector::save( gzFile backup_file )
{
  save_genetic_unit(backup_file, _dX);
}

/**
 * \brief    Clear random vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void MutationVector::clear( void )
{
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = NON_CODING;
  _dX->identifier        = 0;
  _dX->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s             = 0;
  _dX->p             = 0;
  _dX->kcat          = 0.0;
  _dX->kcat_km_ratio = 0.0;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = 0;
  _dX->coE_tag        = 0;
  _dX->free_activity  = false;
  _dX->bound_activity = false;
  _dX->binding_window = 0;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = false;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Load a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit* unit
 * \return   \e void
 */
void MutationVector::load_genetic_unit( gzFile backup_file, genetic_unit* unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzread( backup_file, &unit->type,              sizeof(unit->type) );
  gzread( backup_file, &unit->identifier,        sizeof(unit->identifier) );
  gzread( backup_file, &unit->parent_identifier, sizeof(unit->parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzread( backup_file, &unit->s,             sizeof(unit->s) );
  gzread( backup_file, &unit->p,             sizeof(unit->p) );
  gzread( backup_file, &unit->kcat,          sizeof(unit->kcat) );
  gzread( backup_file, &unit->kcat_km_ratio, sizeof(unit->kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzread( backup_file, &unit->BS_tag,         sizeof(unit->BS_tag) );
  gzread( backup_file, &unit->coE_tag,        sizeof(unit->coE_tag) );
  gzread( backup_file, &unit->free_activity,  sizeof(unit->free_activity) );
  gzread( backup_file, &unit->bound_activity, sizeof(unit->bound_activity) );
  gzread( backup_file, &unit->binding_window, sizeof(unit->binding_window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &unit->TF_tag, sizeof(unit->TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzread( backup_file, &unit->basal_expression_level, sizeof(unit->basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzread( backup_file, &unit->functional, sizeof(unit->functional) );
}

/**
 * \brief    Save a genetic unit in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    genetic_unit* unit
 * \return   \e void
 */
void MutationVector::save_genetic_unit( gzFile backup_file, genetic_unit* unit )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzwrite( backup_file, &unit->type,              sizeof(unit->type) );
  gzwrite( backup_file, &unit->identifier,        sizeof(unit->identifier) );
  gzwrite( backup_file, &unit->parent_identifier, sizeof(unit->parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzwrite( backup_file, &unit->s,             sizeof(unit->s) );
  gzwrite( backup_file, &unit->p,             sizeof(unit->p) );
  gzwrite( backup_file, &unit->kcat,          sizeof(unit->kcat) );
  gzwrite( backup_file, &unit->kcat_km_ratio, sizeof(unit->kcat_km_ratio) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzwrite( backup_file, &unit->BS_tag,         sizeof(unit->BS_tag) );
  gzwrite( backup_file, &unit->coE_tag,        sizeof(unit->coE_tag) );
  gzwrite( backup_file, &unit->free_activity,  sizeof(unit->free_activity) );
  gzwrite( backup_file, &unit->bound_activity, sizeof(unit->bound_activity) );
  gzwrite( backup_file, &unit->binding_window, sizeof(unit->binding_window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &unit->TF_tag, sizeof(unit->TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &unit->basal_expression_level, sizeof(unit->basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &unit->functional, sizeof(unit->functional) );
}
