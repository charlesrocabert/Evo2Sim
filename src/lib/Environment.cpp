
/**
 * \file      Environment.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Environment class definition
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

#include "Environment.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
Environment::Environment( Parameters* parameters )
{
  /*------------------------------------------------------------------ parameters */
  
  _parameters = parameters;
  _prng       = parameters->get_environment_prng();
  _width      = parameters->get_width();
  _height     = parameters->get_height();
  _properties = parameters->get_environment_properties();
  
  /*------------------------------------------------------------------ environment structure */
  
  _grid       = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList();
  }
  _X                 = new SpeciesList();
  _size              = 0;
  _total_amount      = 0.0;
  _min_amount        = 0.0;
  _max_amount        = 0.0;
  _inflowing_amount  = 0.0;
  _outflowing_amount = 0.0;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Environment class from backup file
 * \param    Parameters* parameters
 * \param    gzFile backup_file
 * \return   \e void
 */
Environment::Environment( Parameters* parameters, gzFile backup_file )
{
  /*------------------------------------------------------------------ parameters */
  
  _parameters = parameters;
  _prng       = parameters->get_environment_prng();
  _width      = parameters->get_width();
  _height     = parameters->get_height();
  _properties = parameters->get_environment_properties();
  
  /*------------------------------------------------------------------ environment structure */
  
  _grid = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList(backup_file);
  }
  _X = new SpeciesList(backup_file);
  gzread( backup_file, &_size,              sizeof(_size) );
  gzread( backup_file, &_total_amount,      sizeof(_total_amount) );
  gzread( backup_file, &_min_amount,        sizeof(_min_amount) );
  gzread( backup_file, &_max_amount,        sizeof(_max_amount) );
  gzread( backup_file, &_inflowing_amount,  sizeof(_inflowing_amount) );
  gzread( backup_file, &_outflowing_amount, sizeof(_outflowing_amount) );
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Environment& environment 
 * \return   \e void
 */
Environment::Environment( Environment const &environment )
{
  /*------------------------------------------------------------------ parameters */
  
  _parameters = environment._parameters;
  _prng       = environment._prng;
  _width      = environment._width;
  _height     = environment._height;
  _properties = environment._properties;
  
  /*------------------------------------------------------------------ environment structure */
  
  _grid       = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList(*environment._grid[i]);
  }
  _X                 = new SpeciesList(*environment._X);
  _size              = environment._size;
  _total_amount      = environment._total_amount;
  _min_amount        = environment._min_amount;
  _max_amount        = environment._max_amount;
  _inflowing_amount  = environment._inflowing_amount;
  _outflowing_amount = environment._outflowing_amount;
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
Environment::~Environment( void )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    delete _grid[i];
    _grid[i] = NULL;
  }
  delete[] _grid;
  _grid = NULL;
  delete _X;
  _X = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Actualize the state of the environment
 * \details  Compute total amount and X vector
 * \param    void
 * \return   \e void
 */
void Environment::actualize_environment_state( void )
{
  _size         = 0;
  _total_amount = 0.0;
  _min_amount   = 1e+06;
  _max_amount   = 0.0;
  _X->reset(false);
  /*------------------------------------------------*/
  /* 1) Compute the current amount and get max size */
  /*------------------------------------------------*/
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    if (_size < _grid[pos]->get_size())
    {
      _size = _grid[pos]->get_size();
    }
    for (int met = 0; met < (int)_grid[pos]->get_size(); met++)
    {
      _X->add(met+1, _grid[pos]->get(met+1));
    }
    _total_amount  += _grid[pos]->get_amount();
    if (_min_amount > _grid[pos]->get_amount())
    {
      _min_amount = _grid[pos]->get_amount();
    }
    if (_max_amount < _grid[pos]->get_amount())
    {
      _max_amount = _grid[pos]->get_amount();
    }
  }
  
  /*------------------------------------------------*/
  /* 2) Actualize the size of all the species lists */
  /*------------------------------------------------*/
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    if (_grid[pos]->get_size() < _size)
    {
      _grid[pos]->increase_size(_size);
    }
  }
}

/**
 * \brief    Compute diffusion and degradation
 * \details  Compute diffusion and degradation of metabolites in environment
 * \param    void
 * \return   \e void
 */
void Environment::compute_diffusion_and_degradation( void )
{
  if (_properties->diffusion_coefficient <= 0.1)
  {
    compute_classic_diffusion_and_degradation();
  }
  else
  {
    compute_perfect_diffusion_and_degradation();
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Environment::save( gzFile backup_file )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i]->save(backup_file);
  }
  _X->save(backup_file);
  gzwrite( backup_file, &_size,              sizeof(_size) );
  gzwrite( backup_file, &_total_amount,      sizeof(_total_amount) );
  gzwrite( backup_file, &_min_amount,        sizeof(_min_amount) );
  gzwrite( backup_file, &_max_amount,        sizeof(_max_amount) );
  gzwrite( backup_file, &_inflowing_amount,  sizeof(_inflowing_amount) );
  gzwrite( backup_file, &_outflowing_amount, sizeof(_outflowing_amount) );
}

/**
 * \brief    Write the 1D metabolic state vector of the environment
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Environment::write_metabolic_state_vector( std::string filename )
{
  if (_X->get_size() > 0)
  {
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    for (int i = 0; i < (int)_X->get_size(); i++)
    {
      file << i+1 << " " << _X->get(i+1) << "\n";
    }
    file.close();
  }
}

/**
 * \brief    Write the 1D metabolic state vector of one local environment
 * \details  --
 * \param    std::string filename
 * \param    size_t pos
 * \return   \e void
 */
void Environment::write_local_metabolic_state_vector( std::string filename, size_t pos )
{
  assert(pos < _width*_height);
  if (_grid[pos]->get_size() > 0)
  {
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    for (int i = 0; i < (int)_grid[pos]->get_size(); i++)
    {
      file << i+1 << " " << _grid[pos]->get(i+1) << "\n";
    }
    file.close();
  }
}

/**
 * \brief    Update the environment
 * \details  Apply environment variations
 * \param    bool initialize
 * \param    double factor
 * \param    size_t time
 * \return   \e void
 */
void Environment::update( bool initialize, size_t time )
{
  _inflowing_amount  = 0.0;
  _outflowing_amount = 0.0;
  
  /*-------------------------------------------*/
  /* 1) Decide if the environment must vary    */
  /*-------------------------------------------*/
  bool   apply_variation = false;
  double factor          = 1.0;
  if (_properties->variation_scheme == RANDOM_SCHEME)
  {
    apply_variation = (_prng->uniform() < _properties->introduction_rate);
  }
  else if (_properties->variation_scheme == PERIODIC_SCHEME)
  {
    int period      = (int)1.0/_properties->introduction_rate;
    apply_variation = (time%period == 0);
  }
  else if (_properties->variation_scheme == CYCLIC_SCHEME)
  {
    double period   = 1.0/_properties->introduction_rate;
    factor          = (sin(time*2.0*3.14159265359/period)+1.0)/2.0;
    apply_variation = true;
  }
  if (initialize)
  {
    apply_variation = true;
  }
  
  /*-------------------------------------------*/
  /* 2) If the environment must vary, evaluate */
  /*    if it must be rinsed                   */
  /*-------------------------------------------*/
  if (apply_variation && _properties->renewal_scheme == CLEAR_MATTER)
  {
    clear();
  }
  
  /*-------------------------------------------*/
  /* 3) Update the environment                 */
  /*-------------------------------------------*/
  if (apply_variation && _properties->localization_scheme == GLOBAL_LOCALIZATION && (_properties->metabolic_scheme == UNIQUE_METABOLITE || _properties->metabolic_scheme == BOUNDARIES))
  {
    update_global_unique(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == GLOBAL_LOCALIZATION && _properties->metabolic_scheme == MULTIPLE_METABOLITES)
  {
    update_global_multiple(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == RANDOM_LOCALIZATION && (_properties->metabolic_scheme == UNIQUE_METABOLITE || _properties->metabolic_scheme == BOUNDARIES))
  {
    update_random_unique(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == RANDOM_LOCALIZATION && _properties->metabolic_scheme == MULTIPLE_METABOLITES)
  {
    update_random_multiple(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == SPOT_LOCALIZATION && (_properties->metabolic_scheme == UNIQUE_METABOLITE || _properties->metabolic_scheme == BOUNDARIES))
  {
    update_spot_unique(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == SPOT_LOCALIZATION && _properties->metabolic_scheme == MULTIPLE_METABOLITES)
  {
    update_spot_multiple(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == CENTER_LOCALIZATION && (_properties->metabolic_scheme == UNIQUE_METABOLITE || _properties->metabolic_scheme == BOUNDARIES))
  {
    update_center_unique(initialize, factor);
  }
  else if (apply_variation && _properties->localization_scheme == CENTER_LOCALIZATION && _properties->metabolic_scheme == MULTIPLE_METABOLITES)
  {
    update_center_multiple(initialize, factor);
  }
}

/**
 * \brief    Clear the environment
 * \details  Remove all metabolites from the environment
 * \param    void
 * \return   \e void
 */
void Environment::clear( void )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    _outflowing_amount += _grid[i]->get_amount();
    _grid[i]->reset(false);
  }
  _X->reset(false);
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Compute classic diffusion and degradation
 * \details  Compute classic diffusion and degradation of metabolites in environment
 * \param    void
 * \return   \e void
 */
void Environment::compute_classic_diffusion_and_degradation( void )
{
  /*-----------------------------------------------------------------*/
  /* 1) Copy the grid content and compute diffusion first step       */
  /*-----------------------------------------------------------------*/
  double* new_grid  = new double[_width*_height*_size];
  double  diff_coef = _properties->diffusion_coefficient;
  double  degr_rate = _properties->degradation_rate;
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      memcpy(&new_grid[x*_height*_size+y*_size], _grid[x*_height+y]->get_X(), sizeof(double)*_size);
      for (int i = -1; i < 2; i++)
      {
        for (int j = -1; j < 2; j++)
        {
          size_t  cur_x = (x + i + _width)  % _width;
          size_t  cur_y = (y + j + _height) % _height;
          double* cur_v = _grid[cur_x*_height+cur_y]->get_X();
          for (size_t index = 0; index < _size; index++)
          {
            new_grid[x*_height*_size+y*_size+index] += cur_v[index]*diff_coef;
          }
        }
      }
    }
  }
  /*-----------------------------------------------------------------*/
  /* 2) Compute the second step of diffusion and compute degradation */
  /*-----------------------------------------------------------------*/
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      for (int index = 0; index < (int)_size; index++)
      {
        double fresh_amount = new_grid[x*_height*_size+y*_size+index] - 9.0*_grid[x*_height+y]->get_X()[index]*diff_coef;
        double new_amount   = fresh_amount*(1.0-degr_rate);
        _outflowing_amount += fresh_amount*degr_rate;
        _grid[x*_height+y]->set(index+1, new_amount);
      }
    }
  }
  delete[] new_grid;
  new_grid = NULL;
}

/**
 * \brief    Compute perfect diffusion and degradation
 * \details  Compute perfect diffusion and degradation of metabolites in environment
 * \param    void
 * \return   \e void
 */
void Environment::compute_perfect_diffusion_and_degradation( void )
{
  /*---------------------------------------------------*/
  /* 1) Initialize and compute the sum vector          */
  /*---------------------------------------------------*/
  double* sum_vector = new double[_size];
  for (size_t index = 0; index < _size; index++)
  {
    sum_vector[index] = 0.0;
  }
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      double* vector = _grid[x*_height+y]->get_X();
      for (size_t index = 0; index < _size; index++)
      {
        sum_vector[index] += vector[index];
      }
    }
  }
  for (size_t index = 0; index < _size; index++)
  {
    sum_vector[index] /= (double)_width*_height;
  }
  
  /*---------------------------------------------------*/
  /* 2) Compute the perfect diffusion with degradation */
  /*---------------------------------------------------*/
  double degr_rate = _properties->degradation_rate;
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      for (int index = 0; index < (int)_size; index++)
      {
        double fresh_amount = sum_vector[index];
        double new_amount   = fresh_amount*(1.0-degr_rate);
        _outflowing_amount += fresh_amount*degr_rate;
        _grid[x*_height+y]->set(index+1, new_amount);
      }
    }
  }
  delete[] sum_vector;
  sum_vector = NULL;
}

/**
 * \brief    Update the environment with global localization and unique metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_global_unique( bool initialize, double factor )
{
  int species_tag = 0;
  if (_properties->metabolic_scheme == UNIQUE_METABOLITE)
  {
    species_tag = _parameters->draw_environment_species_tag();
  }
  else if (_properties->metabolic_scheme == BOUNDARIES)
  {
    species_tag = (_prng->uniform() < 0.5 ? _properties->species_tag_range.min : _properties->species_tag_range.max);
  }
  double concentration = _parameters->draw_environment_concentration();
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      if (initialize)
      {
        set(x, y, species_tag, concentration*factor);
      }
      else
      {
        add(x, y, species_tag, concentration*factor);
      }
      _inflowing_amount += concentration*factor;
    }
  }
}

/**
 * \brief    Update the environment with global localization and multiple metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_global_multiple( bool initialize, double factor )
{
  size_t  number_of_species = _parameters->draw_environment_number_of_species();
  int*    species_tag       = new int[number_of_species];
  double* concentration     = new double[number_of_species];
  for (size_t nb = 0; nb < number_of_species; nb++)
  {
    species_tag[nb]   = _parameters->draw_environment_species_tag();
    concentration[nb] = _parameters->draw_environment_concentration();
  }
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      for (size_t nb = 0; nb < number_of_species; nb++)
      {
        if (initialize)
        {
          set(x, y, species_tag[nb], concentration[nb]*factor);
        }
        else
        {
          add(x, y, species_tag[nb], concentration[nb]*factor);
        }
        _inflowing_amount += concentration[nb]*factor;
      }
    }
  }
  delete[] species_tag;
  species_tag = NULL;
  delete[] concentration;
  concentration = NULL;
}

/**
 * \brief    Update the environment with random localization and unique metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_random_unique( bool initialize, double factor )
{
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      int species_tag = 0;
      if (_properties->metabolic_scheme == UNIQUE_METABOLITE)
      {
        species_tag = _parameters->draw_environment_species_tag();
      }
      else if (_properties->metabolic_scheme == BOUNDARIES)
      {
        species_tag = (_prng->uniform() < 0.5 ? _properties->species_tag_range.min : _properties->species_tag_range.max);
      }
      double concentration = _parameters->draw_environment_concentration();
      if (initialize)
      {
        set(x, y, species_tag, concentration*factor);
      }
      else
      {
        add(x, y, species_tag, concentration*factor);
      }
      _inflowing_amount += concentration*factor;
    }
  }
}

/**
 * \brief    Update the environment with random localization and multiple metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_random_multiple( bool initialize, double factor )
{
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      size_t number_of_species = _parameters->draw_environment_number_of_species();
      for (size_t nb = 0; nb < number_of_species; nb++)
      {
        int    species_tag   = _parameters->draw_environment_species_tag();
        double concentration = _parameters->draw_environment_concentration();
        if (initialize)
        {
          set(x, y, species_tag, concentration*factor);
        }
        else
        {
          add(x, y, species_tag, concentration*factor);
        }
        _inflowing_amount += concentration*factor;
      }
    }
  }
}

/**
 * \brief    Update the environment with spot localization and unique metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_spot_unique( bool initialize, double factor )
{
  size_t x_draw        = (size_t)_prng->uniform(0, (int)_width);
  size_t y_draw        = (size_t)_prng->uniform(0, (int)_height);
  int species_tag = 0;
  if (_properties->metabolic_scheme == UNIQUE_METABOLITE)
  {
    species_tag = _parameters->draw_environment_species_tag();
  }
  else if (_properties->metabolic_scheme == BOUNDARIES)
  {
    species_tag = (_prng->uniform() < 0.5 ? _properties->species_tag_range.min : _properties->species_tag_range.max);
  }
  double concentration = _parameters->draw_environment_concentration();
  if (initialize)
  {
    set(x_draw, y_draw, species_tag, concentration*factor);
  }
  else
  {
    add(x_draw, y_draw, species_tag, concentration*factor);
  }
  _inflowing_amount += concentration*factor;
}

/**
 * \brief    Update the environment with spot localization and multiple metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_spot_multiple( bool initialize, double factor )
{
  size_t x_draw            = (size_t)_prng->uniform(0, (int)_width);
  size_t y_draw            = (size_t)_prng->uniform(0, (int)_height);
  size_t number_of_species = _parameters->draw_environment_number_of_species();
  for (size_t nb = 0; nb < number_of_species; nb++)
  {
    int    species_tag   = _parameters->draw_environment_species_tag();
    double concentration = _parameters->draw_environment_concentration();
    if (initialize)
    {
      set(x_draw, y_draw, species_tag, concentration*factor);
    }
    else
    {
      add(x_draw, y_draw, species_tag, concentration*factor);
    }
    _inflowing_amount += concentration*factor;
  }
}

/**
 * \brief    Update the environment with center localization and unique metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_center_unique( bool initialize, double factor )
{
  int species_tag = 0;
  if (_properties->metabolic_scheme == UNIQUE_METABOLITE)
  {
    species_tag = _parameters->draw_environment_species_tag();
  }
  else if (_properties->metabolic_scheme == BOUNDARIES)
  {
    species_tag = (_prng->uniform() < 0.5 ? _properties->species_tag_range.min : _properties->species_tag_range.max);
  }
  double concentration = _parameters->draw_environment_concentration();
  if (initialize)
  {
    set(_width/2, _height/2, species_tag, concentration*factor);
  }
  else
  {
    add(_width/2, _height/2, species_tag, concentration*factor);
  }
  _inflowing_amount += concentration*factor;
}

/**
 * \brief    Update the environment with center localization and multiple metabolite scheme
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Environment::update_center_multiple( bool initialize, double factor )
{
  size_t number_of_species = _parameters->draw_environment_number_of_species();
  for (size_t nb = 0; nb < number_of_species; nb++)
  {
    int    species_tag   = _parameters->draw_environment_species_tag();
    double concentration = _parameters->draw_environment_concentration();
    if (initialize)
    {
      set(_width/2, _height/2, species_tag, concentration*factor);
    }
    else
    {
      add(_width/2, _height/2, species_tag, concentration*factor);
    }
    _inflowing_amount += concentration*factor;
  }
}
