
/**
 * \file      Cell.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Cell class definition
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

#include "Cell.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Prng* prng
 * \return   \e void
 */
Cell::Cell( Parameters* parameters, Prng* prng )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ prng */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ main cell classes */
  
  _replication_report = new ReplicationReport();
  _genome             = new Genome(_parameters, _prng, _replication_report);
  _inherited_proteins = NULL;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _inherited_proteins = new InheritedProteins(_parameters);
  }
  _inherited_species_list = new SpeciesList();
  _species_list           = new SpeciesList();
  
  /*------------------------------------------------------------------ ODE system */
  
  _ode = NULL;
  
  /*------------------------------------------------------------------ main cell variables */
  
  _id                  = 0;
  _parent_id           = 0;
  _generation          = 0;
  _energy              = 0.0;
  _active              = false;
  _alive               = false;
  _x                   = 0;
  _y                   = 0;
  _score               = 0.0;
  _number_of_updates   = 0;
  _number_of_divisions = 0;
  _tagged              = false;
  
  /*------------------------------------------------------------------ mutation rates */
  
  initialize_mutation_rates();
  
  /*------------------------------------------------------------------ time variables */
  
  _birth_time = 0;
  _death_time = 0;
  _lifespan   = 0;
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  _toxicity                   = 0.0;
  _inherited_TF_amount        = 0.0;
  _inherited_E_amount         = 0.0;
  _TF_amount                  = 0.0;
  _E_amount                   = 0.0;
  _min_metabolic_amount       = 1e+6;
  _max_metabolic_amount       = 0.0;
  _metabolic_uptake           = 0.0;
  _metabolic_release          = 0.0;
  _previous_metabolic_amount  = 0.0;
  _min_energy                 = 1e+6;
  _mean_energy                = 0.0;
  _max_energy                 = 0.0;
  _min_score                  = 1e+6;
  _mean_score                 = 0.0;
  _max_score                  = 0.0;
  _metabolic_growth_rate      = 0.0;
  _diff_metabolic_growth_rate = 0.0;
  _grn_nb_nodes               = 0;
  _grn_nb_edges               = 0;
  _metabolic_nb_nodes         = 0;
  _metabolic_nb_edges         = 0;
  _inflowing_pumps.clear();
  _outflowing_pumps.clear();
  _trophic_group = 0;
  _trophic_level = NO_LEVEL;
  
  /*------------------------------------------------------------------ cell color */
  
  _red_color   = 0.0;
  _green_color = 0.0;
  _blue_color  = 0.0;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Cell class from backup file
 * \param    Parameters* parameters
 * \param    Prng* prng
 * \param    gzFile backup_file
 * \return   \e void
 */
Cell::Cell( Parameters* parameters, Prng* prng, gzFile backup_file )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ prng */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ main cell classes */
  
  _replication_report = new ReplicationReport(backup_file);
  _genome             = new Genome(_parameters, _prng, _replication_report, backup_file);
  _inherited_proteins = NULL;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _inherited_proteins = new InheritedProteins(_parameters, backup_file);
  }
  _inherited_species_list = new SpeciesList(backup_file);
  _species_list           = new SpeciesList(backup_file);
  
  /*------------------------------------------------------------------ ODE system */
  
  _ode = NULL;
  
  /*------------------------------------------------------------------ main cell variables */
  
  gzread( backup_file, &_id,                  sizeof(_id) );
  gzread( backup_file, &_parent_id,           sizeof(_parent_id) );
  gzread( backup_file, &_generation,          sizeof(_generation) );
  gzread( backup_file, &_energy,              sizeof(_energy) );
  gzread( backup_file, &_active,              sizeof(_active) );
  gzread( backup_file, &_alive,               sizeof(_alive) );
  gzread( backup_file, &_x,                   sizeof(_x) );
  gzread( backup_file, &_y,                   sizeof(_y) );
  gzread( backup_file, &_score,               sizeof(_score) );
  gzread( backup_file, &_number_of_updates,   sizeof(_number_of_updates) );
  gzread( backup_file, &_number_of_divisions, sizeof(_number_of_divisions) );
  gzread( backup_file, &_tagged,              sizeof(_tagged) );
  
  /*------------------------------------------------------------------ mutation rates */
  
  load_mutation_rates(backup_file);
  
  /*------------------------------------------------------------------ time variables */
  
  gzread( backup_file, &_birth_time, sizeof(_birth_time) );
  gzread( backup_file, &_death_time, sizeof(_death_time) );
  gzread( backup_file, &_lifespan,   sizeof(_lifespan) );
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  gzread( backup_file, &_toxicity,                   sizeof(_toxicity) );
  gzread( backup_file, &_inherited_TF_amount,        sizeof(_inherited_TF_amount) );
  gzread( backup_file, &_inherited_E_amount,         sizeof(_inherited_E_amount) );
  gzread( backup_file, &_TF_amount,                  sizeof(_TF_amount) );
  gzread( backup_file, &_E_amount,                   sizeof(_E_amount) );
  gzread( backup_file, &_min_metabolic_amount,       sizeof(_min_metabolic_amount) );
  gzread( backup_file, &_max_metabolic_amount,       sizeof(_max_metabolic_amount) );
  gzread( backup_file, &_metabolic_uptake,           sizeof(_metabolic_uptake) );
  gzread( backup_file, &_metabolic_release,          sizeof(_metabolic_release) );
  gzread( backup_file, &_previous_metabolic_amount,  sizeof(_previous_metabolic_amount) );
  gzread( backup_file, &_min_energy,                 sizeof(_min_energy) );
  gzread( backup_file, &_mean_energy,                sizeof(_mean_energy) );
  gzread( backup_file, &_max_energy,                 sizeof(_max_energy) );
  gzread( backup_file, &_min_score,                  sizeof(_min_score) );
  gzread( backup_file, &_mean_score,                 sizeof(_mean_score) );
  gzread( backup_file, &_max_score,                  sizeof(_max_score) );
  gzread( backup_file, &_metabolic_growth_rate,      sizeof(_metabolic_growth_rate) );
  gzread( backup_file, &_diff_metabolic_growth_rate, sizeof(_diff_metabolic_growth_rate) );
  gzread( backup_file, &_grn_nb_nodes,               sizeof(_grn_nb_nodes) );
  gzread( backup_file, &_grn_nb_edges,               sizeof(_grn_nb_edges) );
  gzread( backup_file, &_metabolic_nb_nodes,         sizeof(_metabolic_nb_nodes) );
  gzread( backup_file, &_metabolic_nb_edges,         sizeof(_metabolic_nb_edges) );
  size_t N_inflowing_pumps  = 0;
  size_t N_outflowing_pumps = 0;
  gzread( backup_file, &N_inflowing_pumps, sizeof(N_inflowing_pumps) );
  gzread( backup_file, &N_outflowing_pumps, sizeof(N_outflowing_pumps) );
  _inflowing_pumps.assign(N_inflowing_pumps, 0);
  _outflowing_pumps.assign(N_outflowing_pumps, 0);
  for (size_t i = 0; i < N_inflowing_pumps; i++)
  {
    gzread( backup_file, &_inflowing_pumps[i], sizeof(_inflowing_pumps[i]) );
  }
  for (size_t i = 0; i < N_outflowing_pumps; i++)
  {
    gzread( backup_file, &_outflowing_pumps[i], sizeof(_outflowing_pumps[i]) );
  }
  gzread( backup_file, &_trophic_group, sizeof(_trophic_group) );
  gzread( backup_file, &_trophic_level, sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ cell color */
  
  gzread( backup_file, &_red_color,   sizeof(_red_color) );
  gzread( backup_file, &_green_color, sizeof(_green_color) );
  gzread( backup_file, &_blue_color,  sizeof(_blue_color) );
  
}

/**
 * \brief    Constructor from parental cell
 * \details  Create an offspring cell from the parent
 * \param    Cell& parent
 * \param    Prng* prng
 * \param    unsigned long long int new_id
 * \param    size_t x
 * \param    size_t y
 * \param    size_t time
 * \return   \e void
 */
Cell::Cell( Cell& parent, Prng* prng, unsigned long long int new_id, size_t x, size_t y, size_t time)
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parent._parameters;
  
  /*------------------------------------------------------------------ prng */
  
  _prng = prng;
  
  /*------------------------------------------------------------------ main cell classes */
  
  _replication_report = new ReplicationReport();
  
  /*------------------------*/
  /* copy genome            */
  /*------------------------*/
  _genome = new Genome(*parent._genome, _prng, _replication_report);
  
  /*------------------------*/
  /* get inherited proteins */
  /*------------------------*/
  _inherited_proteins  = NULL;
  _inherited_TF_amount = 0.0;
  _inherited_E_amount  = 0.0;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    parent.load_genome_in_inherited_proteins();
    _inherited_proteins = new InheritedProteins(*parent._inherited_proteins);
    for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
    {
      if (_inherited_proteins->get_genetic_unit(i)->type == TRANSCRIPTION_FACTOR)
      {
        _inherited_TF_amount += _inherited_proteins->get_concentration_vector()[i];
      }
      else if (_inherited_proteins->get_genetic_unit(i)->type == ENZYME)
      {
        _inherited_E_amount += _inherited_proteins->get_concentration_vector()[i];
      }
    }
  }
  
  /*-------------------*/
  /* copy species list */
  /*-------------------*/
  
  /* copy the parent list */
  _inherited_species_list = new SpeciesList(*parent._species_list);
  _species_list           = new SpeciesList(*parent._species_list);
  
  /* apply metabolic inheritance rules */
  _inherited_species_list->reset(_parameters->get_metabolic_inheritance_scheme());
  _species_list->reset(_parameters->get_metabolic_inheritance_scheme());
  parent._species_list->reset(_parameters->get_metabolic_inheritance_scheme());
  
  /* compute inherited amounts and growth rate */
  _previous_metabolic_amount  = _inherited_species_list->get_amount();
  _metabolic_growth_rate      = _species_list->get_amount()-_previous_metabolic_amount;
  _diff_metabolic_growth_rate = 0.0;
  
  /*------------------------------------------------------------------ ODE system */
  
  _ode = NULL;
  
  /*------------------------------------------------------------------ identifiers */
  
  _id        = new_id;
  _parent_id = parent._id;
  
  /*------------------------------------------------------------------ generation */
  
  parent.set_generation(parent._generation+1);
  _generation = parent._generation+1;
  
  /*------------------------------------------------------------------ energy carrier */
  
  if (_parameters->get_energy_costs_scheme())
  {
    _energy = parent.get_energy()/2.0;
    parent.set_energy(_energy);
    if (_energy < MINIMUM_CONCENTRATION)
    {
      _energy = 0.0;
      parent.set_energy(0.0);
    }
  }
  else
  {
    _energy = 0.0;
  }
  
  /*------------------------------------------------------------------ active and alive states */
  
  _active = true;
  _alive  = true;
  
  /*------------------------------------------------------------------ coordinates */
  
  _x = x;
  _y = y;
  
  /*------------------------------------------------------------------ score, ability to divide and means */
  
  _score               = 0.0;
  _number_of_updates   = 0;
  _number_of_divisions = 0;
  _tagged              = false;
  
  /*------------------------------------------------------------------ mutation rates */
  
  copy_mutation_rates(parent._mutation_rates);
  
  /*------------------------------------------------------------------ time variables */
  
  _birth_time = time;
  _death_time = time;
  _lifespan   = 0;
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  _toxicity             = 0.0;
  _TF_amount            = 0.0;
  _E_amount             = 0.0;
  _min_metabolic_amount = 1e+6;
  _max_metabolic_amount = 0.0;
  _metabolic_uptake     = 0.0;
  _metabolic_release    = 0.0;
  _min_energy           = 1e+6;
  _mean_energy          = 0.0;
  _max_energy           = 0.0;
  _min_score            = 1e+6;
  _mean_score           = 0.0;
  _max_score            = 0.0;
  _grn_nb_nodes         = 0;
  _grn_nb_edges         = 0;
  _metabolic_nb_nodes   = 0;
  _metabolic_nb_edges   = 0;
  _inflowing_pumps.clear();
  _outflowing_pumps.clear();
  _trophic_group = 0;
  _trophic_level = NO_LEVEL;
  
  /*------------------------------------------------------------------ cell color */
  
  _red_color   = 0.0;
  _green_color = 0.0;
  _blue_color  = 0.0;
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
Cell::~Cell( void )
{
  delete _replication_report;
  _replication_report = NULL;
  delete _genome;
  _genome = NULL;
  delete _inherited_proteins;
  _inherited_proteins = NULL;
  delete _inherited_species_list;
  _inherited_species_list = NULL;
  delete _species_list;
  _species_list = NULL;
  delete _ode;
  _ode = NULL;
  delete[] _mutation_rates;
  _mutation_rates = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief   Initialize the inherited species list with the current species list
 * \details --
 * \param   void
 * \return  \e double
 */
void Cell::initialize_inherited_species_list( void )
{
  delete _inherited_species_list;
  _inherited_species_list = new SpeciesList(*_species_list);
}

/**
 * \brief   Load genome in inherited proteins
 * \details At replication, proteins production must be saved before cell's division. First, the inherited protein list is cleaned (empty locations are removed), then the genome content is loaded. All concentrations are divided by 2. The new inherited protein list is copied as is in the daughter cell.
 * \param   void
 * \return  \e void
 */
void Cell::load_genome_in_inherited_proteins( void )
{
  double* new_conc_vector = new double[_genome->get_size()+_inherited_proteins->get_size()];
  size_t  new_N           = 0;
  
  /*----------------------------------------------------------*/
  /* 1) clean inherited proteins list (remove empty content)  */
  /*----------------------------------------------------------*/
  double* conc_vector        = _inherited_proteins->get_concentration_vector();
  genetic_unit* str          = _inherited_proteins->get_genetic_sequence()->x;
  InheritedProteins* new_seq = new InheritedProteins(_parameters);
  for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
  {
    if (conc_vector[i]*0.5 > MINIMUM_HERITABLE_ENZYME_CONCENTRATION)
    {
      new_seq->add_genetic_unit(&str[i]);
      new_conc_vector[new_N] = conc_vector[i]*0.5;
      new_N++;
    }
  }
  
  /*----------------------------------------------------------*/
  /* 2) load genome content                                   */
  /*----------------------------------------------------------*/
  conc_vector = _genome->get_concentration_vector();
  str         = _genome->get_genetic_sequence()->x;
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    if ((str[i].type == ENZYME || str[i].type == TRANSCRIPTION_FACTOR) && conc_vector[i]*0.5 > MINIMUM_HERITABLE_ENZYME_CONCENTRATION)
    {
      new_seq->add_genetic_unit(&str[i]);
      new_conc_vector[new_N] = conc_vector[i]*0.5;
      new_N++;
    }
  }
  
  /*----------------------------------------------------------*/
  /* 3) replace the old inherited protein list by the new one */
  /*----------------------------------------------------------*/
  delete _inherited_proteins;
  _inherited_proteins = new_seq;
  _inherited_proteins->initialize_concentration_vector();
  assert(_inherited_proteins->get_size() == new_N);
  memcpy(_inherited_proteins->get_concentration_vector(), new_conc_vector, sizeof(double)*new_N);
  delete[] new_conc_vector;
  new_conc_vector = NULL;
  
  /*----------------------------------------------------------*/
  /* 4) build E and TF indexes list                           */
  /*----------------------------------------------------------*/
  _inherited_proteins->build_index_list();
}

/**
 * \brief   Load genome in species list
 * \details After replication, some metabolic species may appear. The size of species list is increased consequently.
 * \param   void
 * \return  \e void
 */
void Cell::load_genome_in_species_lists( void )
{
  genetic_unit* str = _genome->get_genetic_sequence()->x;
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    if (str[i].type == ENZYME)
    {
      if (_species_list->get_size() < (size_t)str[i].s)
      {
        _inherited_species_list->increase_size((size_t)str[i].s);
        _species_list->increase_size((size_t)str[i].s);
      }
      if (_species_list->get_size() < (size_t)str[i].p)
      {
        _inherited_species_list->increase_size((size_t)str[i].p);
        _species_list->increase_size((size_t)str[i].p);
      }
    }
  }
}

/**
 * \brief    Load genome in ODE system
 * \details  Load genome reactions rules in ODE system
 * \param    Environment* environment
 * \param    bool from_backup
 * \param    bool new_individual
 * \return   \e void
 */
void Cell::load_genome_in_ODE_system( Environment* environment, bool from_backup, bool new_individual )
{
  delete _ode;
  _ode = NULL;
  _ode = new ODE(_parameters, _replication_report, _genome, _inherited_proteins, _species_list, environment, _x, _y, &_energy, &_metabolic_uptake, &_metabolic_release, &_grn_nb_nodes, &_grn_nb_edges, &_metabolic_nb_nodes, &_metabolic_nb_edges, &_inflowing_pumps, &_outflowing_pumps);
  _ode->load(from_backup, new_individual);
}

/**
 * \brief    Synchronize state vectors
 * \details  Synchronize metabolic state vectors size between the cell and its environment, to avoid issues during membrane diffusion
 * \param    Environment* environment
 * \return   \e void
 */
void Cell::synchronize_state_vectors( Environment* environment )
{
  /*------------------------------------------------------------------*/
  /* 1) Synchronize inherited species list and species list sizes     */
  /*------------------------------------------------------------------*/
  if (_inherited_species_list->get_size() < _species_list->get_size())
  {
    _inherited_species_list->increase_size(_species_list->get_size());
  }
  else if (_inherited_species_list->get_size() > _species_list->get_size())
  {
    _species_list->increase_size(_inherited_species_list->get_size());
  }
  
  /*------------------------------------------------------------------*/
  /* 2) Synchronize cell's species lists and environment species list */
  /*------------------------------------------------------------------*/
  if (_species_list->get_size() < environment->get_size(_x, _y))
  {
    _inherited_species_list->increase_size(environment->get_size(_x, _y));
    _species_list->increase_size(environment->get_size(_x, _y));
  }
  else if (environment->get_size(_x, _y) < _species_list->get_size())
  {
    environment->increase_size(_x, _y, (int)_species_list->get_size());
  }
}


/**
 * \brief    Update the replication report data
 * \details  --
 * \param    void
 * \return   \e void
 */
void Cell::update_replication_report_data( void )
{
  _replication_report->set_id(_id);
  _replication_report->set_parent_id(_parent_id);
  _replication_report->set_generation(_generation);
  _replication_report->set_x(_x);
  _replication_report->set_y(_y);
  _replication_report->set_number_of_updates(_number_of_updates);
  _replication_report->set_number_of_divisions(_number_of_divisions);
  _replication_report->set_birth_time(_birth_time);
  _replication_report->set_death_time(_death_time);
  _replication_report->set_lifespan(_lifespan);
  _replication_report->set_toxicity(_toxicity);
  _replication_report->set_inherited_TF_amount(_inherited_TF_amount);
  _replication_report->set_inherited_E_amount(_inherited_E_amount);
  _replication_report->set_TF_amount(_TF_amount);
  _replication_report->set_E_amount(_E_amount);
  _replication_report->set_inherited_metabolic_amount(_inherited_species_list->get_amount());
  _replication_report->set_min_metabolic_amount(_min_metabolic_amount);
  _replication_report->set_metabolic_amount(_species_list->get_amount());
  _replication_report->set_max_metabolic_amount(_max_metabolic_amount);
  _replication_report->set_metabolic_uptake(_metabolic_uptake);
  _replication_report->set_metabolic_release(_metabolic_release);
  _replication_report->set_min_energy(_min_energy);
  if (_alive)
  {
    _replication_report->set_mean_energy(_energy);
  }
  else
  {
    _replication_report->set_mean_energy(_mean_energy);
  }
  _replication_report->set_mean_energy(_mean_energy);
  _replication_report->set_max_energy(_max_energy);
  _replication_report->set_min_score(_min_score);
  if (_alive)
  {
    _replication_report->set_mean_score(_score);
  }
  else
  {
    _replication_report->set_mean_score(_mean_score);
  }
  _replication_report->set_max_score(_max_score);
  _replication_report->set_metabolic_growth_rate(_metabolic_growth_rate);
  _replication_report->set_Dmetabolic_growth_rate(_diff_metabolic_growth_rate);
  _replication_report->set_grn_nb_nodes(_grn_nb_nodes);
  _replication_report->set_grn_nb_edges(_grn_nb_edges);
  _replication_report->set_metabolic_nb_nodes(_metabolic_nb_nodes);
  _replication_report->set_metabolic_nb_edges(_metabolic_nb_edges);
  _replication_report->set_trophic_group(_trophic_group);
  _replication_report->set_trophic_level(_trophic_level);
}

/**
 * \brief    Update cell's state
 * \details  --
 * \param    size_t time
 * \return   \e void
 */
void Cell::update( size_t time )
{
  /*---------------------------*/
  /* 1) update cell's networks */
  /*---------------------------*/
  if (!_parameters->get_energy_costs_scheme())
  {
    _energy = 0.0;
  }
  _ode->solve();
  _ode->update();
  compute_score();
  
  /*---------------------------*/
  /* 2) update mean statistics */
  /*---------------------------*/
  _inherited_TF_amount = 0.0;
  _inherited_E_amount  = 0.0;
  _TF_amount           = 0.0;
  _E_amount            = 0.0;
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
    {
      if (_inherited_proteins->get_genetic_unit(i)->type == TRANSCRIPTION_FACTOR)
      {
        _inherited_TF_amount += _inherited_proteins->get_concentration_vector()[i];
      }
      else if (_inherited_proteins->get_genetic_unit(i)->type == ENZYME)
      {
        _inherited_E_amount += _inherited_proteins->get_concentration_vector()[i];
      }
    }
  }
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    if (_genome->get_genetic_unit(i)->type == TRANSCRIPTION_FACTOR)
    {
      _TF_amount += _genome->get_concentration_vector()[i];
    }
    else if (_genome->get_genetic_unit(i)->type == ENZYME)
    {
      _E_amount += _genome->get_concentration_vector()[i];
    }
  }
  _min_metabolic_amount = (_min_metabolic_amount < _species_list->get_amount() ? _min_metabolic_amount : _species_list->get_amount());
  _max_metabolic_amount = (_max_metabolic_amount > _species_list->get_amount() ? _max_metabolic_amount : _species_list->get_amount());
  if (_score > 0.0)
  {
    _min_energy   = (_min_energy < _energy ? _min_energy : _energy);
    _max_energy   = (_max_energy > _energy ? _max_energy : _energy);
    _min_score    = (_min_score < _score ? _min_score : _score);
    _max_score    = (_max_score > _score ? _max_score : _score);
    _mean_score  += _score;
    _mean_energy += _energy;
    _number_of_updates++;
  }
  _diff_metabolic_growth_rate = _species_list->get_amount()-_previous_metabolic_amount-_metabolic_growth_rate;
  _metabolic_growth_rate      = _species_list->get_amount()-_previous_metabolic_amount;
  _previous_metabolic_amount  = _species_list->get_amount();
  
  /*-----------------------------------*/
  /* 3) update replication report data */
  /*-----------------------------------*/
  _lifespan = time-_birth_time;
  update_replication_report_data();
}

/**
 * \brief    Compute score
 * \details  Compute cell's score depending on division criteria
 * \param    void
 * \return   \e void
 */
void Cell::compute_score( void )
{
  /*************************************************************************/
  /* A) Evaluate energy constraints                                        */
  /*************************************************************************/
  
  /*-----------------------------------------------*/
  /* A.1) If energy constraints are activated, and */
  /*      energy is below 0, kill the cell.        */
  /*-----------------------------------------------*/
  if (_parameters->get_energy_costs_scheme() && _energy <= 0.0)
  {
    _score = 0.0;
    return;
  }
  
  /*-----------------------------------------------*/
  /* A.2) If energy constraints are activated,     */
  /*      and energy is above the energy toxicity  */
  /*      threshold, kill the cell.                */
  /*-----------------------------------------------*/
  else if (_parameters->get_energy_costs_scheme() && _energy > _parameters->get_energy_toxicity_threshold())
  {
    _score = 0.0;
    return;
  }
  
  /*-----------------------------------------------*/
  /* A.3) If energy constraints are not activated, */
  /*      just set the energy at 0.                */
  /*-----------------------------------------------*/
  else if (!_parameters->get_energy_costs_scheme())
  {
    _energy = 0.0;
  }
  
  /*************************************************************************/
  /* B) If the score is the sum of the essential metabolites               */
  /*************************************************************************/
  if (_parameters->get_score_scheme() == ESSENTIAL_METABOLITES_SUM)
  {
    _score                             = 0.0;
    double  toxicity                   = 0.0;
    int     toxicity_count             = 0;
    bool    toxicity_threshold_reached = false;
    int*    prime_numbers              = _parameters->get_prime_numbers();
    double* X                          = _species_list->get_X();
    double  essential_threshold        = _parameters->get_essential_metabolites_toxicity_threshold();
    double  non_essential_threshold    = _parameters->get_non_essential_metabolites_toxicity_threshold();
    
    /*----------------------------------------*/
    /* B.1) explore the metabolites state     */
    /*      vector                            */
    /*----------------------------------------*/
    for (size_t i = 0; i < _species_list->get_size(); i++)
    {
      /*------------------------------------*/
      /* IF THE METABOLITE IS ESSENTIAL     */
      /*------------------------------------*/
      if (prime_numbers[i] == 1)
      {
        if (X[i] < essential_threshold)
        {
          _score += X[i];
        }
        else if (X[i] >= essential_threshold)
        {
          //toxicity += (X[i]-essential_threshold)*(X[i]-essential_threshold);
          toxicity += X[i]-essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
      /*------------------------------------*/
      /* IF THE METABOLITE IS NON ESSENTIAL */
      /*------------------------------------*/
      else if (prime_numbers[i] == 0)
      {
        if (X[i] >= non_essential_threshold)
        {
          //toxicity += (X[i]-non_essential_threshold)*(X[i]-non_essential_threshold);
          toxicity += X[i]-non_essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
    }
    /*----------------------------------------*/
    /* B.2) Apply toxicity threshold          */
    /*----------------------------------------*/
    if (toxicity_threshold_reached)
    {
      _score = 0.0;
    }
    else
    {
      _score = (_score > 0.0 ? _score : 0.0);
    }
  }
  
  /*************************************************************************/
  /* C) If the score is the essential metabolites mean minus the deviation */
  /*************************************************************************/
  else if (_parameters->get_score_scheme() == ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION)
  {
    std::vector<double> conc_list;
    _score                             = 0.0;
    double  toxicity                   = 0.0;
    int     toxicity_count             = 0;
    bool    toxicity_threshold_reached = false;
    int*    prime_numbers              = _parameters->get_prime_numbers();
    double* X                          = _species_list->get_X();
    double  essential_threshold        = _parameters->get_essential_metabolites_toxicity_threshold();
    double  non_essential_threshold    = _parameters->get_non_essential_metabolites_toxicity_threshold();
    
    /*-------------------------------------------*/
    /* C.1) explore the metabolites state vector */
    /*-------------------------------------------*/
    for (size_t i = 0; i < _species_list->get_size(); i++)
    {
      /*------------------------------------*/
      /* IF THE METABOLITE IS ESSENTIAL     */
      /*------------------------------------*/
      if (prime_numbers[i] == 1)
      {
        if (X[i] > 0.0 && X[i] < essential_threshold)
        {
          _score += X[i];
          conc_list.push_back(X[i]);
        }
        else if (X[i] >= essential_threshold)
        {
          //toxicity += (X[i]-essential_threshold)*(X[i]-essential_threshold);
          toxicity += X[i]-essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
      /*------------------------------------*/
      /* IF THE METABOLITE IS NON ESSENTIAL */
      /*------------------------------------*/
      else if (prime_numbers[i] == 0)
      {
        if (X[i] >= non_essential_threshold)
        {
          //toxicity += (X[i]-non_essential_threshold)*(X[i]-non_essential_threshold);
          toxicity += X[i]-non_essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
    }
    /*-------------------------------------------*/
    /* C.2) compute new score                    */
    /*-------------------------------------------*/
    if (conc_list.size() > 0.0)
    {
      double mean = 0.0;
      for (size_t i = 0; i < conc_list.size(); i++)
      {
        mean += conc_list[i];
      }
      mean /= conc_list.size();
      double quaderror = 0.0;
      for (size_t i = 0; i < conc_list.size(); i++)
      {
        quaderror += (conc_list[i]-mean)*(conc_list[i]-mean);
      }
      quaderror /= conc_list.size();
      _score -= quaderror;
      _score  = (_score > 0.0 ? _score : 0.0);
    }
    /*-------------------------------------------*/
    /* C.3) Apply toxicity threshold             */
    /*-------------------------------------------*/
    if (toxicity_threshold_reached)
    {
      _score = 0.0;
    }
  }
  
  /*************************************************************************/
  /* D) If the score is the essential metabolites "biggest complexes" sum  */
  /*************************************************************************/
  else if (_parameters->get_score_scheme() == ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION)
  {
    std::vector<double> conc_list;
    _score                             = 0.0;
    double  toxicity                   = 0.0;
    int     toxicity_count             = 0;
    bool    toxicity_threshold_reached = false;
    int*    prime_numbers              = _parameters->get_prime_numbers();
    double* X                          = _species_list->get_X();
    double  essential_threshold        = _parameters->get_essential_metabolites_toxicity_threshold();
    double  non_essential_threshold    = _parameters->get_non_essential_metabolites_toxicity_threshold();
    for (size_t i = 0; i < _species_list->get_size(); i++)
    {
      /*------------------------------------*/
      /* IF THE METABOLITE IS ESSENTIAL     */
      /*------------------------------------*/
      if (prime_numbers[i] == 1)
      {
        if (X[i] > 0.0 && X[i] < essential_threshold)
        {
          conc_list.push_back(X[i]);
        }
        else if (X[i] >= essential_threshold)
        {
          //toxicity += (X[i]-essential_threshold)*(X[i]-essential_threshold);
          toxicity += X[i]-essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
      /*------------------------------------*/
      /* IF THE METABOLITE IS NON ESSENTIAL */
      /*------------------------------------*/
      else if (prime_numbers[i] == 0)
      {
        if (X[i] >= non_essential_threshold)
        {
          //toxicity += (X[i]-non_essential_threshold)*(X[i]-non_essential_threshold);
          toxicity += X[i]-non_essential_threshold;
          toxicity_count++;
          toxicity_threshold_reached = true;
        }
      }
    }
    if (conc_list.size() > 0)
    {
      for (size_t i = 0; i < conc_list.size(); i++)
      {
        for (size_t j = i+1; j < conc_list.size(); j++)
        {
          if (conc_list[i] > conc_list[j])
          {
            double tmp = conc_list[i];
            conc_list[i] = conc_list[j];
            conc_list[j] = tmp;
          }
        }
      }
      size_t count = conc_list.size();
      for (size_t i = 0; i < conc_list.size(); i++)
      {
        if (i == 0)
        {
          _score += count*conc_list[i];
          count--;
        }
        else
        {
          _score += count*(conc_list[i]-conc_list[i-1]);
          count--;
        }
      }
      _score = (_score > 0.0 ? _score : 0.0);
    }
    if (toxicity_threshold_reached)
    {
      _score = 0.0;
    }
  }
}

/**
 * \brief    Mutate the cell
 * \details  Mutate cell's genome and mutation rates if required
 * \param    void
 * \return   \e void
 */
void Cell::mutate( void )
{
  /*----------------------------------------------*/
  /* 1) Clear the replication report              */
  /*----------------------------------------------*/
  _replication_report->clear();
  
  /*----------------------------------------------*/
  /* 2) Set the old genome size                   */
  /*----------------------------------------------*/
  _replication_report->set_old_genome_size(_genome->get_size());
  
  /*----------------------------------------------*/
  /* 3) Mutate the genome                         */
  /*----------------------------------------------*/
  if (_genome->get_size() > 0)
  {
    _genome->mutate(_mutation_rates);
  }
  
  /*----------------------------------------------*/
  /* 4) compute replication report mean data      */
  /*----------------------------------------------*/
  _replication_report->compute_mean();
  
  /*----------------------------------------------*/
  /* 5) update genome and inherited proteins data */
  /*----------------------------------------------*/
  _replication_report->set_new_genome_size(_genome->get_size());
  _replication_report->set_genome_nb_NC(_genome->get_nb_NC());
  _replication_report->set_genome_nb_E(_genome->get_nb_E());
  _replication_report->set_genome_nb_TF(_genome->get_nb_TF());
  _replication_report->set_genome_nb_BS(_genome->get_nb_BS());
  _replication_report->set_genome_nb_P(_genome->get_nb_P());
  _replication_report->set_genome_nb_inner_enzymes(_genome->get_nb_inner_enzymes());
  _replication_report->set_genome_nb_inflow_pumps(_genome->get_nb_inflow_pumps());
  _replication_report->set_genome_nb_outflow_pumps(_genome->get_nb_outflow_pumps());
  
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _replication_report->set_inherited_size(_inherited_proteins->get_size());
    _replication_report->set_inherited_nb_E(_inherited_proteins->get_nb_E());
    _replication_report->set_inherited_nb_TF(_inherited_proteins->get_nb_TF());
    _replication_report->set_inherited_nb_inner_enzymes(_inherited_proteins->get_nb_inner_enzymes());
    _replication_report->set_inherited_nb_inflow_pumps(_inherited_proteins->get_nb_inflow_pumps());
    _replication_report->set_inherited_nb_outflow_pumps(_inherited_proteins->get_nb_outflow_pumps());
  }
}

/**
 * \brief    Kill the cell (set alive to false)
 * \details  --
 * \param    size_t death_time
 * \return   \e void
 */
void Cell::kill( size_t death_time )
{
  /*----------------------------------*/
  /* 1) Kill the cell                 */
  /*----------------------------------*/
  assert(death_time >= _birth_time);
  _death_time = death_time;
  _lifespan   = _death_time-_birth_time;
  _active     = false;
  _alive      = false;
  if (_number_of_updates > 0)
  {
    _mean_score  /= _number_of_updates;
    _mean_energy /= _number_of_updates;
  }
  _prng = NULL;
  
  /*----------------------------------*/
  /* 2) Update the replication report */
  /*----------------------------------*/
  update_replication_report_data();
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Cell::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ main cell classes */
  
  _replication_report->save(backup_file);
  _genome->save(backup_file);
  if (_parameters->get_enzymatic_inheritance_scheme())
  {
    _inherited_proteins->save(backup_file);
  }
  _inherited_species_list->save(backup_file);
  _species_list->save(backup_file);
  
  /*------------------------------------------------------------------ main cell variables */
  
  gzwrite( backup_file, &_id,                  sizeof(_id) );
  gzwrite( backup_file, &_parent_id,           sizeof(_parent_id) );
  gzwrite( backup_file, &_generation,          sizeof(_generation) );
  gzwrite( backup_file, &_energy,              sizeof(_energy) );
  gzwrite( backup_file, &_active,              sizeof(_active) );
  gzwrite( backup_file, &_alive,               sizeof(_alive) );
  gzwrite( backup_file, &_x,                   sizeof(_x) );
  gzwrite( backup_file, &_y,                   sizeof(_y) );
  gzwrite( backup_file, &_score,               sizeof(_score) );
  gzwrite( backup_file, &_number_of_updates,   sizeof(_number_of_updates) );
  gzwrite( backup_file, &_number_of_divisions, sizeof(_number_of_divisions) );
  gzwrite( backup_file, &_tagged,              sizeof(_tagged) );
  
  /*------------------------------------------------------------------ mutation rates */
  
  save_mutation_rates(backup_file);
  
  /*------------------------------------------------------------------ time variables */
  
  gzwrite( backup_file, &_birth_time, sizeof(_birth_time) );
  gzwrite( backup_file, &_death_time, sizeof(_death_time) );
  gzwrite( backup_file, &_lifespan,   sizeof(_lifespan) );
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  gzwrite( backup_file, &_toxicity,                   sizeof(_toxicity) );
  gzwrite( backup_file, &_inherited_TF_amount,        sizeof(_inherited_TF_amount) );
  gzwrite( backup_file, &_inherited_E_amount,         sizeof(_inherited_E_amount) );
  gzwrite( backup_file, &_TF_amount,                  sizeof(_TF_amount) );
  gzwrite( backup_file, &_E_amount,                   sizeof(_E_amount) );
  gzwrite( backup_file, &_min_metabolic_amount,       sizeof(_min_metabolic_amount) );
  gzwrite( backup_file, &_max_metabolic_amount,       sizeof(_max_metabolic_amount) );
  gzwrite( backup_file, &_metabolic_uptake,           sizeof(_metabolic_uptake) );
  gzwrite( backup_file, &_metabolic_release,          sizeof(_metabolic_release) );
  gzwrite( backup_file, &_previous_metabolic_amount,  sizeof(_previous_metabolic_amount) );
  gzwrite( backup_file, &_min_energy,                 sizeof(_min_energy) );
  gzwrite( backup_file, &_mean_energy,                sizeof(_mean_energy) );
  gzwrite( backup_file, &_max_energy,                 sizeof(_max_energy) );
  gzwrite( backup_file, &_min_score,                  sizeof(_min_score) );
  gzwrite( backup_file, &_mean_score,                 sizeof(_mean_score) );
  gzwrite( backup_file, &_max_score,                  sizeof(_max_score) );
  gzwrite( backup_file, &_metabolic_growth_rate,      sizeof(_metabolic_growth_rate) );
  gzwrite( backup_file, &_diff_metabolic_growth_rate, sizeof(_diff_metabolic_growth_rate) );
  gzwrite( backup_file, &_grn_nb_nodes,               sizeof(_grn_nb_nodes) );
  gzwrite( backup_file, &_grn_nb_edges,               sizeof(_grn_nb_edges) );
  gzwrite( backup_file, &_metabolic_nb_nodes,         sizeof(_metabolic_nb_nodes) );
  gzwrite( backup_file, &_metabolic_nb_edges,         sizeof(_metabolic_nb_edges) );
  size_t N_inflowing_pumps  = _inflowing_pumps.size();
  size_t N_outflowing_pumps = _outflowing_pumps.size();
  gzwrite( backup_file, &N_inflowing_pumps, sizeof(N_inflowing_pumps) );
  gzwrite( backup_file, &N_outflowing_pumps, sizeof(N_outflowing_pumps) );
  for (size_t i = 0; i < N_inflowing_pumps; i++)
  {
    gzwrite( backup_file, &_inflowing_pumps[i], sizeof(_inflowing_pumps[i]) );
  }
  for (size_t i = 0; i < N_outflowing_pumps; i++)
  {
    gzwrite( backup_file, &_outflowing_pumps[i], sizeof(_outflowing_pumps[i]) );
  }
  gzwrite( backup_file, &_trophic_group, sizeof(_trophic_group) );
  gzwrite( backup_file, &_trophic_level, sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ cell color */
  
  gzwrite( backup_file, &_red_color,   sizeof(_red_color) );
  gzwrite( backup_file, &_green_color, sizeof(_green_color) );
  gzwrite( backup_file, &_blue_color,  sizeof(_blue_color) );
}

/**
 * \brief    Replace the cell's data
 * \details  Replace every structures, excepted the identifiers and coordinates
 * \param    Cell* cell
 * \return   \e void
 */
void Cell::replace_data( Cell* cell )
{
  /*------------------------------------------------------------------ main cell classes */
  
  /*--------------------------------*/
  /* Replace the replication report */
  /*--------------------------------*/
  /*
  delete _replication_report;
  _replication_report = NULL;
  _replication_report = new ReplicationReport(*cell->get_replication_report());
  */
  
  /*--------------------------------*/
  /* Replace the genome             */
  /*--------------------------------*/
  _genome->replace_data(cell->get_genome());
  
  /*--------------------------------*/
  /* Replace the inherited proteins */
  /*--------------------------------*/
  _inherited_proteins->replace_data(cell->get_inherited_proteins());
  
  /*--------------------------------*/
  /* Replace species lists          */
  /*--------------------------------*/
  delete _inherited_species_list;
  _inherited_species_list = NULL;
  _inherited_species_list = new SpeciesList(*cell->get_inherited_species_list());
  delete _species_list;
  _species_list = NULL;
  _species_list = new SpeciesList(*cell->get_species_list());
  
  _trophic_level = cell->get_trophic_level();
  
  /*------------------------------------------------------------------ main cell variables */
  
  /*
  _energy = cell->get_energy();
  _active = cell->get_active();
  _alive  = cell->get_alive();
  _score  = cell->get_score();
  */
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  /*
  _toxicity                   = cell->get_toxicity();
  _inherited_TF_amount        = cell->get_inherited_TF_amount();
  _inherited_E_amount         = cell->get_inherited_E_amount();
  _TF_amount                  = cell->get_TF_amount();
  _E_amount                   = cell->get_E_amount();
  _min_metabolic_amount       = cell->get_min_metabolic_amount();
  _max_metabolic_amount       = cell->get_max_metabolic_amount();
  _metabolic_uptake           = cell->get_metabolic_uptake();
  _metabolic_release          = cell->get_metabolic_release();
  _min_energy                 = cell->get_min_energy();
  _mean_energy                = cell->get_mean_energy();
  _max_energy                 = cell->get_max_energy();
  _min_score                  = cell->get_min_score();
  _mean_score                 = cell->get_mean_score();
  _max_score                  = cell->get_max_score();
  _metabolic_growth_rate      = cell->get_metabolic_growth_rate();
  _diff_metabolic_growth_rate = cell->get_Dmetabolic_growth_rate();
  _grn_nb_nodes               = cell->get_grn_nb_nodes();
  _grn_nb_edges               = cell->get_grn_nb_edges();
  _metabolic_nb_nodes         = cell->get_metabolic_nb_nodes();
  _metabolic_nb_edges         = cell->get_metabolic_nb_edges();
   */
}

/**
 * \brief    Write the genome structure in a file
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void Cell::write_genome( std::ofstream& filestream )
{
  filestream << "type functional id parent_id s p km kcat BStag CoEtag free_activity bound_activity window TFtag beta\n";
  
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    genetic_unit* p = _genome->get_genetic_unit(i);
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (p->type == NON_CODING)
    {
      filestream << "NON_CODING " << p->functional << " ";
    }
    else if (p->type == ENZYME)
    {
      filestream << "ENZYME " << p->functional << " ";
    }
    else if (p->type == TRANSCRIPTION_FACTOR)
    {
      filestream << "TRANSCRIPTION_FACTOR " << p->functional << " ";
    }
    else if (p->type == BINDING_SITE)
    {
      filestream << "BINDING_SITE " << p->functional << " ";
    }
    else if (p->type == PROMOTER)
    {
      filestream << "PROMOTER " << p->functional << " ";
    }
    filestream << p->identifier << " " << p->parent_identifier << " ";
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    filestream << p->s << " " << p->p << " " << p->kcat/p->kcat_km_ratio << " " << p->kcat << " ";
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    filestream << p->BS_tag << " " << p->coE_tag << " ";
    filestream << p->free_activity << " " << p->bound_activity << " " << p->binding_window << " ";
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    filestream << p->TF_tag << " ";
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    filestream << p->basal_expression_level << "\n";
  }
}

/**
 * \brief    Write inherited proteins in a file
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void Cell::write_inherited_proteins( std::ofstream& filestream )
{
  filestream << "type functional id parent_id s p km kcat BStag CoEtag free_activity bound_activity window TFtag beta\n";
  
  for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
  {
    genetic_unit* unit = _inherited_proteins->get_genetic_unit(i);
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (unit->type == ENZYME)
    {
      filestream << "ENZYME 1 ";
    }
    else if (unit->type == TRANSCRIPTION_FACTOR)
    {
      filestream << "TRANSCRIPTION_FACTOR 1 ";
    }
    filestream << unit->identifier << " " << unit->parent_identifier << " ";
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    filestream << unit->s << " " << unit->p << " " << unit->kcat/unit->kcat_km_ratio << " " << unit->kcat << " ";
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    filestream << unit->BS_tag << " " << unit->coE_tag << " ";
    filestream << unit->free_activity << " " << unit->bound_activity << " " << unit->binding_window << " ";
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    filestream << unit->TF_tag << " ";
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    filestream << unit->basal_expression_level << "\n";
  }
}

/**
 * \brief    Write cell's genetic regulation network
 * \details  --
 * \param    std::ofstream& nodeslist
 * \param    std::ofstream& edgeslist
 * \return   \e void
 */
void Cell::write_genetic_regulation_network( std::ofstream& nodeslist, std::ofstream& edgeslist )
{
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) declare the nodes map                         */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  std::unordered_map<size_t, grn_node_type> nodes_map;
  nodes_map.clear();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) write headers                                 */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  nodeslist << "id conc type s p km kcat inherited\n";
  edgeslist << "source target weight\n";
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) explore the reactions list and write edgelist */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  reaction_list*               list                   = _ode->get_reaction_list();
  size_t                       Ngenome                = list->Ngenome;
  size_t                       Ninherited             = list->Ninherited;
  std::vector<size_t>&         Nenhancer              = list->grn_Nenhancer;
  std::vector<size_t>&         Noperator              = list->grn_Noperator;
  std::vector<size_t>&         Ngenes                 = list->grn_Ngenes;
  std::vector<size_t>&         enhancer_TF_list       = list->grn_enhancer_TF_list;
  std::vector<double>&         enhancer_affinity_list = list->grn_enhancer_affinity_list;
  std::vector<int>&            enhancer_coe_list      = list->grn_enhancer_coe_list;
  std::vector<co_enzyme_type>& enhancer_coe_type      = list->grn_enhancer_coe_type;
  std::vector<size_t>&         operator_TF_list       = list->grn_operator_TF_list;
  std::vector<double>&         operator_affinity_list = list->grn_operator_affinity_list;
  std::vector<int>&            operator_coe_list      = list->grn_operator_coe_list;
  std::vector<co_enzyme_type>& operator_coe_type      = list->grn_operator_coe_type;
  std::vector<size_t>&         regulated_genes        = list->grn_regulated_genes;
  
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < list->grn_N; region++)
  {
    /*---------------------------------------------*/
    /* 1) For each TF binding to the enhancer site */
    /*---------------------------------------------*/
    for (size_t pos = 0; pos < Nenhancer[region]; pos++)
    {
      size_t yTF          = enhancer_TF_list[enhancer_TF_index];
      size_t yCoE         = enhancer_coe_list[enhancer_TF_index]+Ngenome+Ninherited-1;
      co_enzyme_type tCoE = enhancer_coe_type[enhancer_TF_index];
      
      /* 1.1) Add TF */
      if (nodes_map.find(yTF) == nodes_map.end())
      {
        nodes_map[yTF] = GRN_TF;
      }
      /* 1.2) Add CoE */
      if (nodes_map.find(yCoE) == nodes_map.end())
      {
        nodes_map[yCoE] = GRN_COE;
      }
      /* 1.3) Add link */
      edgeslist << yCoE << " " << yTF << " " << tCoE << "\n";
      
      /* 1.4) Explore regulated genes */
      for (size_t gene_index = 0; gene_index < Ngenes[region]; gene_index++)
      {
        size_t yRGene           = regulated_genes[regulated_genes_index+gene_index];
        genetic_unit_type TGene = _genome->get_genetic_unit(regulated_genes[regulated_genes_index+gene_index])->type;
        /* Add regulated gene */
        if (nodes_map.find(yRGene) == nodes_map.end())
        {
          if (TGene == ENZYME)
          {
            nodes_map[yRGene] = GRN_E;
          }
          else if (TGene == TRANSCRIPTION_FACTOR)
          {
            nodes_map[yRGene] = GRN_TF;
          }
        }
        /* Add link */
        edgeslist << yTF << " " << yRGene << " " << enhancer_affinity_list[enhancer_TF_index] << "\n";
      }
      /* increment enhancer TF index */
      enhancer_TF_index++;
    }
    
    /*---------------------------------------------*/
    /* 2) For each TF binding to the operator site */
    /*---------------------------------------------*/
    for (size_t pos = 0; pos < Noperator[region]; pos++)
    {
      size_t yTF          = operator_TF_list[operator_TF_index];
      size_t yCoE         = operator_coe_list[operator_TF_index]+Ngenome+Ninherited-1;
      co_enzyme_type tCoE = operator_coe_type[operator_TF_index];
      
      /* 2.1) Add TF */
      if (nodes_map.find(yTF) == nodes_map.end())
      {
        nodes_map[yTF] = GRN_TF;
      }
      /* 2.2) Add CoE */
      if (nodes_map.find(yCoE) == nodes_map.end())
      {
        nodes_map[yCoE] = GRN_COE;
      }
      /* 2.3) Add link */
      edgeslist << yCoE << " " << yTF << " " << tCoE << "\n";
      
      /* 2.4) Explore regulated genes */
      for (size_t gene_index = 0; gene_index < Ngenes[region]; gene_index++)
      {
        size_t yRGene           = regulated_genes[regulated_genes_index+gene_index];
        genetic_unit_type TGene = _genome->get_genetic_unit(regulated_genes[regulated_genes_index+gene_index])->type;
        /* Add regulated gene */
        if (nodes_map.find(yRGene) == nodes_map.end())
        {
          if (TGene == ENZYME)
          {
            nodes_map[yRGene] = GRN_E;
          }
          else if (TGene == TRANSCRIPTION_FACTOR)
          {
            nodes_map[yRGene] = GRN_TF;
          }
        }
        /* Add link */
        edgeslist << yTF << " " << yRGene << " " << -operator_affinity_list[operator_TF_index] << "\n";
      }
      /* Increment operator TF index */
      operator_TF_index++;
    }
    
    /*---------------------------------------------*/
    /* 3) Update regulated genes index             */
    /*---------------------------------------------*/
    regulated_genes_index += Ngenes[region];
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) write the nodelist                            */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  for (std::unordered_map<size_t, grn_node_type>::iterator it = nodes_map.begin(); it != nodes_map.end(); ++it)
  {
    /*------------------------------------*/
    /* 4.1) First case: the node is a CoE */
    /*------------------------------------*/
    if (it->second == GRN_COE)
    {
      nodeslist << it->first << " " << _species_list->get((int)it->first+1) << " COE " << it->first-Ngenome-Ninherited << " 0 0 0 0\n";
    }
    /*------------------------------------*/
    /* 4.2) Second case: the node is a E  */
    /*------------------------------------*/
    else if (it->second == GRN_E)
    {
      genetic_unit* unit = _genome->get_genetic_unit(it->first);
      nodeslist << it->first << " " << _genome->get_concentration_vector()[it->first] << " E " << unit->s << " " << unit->p << " " << unit->kcat/unit->kcat_km_ratio << " " << unit->kcat << " 0\n";
    }
    /*------------------------------------*/
    /* 4.3) Third case: the node is a TF  */
    /*------------------------------------*/
    else if (it->second == GRN_TF)
    {
      bool inherited = false;
      if (it->first >= Ngenome)
      {
        inherited = true;
      }
      if (!inherited)
      {
        nodeslist << it->first << " " << _genome->get_concentration_vector()[it->first] << " TF 0 0 0 0 0\n";
      }
      else
      {
        nodeslist << it->first << " " << _genome->get_concentration_vector()[it->first] << " TF 0 0 0 0 1\n";
      }
    }
  }
}

/**
 * \brief    Write cell's metabolic network
 * \details  --
 * \param    std::ofstream& nodeslist
 * \param    std::ofstream& edgeslist
 * \return   \e void
 */
void Cell::write_metabolic_network( std::ofstream& nodeslist, std::ofstream& edgeslist )
{
  /****************************************************/
  /* 1) Declare the nodes map                         */
  /****************************************************/
  std::unordered_map<int, unsigned char> nodes_map;
  nodes_map.clear();
  
  /****************************************************/
  /* 2) Write headers                                 */
  /****************************************************/
  edgeslist << "s p km kcat delta_g e flux\n";
  nodeslist << "tag conc\n";
  
  /****************************************************/
  /* 3) Explore the reactions list and write edgelist */
  /****************************************************/
  reaction_list* list = _ode->get_reaction_list();
  
  for (size_t i = 0; i < list->metabolic_N; i++)
  {
    /*--------------------------------------------------------*/
    /* 3.1) Get the reaction                                  */
    /*--------------------------------------------------------*/
    reaction_type type = list->metabolic_type[i];
    int s              = list->metabolic_s[i];
    int p              = list->metabolic_p[i];
    double km          = list->metabolic_km[i];
    double kcat        = list->metabolic_kcat[i];
    double delta_g     = list->metabolic_delta_g[i];
    double e           = _genome->get_concentration_vector()[list->metabolic_e[i]];
    
    /*--------------------------------------------------------*/
    /* 3.2) If reaction is metabolic inflowing pump activity  */
    /*--------------------------------------------------------*/
    if (type == INFLOWING_PUMP_ACTIVITY)
    {
      if (nodes_map.find(-s) == nodes_map.end())
      {
        nodes_map[-s] = 1;
      }
      if (nodes_map.find(p) == nodes_map.end())
      {
        nodes_map[p] = 1;
      }
      edgeslist << -s << " " << p << " " << km << " " << kcat << " " << delta_g << " " << e << " " << "1.0" << "\n";
    }
    /*--------------------------------------------------------*/
    /* 3.3) If reaction is metabolic outflowing pump activity */
    /*--------------------------------------------------------*/
    else if (type == OUTFLOWING_PUMP_ACTIVITY)
    {
      if (nodes_map.find(s) == nodes_map.end())
      {
        nodes_map[s] = 1;
      }
      if (nodes_map.find(-p) == nodes_map.end())
      {
        nodes_map[-p] = 1;
      }
      edgeslist << s << " " << -p << " " << km << " " << kcat << " " << delta_g << " " << e << " " << "1.0" << "\n";
    }
    /*--------------------------------------------------------*/
    /* 3.4) If reaction is metabolic catalytic activity       */
    /*--------------------------------------------------------*/
    else if (type == CATALYTIC_CONSUMING_ACTIVITY || type == CATALYTIC_REWARDING_ACTIVITY)
    {
      if (nodes_map.find(s) == nodes_map.end())
      {
        nodes_map[s] = 1;
      }
      if (nodes_map.find(p) == nodes_map.end())
      {
        nodes_map[p] = 1;
      }
      edgeslist << s << " " << p << " " << km << " " << kcat << " " << delta_g << " " << e << " " << "1.0" << "\n";
    }
  }
  /****************************************************/
  /* 4) Write node list with concentrations           */
  /****************************************************/
  for (std::unordered_map<int, unsigned char>::iterator it = nodes_map.begin(); it != nodes_map.end(); ++it)
  {
    if (it->first > 0)
    {
      nodeslist << it->first << " " << _species_list->get(it->first) << "\n";
    }
  }
}

/**
 * \brief    Write metabolic amounts
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void Cell::write_metabolic_amounts( std::ofstream& filestream )
{
  assert(_species_list->get_size() == _inherited_species_list->get_size());
  size_t size = _species_list->get_size();
  for (size_t i = 0; i < size; i++)
  {
    if (_species_list->get_size() > i)
    {
      filestream << i+1 << " " << _species_list->get((int)i+1) << " ";
    }
    else
    {
      filestream << "0.0 ";
    }
    if (_inherited_species_list->get_size() > i)
    {
      filestream << _inherited_species_list->get((int)i+1) << "\n";
    }
    else
    {
      filestream << "0.0\n";
    }
  }
}

/**
 * \brief    DEPRECATED METHOD - Write the current ODE system
 * \details  The reaction_list and the X vector must be built
 * \param    void
 * \return   \e void
 */
/*
void Cell::write_current_ODE_system( void )
{
  std::cout << "----------------------------------------------------------\n";
  std::cout << "WARNING !!! Cell::write_current_ODE_system() is deprecated\n";
  std::cout << "----------------------------------------------------------\n";
  std::ofstream file("ODEsystem.txt", std::ios::out | std::ios::trunc);
  file << "----------------------------------------------------------\n";
  file << "WARNING !!! Cell::write_current_ODE_system() is deprecated\n";
  file << "----------------------------------------------------------\n\n";
  double* y           = _ode->get_X();
  reaction_list* list = _ode->get_reaction_list();
  
  //-----------------------------------
  // A) Write reactions
  //-----------------------------------
  
  size_t N          = list->N;
  size_t Ngenome    = list->Ngenome;
  size_t Ninherited = list->Ninherited;
  size_t Ncell      = list->Ncell;
  
  file << "Current cell state for cell " << _id << " (generation=" << _generation << ", birth time=" << _birth_time << ", nb updates=" << _number_of_updates << ")\n";
  file << "\n";
  file << "Genome size = " << Ngenome << "\n";
  file << "Nb inherited proteins = " << Ninherited << "\n";
  file << "Cell metabolite vector length = " << Ncell << " (=environment vector length)\n";
  file << "State vector length = " << N << " (=" << Ngenome << "+" << Ninherited << "+2*" << Ncell << ")\n\n";
  
  //-----------------------------------
  // B) Write genome vector
  //-----------------------------------
  
  file << "GENOME STRUCTURE:\n";
  file << "-----------------\n";
  file << "For enzyme units: if kcat < 0.0, the reaction is reverted (e.g. 4->4 with kcat > 0.0 is an uptake pump ; the same reaction with kcat < 0.0 is a release pump).\n";
  file << "\n";
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    genetic_unit* unit = _genome->get_genetic_unit(i);
    if (unit->type == NON_CODING)
    {
      file << "y[" << i << "] (NONCODING)\n";
    }
    else if (unit->type == ENZYME)
    {
      file << "y[" << i << "] (ENZYME) s=" << unit->s << ", p=" << unit->p << ", km=" << unit->km << ", kcat=" << unit->kcat << "\n";
    }
    else if (unit->type == TRANSCRIPTION_FACTOR)
    {
      file << "y[" << i << "] (TRANSCRIPTION_FACTOR) BS_tag=" << unit->BS_tag << ", coE=" << unit->coE_tag << ", bound_activity=" << unit->bound_activity << ", free_activity=" << unit->free_activity << "\n";
    }
    else if (unit->type == BINDING_SITE)
    {
      file << "y[" << i << "] (BINDING_SITE) TF_tag=" << unit->TF_tag << "\n";
    }
    else if (unit->type == PROMOTER)
    {
      file << "y[" << i << "] (PROMOTER) basal_expression_level=" << unit->basal_expression_level << "\n";
    }
  }
  file << "\n";
  
  //-----------------------------------
  // C) Write inherited vector
  //-----------------------------------
  
  file << "INHERITED PROTEINS:\n";
  file << "-------------------\n";
  file << "For enzyme units: if kcat < 0.0, the reaction is reverted (e.g. 4->4 with kcat > 0.0 is an uptake pump ; the same reaction with kcat < 0.0 is a release pump).\n";
  file << "\n";
  for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
  {
    genetic_unit* unit = _inherited_proteins->get_genetic_unit(i);
    if (unit->type == ENZYME)
    {
      file << "y[" << i << "] (ENZYME) s=" << unit->s << ", p=" << unit->p << ", km=" << unit->km << ", kcat=" << unit->kcat << "\n";
    }
    else if (unit->type == TRANSCRIPTION_FACTOR)
    {
      file << "y[" << i << "] (TRANSCRIPTION_FACTOR) BS_tag=" << unit->BS_tag << ", coE=" << unit->coE_tag << ", bound_activity=" << unit->bound_activity << ", free_activity=" << unit->free_activity << "\n";
    }
  }
  file << "\n";
  
  //-----------------------------------
  // D) Write the current state vector
  //-----------------------------------
  
  file << "STATE VECTOR:\n";
  file << "-------------\n";
  file << "\n";
  file << "double y[" << list->N << "] = {";
  for (size_t i = 0; i < list->N; i++)
  {
    if (i < list->N-1)
    {
      file << y[i] << ", ";
    }
    else
    {
      file << y[i] << "};\n\n";
    }
  }
  file << "\n";
  file << "y = [";
  for (size_t i = 0; i < list->N; i++)
  {
    if (i < list->N-1)
    {
      file << y[i] << ", ";
    }
    else
    {
      file << y[i] << "]\n\n";
    }
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // WRITE GENETIC REGULATION REACTIONS
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  file << "ODE SYSTEM:\n";
  file << "-----------\n";
  file << "\n";
  std::vector<size_t>&         Nenhancer                 = list->grn_Nenhancer;
  std::vector<size_t>&         Noperator                 = list->grn_Noperator;
  std::vector<size_t>&         Ngenes                    = list->grn_Ngenes;
  std::vector<size_t>&         enhancer_TF_list          = list->grn_enhancer_TF_list;
  std::vector<double>&         enhancer_affinity_list    = list->grn_enhancer_affinity_list;
  std::vector<int>&            enhancer_coe_list         = list->grn_enhancer_coe_list;
  std::vector<co_enzyme_type>& enhancer_coe_type         = list->grn_enhancer_coe_type;
  std::vector<double>&         enhancer_coe_km           = list->grn_enhancer_coe_km;
  std::vector<size_t>&         operator_TF_list          = list->grn_operator_TF_list;
  std::vector<double>&         operator_affinity_list    = list->grn_operator_affinity_list;
  std::vector<int>&            operator_coe_list         = list->grn_operator_coe_list;
  std::vector<co_enzyme_type>& operator_coe_type         = list->grn_operator_coe_type;
  std::vector<double>&         operator_coe_km           = list->grn_operator_coe_km;
  std::vector<size_t>&         regulated_genes           = list->grn_regulated_genes;
  std::vector<double>&         beta                      = list->grn_beta;
  double                       hill_n                    = list->grn_hill_n;
  double                       theta_n                   = list->grn_theta_n;
  double                       protein_degradation_rate  = list->grn_protein_degradation_rate;
  double                       timestep_ratio            = list->timestep_ratio;
  bool                         coe_activity              = list->co_enzyme_activity;
  
  size_t enhancer_TF_index     = 0;
  size_t operator_TF_index     = 0;
  size_t regulated_genes_index = 0;
  for (size_t region = 0; region < list->grn_N; region++)
  {
    //------------------------------------------------------------------
    // Write Ei(t) and Oi(t) equations
    //------------------------------------------------------------------
    
    // A) Write Ei(t) equation -> enhancer contribution to the transcriptional activity
    std::stringstream ei_stream;
    std::stringstream coe_stream;
    
    // For each TF binding to the enhancer site
    for (size_t pos = 0; pos < Nenhancer[region]; pos++)
    {
      coe_stream.str("");
      
      // if co-enzyme activity is activated, compute its contribution
      if (coe_activity)
      {
        //-------------------------
        // If the co-enzyme exists
        //-------------------------
        if (enhancer_coe_list[enhancer_TF_index] <= (int)Ncell)
        {
          size_t coE_index = enhancer_coe_list[enhancer_TF_index]+Ngenome+Ninherited-1;
          
          if (enhancer_coe_type[enhancer_TF_index] == REPRESSOR)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "]-y[" << coE_index << "] > 0.0 ? y[" << enhancer_TF_list[enhancer_TF_index] << "]-y[" << coE_index << "] : 0.0)";
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ACTIVATOR)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "] < y[" << coE_index << "] ? y[" << enhancer_TF_list[enhancer_TF_index] << "] : y[" << coE_index << "])";
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ALWAYS_REPRESSED)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "]*0.0)";
          }
          else
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "])";
          }
        }
        // Else the co-enzyme cannot exist
        else
        {
          if (enhancer_coe_type[enhancer_TF_index] == REPRESSOR)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "] > 0.0 ? y[" << enhancer_TF_list[enhancer_TF_index] << "] : 0.0)";
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ACTIVATOR)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "] < 0.0 ? y[" << enhancer_TF_list[enhancer_TF_index] << "] : 0.0)";
          }
          else if (enhancer_coe_type[enhancer_TF_index] == ALWAYS_REPRESSED)
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "]*0.0)";
          }
          else
          {
            coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "])";
          }
        }
      }
      // Else if co-enzyme activity is not activated
      else
      {
        coe_stream << "(y[" << enhancer_TF_list[enhancer_TF_index] << "])";
      }
      
      // Then write the local Ei equation
      if (Nenhancer[region] > 1)
      {
        if (pos == 0)
        {
          ei_stream << "(" << coe_stream.str() << "*" << enhancer_affinity_list[enhancer_TF_index] << " + ";
        }
        else if (pos > 0 && pos < Nenhancer[region]-1)
        {
          ei_stream << coe_stream.str() << "*" << enhancer_affinity_list[enhancer_TF_index] << " + ";
        }
        else
        {
          ei_stream << coe_stream.str() << "*" << enhancer_affinity_list[enhancer_TF_index] << ")";
        }
      }
      else if (Nenhancer[region] == 1)
      {
        ei_stream << "(" << coe_stream.str() << "*" << enhancer_affinity_list[enhancer_TF_index] << ")";
      }
      
      // Increment enhancer TF index
      enhancer_TF_index++;
    }
    if (Nenhancer[region] == 0)
    {
      ei_stream << "0.0";
    }
    
    // B) Compute Oi(t) -> operator contribution to the transcriptional activity
    std::stringstream oi_stream;
    
    // For each TF binding to the operator site
    for (size_t pos = 0; pos < Noperator[region]; pos++)
    {
      coe_stream.str("");
      
      // if co-enzyme activity is activated, compute its contribution
      if (coe_activity)
      {
        // If the co-enzyme exists
        if (operator_coe_list[operator_TF_index] <= (int)Ncell)
        {
          size_t coE_index = operator_coe_list[operator_TF_index]+Ngenome+Ninherited-1;
          
          if (operator_coe_type[operator_TF_index] == REPRESSOR)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "]-y[" << coE_index << "] > 0.0 ? y[" << operator_TF_list[operator_TF_index] << "]-y[" << coE_index << "] : 0.0)";
          }
          else if (operator_coe_type[operator_TF_index] == ACTIVATOR)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "] < y[" << coE_index << "] ? y[" << operator_TF_list[operator_TF_index] << "] : y[" << coE_index << "])";
          }
          else if (operator_coe_type[operator_TF_index] == ALWAYS_REPRESSED)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "]*0.0)";
          }
          else
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "])";
          }
        }
        // Else the co-enzyme cannot exist
        else
        {
          if (operator_coe_type[operator_TF_index] == REPRESSOR)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "] > 0.0 ? y[" << operator_TF_list[operator_TF_index] << "] : 0.0)";
          }
          else if (operator_coe_type[operator_TF_index] == ACTIVATOR)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "] < 0.0 ? y[" << operator_TF_list[operator_TF_index] << "] : 0.0)";
          }
          else if (operator_coe_type[operator_TF_index] == ALWAYS_REPRESSED)
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "]*0.0)";
          }
          else
          {
            coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "])";
          }
        }
      }
      // Else if co-enzyme activity is not activated
      else
      {
        coe_stream << "(y[" << operator_TF_list[operator_TF_index] << "])";
      }
      
      // Then write the local Oi equation
      if (Noperator[region] > 1)
      {
        if (pos == 0)
        {
          oi_stream << "(" << coe_stream.str() << "*" << operator_affinity_list[operator_TF_index] << " + ";
        }
        else if (pos > 0 && pos < Noperator[region]-1)
        {
          oi_stream << coe_stream.str() << "*" << operator_affinity_list[operator_TF_index] << " + ";
        }
        else
        {
          oi_stream << coe_stream.str() << "*" << operator_affinity_list[operator_TF_index] << ")";
        }
      }
      else if (Noperator[region] == 1)
      {
        oi_stream << "(" << coe_stream.str() << "*" << operator_affinity_list[operator_TF_index] << ")";
      }
      
      // Increment operator TF index
      operator_TF_index++;
    }
    if (Noperator[region] == 0)
    {
      oi_stream << "0.0";
    }
    
    //------------------------------------------------------------------
    // Write the transcription rate ei(t) equation
    //------------------------------------------------------------------
    std::stringstream express_stream;
    if (beta[region] > 0.0)
    {
      express_stream << "( " << beta[region] << " * " << theta_n << " / (pow(" << oi_stream.str() << ", " << hill_n << ")+" << theta_n << ") * (1.0 + (1.0/" << beta[region] << " - 1.0) * (pow(" << ei_stream.str() << ", " << hill_n << ")/(pow(" << ei_stream.str() << ", " << hill_n << ")+" << theta_n << "))) )";
    }
    else
    {
      express_stream << "0.0";
    }
    
    //------------------------------------------------------------------
    // For each regulated gene, write the equation ei(t)*timestep_ratio
    //------------------------------------------------------------------
    for (size_t gene_index = 0; gene_index < Ngenes[region]; gene_index++)
    {
      file << "dydt[" << regulated_genes[regulated_genes_index] << "] = (" << express_stream.str() << "-y[" << regulated_genes[regulated_genes_index] << "]*" << protein_degradation_rate << ") * " << timestep_ratio << ";\n";
      regulated_genes_index++;
    }
  }
  file << "\n";
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // INHERITED PROTEINS DEGRADATION
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  for (size_t i = 0; i < Ninherited; i++)
  {
    size_t Cindex = Ngenome+i;
    file << "dydt[" << Cindex << "] -= y[" << Cindex << "]*" << protein_degradation_rate << "*" << timestep_ratio << ";\n";
  }
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // APPLY METABOLIC REACTIONS
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  std::vector<reaction_type>& type = list->metabolic_type;
  std::vector<int>&           s    = list->metabolic_s;
  std::vector<int>&           p    = list->metabolic_p;
  std::vector<double>&        km   = list->metabolic_km;
  std::vector<double>&        kcat = list->metabolic_kcat;
  std::vector<size_t>&        e    = list->metabolic_e;
  for (size_t i = 0; i < list->metabolic_N; i++)
  {
    //-----------------------------------------------------
    // 4.2.1) if reaction is inflowing pump activity
    //-----------------------------------------------------
    if (type[i] == INFLOWING_PUMP_ACTIVITY)
    {
      size_t Cindex1 = Ngenome+Ninherited+Ncell+s[i]-1;
      size_t Cindex2 = Ngenome+Ninherited+p[i]-1;
      file << "\n";
      file << "dydt[" << Cindex1 << "] -= (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
      file << "dydt[" << Cindex2 << "] += (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
    }
    
    //-----------------------------------------------------
    // 4.2.2) if reaction is outflowing pump activity
    //-----------------------------------------------------
    else if (type[i] == OUTFLOWING_PUMP_ACTIVITY)
    {
      size_t Cindex1 = Ngenome+Ninherited+s[i]-1;
      size_t Cindex2 = Ngenome+Ninherited+Ncell+p[i]-1;
      file << "\n";
      file << "dydt[" << Cindex1 << "] -= (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
      file << "dydt[" << Cindex2 << "] += (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
    }
    
    //-----------------------------------------------------
    // 4.2.3) if reaction is inner cell catalytic activity
    //-----------------------------------------------------
    else if (type[i] == CATALYTIC_CONSUMING_ACTIVITY || type[i] == CATALYTIC_REWARDING_ACTIVITY)
    {
      size_t Cindex1 = Ngenome+Ninherited+s[i]-1;
      size_t Cindex2 = Ngenome+Ninherited+p[i]-1;
      file << "\n";
      file << "dydt[" << Cindex1 << "] -= (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
      file << "dydt[" << Cindex2 << "] += (" << kcat[i] << "*" << "y[" << Cindex1 << "]*y[" << e[i] << "])/(" << km[i] << "+" << "y[" << Cindex1 << "]);\n";
    }
  }
  
  file.close();
}
*/

/**
 * \brief    Write the current cell's state header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void Cell::write_state_header( std::ofstream& filestream )
{
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    filestream << "g" << i << " ";
  }
  for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
  {
    filestream << "i" << i << " ";
  }
  for (size_t i = 0; i < _species_list->get_size(); i++)
  {
    filestream << "in" << (i+1) << " ";
  }
  for (size_t i = 0; i < _species_list->get_size(); i++)
  {
    filestream << "out" << (i+1) << " ";
  }
  filestream << "energy score t\n";
}

/**
 * \brief    Write the current cell's state
 * \details  --
 * \param    std::ofstream& filestream
 * \param    size_t time
 * \param    Environment* environment
 * \return   \e void
 */
void Cell::write_current_state( std::ofstream& filestream, size_t time, Environment* environment )
{
  for (size_t i = 0; i < _genome->get_size(); i++)
  {
    filestream << _genome->get_concentration_vector()[i] << " ";
  }
  for (size_t i = 0; i < _inherited_proteins->get_size(); i++)
  {
    filestream << _inherited_proteins->get_concentration_vector()[i] << " ";
  }
  for (size_t i = 0; i < _species_list->get_size(); i++)
  {
    filestream << _species_list->get((int)i+1) << " ";
  }
  for (size_t i = 0; i < environment->get_species_list(_x, _y)->get_size(); i++)
  {
    filestream << environment->get_species_list(_x, _y)->get((int)i+1) << " ";
  }
  filestream << _energy << " " << _score << " " << time << "\n";
}

/**
 * \brief    Get pumped and produced metabolites
 * \details  --
 * \param    std::vector<int>& pumped
 * \param    std::vector<int>& produced
 * \param    std::vector<int>& pm_pumped
 * \param    std::vector<int>& pm_produced
 * \return   \e void
 */
void Cell::get_pumped_and_produced_metabolites( std::vector<int>& pumped, std::vector<int>& produced, std::vector<int>& pm_pumped, std::vector<int>& pm_produced )
{
  pumped.clear();
  produced.clear();
  pm_pumped.clear();
  pm_produced.clear();
  
  /*********************************/
  /* 1) Declare the nodes map      */
  /*********************************/
  std::unordered_map<int, unsigned char> pumped_map;
  pumped_map.clear();
  std::unordered_map<int, unsigned char> produced_map;
  produced_map.clear();
  
  /*********************************/
  /* 2) Explore the reactions list */
  /*********************************/
  reaction_list* list = _ode->get_reaction_list();
  
  for (size_t i = 0; i < list->metabolic_N; i++)
  {
    /*-------------------------------------------------------*/
    /* 3.1) Get the reaction                                 */
    /*-------------------------------------------------------*/
    reaction_type type = list->metabolic_type[i];
    int           p    = list->metabolic_p[i];
    //double      conc = _species_list->get(p);
    
    /*-------------------------------------------------------*/
    /* 3.2) If reaction is metabolic inflowing pump activity */
    /*-------------------------------------------------------*/
    if (type == INFLOWING_PUMP_ACTIVITY)
    {
      pumped_map[abs(p)] = 1;
    }
    
    /*-------------------------------------------------------*/
    /* 3.3) If reaction is metabolic catalytic activity      */
    /*-------------------------------------------------------*/
    else if (type == CATALYTIC_CONSUMING_ACTIVITY || type == CATALYTIC_REWARDING_ACTIVITY)
    {
      produced_map[abs(p)] = 1;
    }
  }
  
  /*********************************/
  /* 3) Save the metabolites       */
  /*********************************/
  for (std::unordered_map<int, unsigned char>::iterator it = pumped_map.begin(); it != pumped_map.end(); ++it)
  {
    pumped.push_back(it->first);
  }
  for (std::unordered_map<int, unsigned char>::iterator it = produced_map.begin(); it != produced_map.end(); ++it)
  {
    produced.push_back(it->first);
  }
  pumped_map.clear();
  produced_map.clear();
  
  /*********************************/
  /* 4) Save prime numbers         */
  /*********************************/
  int* prime_numbers = _parameters->get_prime_numbers();
  for (size_t i = 0; i < pumped.size(); i++)
  {
    if (prime_numbers[pumped[i]-1] == 1)
    {
      pm_pumped.push_back(pumped[i]);
    }
  }
  for (size_t i = 0; i < produced.size(); i++)
  {
    if (prime_numbers[produced[i]-1] == 1)
    {
      pm_produced.push_back(produced[i]);
    }
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Initialize the mutation rates vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void Cell::initialize_mutation_rates( void )
{
  _mutation_rates = new double[NUMBER_OF_MUTATION_RATES];
  _mutation_rates[POINT_MUTATION_RATE]                    = _parameters->get_point_mutation_rate();
  _mutation_rates[DUPLICATION_RATE]                       = _parameters->get_duplication_rate();
  _mutation_rates[DELETION_RATE]                          = _parameters->get_deletion_rate();
  _mutation_rates[TRANSLOCATION_RATE]                     = _parameters->get_translocation_rate();
  _mutation_rates[INVERSION_RATE]                         = _parameters->get_inversion_rate();
  _mutation_rates[TRANSITION_RATE]                        = _parameters->get_transition_rate();
  _mutation_rates[BREAKPOINT_RATE]                        = _parameters->get_breakpoint_rate();
  _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]            = _parameters->get_substrate_tag_mutation_size();
  _mutation_rates[PRODUCT_TAG_MUTATION_SIZE]              = _parameters->get_product_tag_mutation_size();
  _mutation_rates[KCAT_MUTATION_SIZE]                     = _parameters->get_kcat_mutation_size();
  _mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]            = _parameters->get_kcat_km_ratio_mutation_size();
  _mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]         = _parameters->get_binding_site_tag_mutation_size();
  _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]            = _parameters->get_co_enzyme_tag_mutation_size();
  _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE] = _parameters->get_transcription_factor_tag_mutation_size();
  _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]   = _parameters->get_basal_expression_level_mutation_size();
}

/**
 * \brief    Save the mutation rates vector in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Cell::save_mutation_rates( gzFile backup_file )
{
  gzwrite( backup_file, &_mutation_rates[POINT_MUTATION_RATE],                    sizeof(_mutation_rates[POINT_MUTATION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[DUPLICATION_RATE],                       sizeof(_mutation_rates[DUPLICATION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[DELETION_RATE],                          sizeof(_mutation_rates[DELETION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[TRANSLOCATION_RATE],                     sizeof(_mutation_rates[TRANSLOCATION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[INVERSION_RATE],                         sizeof(_mutation_rates[INVERSION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[TRANSITION_RATE],                        sizeof(_mutation_rates[TRANSITION_RATE]) );
  gzwrite( backup_file, &_mutation_rates[BREAKPOINT_RATE],                        sizeof(_mutation_rates[BREAKPOINT_RATE]) );
  gzwrite( backup_file, &_mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE],            sizeof(_mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[PRODUCT_TAG_MUTATION_SIZE],              sizeof(_mutation_rates[PRODUCT_TAG_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[KCAT_MUTATION_SIZE],                     sizeof(_mutation_rates[KCAT_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE],            sizeof(_mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE],         sizeof(_mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE],            sizeof(_mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], sizeof(_mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]) );
  gzwrite( backup_file, &_mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE],   sizeof(_mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]) );
}

/**
 * \brief    Load the mutation rates vector from backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Cell::load_mutation_rates( gzFile backup_file )
{
  _mutation_rates = new double[NUMBER_OF_MUTATION_RATES];
  gzread( backup_file, &_mutation_rates[POINT_MUTATION_RATE],                    sizeof(_mutation_rates[POINT_MUTATION_RATE]) );
  gzread( backup_file, &_mutation_rates[DUPLICATION_RATE],                       sizeof(_mutation_rates[DUPLICATION_RATE]) );
  gzread( backup_file, &_mutation_rates[DELETION_RATE],                          sizeof(_mutation_rates[DELETION_RATE]) );
  gzread( backup_file, &_mutation_rates[TRANSLOCATION_RATE],                     sizeof(_mutation_rates[TRANSLOCATION_RATE]) );
  gzread( backup_file, &_mutation_rates[INVERSION_RATE],                         sizeof(_mutation_rates[INVERSION_RATE]) );
  gzread( backup_file, &_mutation_rates[TRANSITION_RATE],                        sizeof(_mutation_rates[TRANSITION_RATE]) );
  gzread( backup_file, &_mutation_rates[BREAKPOINT_RATE],                        sizeof(_mutation_rates[BREAKPOINT_RATE]) );
  gzread( backup_file, &_mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE],            sizeof(_mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[PRODUCT_TAG_MUTATION_SIZE],              sizeof(_mutation_rates[PRODUCT_TAG_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[KCAT_MUTATION_SIZE],                     sizeof(_mutation_rates[KCAT_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE],            sizeof(_mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE],         sizeof(_mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE],            sizeof(_mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], sizeof(_mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]) );
  gzread( backup_file, &_mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE],   sizeof(_mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]) );
}

/**
 * \brief    Copy a mutation rates vector
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Cell::copy_mutation_rates( double* mutation_rates )
{
  _mutation_rates = new double[NUMBER_OF_MUTATION_RATES];
  _mutation_rates[POINT_MUTATION_RATE]                    = mutation_rates[POINT_MUTATION_RATE];
  _mutation_rates[DUPLICATION_RATE]                       = mutation_rates[DUPLICATION_RATE];
  _mutation_rates[DELETION_RATE]                          = mutation_rates[DELETION_RATE];
  _mutation_rates[TRANSLOCATION_RATE]                     = mutation_rates[TRANSLOCATION_RATE];
  _mutation_rates[INVERSION_RATE]                         = mutation_rates[INVERSION_RATE];
  _mutation_rates[TRANSITION_RATE]                        = mutation_rates[TRANSITION_RATE];
  _mutation_rates[BREAKPOINT_RATE]                        = mutation_rates[BREAKPOINT_RATE];
  _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]            = mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE];
  _mutation_rates[PRODUCT_TAG_MUTATION_SIZE]              = mutation_rates[PRODUCT_TAG_MUTATION_SIZE];
  _mutation_rates[KCAT_MUTATION_SIZE]                     = mutation_rates[KCAT_MUTATION_SIZE];
  _mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE]            = mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE];
  _mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE]         = mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE];
  _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]            = mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE];
  _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE] = mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE];
  _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]   = mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE];
}

