
/**
 * \file      Simulation.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Simulation class definition
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

#include "Simulation.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
Simulation::Simulation( Parameters* parameters )
{
  _parameters        = parameters;
  _population        = new Population(_parameters);
  _environment       = new Environment(_parameters);
  _trophic_network   = new TrophicNetwork(_parameters, _population, _environment);
  _lineage_tree      = new Tree(_parameters);
  _phylogenetic_tree = new Tree(_parameters);
  _statistics        = new Statistics(_parameters, _population, _environment, _trophic_network, _lineage_tree, _phylogenetic_tree, true);
  _min_score         = 0.0;
  _max_score         = 0.0;
  _width             = _parameters->get_width();
  _height            = _parameters->get_height();
}

/**
 * \brief    Constructor from backup files
 * \details  Load the simulation from a specified backup file
 * \param    Parameters* parameters
 * \param    size_t backup_time
 * \param    bool clean_statistic_files
 * \return   \e void
 */
Simulation::Simulation( Parameters* parameters, size_t backup_time, bool clean_statistic_files )
{
  /*---------------------------*/
  /* 1) Set parameters         */
  /*---------------------------*/
  _parameters = parameters;
  
  /*---------------------------*/
  /* 2) Load population        */
  /*---------------------------*/
  std::stringstream pop_file_name;
  pop_file_name << "./population/population_" << backup_time;
  gzFile pop_file = gzopen(pop_file_name.str().c_str(), "r");
  _population = new Population(_parameters, pop_file);
  gzclose(pop_file);
  
  /*---------------------------*/
  /* 3) Load environment       */
  /*---------------------------*/
  std::stringstream env_file_name;
  env_file_name << "./environment/environment_" << backup_time;
  gzFile env_file = gzopen(env_file_name.str().c_str(), "r");
  _environment = new Environment(_parameters, env_file);
  gzclose(env_file);
  
  /*---------------------------*/
  /* 4) Load trophic network   */
  /*---------------------------*/
  std::stringstream trophic_file_name;
  trophic_file_name << "./trophic_network/trophic_network_" << backup_time;
  gzFile trophic_file = gzopen(trophic_file_name.str().c_str(), "r");
  _trophic_network = new TrophicNetwork(_parameters, _population, _environment, trophic_file);
  gzclose(trophic_file);
  
  /*---------------------------*/
  /* 5) Load lineage tree      */
  /*---------------------------*/
  std::stringstream tree_file_name;
  tree_file_name << "./tree/lineage_tree_" << backup_time;
  gzFile tree_file = gzopen(tree_file_name.str().c_str(), "r");
  _lineage_tree = new Tree(_parameters, _population, tree_file);
  gzclose(tree_file);
  
  /*---------------------------*/
  /* 6) Load phylogenetic tree */
  /*---------------------------*/
  tree_file_name.str("");
  tree_file_name << "./tree/phylogenetic_tree_" << backup_time;
  tree_file = gzopen(tree_file_name.str().c_str(), "r");
  _phylogenetic_tree = new Tree(_parameters, _population, tree_file);
  gzclose(tree_file);
  
  /*---------------------------*/
  /* 7) Initialize statistics  */
  /*---------------------------*/
  _statistics = new Statistics(_parameters, _population, _environment, _trophic_network, _lineage_tree, _phylogenetic_tree, clean_statistic_files);
  _statistics->flush_files();
  
  /*---------------------------*/
  /* 8) Initialize others      */
  /*---------------------------*/
  _min_score = 0.0;
  _max_score = 0.0;
  _width     = _parameters->get_width();
  _height    = _parameters->get_height();
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
Simulation::~Simulation( void )
{
  delete _population;
  _population = NULL;
  delete _environment;
  _environment = NULL;
  delete _trophic_network;
  _trophic_network = NULL;
  delete _lineage_tree;
  _lineage_tree = NULL;
  delete _phylogenetic_tree;
  _phylogenetic_tree = NULL;
  delete _statistics;
  _statistics = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Initialize the simulation
 * \details  --
 * \param    simulation_state state
 * \return   \e void
 */
void Simulation::initialize( simulation_state state )
{
  /*-----------------------------------*/
  /* 1) Initialize the environment     */
  /*-----------------------------------*/
  if (state == NEW_SIMULATION)
  {
    initialize_environment();
  }
  
  /*-----------------------------------*/
  /* 2) Initialize the population      */
  /*-----------------------------------*/
  initialize_population(state);
  mix();
  
  /*-----------------------------------*/
  /* 3) Actualize environment state    */
  /*-----------------------------------*/
  _environment->actualize_environment_state();
  
  /*-----------------------------------*/
  /* 4) Initialize the trophic network */
  /*-----------------------------------*/
  if (state == NEW_SIMULATION)
  {
    _trophic_network->initialize_trophic_network();
    _trophic_network->load_population();
  }
  
  /*-----------------------------------*/
  /* 5) Apply tests                    */
  /*-----------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
}

/**
 * \brief    Initialize the simulation from an evolved one (cells have already been modified)
 * \details  This method should only be called at simulation creation
 * \return   \e void
 */
void Simulation::initialize_from_evolved_population( void )
{
  /*-----------------------------------*/
  /* 1) Initialize the environment     */
  /*-----------------------------------*/
  initialize_environment();
  
  /*-----------------------------------*/
  /* 2) Initialize the population      */
  /*-----------------------------------*/
  initialize_population_from_evolved_population();
  mix();
  
  /*-----------------------------------*/
  /* 3) Actualize environment state    */
  /*-----------------------------------*/
  _environment->actualize_environment_state();
  
  /*-----------------------------------*/
  /* 4) Initialize the trophic network */
  /*-----------------------------------*/
  _trophic_network->initialize_trophic_network();
  _trophic_network->load_population();
  
  /*-----------------------------------*/
  /* 5) Apply tests                    */
  /*-----------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
}

/**
 * \brief    Update the simulation
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update( void )
{
  /*----------------------------------------------*/
  /* 1) Update the environment                    */
  /*----------------------------------------------*/
  if (_population->get_time() > 0)
  {
    update_environment();
  }
  
  /*----------------------------------------------*/
  /* 2) Update (an mix) the population            */
  /*----------------------------------------------*/
  update_population();
  mix();
  
  /*----------------------------------------------*/
  /* 3) Actualize environment state               */
  /*----------------------------------------------*/
  _environment->actualize_environment_state();
  
  /*----------------------------------------------*/
  /* 4) Update the trees, the trophic network and */
  /*    the environment                           */
  /*----------------------------------------------*/
  if (!_parameters->get_parallel_computing())
  {
    update_trophic_network();
    update_trees();
  }
  else if (_parameters->get_parallel_computing())
  {
    tbb::task_group tasks;
    tasks.run([=]{update_trophic_network();});
    tasks.run([=]{update_trees();});
    tasks.wait();
  }
  
  /*----------------------------------------------*/
  /* 4) Apply tests                               */
  /*----------------------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
  
}

/**
 * \brief    Save simulation in backup
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_simulation( void )
{
  if (!_parameters->get_parallel_computing())
  {
    save_parameters();
    save_population();
    save_environment();
    save_trophic_network();
    save_lineage_tree();
    save_phylogenetic_tree();
  }
  else if (_parameters->get_parallel_computing())
  {
    tbb::task_group tasks;
    tasks.run([=]{save_parameters();});
    tasks.run([=]{save_population();});
    tasks.run([=]{save_environment();});
    tasks.run([=]{save_trophic_network();});
    tasks.run([=]{save_lineage_tree();});
    tasks.run([=]{save_phylogenetic_tree();});
    tasks.wait();
  }
}

/**
 * \brief    Add random genetic units to each individual, and shuffle it if shuffle is true
 * \details  --
 * \param    size_t cell_pos
 * \param    size_t NC_type
 * \param    size_t E_type
 * \param    size_t TF_type
 * \param    size_t BS_type
 * \param    size_t P_type
 * \param    bool shuffle
 * \return   \e void
 */
void Simulation::add_random_genetic_units( size_t cell_pos, size_t NC_type, size_t E_type, size_t TF_type, size_t BS_type, size_t P_type, bool shuffle )
{
  Cell* cell = _population->get_cell(cell_pos);
  for (size_t i = 0; i < NC_type; i++)
  {
    genetic_unit unit;
    draw_random_genetic_unit(unit, NON_CODING);
    cell->get_genome()->add_genetic_unit(&unit);
  }
  for (size_t i = 0; i < E_type; i++)
  {
    genetic_unit unit;
    draw_random_genetic_unit(unit, ENZYME);
    cell->get_genome()->add_genetic_unit(&unit);
  }
  for (size_t i = 0; i < TF_type; i++)
  {
    genetic_unit unit;
    draw_random_genetic_unit(unit, TRANSCRIPTION_FACTOR);
    cell->get_genome()->add_genetic_unit(&unit);
  }
  for (size_t i = 0; i < BS_type; i++)
  {
    genetic_unit unit;
    draw_random_genetic_unit(unit, BINDING_SITE);
    cell->get_genome()->add_genetic_unit(&unit);
  }
  for (size_t i = 0; i < P_type; i++)
  {
    genetic_unit unit;
    draw_random_genetic_unit(unit, PROMOTER);
    cell->get_genome()->add_genetic_unit(&unit);
  }
  if (shuffle)
  {
    cell->get_genome()->shuffle();
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Initialize environment
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::initialize_environment( void )
{
  _environment->update(true, _population->get_time());
  _environment->compute_diffusion_and_degradation();
  for (size_t i = 1; i < _parameters->get_environment_properties()->number_of_init_cycles; i++)
  {
    _environment->update(false, _population->get_time());
    _environment->compute_diffusion_and_degradation();
  }
}

/**
 * \brief    Initialize population
 * \details  --
 * \param    simulation_state state
 * \return   \e void
 */
void Simulation::initialize_population( simulation_state state )
{
  _min_score = 1e+6;
  _max_score = 0.0;
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    
    /**************************************************/
    /* A) If simulation is freshly created            */
    /**************************************************/
    
    if (state == NEW_SIMULATION)
    {
      /*--------------------------------------------------*/
      /* A.1) Set coordinates                             */
      /*--------------------------------------------------*/
      size_t y = pos%_height;
      size_t x = (pos-y)/_height;
      cell->set_x(x);
      cell->set_y(y);
      
      /*--------------------------------------------------*/
      /* A.2) Build the genome                            */
      /*--------------------------------------------------*/
      if (!_parameters->get_load_genome_from_file())
      {
        add_random_genetic_units(pos, _parameters->get_initial_number_of_NC_units(), _parameters->get_initial_number_of_E_units(), _parameters->get_initial_number_of_TF_units(), _parameters->get_initial_number_of_BS_units(), _parameters->get_initial_number_of_P_units(), true);
      }
      else
      {
        std::vector<genetic_unit> sequence;
        _parameters->load_genome_from_file(sequence);
        for (size_t i = 0; i < sequence.size(); i++)
        {
          cell->get_genome()->add_genetic_unit(&sequence[i]);
        }
        sequence.clear();
      }
      
      /*--------------------------------------------------*/
      /* A.3) Mutate the genome                           */
      /*--------------------------------------------------*/
      cell->mutate();
      
      /*--------------------------------------------------*/
      /* A.4) Initialize metabolic amounts                */
      /*--------------------------------------------------*/
      int    tag  = _parameters->draw_initial_metabolite_tag();
      double conc = _parameters->get_initial_metabolites_amount_in_cells();
      cell->get_species_list()->add(tag, conc);
      cell->initialize_inherited_species_list();
      
      /*--------------------------------------------------*/
      /* A.5) Load genome in species list and environment */
      /*--------------------------------------------------*/
      cell->load_genome_in_species_lists();
      
      /*--------------------------------------------------*/
      /* A.6) Initialize energy                           */
      /*--------------------------------------------------*/
      cell->set_energy(0.0);
      if (_parameters->get_energy_costs_scheme())
      {
        cell->set_energy(_parameters->get_initial_energy_amount_in_cells());
      }
      
      /*--------------------------------------------------*/
      /* A.7) Activate the cell                           */
      /*--------------------------------------------------*/
      cell->activate();
      cell->set_alive(true);
      
      /*--------------------------------------------------*/
      /* A.8) Initialize time variables                   */
      /*--------------------------------------------------*/
      cell->set_birth_time(_population->get_time());
      
      /*--------------------------------------------------*/
      /* A.9) Set cell id                                 */
      /*--------------------------------------------------*/
      cell->set_parent_id(0);
      cell->set_id(_population->get_new_id());
      
      /*--------------------------------------------------*/
      /* A.10) Compute score                              */
      /*--------------------------------------------------*/
      cell->synchronize_state_vectors(_environment);
      cell->load_genome_in_ODE_system(_environment, false, true);
      cell->update(_population->get_time());
      if (_min_score > cell->get_score())
      {
        _min_score = cell->get_score();
      }
      if (_max_score < cell->get_score())
      {
        _max_score = cell->get_score();
        _population->set_best(cell->get_id(), pos);
      }
      
      /*--------------------------------------------------*/
      /* A.11) Add the cell as a root in trees            */
      /*--------------------------------------------------*/
      _lineage_tree->add_root(cell);
      _phylogenetic_tree->add_root(cell);
    }
    
    /**************************************************/
    /* B) Else if simulation is loaded from backup    */
    /**************************************************/
    
    else if (state == FROM_BACKUP)
    {
      /*-----------------------------------------------------*/
      /* B.1) Synchronize cell and environment state vectors */
      /*-----------------------------------------------------*/
      cell->synchronize_state_vectors(_environment);
      
      /*-----------------------------------------------------*/
      /* B.2) Load genome in ODE system                      */
      /*-----------------------------------------------------*/
      cell->load_genome_in_ODE_system(_environment, true, false);
    }
  }
}

/**
 * \brief    Initialize the population from an evolved population (cells have already been modified)
 * \details  This method should only be called at simulation creation
 * \return   \e void
 */
void Simulation::initialize_population_from_evolved_population( void )
{
  _min_score = 1e+6;
  _max_score = 0.0;
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    
    /*------------------------------------------------*/
    /* 1) Set coordinates                             */
    /*------------------------------------------------*/
    size_t y = pos%_height;
    size_t x = (pos-y)/_height;
    cell->set_x(x);
    cell->set_y(y);
    
    /*------------------------------------------------*/
    /* 2) Load genome in species list and environment */
    /*------------------------------------------------*/
    cell->load_genome_in_species_lists();
    
    /*------------------------------------------------*/
    /* 3) Activate the cell                           */
    /*------------------------------------------------*/
    cell->activate();
    cell->set_alive(true);
    
    /*------------------------------------------------*/
    /* 4) Initialize time variables                   */
    /*------------------------------------------------*/
    cell->set_birth_time(_population->get_time());
    
    /*------------------------------------------------*/
    /* 5) Set cell id                                 */
    /*------------------------------------------------*/
    cell->set_parent_id(0);
    cell->set_id(_population->get_new_id());
    
    /*------------------------------------------------*/
    /* 6) Compute score                               */
    /*------------------------------------------------*/
    cell->synchronize_state_vectors(_environment);
    cell->load_genome_in_ODE_system(_environment, false, true);
    cell->update(_population->get_time());
    if (_min_score > cell->get_score())
    {
      _min_score = cell->get_score();
    }
    if (_max_score < cell->get_score())
    {
      _max_score = cell->get_score();
      _population->set_best(cell->get_id(), pos);
    }
    
    /*------------------------------------------------*/
    /* 7) Add the cell as a root in trees             */
    /*------------------------------------------------*/
    if (_parameters->get_simulation_backup_step() > 0)
    {
      _lineage_tree->add_root(cell);
      _phylogenetic_tree->add_root(cell);
    }
  }
}

/**
 * \brief    Draw genetic unit's attributes at random
 * \details  --
 * \param    genetic_unit& unit
 * \param    genetic_unit_type type
 * \return   \e void
 */
void Simulation::draw_random_genetic_unit( genetic_unit& unit, genetic_unit_type type )
{
  /*------------------------------------------------------------------ Global attributes */
  
  unit.type              = type;
  unit.identifier        = 0;
  unit.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  unit.s             = _parameters->draw_initial_metabolite_tag();
  unit.s             = (unit.s > 0 ? unit.s : 1);
  unit.p             = _parameters->draw_initial_metabolite_tag();
  unit.p             = (unit.p > 0 ? unit.p : 1);
  unit.kcat          = pow(10.0, _parameters->get_simulation_prng()->uniform()*(KCAT_MAX_LOG-KCAT_MIN_LOG)+KCAT_MIN_LOG);
  unit.kcat          = unit.kcat*(_parameters->get_simulation_prng()->uniform() < 0.5 ? -1.0 : 1.0);
  unit.kcat_km_ratio = pow(10.0, _parameters->get_simulation_prng()->uniform()*(KCAT_KM_RATIO_MAX_LOG-KCAT_KM_RATIO_MIN_LOG)+KCAT_KM_RATIO_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  unit.BS_tag         = _parameters->draw_initial_binding_site_tag();
  unit.coE_tag        = _parameters->draw_initial_co_enzyme_tag();
  unit.coE_tag        = (unit.coE_tag > 0 ? unit.coE_tag : 1);
  unit.free_activity  = (_parameters->get_simulation_prng()->uniform() < 0.5 ? true : false);
  unit.bound_activity = (_parameters->get_simulation_prng()->uniform() < 0.5 ? true : false);
  unit.binding_window = _parameters->get_transcription_factor_binding_window();
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  unit.TF_tag = _parameters->draw_initial_transcription_factor_tag();
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  unit.basal_expression_level = _parameters->get_simulation_prng()->uniform();
}

/**
 * \brief    Update the environment
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_environment( void )
{
  _environment->update(false, _population->get_time());
  _environment->compute_diffusion_and_degradation();
}

/**
 * \brief    Update the population
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_population( void )
{
  _population->set_previous_size();
  
  /***************************************************************************/
  /* 1) Initialize population size and statistics                            */
  /***************************************************************************/
  _population->set_population_size(0);
  _statistics->init_variables();
  _min_score              = 1e+6;
  _max_score              = 0.0;
  cell_state* cell_states = new cell_state[_width*_height];
  
  /***************************************************************************/
  /* 2) Explore population: save gaps, kill or update cells unable to divide */
  /***************************************************************************/
  std::vector<size_t> gaps;
  
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    cell->untag();
    cell_states[pos] = NOTHING;
    
    /*---------------------------------------------------------*/
    /* 2.1) If the cell is alive                               */
    /*---------------------------------------------------------*/
    if (cell->isAlive())
    {
      /* A) Increase population size --------*/
      _population->increase_population_size(1);
      
      /* B) Compute statistics --------------*/
      _statistics->add_individual(pos);
      if (_min_score > cell->get_score())
      {
        _min_score = cell->get_score();
      }
      if (_max_score < cell->get_score())
      {
        _max_score = cell->get_score();
        _population->set_best(cell->get_id(), pos);
      }
      
      /* C) Evaluate the individual ---------*/
      size_t child_position = 0;
      bool   isActive       = cell->isActive();
      bool   amount         = (cell->get_amount() > MINIMUM_CONCENTRATION);
      bool   score          = (cell->get_score() > MINIMUM_SCORE);
      bool   energy         = false;
      if ((_parameters->get_energy_costs_scheme() && cell->get_energy() > MINIMUM_CONCENTRATION) || !_parameters->get_energy_costs_scheme())
      {
        energy = true;
      }
      bool dying    = (_parameters->get_simulation_prng()->uniform() < _parameters->get_death_probability() || !energy || !amount || !isActive || !score || cell->get_replication_report()->get_genome_functional_size() == 0);
      bool updating = !find_empty_cases(pos, child_position);
      
      /* D) Update or kill individual -------*/
      if (dying)
      {
        cell_states[pos] = TO_KILL;
        cell->tag();
      }
      else if (!dying && updating)
      {
        cell_states[pos] = TO_UPDATE;
        cell->tag();
      }
      else if (!dying && !updating)
      {
        cell_states[pos] = TO_UPDATE;
      }
    }
    /*---------------------------------------------------------*/
    /* 2.2) Else if the cell is dead, add gap to the gaps list */
    /*---------------------------------------------------------*/
    else
    {
      gaps.push_back(pos);
      cell_states[pos] = GAP;
    }
  }
  
  /***************************************************************************/
  /* 3) Compute statistics                                                   */
  /***************************************************************************/
  _statistics->compute_mean_and_var();
  
  /***************************************************************************/
  /* 4) Shuffle gaps vector                                                  */
  /***************************************************************************/
  int size = (int)gaps.size();
  for (int i = 0; i < size; i++)
  {
    size_t pos1 = _parameters->get_simulation_prng()->uniform(0, size-1);
    size_t pos2 = _parameters->get_simulation_prng()->uniform(0, size-1);
    size_t tmp  = gaps[pos1];
    gaps[pos1]  = gaps[pos2];
    gaps[pos2]  = tmp;
  }
  
  /***************************************************************************/
  /* 5) Explore gaps and compute cell divisions                              */
  /***************************************************************************/
  double threshold = _parameters->get_selection_threshold();
  
  /* A) For each gap -------------------------------------------------------*/
  for (size_t gap_index = 0; gap_index < gaps.size(); gap_index++)
  {
    size_t y_gap = gaps[gap_index]%_height;
    size_t x_gap = (gaps[gap_index]-y_gap)/_height;
    size_t selected_position = 0;
    double propensity        = 0.0;
    double best_propensity   = 0.0;
    
    /* B) Explore gap moore neighbourhood ----------------------------------*/
    for (int i = -1; i < 2; i++)
    {
      for (int j = -1; j < 2; j++)
      {
        if ( !(i == 0 && j == 0) )
        {
          size_t x     = (x_gap + i + _width) % _width;
          size_t y     = (y_gap + j + _height) % _height;
          Cell*  cell  = _population->get_cell(x, y);
          double score = cell->get_score();
          
          /* C) If the cell is able to divide itself, save its coordinates -*/
          if (cell->isActive() && !cell->isTagged() && score > threshold*_max_score)
          {
            propensity = score;
            if (best_propensity < propensity)
            {
              best_propensity   = propensity;
              selected_position = x*_height+y;
            }
          }
        }
      }
    }
    
    /* C) if best propensity is not null, divide the cell ------------------*/
    if (best_propensity > 0.0)
    {
      /* Compute cell replication */
      divide_cell(selected_position, gaps[gap_index]);
      
      /* Update parent cell ------*/
      _population->get_cell(selected_position)->update_number_of_divisions();
      cell_states[selected_position] = TO_MUTATE;
      _population->get_cell(selected_position)->tag();
      
      /* Update daughter cell ----*/
      cell_states[gaps[gap_index]] = TO_MUTATE;
      _population->get_cell(gaps[gap_index])->tag();
    }
  }
  
  /***************************************************************************/
  /* 6) Update cell state                                                    */
  /***************************************************************************/
  if (!_parameters->get_parallel_computing())
  {
    for (size_t pos = 0; pos < _parameters->get_width()*_parameters->get_height(); pos++)
    {
      if (cell_states[pos] == TO_KILL)
      {
        kill_cell(pos);
      }
      else if (cell_states[pos] == TO_UPDATE)
      {
        update_cell(pos);
      }
      else if (cell_states[pos] == TO_MUTATE)
      {
        mutate_cell(pos);
      }
    }
  }
  else if (_parameters->get_parallel_computing())
  {
    tbb::task_group tasks;
    for (size_t pos = 0; pos < _width*_height; pos++)
    {
      if (cell_states[pos] == TO_KILL)
      {
        tasks.run([=]{kill_cell(pos);});
      }
      else if (cell_states[pos] == TO_UPDATE)
      {
        tasks.run([=]{update_cell(pos);});
      }
      else if (cell_states[pos] == TO_MUTATE)
      {
        tasks.run([=]{mutate_cell(pos);});
      }
    }
    tasks.wait();
  }
  
  /***************************************************************************/
  /* 7) Freeze dead cell's nodes outside of the loop                         */
  /***************************************************************************/
  if (_parameters->get_simulation_backup_step() > 0)
  {
    for (size_t pos = 0; pos < _width*_height; pos++)
    {
      if (cell_states[pos] == TO_KILL)
      {
        Cell* cell = _population->get_cell(pos);
        _lineage_tree->freeze_node(cell->get_id(), _population->get_time());
        _phylogenetic_tree->freeze_node(cell->get_id(), _population->get_time());
      }
    }
  }
  delete[] cell_states;
  cell_states = NULL;
  
  /***************************************************************************/
  /* 8) Compute population growth rate                                       */
  /***************************************************************************/
  _population->compute_growth_rate();
  
  /***************************************************************************/
  /* 9) Update population time                                               */
  /***************************************************************************/
  _population->update_time();
  
}

/**
 * \brief    Save parameters in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_parameters( void )
{
  _parameters->save(_population->get_time());
}

/**
 * \brief    Save the population in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_population( void )
{
  std::stringstream pop_file_name;
  pop_file_name << "./population/population_" << _population->get_time();
  gzFile pop_file = gzopen(pop_file_name.str().c_str(), "w");
  _population->save(pop_file);
  gzclose(pop_file);
}

/**
 * \brief    Save the environment in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_environment( void )
{
  std::stringstream env_file_name;
  env_file_name << "./environment/environment_" << _population->get_time();
  gzFile env_file = gzopen(env_file_name.str().c_str(), "w");
  _environment->save(env_file);
  gzclose(env_file);
}

/**
 * \brief    Save the trophic network in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_trophic_network( void )
{
  std::stringstream trophic_file_name;
  trophic_file_name << "./trophic_network/trophic_network_" << _population->get_time();
  gzFile trophic_file = gzopen(trophic_file_name.str().c_str(), "w");
  _trophic_network->save(trophic_file);
  gzclose(trophic_file);
}

/**
 * \brief    Save the lineage tree in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_lineage_tree( void )
{
  std::stringstream tree_file_name;
  tree_file_name << "./tree/lineage_tree_" << _population->get_time();
  gzFile tree_file = gzopen(tree_file_name.str().c_str(), "w");
  _lineage_tree->save(tree_file);
  gzclose(tree_file);
}

/**
 * \brief    Save the phylogenetic tree in backup file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_phylogenetic_tree( void )
{
  std::stringstream tree_file_name;
  tree_file_name << "./tree/phylogenetic_tree_" << _population->get_time();
  gzFile tree_file = gzopen(tree_file_name.str().c_str(), "w");
  _phylogenetic_tree->save(tree_file);
  gzclose(tree_file);
}

/**
 * \brief    Update the trophic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_trophic_network( void )
{
  _trophic_network->load_population();
}

/**
 * \brief    Update trees
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_trees( void )
{
  if (_parameters->get_simulation_backup_step() > 0)
  {
    _lineage_tree->clean_cell_map();
    _lineage_tree->prune();
    _phylogenetic_tree->clean_cell_map();
    _phylogenetic_tree->prune();
    _phylogenetic_tree->shorten();
  }
}

/**
 * \brief    Mix the population
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::mix( void )
{
  if (_parameters->get_migration_rate() <= 0.0)
  {
    return;
  }
  size_t n = _parameters->get_simulation_prng()->binomial(_width*_height, _parameters->get_migration_rate());
  while (n > 0)
  {
    size_t x1 = _parameters->get_simulation_prng()->uniform(0, (int)_width-1);
    size_t y1 = _parameters->get_simulation_prng()->uniform(0, (int)_height-1);
    size_t x2 = _parameters->get_simulation_prng()->uniform(0, (int)_width-1);
    size_t y2 = _parameters->get_simulation_prng()->uniform(0, (int)_height-1);
    Cell* tmp1 = _population->get_cell(x1, y1);
    tmp1->set_x(x2);
    tmp1->set_y(y2);
    Cell* tmp2 = _population->get_cell(x2, y2);
    tmp2->set_x(x1);
    tmp2->set_y(y1);
    _population->set_cell(x1*_height+y1, tmp2);
    _population->set_cell(x2*_height+y2, tmp1);
    n--;
  }
}

/**
 * \brief    Kill the cell at position i
 * \details  --
 * \param    size_t i
 * \return   \e void
 */
void Simulation::kill_cell( size_t i )
{
  assert(i < _width*_height);
  Cell* cell = _population->get_cell(i);
  size_t y = i%_height;
  size_t x = (i-y)/_height;
  
  /*--------------------------------------------------*/
  /* 1) Copy metabolite concentrations in environment */
  /*--------------------------------------------------*/
  if (_parameters->get_environment_properties()->interaction_scheme == INTERACTION)
  {
    SpeciesList *cell_list = cell->get_species_list();
    for (int tag = 1; tag <= (int)cell_list->get_size(); tag++)
    {
      _environment->add(x, y, tag, cell_list->get(tag));
    }
  }
  
  /*--------------------------------------------------*/
  /* 2) Kill the cell                                 */
  /*--------------------------------------------------*/
  /*
  if (_parameters->get_simulation_backup_step() > 0)
  {
    _lineage_tree->freeze_node(cell->get_id(), _population->get_time());
    _phylogenetic_tree->freeze_node(cell->get_id(), _population->get_time());
  }
   */
  cell->kill(_population->get_time());
}

/**
 * \brief    Mutate the cell at position i
 * \details  --
 * \param    size_t i
 * \return   \e void
 */
void Simulation::mutate_cell( size_t i )
{
  /*-------------------------------------------------------*/
  /* 1) check cell position                                */
  /*-------------------------------------------------------*/
  assert(i < _width*_height);
  
  /*-------------------------------------------------------*/
  /* 2) get cell                                           */
  /*-------------------------------------------------------*/
  Cell* cell = _population->get_cell(i);
  
  /*-------------------------------------------------------*/
  /* 3) mutate the cell                                    */
  /*-------------------------------------------------------*/
  cell->mutate();
  cell->load_genome_in_species_lists();
  cell->synchronize_state_vectors(_environment);
  cell->load_genome_in_ODE_system(_environment, false, false);
  
  /*-------------------------------------------------------*/
  /* 4) solve ODEs, update species lists and compute score */
  /*-------------------------------------------------------*/
  cell->update(_population->get_time());
}

/**
 * \brief    Update the cell at position i
 * \details  --
 * \param    size_t i
 * \return   \e void
 */
void Simulation::update_cell( size_t i )
{
  /*----------------------------------------------------------------*/
  /* 1) check cell position                                         */
  /*----------------------------------------------------------------*/
  assert(i < _width*_height);
  
  /*----------------------------------------------------------------*/
  /* 2) get cell                                                    */
  /*----------------------------------------------------------------*/
  Cell* cell = _population->get_cell(i);
  
  /*----------------------------------------------------------------*/
  /* 3) synchronize species lists size between cell and environment */
  /*----------------------------------------------------------------*/
  cell->synchronize_state_vectors(_environment);
  
  /*----------------------------------------------------------------*/
  /* 4) solve ODEs, update species lists and compute score          */
  /*----------------------------------------------------------------*/
  cell->update(_population->get_time());
}

/**
 * \brief    Divide cell at position i in position child_position
 * \details  --
 * \param    size_t i
 * \param    size_t child_position
 * \return   \e void
 */
void Simulation::divide_cell( size_t i, size_t child_position )
{
  assert(i < _width*_height);
  assert(child_position < _width*_height);
  Cell* parent = _population->get_cell(i);
  if (_parameters->get_simulation_backup_step() > 0)
  {
    _lineage_tree->freeze_node(parent->get_id(), _population->get_time());
    _phylogenetic_tree->freeze_node(parent->get_id(), _population->get_time());
  }
  
  /*-----------------*/
  /* 1) Create child */
  /*-----------------*/
  _population->new_cell(i, child_position);
  Cell* child = _population->get_cell(child_position);
  
  /*-----------------*/
  /* 2) Update tree  */
  /*-----------------*/
  if (_parameters->get_simulation_backup_step() > 0)
  {
    _lineage_tree->add_division(parent, parent, child);
    _phylogenetic_tree->add_division(parent, parent, child);
  }
}

/**
 * \brief    Find empty cases around cell at position i
 * \details  --
 * \param    size_t i
 * \param    size_t& child_position
 * \return   \e bool
 */
bool Simulation::find_empty_cases( size_t i, size_t& child_position )
{
  assert(i < _width*_height);
  size_t x_coord[8];
  size_t y_coord[8];
  int    N = 0;
  size_t y = i%_height;
  size_t x = (i-y)/_height;
  for (int i = -1; i < 2; i++)
  {
    for (int j = -1; j < 2; j++)
    {
      if ( !(i == 0 && j == 0) )
      {
        size_t new_x = (x + i + _width) % _width;
        size_t new_y = (y + j + _height) % _height;
        if ( !_population->get_cell(new_x, new_y)->isAlive() )
        {
          x_coord[N] = new_x;
          y_coord[N] = new_y;
          N++;
        }
      }
    }
  }
  if (N == 0)
  {
    return false;
  }
  else
  {
    size_t draw = _parameters->get_simulation_prng()->uniform(0, N-1);
    child_position = x_coord[draw]*_height + y_coord[draw];
    return true;
  }
}

/**
 * \brief    Test the structure of the lineage tree
 * \details  --
 * \param    void
 * \return   \e void
 */
#ifdef DEBUG
void Simulation::test_lineage_tree_structure( void )
{
  /*--------------------------*/
  /* 1) Check tree structure  */
  /*--------------------------*/
  size_t master_root_count = 0;
  size_t root_count        = 0;
  size_t normal_count      = 0;
  size_t dead_count        = 0;
  size_t alive_count       = 0;
  Node* node = _lineage_tree->get_first_node();
  while (node != NULL)
  {
    /*----------------------------*/
    /* 1.1) Test master root node */
    /*----------------------------*/
    if (node->isMasterRoot())
    {
      master_root_count++;
      dead_count++;
      assert(node->get_parent() == NULL);
      assert(!node->isRoot());
      assert(!node->isNormal());
      assert(node->isDead());
      assert(!node->isAlive());
      assert(node->get_alive_cell() == NULL);
      assert(node->get_replication_report() == NULL);
    }
    /*----------------------------*/
    /* 1.2) Test root nodes       */
    /*----------------------------*/
    else if (node->isRoot())
    {
      root_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isNormal());
      assert(node->get_parent()->isMasterRoot());
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
      }
    }
    /*----------------------------*/
    /* 1.3) Test normal nodes     */
    /*----------------------------*/
    else if (node->isNormal())
    {
      normal_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isRoot());
      assert(node->get_parent()->isRoot() || node->get_parent()->isNormal());
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
      }
    }
    
    node = _lineage_tree->get_next_node();
  }
  
  /*--------------------------*/
  /* 2) Test nodes counts     */
  /*--------------------------*/
  assert(master_root_count == 1);
  assert(master_root_count+root_count+normal_count == _lineage_tree->get_number_of_nodes());
  assert(dead_count+alive_count == _lineage_tree->get_number_of_nodes());
  std::vector<unsigned long long int> alive_nodes;
  _lineage_tree->get_alive_nodes(&alive_nodes);
  assert(alive_nodes.size() == alive_count);
  for (size_t i = 0; i < alive_nodes.size(); i++)
  {
    Node* node = _lineage_tree->get_node(alive_nodes[i]);
    assert(!node->isMasterRoot());
    assert(node->isRoot() || node->isNormal());
    assert(node->isAlive());
    assert(node->get_alive_cell()->isAlive());
    Cell* cell = _population->get_cell_by_id(node->get_alive_cell()->get_id());
    assert(cell->isAlive());
    assert(node->get_replication_report()->get_number_of_events() == cell->get_replication_report()->get_number_of_events());
    (void)cell;
  }
  
  /*--------------------------*/
  /* 3) Check the master root */
  /*--------------------------*/
  Node* master_root = _lineage_tree->get_node(0);
  assert(master_root->isMasterRoot());
  assert(master_root->get_number_of_children() == root_count);
  (void)master_root;
}
#endif

/**
 * \brief    Test the structure of the phylogenetic tree
 * \details  --
 * \param    void
 * \return   \e void
 */
#ifdef DEBUG
void Simulation::test_phylogenetic_tree_structure( void )
{
  /*--------------------------*/
  /* 1) Check tree structure  */
  /*--------------------------*/
  size_t master_root_count = 0;
  size_t root_count        = 0;
  size_t normal_count      = 0;
  size_t dead_count        = 0;
  size_t alive_count       = 0;
  Node* node = _phylogenetic_tree->get_first_node();
  while (node != NULL)
  {
    /*----------------------------*/
    /* 1.1) Test master root node */
    /*----------------------------*/
    if (node->isMasterRoot())
    {
      master_root_count++;
      dead_count++;
      assert(node->get_parent() == NULL);
      assert(!node->isRoot());
      assert(!node->isNormal());
      assert(node->isDead());
      assert(!node->isAlive());
      assert(node->get_alive_cell() == NULL);
      assert(node->get_replication_report() == NULL);
    }
    /*----------------------------*/
    /* 1.2) Test root nodes       */
    /*----------------------------*/
    else if (node->isRoot())
    {
      root_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isNormal());
      assert(node->get_parent()->isMasterRoot());
      bool alive_child = false;
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
        if (node->get_child(i)->isAlive())
        {
          assert(_population->get_cell_by_id(node->get_child(i)->get_alive_cell()->get_id()) != NULL);
          alive_child = true;
        }
      }
      if (node->isDead() && !alive_child)
      {
        assert(node->get_number_of_children() == 2);
      }
    }
    /*----------------------------*/
    /* 1.3) Test normal nodes     */
    /*----------------------------*/
    else if (node->isNormal())
    {
      normal_count++;
      if (node->isDead())
      {
        assert(node->get_number_of_children() == 2);
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isRoot());
      assert(node->get_parent()->isRoot() || node->get_parent()->isNormal());
      bool alive_child = false;
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
        if (node->get_child(i)->isAlive())
        {
          assert(_population->get_cell_by_id(node->get_child(i)->get_alive_cell()->get_id()) != NULL);
          alive_child = true;
        }
      }
      if (node->isDead() && !alive_child)
      {
        assert(node->get_number_of_children() == 2);
      }
    }
    
    node = _phylogenetic_tree->get_next_node();
  }
  
  /*--------------------------*/
  /* 2) Test nodes counts     */
  /*--------------------------*/
  assert(master_root_count == 1);
  assert(master_root_count+root_count+normal_count == _phylogenetic_tree->get_number_of_nodes());
  assert(dead_count+alive_count == _phylogenetic_tree->get_number_of_nodes());
  std::vector<unsigned long long int> alive_nodes;
  _phylogenetic_tree->get_alive_nodes(&alive_nodes);
  assert(alive_nodes.size() == alive_count);
  for (size_t i = 0; i < alive_nodes.size(); i++)
  {
    Node* node = _phylogenetic_tree->get_node(alive_nodes[i]);
    assert(!node->isMasterRoot());
    assert(node->isRoot() || node->isNormal());
    assert(node->isAlive());
    assert(node->get_alive_cell()->isAlive());
    Cell* cell = _population->get_cell_by_id(node->get_alive_cell()->get_id());
    assert(cell->isAlive());
    assert(node->get_replication_report()->get_number_of_events() == cell->get_replication_report()->get_number_of_events());
    (void)cell;
  }
  
  /*--------------------------*/
  /* 3) Check the master root */
  /*--------------------------*/
  Node* master_root = _phylogenetic_tree->get_node(0);
  assert(master_root->isMasterRoot());
  assert(master_root->get_number_of_children() == root_count);
  (void)master_root;
}
#endif
