
/**
 * \file      TrophicNetwork.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicNetwork class definition
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

#include "TrophicNetwork.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* population
 * \param    Environment* environment
 * \return   \e void
 */
TrophicNetwork::TrophicNetwork( Parameters* parameters, Population* population, Environment* environment )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  _parameters  = parameters;
  _population  = population;
  _environment = environment;
  _current_id  = 0;
  _group_map.clear();
  _iterator = _group_map.begin();
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  _nb_level_0_groups  = 0;
  _nb_level_1_groups  = 0;
  _nb_level_2_groups  = 0;
  _nb_no_level_groups = 0;
  
  _nb_level_0_cells  = 0;
  _nb_level_1_cells  = 0;
  _nb_level_2_cells  = 0;
  _nb_no_level_cells = 0;
  
  _nb_group_appearances = 0;
  _nb_group_extinctions = 0;
  
  _mean_group_lifespan = 0.0;
  
}

/**
 * \brief    Constructor from backup file
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* population
 * \param    Environment* environment
 * \param    gzFile backup_file
 * \return   \e void
 */
TrophicNetwork::TrophicNetwork( Parameters* parameters, Population* population, Environment* environment, gzFile backup_file )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  _parameters  = parameters;
  _population  = population;
  _environment = environment;
  gzread( backup_file, &_current_id, sizeof(_current_id) );
  size_t n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _group_map.clear();
  for (size_t i = 0; i < n; i++)
  {
    TrophicGroup* current_group = new TrophicGroup(backup_file);
    _group_map[current_group->get_identifier()] = current_group;
  }
  _iterator = _group_map.begin();
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  gzread( backup_file, &_nb_level_0_groups,  sizeof(_nb_level_0_groups) );
  gzread( backup_file, &_nb_level_1_groups,  sizeof(_nb_level_1_groups) );
  gzread( backup_file, &_nb_level_2_groups,  sizeof(_nb_level_2_groups) );
  gzread( backup_file, &_nb_no_level_groups, sizeof(_nb_no_level_groups) );
  
  gzread( backup_file, &_nb_level_0_cells,  sizeof(_nb_level_0_cells) );
  gzread( backup_file, &_nb_level_1_cells,  sizeof(_nb_level_1_cells) );
  gzread( backup_file, &_nb_level_2_cells,  sizeof(_nb_level_2_cells) );
  gzread( backup_file, &_nb_no_level_cells, sizeof(_nb_no_level_cells) );
  
  gzread( backup_file, &_nb_group_appearances, sizeof(_nb_group_appearances) );
  gzread( backup_file, &_nb_group_extinctions, sizeof(_nb_group_extinctions) );
  
  gzread( backup_file, &_mean_group_lifespan, sizeof(_mean_group_lifespan) );
  
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
TrophicNetwork::~TrophicNetwork( void )
{
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    delete _iterator->second;
    _iterator->second = NULL;
  }
  _group_map.clear();
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
void TrophicNetwork::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  gzwrite( backup_file, &_current_id, sizeof(_current_id) );
  size_t n = _group_map.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    _iterator->second->save(backup_file);
  }
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  gzwrite( backup_file, &_nb_level_0_groups,  sizeof(_nb_level_0_groups) );
  gzwrite( backup_file, &_nb_level_1_groups,  sizeof(_nb_level_1_groups) );
  gzwrite( backup_file, &_nb_level_2_groups,  sizeof(_nb_level_2_groups) );
  gzwrite( backup_file, &_nb_no_level_groups, sizeof(_nb_no_level_groups) );
  
  gzwrite( backup_file, &_nb_level_0_cells,  sizeof(_nb_level_0_cells) );
  gzwrite( backup_file, &_nb_level_1_cells,  sizeof(_nb_level_1_cells) );
  gzwrite( backup_file, &_nb_level_2_cells,  sizeof(_nb_level_2_cells) );
  gzwrite( backup_file, &_nb_no_level_cells, sizeof(_nb_no_level_cells) );
  
  gzwrite( backup_file, &_nb_group_appearances, sizeof(_nb_group_appearances) );
  gzwrite( backup_file, &_nb_group_extinctions, sizeof(_nb_group_extinctions) );
  
  gzwrite( backup_file, &_mean_group_lifespan, sizeof(_mean_group_lifespan) );
  
}

/**
 * \brief    Initialize the trophic network
 * \details  This method basically creates the environment group
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::initialize_trophic_network( void )
{
  size_t N    = _environment->get_size();
  _current_id = 0;
  
  /*-----------------------------------------------------------------*/
  /* 1) Initialize the environment trophic group and profile strings */
  /*-----------------------------------------------------------------*/
  TrophicGroup* env_group = new TrophicGroup(get_new_id(), _population->get_time(), 0.0, 0.0, 0.0);
  env_group->_trophic_profile    = "";
  env_group->_production_profile = "";
  env_group->_uptake_profile     = "";
  env_group->_release_profile    = "";
  
  /*-----------------------------------------------------------------*/
  /* 2) Build the environment production profile                     */
  /*-----------------------------------------------------------------*/
  for (size_t i = 0; i < N; i++)
  {
    if ((i+1) >= _parameters->get_environment_properties()->species_tag_range.min && (i+1) <= _parameters->get_environment_properties()->species_tag_range.max)
    {
      env_group->_production_profile += "1";
    }
    else
    {
      env_group->_production_profile += "0";
    }
  }
  
  /*-----------------------------------------------------------------*/
  /* 3) Build the environment uptake and release profiles            */
  /*-----------------------------------------------------------------*/
  for (size_t i = 0; i < N; i++)
  {
    env_group->_uptake_profile  += "0";
    env_group->_release_profile += "0";
  }
  
  /*-----------------------------------------------------------------*/
  /* 4) Build environment trophic profile                            */
  /*-----------------------------------------------------------------*/
  env_group->_trophic_profile = env_group->_production_profile+env_group->_uptake_profile+env_group->_release_profile;
  
  /*-----------------------------------------------------------------*/
  /* 5) Add the environment node to the trophic network              */
  /*-----------------------------------------------------------------*/
  _group_map.clear();
  _group_map[env_group->get_identifier()] = env_group;
  _nb_no_level_groups++;
  _nb_group_appearances++;
}

/**
 * \brief    Load the population in the trophic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::load_population( void )
{
  size_t N = _environment->get_size();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Clear statistics and relationships */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    _iterator->second->update_profile(N);
    _iterator->second->clear();
  }
  _nb_level_0_groups    = 0;
  _nb_level_1_groups    = 0;
  _nb_level_2_groups    = 0;
  _nb_no_level_groups   = 0;
  _nb_level_0_cells     = 0;
  _nb_level_1_cells     = 0;
  _nb_level_2_cells     = 0;
  _nb_no_level_cells    = 0;
  _nb_group_appearances = 0;
  _nb_group_extinctions = 0;
  _mean_group_lifespan  = 0.0;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Explore cell profiles              */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (size_t pos = 0; pos < _population->get_width()*_population->get_height(); pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      double*        Xcell = cell->get_species_list()->get_X();
      double*        Xenv  = _environment->get_X()->get_X();
      reaction_list* rlist = cell->get_ode()->get_reaction_list();
      
      std::string trophic_profile    = "";
      std::string production_profile = "";
      std::string uptake_profile     = "";
      std::string release_profile    = "";
      
      /*-----------------------------------------------------*/
      /* 2.1) Build the uptake profile                       */
      /*-----------------------------------------------------*/
      bool* profile = new bool[N];
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int s = rlist->metabolic_s[i];
        /* If there is by-products in the environment (Xenv) */
        /* or exogeneous nutrient (env_profile), the uptake  */
        /* profile is 1. Else it is 0.                       */
        if (rlist->metabolic_type[i] == INFLOWING_PUMP_ACTIVITY && (Xenv[s-1] > MINIMUM_CONCENTRATION || _group_map[0]->_production_profile[s-1] == '1'))
        {
          profile[s-1] = true;
        }
      }
      /* Build the uptake profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          uptake_profile += "1";
        }
        else
        {
          uptake_profile += "0";
        }
      }
      
      /*-----------------------------------------------------*/
      /* 2.2) Build the release profile                      */
      /*-----------------------------------------------------*/
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int s = rlist->metabolic_s[i];
        /* If the release pump is active, set the profile at 1 */
        if (rlist->metabolic_type[i] == OUTFLOWING_PUMP_ACTIVITY && Xcell[s-1] > MINIMUM_CONCENTRATION)
        {
          profile[s-1] = true;
        }
      }
      /* Build the release profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          release_profile += "1";
        }
        else
        {
          release_profile += "0";
        }
      }
      
      /*-----------------------------------------------------*/
      /* 2.3) Build the production profile                   */
      /*-----------------------------------------------------*/
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int p = rlist->metabolic_p[i];
        /* If the reaction produces a not uptaken metabolite, */
        /* set the production profile to 1                    */
        if (rlist->metabolic_type[i] == CATALYTIC_CONSUMING_ACTIVITY && uptake_profile[p-1] == '0' && Xcell[p-1] > MINIMUM_CONCENTRATION )
        {
          profile[p-1] = true;
        }
        else if (rlist->metabolic_type[i] == CATALYTIC_REWARDING_ACTIVITY && uptake_profile[p-1] == '0' && Xcell[p-1] > MINIMUM_CONCENTRATION )
        {
          profile[p-1] = true;
        }
      }
      /* Build the production profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          production_profile += "1";
        }
        else
        {
          production_profile += "0";
        }
      }
      delete[] profile;
      profile = NULL;
      
      /*-----------------------------------------------------*/
      /* 2.4) Build node profile string                      */
      /*-----------------------------------------------------*/
      trophic_profile = production_profile+uptake_profile+release_profile;
      
      /*-----------------------------------------------------*/
      /* 2.5) Evaluate if the node pre-exists in the network */
      /*-----------------------------------------------------*/
      unsigned long long int group_id;
      
      /* If the group already exists, update its statistics */
      
      if (groupExists(trophic_profile, group_id))
      {
        TrophicGroup* group = _group_map[group_id];
        group->add_alive_cell(_parameters, cell);
        cell->set_trophic_group(group_id);
        cell->set_red_color(group->get_red_color());
        cell->set_green_color(group->get_green_color());
        cell->set_blue_color(group->get_blue_color());
      }
      
      /* If the group does not exist, create it */
      
      else
      {
        double red   = _parameters->get_simulation_prng()->uniform()*(255.0-50.0)+50.0;
        double green = _parameters->get_simulation_prng()->uniform()*(255.0-50.0)+50.0;
        double blue  = _parameters->get_simulation_prng()->uniform()*(255.0-50.0)+50.0;
        TrophicGroup* group = new TrophicGroup(get_new_id(), _population->get_time(), red, green, blue);
        group->_trophic_profile    = trophic_profile;
        group->_production_profile = production_profile;
        group->_uptake_profile     = uptake_profile;
        group->_release_profile    = release_profile;
        group->add_alive_cell(_parameters, cell);
        _group_map[group->get_identifier()] = group;
        _nb_group_appearances++;
        cell->set_trophic_group(group->get_identifier());
        cell->set_red_color(group->get_red_color());
        cell->set_green_color(group->get_green_color());
        cell->set_blue_color(group->get_blue_color());
      }
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Build trophic network links        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator it1;
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator it2;
  for (it1 = _group_map.begin(); it1 != _group_map.end(); ++it1)
  {
    TrophicGroup* group_i = it1->second;
    bool level_0 = false;
    bool level_2 = false;
    for (it2 = _group_map.begin(); it2 != _group_map.end(); ++it2)
    {
      if (it1->first != it2->first)
      {
        TrophicGroup* group_j = it2->second;
        /*----------------------------------------------------------------------*/
        /* For each pump in the uptake profile of the group i, if the group j   */
        /* produces or releases the same metabolite, save it in necrophagy link */
        /* or active release link of i                                          */
        /*----------------------------------------------------------------------*/
        bool necrophagy     = false;
        bool active_release = false;
        for (size_t met = 0; met < N; met++)
        {
          /*------------------------------------------------*/
          /* evaluate the production profile of the group j */
          /*------------------------------------------------*/
          if (group_i->_uptake_profile[met] == '1' && group_j->_production_profile[met] == '1')
          {
            necrophagy = true;
            if (group_j->get_identifier() == 0)
            {
              level_0 = true;
            }
            else
            {
              level_2 = true;
            }
          }
          /*------------------------------------------------*/
          /* evaluate the release profile of the group j    */
          /*------------------------------------------------*/
          if (group_i->_uptake_profile[met] == '1' && group_j->_release_profile[met] == '1')
          {
            active_release = true;
            level_2 = true;
          }
        }
        if (necrophagy && !active_release)
        {
          group_i->add_necrophagy_link(it2->first);
        }
        else if (necrophagy && active_release)
        {
          group_i->add_active_release_link(it2->first);
        }
      }
    }
    /*----------------------------------------*/
    /* Evaluate the trophic level the group i */
    /*----------------------------------------*/
    if (!level_0 && !level_2)
    {
      group_i->set_trophic_level(NO_LEVEL);
    }
    else if (level_0 && !level_2)
    {
      group_i->set_trophic_level(LEVEL_0);
    }
    else if (level_0 && level_2)
    {
      group_i->set_trophic_level(LEVEL_1);
    }
    else if (!level_0 && level_2)
    {
      group_i->set_trophic_level(LEVEL_2);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Delete empty groups                */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  std::vector<unsigned long long int> to_remove;
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->second->get_number_of_cells() == 0 && _iterator->first != 0)
    {
      to_remove.push_back(_iterator->first);
    }
  }
  for (size_t i = 0; i < to_remove.size(); i++)
  {
    delete _group_map[to_remove[i]];
    _group_map[to_remove[i]] = NULL;
    _group_map.erase(to_remove[i]);
  }
  _nb_group_extinctions += to_remove.size();
  to_remove.clear();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 5) Update trophic network statistics  */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->first != 0)
    {
      _iterator->second->compute_mean();
      _iterator->second->update_lifespan(_population->get_time());
      _mean_group_lifespan += _iterator->second->get_lifespan();
      if (_iterator->second->get_trophic_level() == LEVEL_0)
      {
        _nb_level_0_groups++;
        _nb_level_0_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == LEVEL_1)
      {
        _nb_level_1_groups++;
        _nb_level_1_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == LEVEL_2)
      {
        _nb_level_2_groups++;
        _nb_level_2_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == NO_LEVEL)
      {
        _nb_no_level_groups++;
        _nb_no_level_cells += _iterator->second->get_number_of_cells();
      }
    }
  }
  if (_group_map.size() > 1)
  {
    _mean_group_lifespan /= (double)(_group_map.size()-1);
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 6) Update cells trophic level         */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (size_t pos = 0; pos < _population->get_width()*_population->get_height(); pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      cell->set_trophic_level(_group_map[cell->get_trophic_group()]->get_trophic_level());
    }
  }
}

/**
 * \brief    Write the trophic network in files
 * \details  --
 * \param    std::string node_filename
 * \param    std::string edge_filename
 * \return   \e void
 */
void TrophicNetwork::write_trophic_network( std::string node_filename, std::string edge_filename )
{
  /*----------------------------*/
  /* 1) Write the list of nodes */
  /*----------------------------*/
  std::ofstream nodes_file(node_filename.c_str(), std::ios::out | std::ios::trunc);
  nodes_file << "id production uptake release level count appearance lifespan\n";
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    TrophicGroup* current_group = _iterator->second;
    nodes_file << current_group->get_identifier() << " " << current_group->_production_profile << " " << current_group->_uptake_profile << " " << current_group->_release_profile << " " << current_group->get_trophic_level() << " " << current_group->get_number_of_cells() << " " << current_group->get_appearance_time() << " " << current_group->get_lifespan() << "\n";
  }
  nodes_file.close();
  
  /*----------------------------*/
  /* 2) Write the list of edges */
  /*----------------------------*/
  std::ofstream edges_file(edge_filename.c_str(), std::ios::out | std::ios::trunc);
  edges_file << "id1 id2 link_type\n";
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    TrophicGroup* current_group = _iterator->second;
    for (size_t j = 0; j < current_group->get_necrophagy_links()->size(); j++)
    {
      edges_file << current_group->get_identifier() << " " << current_group->get_necrophagy_links()->at(j) << " " << "0\n";
    }
    for (size_t j = 0; j < current_group->get_active_release_links()->size(); j++)
    {
      edges_file << current_group->get_identifier() << " " << current_group->get_active_release_links()->at(j) << " " << "1\n";
    }
  }
  edges_file.close();
}

/**
 * \brief    Write the trophic network in one file
 * \details  --
 * \param    double time
 * \param    std::ofstream& file
 * \return   \e void
 */
void TrophicNetwork::write_trophic_network( double time, std::ofstream& file )
{
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    TrophicGroup* current_group = _iterator->second;
    file << time << " " << current_group->get_identifier() << " " << current_group->_trophic_profile << " " << current_group->get_trophic_level() << " " << current_group->get_number_of_cells() << "\n";
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Check if the trophic profile is already in the trophic network
 * \details  If the group already exists, its identifier is copied in group_id
 * \param    std::string trophic_profile
 * \return   \e bool
 */
bool TrophicNetwork::groupExists( std::string trophic_profile, unsigned long long int &group_id )
{
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->second->_trophic_profile.compare(trophic_profile) == 0 && _iterator->first != 0)
    {
      group_id = _iterator->second->get_identifier();
      return true;
    }
  }
  return false;
}
