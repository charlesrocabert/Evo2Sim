
/**
 * \file      TrophicGroup.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      20-10-2015
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicGroup class definition
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

#include "TrophicGroup.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    unsigned long long int identifier
 * \param    size_t time
 * \param    double red
 * \param    double green
 * \param    double blue
 * \return   \e void
 */
TrophicGroup::TrophicGroup( unsigned long long int identifier, size_t time, double red, double green, double blue )
{
  /*------------------------------------------------------------------ trophic profiles */
  
  _trophic_profile    = "";
  _production_profile = "";
  _uptake_profile     = "";
  _release_profile    = "";
  
  /*------------------------------------------------------------------ trophic group properties */
  
  _identifier      = identifier;
  _necrophagy_links.clear();
  _active_release_links.clear();
  _appearance_time = time;
  _lifespan        = 0;
  _number_of_cells = 0;
  _trophic_level   = NO_LEVEL;
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  _mean_generations                = 0.0;
  _mean_inherited_TF_amount        = 0.0;
  _mean_inherited_E_amount         = 0.0;
  _mean_TF_amount                  = 0.0;
  _mean_E_amount                   = 0.0;
  _mean_inherited_metabolic_amount = 0.0;
  _mean_metabolic_amount           = 0.0;
  _mean_energy                     = 0.0;
  _mean_score                      = 0.0;
  _mean_lifespan                   = 0.0;
  _mean_number_of_divisions        = 0.0;
  _mean_toxicity                   = 0.0;
  _mean_metabolic_uptake           = 0.0;
  _mean_metabolic_release          = 0.0;
  _mean_metabolic_growth_rate      = 0.0;
  _mean_Dmetabolic_growth_rate     = 0.0;
  _mean_grn_nb_nodes               = 0.0;
  _mean_grn_nb_edges               = 0.0;
  _mean_metabolic_nb_nodes         = 0.0;
  _mean_metabolic_nb_edges         = 0.0;
  _mean_regulation_redundancy      = 0.0;
  _mean_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   = 0.0;
  _mean_functional_size               = 0.0;
  _mean_genome_nb_NC                  = 0.0;
  _mean_genome_nb_E                   = 0.0;
  _mean_genome_nb_TF                  = 0.0;
  _mean_genome_nb_BS                  = 0.0;
  _mean_genome_nb_P                   = 0.0;
  _mean_genome_nb_inner_enzymes       = 0.0;
  _mean_genome_nb_inflow_pumps        = 0.0;
  _mean_genome_nb_outflow_pumps       = 0.0;
  _mean_genome_nb_functional_regions  = 0.0;
  _mean_genome_nb_enhancers           = 0.0;
  _mean_genome_nb_operators           = 0.0;
  _mean_genome_nb_E_regions           = 0.0;
  _mean_genome_nb_TF_regions          = 0.0;
  _mean_genome_nb_mixed_regions       = 0.0;
  _mean_genome_functional_region_size = 0.0;
  _mean_genome_E_region_size          = 0.0;
  _mean_genome_TF_region_size         = 0.0;
  _mean_genome_mixed_region_size      = 0.0;
  _mean_genome_enhancer_size          = 0.0;
  _mean_genome_operator_size          = 0.0;
  _mean_genome_operon_size            = 0.0;
  _mean_genome_E_operon_size          = 0.0;
  _mean_genome_TF_operon_size         = 0.0;
  _mean_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _mean_inherited_size             = 0.0;
  _mean_inherited_nb_E             = 0.0;
  _mean_inherited_nb_TF            = 0.0;
  _mean_inherited_nb_inner_enzymes = 0.0;
  _mean_inherited_nb_inflow_pumps  = 0.0;
  _mean_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ trophic group color */
  
  _red_color   = red;
  _green_color = green;
  _blue_color  = blue;
  
}

/**
 * \brief    Constructor from backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
TrophicGroup::TrophicGroup( gzFile backup_file )
{
  /*------------------------------------------------------------------ trophic profiles */
  
  size_t n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _trophic_profile = "";
  for (size_t i = 0; i < n; i++)
  {
    char val;
    gzread( backup_file, &val, sizeof(val) );
    _trophic_profile += val;
  }
  
  n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _production_profile = "";
  _uptake_profile     = "";
  _release_profile    = "";
  for (size_t i = 0; i < n; i++)
  {
    char val;
    gzread( backup_file, &val, sizeof(val) );
    _production_profile += val;
  }
  for (size_t i = 0; i < n; i++)
  {
    char val;
    gzread( backup_file, &val, sizeof(val) );
    _uptake_profile += val;
  }
  for (size_t i = 0; i < n; i++)
  {
    char val;
    gzread( backup_file, &val, sizeof(val) );
    _release_profile += val;
  }
  
  /*------------------------------------------------------------------ trophic group properties */
  
  gzread( backup_file, &_identifier, sizeof(_identifier) );
  
  n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _necrophagy_links.clear();
  for (size_t i = 0; i < n; i++)
  {
    unsigned long long int current_id;
    gzread( backup_file, &current_id, sizeof(current_id) );
    _necrophagy_links.push_back(current_id);
  }
  n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _active_release_links.clear();
  for (size_t i = 0; i < n; i++)
  {
    unsigned long long int current_id;
    gzread( backup_file, &current_id, sizeof(current_id) );
    _active_release_links.push_back(current_id);
  }
  
  gzread( backup_file, &_appearance_time, sizeof(_appearance_time) );
  gzread( backup_file, &_lifespan,        sizeof(_lifespan) );
  gzread( backup_file, &_number_of_cells, sizeof(_number_of_cells) );
  gzread( backup_file, &_trophic_level,   sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  gzread( backup_file, &_mean_generations,                sizeof(_mean_generations) );
  gzread( backup_file, &_mean_inherited_TF_amount,        sizeof(_mean_inherited_TF_amount) );
  gzread( backup_file, &_mean_inherited_E_amount,         sizeof(_mean_inherited_E_amount) );
  gzread( backup_file, &_mean_TF_amount,                  sizeof(_mean_TF_amount) );
  gzread( backup_file, &_mean_E_amount,                   sizeof(_mean_E_amount) );
  gzread( backup_file, &_mean_inherited_metabolic_amount, sizeof(_mean_inherited_metabolic_amount) );
  gzread( backup_file, &_mean_metabolic_amount,           sizeof(_mean_metabolic_amount) );
  gzread( backup_file, &_mean_energy,                     sizeof(_mean_energy) );
  gzread( backup_file, &_mean_score,                      sizeof(_mean_score) );
  gzread( backup_file, &_mean_lifespan,                   sizeof(_mean_lifespan) );
  gzread( backup_file, &_mean_number_of_divisions,        sizeof(_mean_number_of_divisions) );
  gzread( backup_file, &_mean_toxicity,                   sizeof(_mean_toxicity) );
  gzread( backup_file, &_mean_metabolic_uptake,           sizeof(_mean_metabolic_uptake) );
  gzread( backup_file, &_mean_metabolic_release,          sizeof(_mean_metabolic_release) );
  gzread( backup_file, &_mean_metabolic_growth_rate,      sizeof(_mean_metabolic_growth_rate) );
  gzread( backup_file, &_mean_Dmetabolic_growth_rate,     sizeof(_mean_Dmetabolic_growth_rate) );
  gzread( backup_file, &_mean_grn_nb_nodes,               sizeof(_mean_grn_nb_nodes) );
  gzread( backup_file, &_mean_grn_nb_edges,               sizeof(_mean_grn_nb_edges) );
  gzread( backup_file, &_mean_metabolic_nb_nodes,         sizeof(_mean_metabolic_nb_nodes) );
  gzread( backup_file, &_mean_metabolic_nb_edges,         sizeof(_mean_metabolic_nb_edges) );
  gzread( backup_file, &_mean_regulation_redundancy,      sizeof(_mean_regulation_redundancy) );
  gzread( backup_file, &_mean_metabolic_redundancy,       sizeof(_mean_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzread( backup_file, &_mean_genome_size,                   sizeof(_mean_genome_size) );
  gzread( backup_file, &_mean_functional_size,               sizeof(_mean_functional_size) );
  gzread( backup_file, &_mean_genome_nb_NC,                  sizeof(_mean_genome_nb_NC) );
  gzread( backup_file, &_mean_genome_nb_E,                   sizeof(_mean_genome_nb_E) );
  gzread( backup_file, &_mean_genome_nb_TF,                  sizeof(_mean_genome_nb_TF) );
  gzread( backup_file, &_mean_genome_nb_BS,                  sizeof(_mean_genome_nb_BS) );
  gzread( backup_file, &_mean_genome_nb_P,                   sizeof(_mean_genome_nb_P) );
  gzread( backup_file, &_mean_genome_nb_inner_enzymes,       sizeof(_mean_genome_nb_inner_enzymes) );
  gzread( backup_file, &_mean_genome_nb_inflow_pumps,        sizeof(_mean_genome_nb_inflow_pumps) );
  gzread( backup_file, &_mean_genome_nb_outflow_pumps,       sizeof(_mean_genome_nb_outflow_pumps) );
  gzread( backup_file, &_mean_genome_nb_functional_regions,  sizeof(_mean_genome_nb_functional_regions) );
  gzread( backup_file, &_mean_genome_nb_enhancers,           sizeof(_mean_genome_nb_enhancers) );
  gzread( backup_file, &_mean_genome_nb_operators,           sizeof(_mean_genome_nb_operators) );
  gzread( backup_file, &_mean_genome_nb_E_regions,           sizeof(_mean_genome_nb_E_regions) );
  gzread( backup_file, &_mean_genome_nb_TF_regions,          sizeof(_mean_genome_nb_TF_regions) );
  gzread( backup_file, &_mean_genome_nb_mixed_regions,       sizeof(_mean_genome_nb_mixed_regions) );
  gzread( backup_file, &_mean_genome_functional_region_size, sizeof(_mean_genome_functional_region_size) );
  gzread( backup_file, &_mean_genome_E_region_size,          sizeof(_mean_genome_E_region_size) );
  gzread( backup_file, &_mean_genome_TF_region_size,         sizeof(_mean_genome_TF_region_size) );
  gzread( backup_file, &_mean_genome_mixed_region_size,      sizeof(_mean_genome_mixed_region_size) );
  gzread( backup_file, &_mean_genome_enhancer_size,          sizeof(_mean_genome_enhancer_size) );
  gzread( backup_file, &_mean_genome_operator_size,          sizeof(_mean_genome_operator_size) );
  gzread( backup_file, &_mean_genome_operon_size,            sizeof(_mean_genome_operon_size) );
  gzread( backup_file, &_mean_genome_E_operon_size,          sizeof(_mean_genome_E_operon_size) );
  gzread( backup_file, &_mean_genome_TF_operon_size,         sizeof(_mean_genome_TF_operon_size) );
  gzread( backup_file, &_mean_genome_mixed_operon_size,      sizeof(_mean_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzread( backup_file, &_mean_inherited_size,             sizeof(_mean_inherited_size) );
  gzread( backup_file, &_mean_inherited_nb_E,             sizeof(_mean_inherited_nb_E) );
  gzread( backup_file, &_mean_inherited_nb_TF,            sizeof(_mean_inherited_nb_TF) );
  gzread( backup_file, &_mean_inherited_nb_inner_enzymes, sizeof(_mean_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_mean_inherited_nb_inflow_pumps,  sizeof(_mean_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_mean_inherited_nb_outflow_pumps, sizeof(_mean_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ trophic group color */
  
  gzread( backup_file, &_red_color,   sizeof(_red_color) );
  gzread( backup_file, &_green_color, sizeof(_green_color) );
  gzread( backup_file, &_blue_color,  sizeof(_blue_color) );
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const TrophicGroup& trophic_group
 * \return   \e void
 */
TrophicGroup::TrophicGroup( const TrophicGroup& trophic_group )
{
  /*------------------------------------------------------------------ trophic profiles */
  
  _trophic_profile    = trophic_group._trophic_profile;
  _production_profile = trophic_group._production_profile;
  _uptake_profile     = trophic_group._uptake_profile;
  _release_profile    = trophic_group._release_profile;
  
  /*------------------------------------------------------------------ trophic group properties */
  
  _identifier = trophic_group._identifier;
  _necrophagy_links.clear();
  for (size_t i = 0; i < trophic_group._necrophagy_links.size(); i++)
  {
    _necrophagy_links.push_back(trophic_group._necrophagy_links[i]);
  }
  _active_release_links.clear();
  for (size_t i = 0; i < trophic_group._active_release_links.size(); i++)
  {
    _active_release_links.push_back(trophic_group._active_release_links[i]);
  }
  _appearance_time = trophic_group._appearance_time;
  _lifespan        = trophic_group._lifespan;
  _number_of_cells = trophic_group._number_of_cells;
  _trophic_level   = trophic_group._trophic_level;
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  _mean_generations                = trophic_group._mean_generations;
  _mean_inherited_TF_amount        = trophic_group._mean_inherited_TF_amount;
  _mean_inherited_E_amount         = trophic_group._mean_inherited_E_amount;
  _mean_TF_amount                  = trophic_group._mean_TF_amount;
  _mean_E_amount                   = trophic_group._mean_E_amount;
  _mean_inherited_metabolic_amount = trophic_group._mean_inherited_metabolic_amount;
  _mean_metabolic_amount           = trophic_group._mean_metabolic_amount;
  _mean_energy                     = trophic_group._mean_energy;
  _mean_score                      = trophic_group._mean_score;
  _mean_lifespan                   = trophic_group._mean_lifespan;
  _mean_number_of_divisions        = trophic_group._mean_number_of_divisions;
  _mean_toxicity                   = trophic_group._mean_toxicity;
  _mean_metabolic_uptake           = trophic_group._mean_metabolic_uptake;
  _mean_metabolic_release          = trophic_group._mean_metabolic_release;
  _mean_metabolic_growth_rate      = trophic_group._mean_metabolic_growth_rate;
  _mean_Dmetabolic_growth_rate     = trophic_group._mean_Dmetabolic_growth_rate;
  _mean_grn_nb_nodes               = trophic_group._mean_grn_nb_nodes;
  _mean_grn_nb_edges               = trophic_group._mean_grn_nb_edges;
  _mean_metabolic_nb_nodes         = trophic_group._mean_metabolic_nb_nodes;
  _mean_metabolic_nb_edges         = trophic_group._mean_metabolic_nb_edges;
  _mean_regulation_redundancy      = trophic_group._mean_regulation_redundancy;
  _mean_metabolic_redundancy       = trophic_group._mean_metabolic_redundancy;
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   = trophic_group._mean_genome_size;
  _mean_functional_size               = trophic_group._mean_functional_size;
  _mean_genome_nb_NC                  = trophic_group._mean_genome_nb_NC;
  _mean_genome_nb_E                   = trophic_group._mean_genome_nb_E;
  _mean_genome_nb_TF                  = trophic_group._mean_genome_nb_TF;
  _mean_genome_nb_BS                  = trophic_group._mean_genome_nb_BS;
  _mean_genome_nb_P                   = trophic_group._mean_genome_nb_P;
  _mean_genome_nb_inner_enzymes       = trophic_group._mean_genome_nb_inner_enzymes;
  _mean_genome_nb_inflow_pumps        = trophic_group._mean_genome_nb_inflow_pumps;
  _mean_genome_nb_outflow_pumps       = trophic_group._mean_genome_nb_outflow_pumps;
  _mean_genome_nb_functional_regions  = trophic_group._mean_genome_nb_functional_regions;
  _mean_genome_nb_enhancers           = trophic_group._mean_genome_nb_enhancers;
  _mean_genome_nb_operators           = trophic_group._mean_genome_nb_operators;
  _mean_genome_nb_E_regions           = trophic_group._mean_genome_nb_E_regions;
  _mean_genome_nb_TF_regions          = trophic_group._mean_genome_nb_TF_regions;
  _mean_genome_nb_mixed_regions       = trophic_group._mean_genome_nb_mixed_regions;
  _mean_genome_functional_region_size = trophic_group._mean_genome_functional_region_size;
  _mean_genome_E_region_size          = trophic_group._mean_genome_E_region_size;
  _mean_genome_TF_region_size         = trophic_group._mean_genome_TF_region_size;
  _mean_genome_mixed_region_size      = trophic_group._mean_genome_mixed_region_size;
  _mean_genome_enhancer_size          = trophic_group._mean_genome_enhancer_size;
  _mean_genome_operator_size          = trophic_group._mean_genome_operator_size;
  _mean_genome_operon_size            = trophic_group._mean_genome_operon_size;
  _mean_genome_E_operon_size          = trophic_group._mean_genome_E_operon_size;
  _mean_genome_TF_operon_size         = trophic_group._mean_genome_TF_operon_size;
  _mean_genome_mixed_operon_size      = trophic_group._mean_genome_mixed_operon_size;
  
  /* INHERITED STRUCTURE */
  _mean_inherited_size             = trophic_group._mean_inherited_size;
  _mean_inherited_nb_E             = trophic_group._mean_inherited_nb_E;
  _mean_inherited_nb_TF            = trophic_group._mean_inherited_nb_TF;
  _mean_inherited_nb_inner_enzymes = trophic_group._mean_inherited_nb_inner_enzymes;
  _mean_inherited_nb_inflow_pumps  = trophic_group._mean_inherited_nb_inflow_pumps;
  _mean_inherited_nb_outflow_pumps = trophic_group._mean_inherited_nb_outflow_pumps;
  
  /*------------------------------------------------------------------ trophic group color */
  
  _red_color   = trophic_group._red_color;
  _green_color = trophic_group._green_color;
  _blue_color  = trophic_group._blue_color;
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
TrophicGroup::~TrophicGroup( void )
{
  _necrophagy_links.clear();
  _active_release_links.clear();
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
void TrophicGroup::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ trophic profiles */
  
  size_t n = _trophic_profile.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    char val = _trophic_profile[i];
    gzwrite( backup_file, &val, sizeof(val) );
  }
  
  n = _production_profile.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    char val = _production_profile[i];
    gzwrite( backup_file, &val, sizeof(val) );
  }
  for (size_t i = 0; i < n; i++)
  {
    char val = _uptake_profile[i];
    gzwrite( backup_file, &val, sizeof(val) );
  }
  for (size_t i = 0; i < n; i++)
  {
    char val = _release_profile[i];
    gzwrite( backup_file, &val, sizeof(val) );
  }
  
  /*------------------------------------------------------------------ trophic group properties */
  
  gzwrite( backup_file, &_identifier, sizeof(_identifier) );
  
  n = _necrophagy_links.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    gzwrite( backup_file, &_necrophagy_links[i], sizeof(_necrophagy_links[i]) );
  }
  n = _active_release_links.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    gzwrite( backup_file, &_active_release_links[i], sizeof(_active_release_links[i]) );
  }
  
  gzwrite( backup_file, &_appearance_time, sizeof(_appearance_time) );
  gzwrite( backup_file, &_lifespan,        sizeof(_lifespan) );
  gzwrite( backup_file, &_number_of_cells, sizeof(_number_of_cells) );
  gzwrite( backup_file, &_trophic_level,   sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  gzwrite( backup_file, &_mean_generations,                sizeof(_mean_generations) );
  gzwrite( backup_file, &_mean_inherited_TF_amount,        sizeof(_mean_inherited_TF_amount) );
  gzwrite( backup_file, &_mean_inherited_E_amount,         sizeof(_mean_inherited_E_amount) );
  gzwrite( backup_file, &_mean_TF_amount,                  sizeof(_mean_TF_amount) );
  gzwrite( backup_file, &_mean_E_amount,                   sizeof(_mean_E_amount) );
  gzwrite( backup_file, &_mean_inherited_metabolic_amount, sizeof(_mean_inherited_metabolic_amount) );
  gzwrite( backup_file, &_mean_metabolic_amount,           sizeof(_mean_metabolic_amount) );
  gzwrite( backup_file, &_mean_energy,                     sizeof(_mean_energy) );
  gzwrite( backup_file, &_mean_score,                      sizeof(_mean_score) );
  gzwrite( backup_file, &_mean_lifespan,                   sizeof(_mean_lifespan) );
  gzwrite( backup_file, &_mean_number_of_divisions,        sizeof(_mean_number_of_divisions) );
  gzwrite( backup_file, &_mean_toxicity,                   sizeof(_mean_toxicity) );
  gzwrite( backup_file, &_mean_metabolic_uptake,           sizeof(_mean_metabolic_uptake) );
  gzwrite( backup_file, &_mean_metabolic_release,          sizeof(_mean_metabolic_release) );
  gzwrite( backup_file, &_mean_metabolic_growth_rate,      sizeof(_mean_metabolic_growth_rate) );
  gzwrite( backup_file, &_mean_Dmetabolic_growth_rate,     sizeof(_mean_Dmetabolic_growth_rate) );
  gzwrite( backup_file, &_mean_grn_nb_nodes,               sizeof(_mean_grn_nb_nodes) );
  gzwrite( backup_file, &_mean_grn_nb_edges,               sizeof(_mean_grn_nb_edges) );
  gzwrite( backup_file, &_mean_metabolic_nb_nodes,         sizeof(_mean_metabolic_nb_nodes) );
  gzwrite( backup_file, &_mean_metabolic_nb_edges,         sizeof(_mean_metabolic_nb_edges) );
  gzwrite( backup_file, &_mean_regulation_redundancy,      sizeof(_mean_regulation_redundancy) );
  gzwrite( backup_file, &_mean_metabolic_redundancy,       sizeof(_mean_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzwrite( backup_file, &_mean_genome_size,                   sizeof(_mean_genome_size) );
  gzwrite( backup_file, &_mean_functional_size,               sizeof(_mean_functional_size) );
  gzwrite( backup_file, &_mean_genome_nb_NC,                  sizeof(_mean_genome_nb_NC) );
  gzwrite( backup_file, &_mean_genome_nb_E,                   sizeof(_mean_genome_nb_E) );
  gzwrite( backup_file, &_mean_genome_nb_TF,                  sizeof(_mean_genome_nb_TF) );
  gzwrite( backup_file, &_mean_genome_nb_BS,                  sizeof(_mean_genome_nb_BS) );
  gzwrite( backup_file, &_mean_genome_nb_P,                   sizeof(_mean_genome_nb_P) );
  gzwrite( backup_file, &_mean_genome_nb_inner_enzymes,       sizeof(_mean_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_mean_genome_nb_inflow_pumps,        sizeof(_mean_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_mean_genome_nb_outflow_pumps,       sizeof(_mean_genome_nb_outflow_pumps) );
  gzwrite( backup_file, &_mean_genome_nb_functional_regions,  sizeof(_mean_genome_nb_functional_regions) );
  gzwrite( backup_file, &_mean_genome_nb_enhancers,           sizeof(_mean_genome_nb_enhancers) );
  gzwrite( backup_file, &_mean_genome_nb_operators,           sizeof(_mean_genome_nb_operators) );
  gzwrite( backup_file, &_mean_genome_nb_E_regions,           sizeof(_mean_genome_nb_E_regions) );
  gzwrite( backup_file, &_mean_genome_nb_TF_regions,          sizeof(_mean_genome_nb_TF_regions) );
  gzwrite( backup_file, &_mean_genome_nb_mixed_regions,       sizeof(_mean_genome_nb_mixed_regions) );
  gzwrite( backup_file, &_mean_genome_functional_region_size, sizeof(_mean_genome_functional_region_size) );
  gzwrite( backup_file, &_mean_genome_E_region_size,          sizeof(_mean_genome_E_region_size) );
  gzwrite( backup_file, &_mean_genome_TF_region_size,         sizeof(_mean_genome_TF_region_size) );
  gzwrite( backup_file, &_mean_genome_mixed_region_size,      sizeof(_mean_genome_mixed_region_size) );
  gzwrite( backup_file, &_mean_genome_enhancer_size,          sizeof(_mean_genome_enhancer_size) );
  gzwrite( backup_file, &_mean_genome_operator_size,          sizeof(_mean_genome_operator_size) );
  gzwrite( backup_file, &_mean_genome_operon_size,            sizeof(_mean_genome_operon_size) );
  gzwrite( backup_file, &_mean_genome_E_operon_size,          sizeof(_mean_genome_E_operon_size) );
  gzwrite( backup_file, &_mean_genome_TF_operon_size,         sizeof(_mean_genome_TF_operon_size) );
  gzwrite( backup_file, &_mean_genome_mixed_operon_size,      sizeof(_mean_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzwrite( backup_file, &_mean_inherited_size,             sizeof(_mean_inherited_size) );
  gzwrite( backup_file, &_mean_inherited_nb_E,             sizeof(_mean_inherited_nb_E) );
  gzwrite( backup_file, &_mean_inherited_nb_TF,            sizeof(_mean_inherited_nb_TF) );
  gzwrite( backup_file, &_mean_inherited_nb_inner_enzymes, sizeof(_mean_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_mean_inherited_nb_inflow_pumps,  sizeof(_mean_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_mean_inherited_nb_outflow_pumps, sizeof(_mean_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ trophic group color */
  
  gzwrite( backup_file, &_red_color,   sizeof(_red_color) );
  gzwrite( backup_file, &_green_color, sizeof(_green_color) );
  gzwrite( backup_file, &_blue_color,  sizeof(_blue_color) );
}

/**
 * \brief    Add an alive cell to the statistics
 * \details  --
 * \param    Parameters* parameters
 * \param    Cell* cell
 * \return   \e void
 */
void TrophicGroup::add_alive_cell( Parameters* parameters, Cell* cell )
{
  /*------------------------------------------------------------------ trophic group properties */
  
  _number_of_cells++;
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  _mean_generations                += cell->get_generation();
  _mean_inherited_TF_amount        += cell->get_inherited_TF_amount();
  _mean_inherited_E_amount         += cell->get_inherited_E_amount();
  _mean_TF_amount                  += cell->get_TF_amount();
  _mean_E_amount                   += cell->get_E_amount();
  _mean_inherited_metabolic_amount += cell->get_inherited_metabolic_amount();
  _mean_metabolic_amount           += cell->get_species_list()->get_amount();
  _mean_energy                     += cell->get_energy();
  _mean_score                      += cell->get_score();
  _mean_lifespan                   += cell->get_lifespan();
  _mean_number_of_divisions        += cell->get_number_of_divisions();
  _mean_toxicity                   += cell->get_toxicity();
  _mean_metabolic_uptake           += cell->get_metabolic_uptake();
  _mean_metabolic_release          += cell->get_metabolic_release();
  _mean_metabolic_growth_rate      += cell->get_metabolic_growth_rate();
  _mean_Dmetabolic_growth_rate     += cell->get_Dmetabolic_growth_rate();
  _mean_grn_nb_nodes               += cell->get_grn_nb_nodes();
  _mean_grn_nb_edges               += cell->get_grn_nb_edges();
  _mean_metabolic_nb_nodes         += cell->get_metabolic_nb_nodes();
  _mean_metabolic_nb_edges         += cell->get_metabolic_nb_edges();
  _mean_regulation_redundancy      += cell->get_replication_report()->get_mean_regulation_redundancy();
  _mean_metabolic_redundancy       += cell->get_replication_report()->get_mean_metabolic_redundancy();
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   += cell->get_genome()->get_size();
  _mean_functional_size               += cell->get_replication_report()->get_genome_functional_size();
  _mean_genome_nb_NC                  += cell->get_genome()->get_nb_NC();
  _mean_genome_nb_E                   += cell->get_genome()->get_nb_E();
  _mean_genome_nb_TF                  += cell->get_genome()->get_nb_TF();
  _mean_genome_nb_BS                  += cell->get_genome()->get_nb_BS();
  _mean_genome_nb_P                   += cell->get_genome()->get_nb_P();
  _mean_genome_nb_inner_enzymes       += cell->get_genome()->get_nb_inner_enzymes();
  _mean_genome_nb_inflow_pumps        += cell->get_genome()->get_nb_inflow_pumps();
  _mean_genome_nb_outflow_pumps       += cell->get_genome()->get_nb_outflow_pumps();
  _mean_genome_nb_functional_regions  += cell->get_replication_report()->get_nb_functional_regions();
  _mean_genome_nb_enhancers           += cell->get_replication_report()->get_nb_enhancers();
  _mean_genome_nb_operators           += cell->get_replication_report()->get_nb_operators();
  _mean_genome_nb_E_regions           += cell->get_replication_report()->get_nb_E_regions();
  _mean_genome_nb_TF_regions          += cell->get_replication_report()->get_nb_TF_regions();
  _mean_genome_nb_mixed_regions       += cell->get_replication_report()->get_nb_mixed_regions();
  _mean_genome_functional_region_size += cell->get_replication_report()->get_mean_functional_region_size();
  _mean_genome_E_region_size          += cell->get_replication_report()->get_mean_E_region_size();
  _mean_genome_TF_region_size         += cell->get_replication_report()->get_mean_TF_region_size();
  _mean_genome_mixed_region_size      += cell->get_replication_report()->get_mean_mixed_region_size();
  _mean_genome_enhancer_size          += cell->get_replication_report()->get_mean_enhancer_size();
  _mean_genome_operator_size          += cell->get_replication_report()->get_mean_operator_size();
  _mean_genome_operon_size            += cell->get_replication_report()->get_mean_operon_size();
  _mean_genome_E_operon_size          += cell->get_replication_report()->get_mean_E_operon_size();
  _mean_genome_TF_operon_size         += cell->get_replication_report()->get_mean_TF_operon_size();
  _mean_genome_mixed_operon_size      += cell->get_replication_report()->get_mean_mixed_operon_size();
  
  /* INHERITED STRUCTURE */
  if (parameters->get_enzymatic_inheritance_scheme())
  {
    _mean_inherited_size             += cell->get_inherited_proteins()->get_size();
    _mean_inherited_nb_E             += cell->get_inherited_proteins()->get_nb_E();
    _mean_inherited_nb_TF            += cell->get_inherited_proteins()->get_nb_TF();
    _mean_inherited_nb_inner_enzymes += cell->get_inherited_proteins()->get_nb_inner_enzymes();
    _mean_inherited_nb_inflow_pumps  += cell->get_inherited_proteins()->get_nb_inflow_pumps();
    _mean_inherited_nb_outflow_pumps += cell->get_inherited_proteins()->get_nb_outflow_pumps();
  }
}

/**
 * \brief    Compute statistics mean
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicGroup::compute_mean( void )
{
  if (_number_of_cells > 0)
  {
    /* PHENOTYPE */
    _mean_generations                /= _number_of_cells;
    _mean_inherited_TF_amount        /= _number_of_cells;
    _mean_inherited_E_amount         /= _number_of_cells;
    _mean_TF_amount                  /= _number_of_cells;
    _mean_E_amount                   /= _number_of_cells;
    _mean_inherited_metabolic_amount /= _number_of_cells;
    _mean_metabolic_amount           /= _number_of_cells;
    _mean_energy                     /= _number_of_cells;
    _mean_score                      /= _number_of_cells;
    _mean_lifespan                   /= _number_of_cells;
    _mean_number_of_divisions        /= _number_of_cells;
    _mean_toxicity                   /= _number_of_cells;
    _mean_metabolic_uptake           /= _number_of_cells;
    _mean_metabolic_release          /= _number_of_cells;
    _mean_metabolic_growth_rate      /= _number_of_cells;
    _mean_Dmetabolic_growth_rate     /= _number_of_cells;
    _mean_grn_nb_nodes               /= _number_of_cells;
    _mean_grn_nb_edges               /= _number_of_cells;
    _mean_metabolic_nb_nodes         /= _number_of_cells;
    _mean_metabolic_nb_edges         /= _number_of_cells;
    _mean_regulation_redundancy      /= _number_of_cells;
    _mean_metabolic_redundancy       /= _number_of_cells;
    
    /* GENOME STRUCTURE */
    _mean_genome_size                   /= _number_of_cells;
    _mean_functional_size               /= _number_of_cells;
    _mean_genome_nb_NC                  /= _number_of_cells;
    _mean_genome_nb_E                   /= _number_of_cells;
    _mean_genome_nb_TF                  /= _number_of_cells;
    _mean_genome_nb_BS                  /= _number_of_cells;
    _mean_genome_nb_P                   /= _number_of_cells;
    _mean_genome_nb_inner_enzymes       /= _number_of_cells;
    _mean_genome_nb_inflow_pumps        /= _number_of_cells;
    _mean_genome_nb_outflow_pumps       /= _number_of_cells;
    _mean_genome_nb_functional_regions  /= _number_of_cells;
    _mean_genome_nb_enhancers           /= _number_of_cells;
    _mean_genome_nb_operators           /= _number_of_cells;
    _mean_genome_nb_E_regions           /= _number_of_cells;
    _mean_genome_nb_TF_regions          /= _number_of_cells;
    _mean_genome_nb_mixed_regions       /= _number_of_cells;
    _mean_genome_functional_region_size /= _number_of_cells;
    _mean_genome_E_region_size          /= _number_of_cells;
    _mean_genome_TF_region_size         /= _number_of_cells;
    _mean_genome_mixed_region_size      /= _number_of_cells;
    _mean_genome_enhancer_size          /= _number_of_cells;
    _mean_genome_operator_size          /= _number_of_cells;
    _mean_genome_operon_size            /= _number_of_cells;
    _mean_genome_E_operon_size          /= _number_of_cells;
    _mean_genome_TF_operon_size         /= _number_of_cells;
    _mean_genome_mixed_operon_size      /= _number_of_cells;
    
    /* INHERITED STRUCTURE */
    _mean_inherited_size             /= _number_of_cells;
    _mean_inherited_nb_E             /= _number_of_cells;
    _mean_inherited_nb_TF            /= _number_of_cells;
    _mean_inherited_nb_inner_enzymes /= _number_of_cells;
    _mean_inherited_nb_inflow_pumps  /= _number_of_cells;
    _mean_inherited_nb_outflow_pumps /= _number_of_cells;
  }
}

/**
 * \brief    Clear statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicGroup::clear( void )
{
  /*------------------------------------------------------------------ trophic group properties */
  
  _necrophagy_links.clear();
  _active_release_links.clear();
  _number_of_cells = 0;
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  _mean_generations                = 0.0;
  _mean_inherited_TF_amount        = 0.0;
  _mean_inherited_E_amount         = 0.0;
  _mean_TF_amount                  = 0.0;
  _mean_E_amount                   = 0.0;
  _mean_inherited_metabolic_amount = 0.0;
  _mean_metabolic_amount           = 0.0;
  _mean_energy                     = 0.0;
  _mean_score                      = 0.0;
  _mean_lifespan                   = 0.0;
  _mean_number_of_divisions        = 0.0;
  _mean_toxicity                   = 0.0;
  _mean_metabolic_uptake           = 0.0;
  _mean_metabolic_release          = 0.0;
  _mean_metabolic_growth_rate      = 0.0;
  _mean_Dmetabolic_growth_rate     = 0.0;
  _mean_grn_nb_nodes               = 0.0;
  _mean_grn_nb_edges               = 0.0;
  _mean_metabolic_nb_nodes         = 0.0;
  _mean_metabolic_nb_edges         = 0.0;
  _mean_regulation_redundancy      = 0.0;
  _mean_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   = 0.0;
  _mean_functional_size               = 0.0;
  _mean_genome_nb_NC                  = 0.0;
  _mean_genome_nb_E                   = 0.0;
  _mean_genome_nb_TF                  = 0.0;
  _mean_genome_nb_BS                  = 0.0;
  _mean_genome_nb_P                   = 0.0;
  _mean_genome_nb_inner_enzymes       = 0.0;
  _mean_genome_nb_inflow_pumps        = 0.0;
  _mean_genome_nb_outflow_pumps       = 0.0;
  _mean_genome_nb_functional_regions  = 0.0;
  _mean_genome_nb_enhancers           = 0.0;
  _mean_genome_nb_operators           = 0.0;
  _mean_genome_nb_E_regions           = 0.0;
  _mean_genome_nb_TF_regions          = 0.0;
  _mean_genome_nb_mixed_regions       = 0.0;
  _mean_genome_functional_region_size = 0.0;
  _mean_genome_E_region_size          = 0.0;
  _mean_genome_TF_region_size         = 0.0;
  _mean_genome_mixed_region_size      = 0.0;
  _mean_genome_enhancer_size          = 0.0;
  _mean_genome_operator_size          = 0.0;
  _mean_genome_operon_size            = 0.0;
  _mean_genome_E_operon_size          = 0.0;
  _mean_genome_TF_operon_size         = 0.0;
  _mean_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _mean_inherited_size             = 0.0;
  _mean_inherited_nb_E             = 0.0;
  _mean_inherited_nb_TF            = 0.0;
  _mean_inherited_nb_inner_enzymes = 0.0;
  _mean_inherited_nb_inflow_pumps  = 0.0;
  _mean_inherited_nb_outflow_pumps = 0.0;
}

/**
 * \brief    Update profile thanks to the new profile length
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicGroup::update_profile( size_t N )
{
  assert(_production_profile.size() == _uptake_profile.size());
  assert(_uptake_profile.size() == _release_profile.size());
  for (size_t i = _production_profile.size(); i < N; i++)
  {
    _production_profile += "0";
    _uptake_profile += "0";
    _release_profile += "0";
  }
  assert(_production_profile.size() == N);
  assert(_uptake_profile.size() == N);
  assert(_release_profile.size() == N);
  _trophic_profile = _production_profile+_uptake_profile+_release_profile;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/




