
/**
 * \file      ReplicationReport.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     ReplicationReport class definition
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

#include "ReplicationReport.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
ReplicationReport::ReplicationReport( void )
{
  /*------------------------------------------------------------------ Genome structure */
  
  _old_genome_size         = 0;
  _new_genome_size         = 0;
  _genome_functional_size  = 0;
  _genome_nb_NC            = 0;
  _genome_nb_E             = 0;
  _genome_nb_TF            = 0;
  _genome_nb_BS            = 0;
  _genome_nb_P             = 0;
  _genome_nb_inner_enzymes = 0;
  _genome_nb_inflow_pumps  = 0;
  _genome_nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ GRN data */
  
  _nb_functional_regions       = 0;
  _nb_enhancers                = 0;
  _nb_operators                = 0;
  _nb_E_regions                = 0;
  _nb_TF_regions               = 0;
  _nb_mixed_regions            = 0;
  _mean_functional_region_size = 0.0;
  _mean_E_region_size          = 0.0;
  _mean_TF_region_size         = 0.0;
  _mean_mixed_region_size      = 0.0;
  _mean_enhancer_size          = 0.0;
  _mean_operator_size          = 0.0;
  _mean_operon_size            = 0.0;
  _mean_E_operon_size          = 0.0;
  _mean_TF_operon_size         = 0.0;
  _mean_mixed_operon_size      = 0.0;
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  _mean_regulation_redundancy = 0.0;
  _mean_metabolic_redundancy  = 0.0;
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  _inherited_size             = 0;
  _inherited_nb_E             = 0;
  _inherited_nb_TF            = 0;
  _inherited_nb_inner_enzymes = 0;
  _inherited_nb_inflow_pumps  = 0;
  _inherited_nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ Phenotype */
  
  _id                         = 0;
  _parent_id                  = 0;
  _generation                 = 0;
  _x                          = 0;
  _y                          = 0;
  _number_of_updates          = 0;
  _number_of_divisions        = 0;
  _birth_time                 = 0.0;
  _death_time                 = 0.0;
  _lifespan                   = 0.0;
  _toxicity                   = 0.0;
  _inherited_TF_amount        = 0.0;
  _inherited_E_amount         = 0.0;
  _TF_amount                  = 0.0;
  _E_amount                   = 0.0;
  _inherited_metabolic_amount = 0.0;
  _min_metabolic_amount       = 0.0;
  _metabolic_amount           = 0.0;
  _max_metabolic_amount       = 0.0;
  _metabolic_uptake           = 0.0;
  _metabolic_release          = 0.0;
  _previous_metabolic_amount  = 0.0;
  _min_energy                 = 0.0;
  _mean_energy                = 0.0;
  _max_energy                 = 0.0;
  _min_score                  = 0.0;
  _mean_score                 = 0.0;
  _max_score                  = 0.0;
  _metabolic_growth_rate      = 0.0;
  _diff_metabolic_growth_rate = 0.0;
  _grn_nb_nodes               = 0;
  _grn_nb_edges               = 0;
  _metabolic_nb_nodes         = 0;
  _metabolic_nb_edges         = 0;
  _trophic_group              = 0;
  _trophic_level              = NO_LEVEL;
  
  /*------------------------------------------------------------------ List of mutation events */
  
  _list_of_events.clear();
  
  /*------------------------------------------------------------------ Point mutations data */
  
  _nb_point_mutations = 0;
  
  _nb_NC_point_mutations = 0;
  _nb_E_point_mutations  = 0;
  _nb_TF_point_mutations = 0;
  _nb_BS_point_mutations = 0;
  _nb_P_point_mutations  = 0;
  
  _nb_NC_to_E_transitions  = 0;
  _nb_NC_to_TF_transitions = 0;
  _nb_NC_to_BS_transitions = 0;
  _nb_NC_to_P_transitions  = 0;
  
  _nb_E_to_NC_transitions = 0;
  _nb_E_to_TF_transitions = 0;
  _nb_E_to_BS_transitions = 0;
  _nb_E_to_P_transitions  = 0;
  
  _nb_TF_to_NC_transitions = 0;
  _nb_TF_to_E_transitions  = 0;
  _nb_TF_to_BS_transitions = 0;
  _nb_TF_to_P_transitions  = 0;
  
  _nb_BS_to_NC_transitions = 0;
  _nb_BS_to_E_transitions  = 0;
  _nb_BS_to_TF_transitions = 0;
  _nb_BS_to_P_transitions  = 0;
  
  _nb_P_to_NC_transitions = 0;
  _nb_P_to_E_transitions  = 0;
  _nb_P_to_TF_transitions = 0;
  _nb_P_to_BS_transitions = 0;
  
  _mean_s_mutation_size             = 0.0;
  _mean_p_mutation_size             = 0.0;
  _mean_kcat_mutation_size          = 0.0;
  _mean_kcat_km_ratio_mutation_size = 0.0;
  
  _mean_BS_tag_mutation_size  = 0.0;
  _mean_coE_tag_mutation_size = 0.0;
  
  _mean_TF_tag_mutation_size = 0.0;
  
  _mean_basal_expression_level_mutation_size = 0.0;
  
  /*------------------------------------------------------------------ HGT data */
  
  _nb_HGT        = 0;
  _mean_HGT_size = 0.0;
  _nb_NC_HGT     = 0;
  _nb_E_HGT      = 0;
  _nb_TF_HGT     = 0;
  _nb_BS_HGT     = 0;
  _nb_P_HGT      = 0;
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  _nb_rearrangements = 0;
  
  _nb_duplicated_NC = 0;
  _nb_duplicated_E  = 0;
  _nb_duplicated_TF = 0;
  _nb_duplicated_BS = 0;
  _nb_duplicated_P  = 0;
  
  _nb_deleted_NC = 0;
  _nb_deleted_E  = 0;
  _nb_deleted_TF = 0;
  _nb_deleted_BS = 0;
  _nb_deleted_P  = 0;
  
  _nb_duplications   = 0;
  _nb_deletions      = 0;
  _nb_translocations = 0;
  _nb_inversions     = 0;
  
  _mean_rearrangement_size = 0.0;
  _mean_duplication_size   = 0.0;
  _mean_deletion_size      = 0.0;
  _mean_translocation_size = 0.0;
  _mean_inversion_size     = 0.0;
  
}

/**
 * \brief    Constructor from backup file
 * \details  Load ReplicationReport class from backup file
 * \param    gzFile backup_file
 * \return   \e void
 */
ReplicationReport::ReplicationReport( gzFile backup_file )
{
  /*------------------------------------------------------------------ Genome structure */
  
  gzread( backup_file, &_old_genome_size,         sizeof(_old_genome_size) );
  gzread( backup_file, &_new_genome_size,         sizeof(_new_genome_size) );
  gzread( backup_file, &_genome_functional_size,  sizeof(_genome_functional_size) );
  gzread( backup_file, &_genome_nb_NC,            sizeof(_genome_nb_NC) );
  gzread( backup_file, &_genome_nb_E,             sizeof(_genome_nb_E) );
  gzread( backup_file, &_genome_nb_TF,            sizeof(_genome_nb_TF) );
  gzread( backup_file, &_genome_nb_BS,            sizeof(_genome_nb_BS) );
  gzread( backup_file, &_genome_nb_P,             sizeof(_genome_nb_P) );
  gzread( backup_file, &_genome_nb_inner_enzymes, sizeof(_genome_nb_inner_enzymes) );
  gzread( backup_file, &_genome_nb_inflow_pumps,  sizeof(_genome_nb_inflow_pumps) );
  gzread( backup_file, &_genome_nb_outflow_pumps, sizeof(_genome_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ GRN data */
  
  gzread( backup_file, &_nb_functional_regions,       sizeof(_nb_functional_regions) );
  gzread( backup_file, &_nb_enhancers,                sizeof(_nb_enhancers) );
  gzread( backup_file, &_nb_operators,                sizeof(_nb_operators) );
  gzread( backup_file, &_nb_E_regions,                sizeof(_nb_E_regions) );
  gzread( backup_file, &_nb_TF_regions,               sizeof(_nb_TF_regions) );
  gzread( backup_file, &_nb_mixed_regions,            sizeof(_nb_mixed_regions) );
  gzread( backup_file, &_mean_functional_region_size, sizeof(_mean_functional_region_size) );
  gzread( backup_file, &_mean_E_region_size,          sizeof(_mean_E_region_size) );
  gzread( backup_file, &_mean_TF_region_size,         sizeof(_mean_TF_region_size) );
  gzread( backup_file, &_mean_mixed_region_size,      sizeof(_mean_mixed_region_size) );
  gzread( backup_file, &_mean_enhancer_size,          sizeof(_mean_enhancer_size) );
  gzread( backup_file, &_mean_operator_size,          sizeof(_mean_operator_size) );
  gzread( backup_file, &_mean_operon_size,            sizeof(_mean_operon_size) );
  gzread( backup_file, &_mean_E_operon_size,          sizeof(_mean_E_operon_size) );
  gzread( backup_file, &_mean_TF_operon_size,         sizeof(_mean_TF_operon_size) );
  gzread( backup_file, &_mean_mixed_operon_size,      sizeof(_mean_mixed_operon_size) );
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  gzread( backup_file, &_mean_regulation_redundancy, sizeof(_mean_regulation_redundancy) );
  gzread( backup_file, &_mean_metabolic_redundancy,  sizeof(_mean_metabolic_redundancy) );
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  gzread( backup_file, &_inherited_size,             sizeof(_inherited_size) );
  gzread( backup_file, &_inherited_nb_E,             sizeof(_inherited_nb_E) );
  gzread( backup_file, &_inherited_nb_TF,            sizeof(_inherited_nb_TF) );
  gzread( backup_file, &_inherited_nb_inner_enzymes, sizeof(_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_inherited_nb_inflow_pumps,  sizeof(_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_inherited_nb_outflow_pumps, sizeof(_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ Phenotype */
  
  gzread( backup_file, &_id,                         sizeof(_id) );
  gzread( backup_file, &_parent_id,                  sizeof(_parent_id) );
  gzread( backup_file, &_generation,                 sizeof(_generation) );
  gzread( backup_file, &_x,                          sizeof(_x) );
  gzread( backup_file, &_y,                          sizeof(_y) );
  gzread( backup_file, &_number_of_updates,          sizeof(_number_of_updates) );
  gzread( backup_file, &_number_of_divisions,        sizeof(_number_of_divisions) );
  gzread( backup_file, &_birth_time,                 sizeof(_birth_time) );
  gzread( backup_file, &_death_time,                 sizeof(_death_time) );
  gzread( backup_file, &_lifespan,                   sizeof(_lifespan) );
  gzread( backup_file, &_toxicity,                   sizeof(_toxicity) );
  gzread( backup_file, &_inherited_TF_amount,        sizeof(_inherited_TF_amount) );
  gzread( backup_file, &_inherited_E_amount,         sizeof(_inherited_E_amount) );
  gzread( backup_file, &_TF_amount,                  sizeof(_TF_amount) );
  gzread( backup_file, &_E_amount,                   sizeof(_E_amount) );
  gzread( backup_file, &_inherited_metabolic_amount, sizeof(_inherited_metabolic_amount) );
  gzread( backup_file, &_min_metabolic_amount,       sizeof(_min_metabolic_amount) );
  gzread( backup_file, &_metabolic_amount,           sizeof(_metabolic_amount) );
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
  gzread( backup_file, &_trophic_group,              sizeof(_trophic_group) );
  gzread( backup_file, &_trophic_level,              sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ List of mutation events */
  
  _list_of_events.clear();
  size_t n = 0;
  gzread( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    _list_of_events.push_back(new MutationEvent(backup_file));
  }
  
  /*------------------------------------------------------------------ Point mutations data */
  
  gzread( backup_file, &_nb_point_mutations, sizeof(_nb_point_mutations) );
  
  gzread( backup_file, &_nb_NC_point_mutations, sizeof(_nb_NC_point_mutations) );
  gzread( backup_file, &_nb_E_point_mutations,  sizeof(_nb_E_point_mutations) );
  gzread( backup_file, &_nb_TF_point_mutations, sizeof(_nb_TF_point_mutations) );
  gzread( backup_file, &_nb_BS_point_mutations, sizeof(_nb_BS_point_mutations) );
  gzread( backup_file, &_nb_P_point_mutations,  sizeof(_nb_P_point_mutations) );
  
  gzread( backup_file, &_nb_NC_to_E_transitions,  sizeof(_nb_NC_to_E_transitions) );
  gzread( backup_file, &_nb_NC_to_TF_transitions, sizeof(_nb_NC_to_TF_transitions) );
  gzread( backup_file, &_nb_NC_to_BS_transitions, sizeof(_nb_NC_to_BS_transitions) );
  gzread( backup_file, &_nb_NC_to_P_transitions,  sizeof(_nb_NC_to_P_transitions) );
  
  gzread( backup_file, &_nb_E_to_NC_transitions, sizeof(_nb_E_to_NC_transitions) );
  gzread( backup_file, &_nb_E_to_TF_transitions, sizeof(_nb_E_to_TF_transitions) );
  gzread( backup_file, &_nb_E_to_BS_transitions, sizeof(_nb_E_to_BS_transitions) );
  gzread( backup_file, &_nb_E_to_P_transitions,  sizeof(_nb_E_to_P_transitions) );
  
  gzread( backup_file, &_nb_TF_to_NC_transitions, sizeof(_nb_TF_to_NC_transitions) );
  gzread( backup_file, &_nb_TF_to_E_transitions,  sizeof(_nb_TF_to_E_transitions) );
  gzread( backup_file, &_nb_TF_to_BS_transitions, sizeof(_nb_TF_to_BS_transitions) );
  gzread( backup_file, &_nb_TF_to_P_transitions,  sizeof(_nb_TF_to_P_transitions) );
  
  gzread( backup_file, &_nb_BS_to_NC_transitions, sizeof(_nb_BS_to_NC_transitions) );
  gzread( backup_file, &_nb_BS_to_E_transitions,  sizeof(_nb_BS_to_E_transitions) );
  gzread( backup_file, &_nb_BS_to_TF_transitions, sizeof(_nb_BS_to_TF_transitions) );
  gzread( backup_file, &_nb_BS_to_P_transitions,  sizeof(_nb_BS_to_P_transitions) );
  
  gzread( backup_file, &_nb_P_to_NC_transitions, sizeof(_nb_P_to_NC_transitions) );
  gzread( backup_file, &_nb_P_to_E_transitions,  sizeof(_nb_P_to_E_transitions) );
  gzread( backup_file, &_nb_P_to_TF_transitions, sizeof(_nb_P_to_TF_transitions) );
  gzread( backup_file, &_nb_P_to_BS_transitions, sizeof(_nb_P_to_BS_transitions) );
  
  gzread( backup_file, &_mean_s_mutation_size,             sizeof(_mean_s_mutation_size) );
  gzread( backup_file, &_mean_p_mutation_size,             sizeof(_mean_p_mutation_size) );
  gzread( backup_file, &_mean_kcat_mutation_size,          sizeof(_mean_kcat_mutation_size) );
  gzread( backup_file, &_mean_kcat_km_ratio_mutation_size, sizeof(_mean_kcat_km_ratio_mutation_size) );
  
  gzread( backup_file, &_mean_BS_tag_mutation_size,  sizeof(_mean_BS_tag_mutation_size) );
  gzread( backup_file, &_mean_coE_tag_mutation_size, sizeof(_mean_coE_tag_mutation_size) );
  
  gzread( backup_file, &_mean_TF_tag_mutation_size, sizeof(_mean_TF_tag_mutation_size) );
  
  gzread( backup_file, &_mean_basal_expression_level_mutation_size, sizeof(_mean_basal_expression_level_mutation_size) );
  
  /*------------------------------------------------------------------ HGT data */
  
  gzread( backup_file, &_nb_HGT,        sizeof(_nb_HGT) );
  gzread( backup_file, &_mean_HGT_size, sizeof(_mean_HGT_size) );
  gzread( backup_file, &_nb_NC_HGT,     sizeof(_nb_NC_HGT) );
  gzread( backup_file, &_nb_E_HGT,      sizeof(_nb_E_HGT) );
  gzread( backup_file, &_nb_TF_HGT,     sizeof(_nb_TF_HGT) );
  gzread( backup_file, &_nb_BS_HGT,     sizeof(_nb_BS_HGT) );
  gzread( backup_file, &_nb_P_HGT,      sizeof(_nb_P_HGT) );
  
  /*------------------------------------------------------------------ Rearrangements data */

  gzread( backup_file, &_nb_rearrangements, sizeof(_nb_rearrangements) );
  
  gzread( backup_file, &_nb_duplicated_NC, sizeof(_nb_duplicated_NC) );
  gzread( backup_file, &_nb_duplicated_E,  sizeof(_nb_duplicated_E) );
  gzread( backup_file, &_nb_duplicated_TF, sizeof(_nb_duplicated_TF) );
  gzread( backup_file, &_nb_duplicated_BS, sizeof(_nb_duplicated_BS) );
  gzread( backup_file, &_nb_duplicated_P,  sizeof(_nb_duplicated_P) );
  
  gzread( backup_file, &_nb_deleted_NC, sizeof(_nb_deleted_NC) );
  gzread( backup_file, &_nb_deleted_E,  sizeof(_nb_deleted_E) );
  gzread( backup_file, &_nb_deleted_TF, sizeof(_nb_deleted_TF) );
  gzread( backup_file, &_nb_deleted_BS, sizeof(_nb_deleted_BS) );
  gzread( backup_file, &_nb_deleted_P,  sizeof(_nb_deleted_P) );
  
  gzread( backup_file, &_nb_duplications,   sizeof(_nb_duplications) );
  gzread( backup_file, &_nb_deletions,      sizeof(_nb_deletions) );
  gzread( backup_file, &_nb_translocations, sizeof(_nb_translocations) );
  gzread( backup_file, &_nb_inversions,     sizeof(_nb_inversions) );
  
  gzread( backup_file, &_mean_rearrangement_size, sizeof(_mean_rearrangement_size) );
  gzread( backup_file, &_mean_duplication_size,   sizeof(_mean_duplication_size) );
  gzread( backup_file, &_mean_deletion_size,      sizeof(_mean_deletion_size) );
  gzread( backup_file, &_mean_translocation_size, sizeof(_mean_translocation_size) );
  gzread( backup_file, &_mean_inversion_size,     sizeof(_mean_inversion_size) );
  
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const ReplicationReport& report
 * \return   \e void
 */
ReplicationReport::ReplicationReport( const ReplicationReport& report )
{
  /*------------------------------------------------------------------ Genome structure */
  
  _old_genome_size         = report._old_genome_size;
  _new_genome_size         = report._new_genome_size;
  _genome_functional_size  = report._genome_functional_size;
  _genome_nb_NC            = report._genome_nb_NC;
  _genome_nb_E             = report._genome_nb_E;
  _genome_nb_TF            = report._genome_nb_TF;
  _genome_nb_BS            = report._genome_nb_BS;
  _genome_nb_P             = report._genome_nb_P;
  _genome_nb_inner_enzymes = report._genome_nb_inner_enzymes;
  _genome_nb_inflow_pumps  = report._genome_nb_inflow_pumps;
  _genome_nb_outflow_pumps = report._genome_nb_outflow_pumps;
  
  /*------------------------------------------------------------------ GRN data */
  
  _nb_functional_regions       = report._nb_functional_regions;
  _nb_enhancers                = report._nb_enhancers;
  _nb_operators                = report._nb_operators;
  _nb_E_regions                = report._nb_E_regions;
  _nb_TF_regions               = report._nb_TF_regions;
  _nb_mixed_regions            = report._nb_mixed_regions;
  _mean_functional_region_size = report._mean_functional_region_size;
  _mean_E_region_size          = report._mean_E_region_size;
  _mean_TF_region_size         = report._mean_TF_region_size;
  _mean_mixed_region_size      = report._mean_mixed_region_size;
  _mean_enhancer_size          = report._mean_enhancer_size;
  _mean_operator_size          = report._mean_operator_size;
  _mean_operon_size            = report._mean_operon_size;
  _mean_E_operon_size          = report._mean_E_operon_size;
  _mean_TF_operon_size         = report._mean_TF_operon_size;
  _mean_mixed_operon_size      = report._mean_mixed_operon_size;
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  _mean_regulation_redundancy = report._mean_regulation_redundancy;
  _mean_metabolic_redundancy  = report._mean_metabolic_redundancy;
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  _inherited_size             = report._inherited_size;
  _inherited_nb_E             = report._inherited_nb_E;
  _inherited_nb_TF            = report._inherited_nb_TF;
  _inherited_nb_inner_enzymes = report._inherited_nb_inner_enzymes;
  _inherited_nb_inflow_pumps  = report._inherited_nb_inflow_pumps;
  _inherited_nb_outflow_pumps = report._inherited_nb_outflow_pumps;
  
  /*------------------------------------------------------------------ Phenotype */
  
  _id                         = report._id;
  _parent_id                  = report._parent_id;
  _generation                 = report._generation;
  _x                          = report._x;
  _y                          = report._y;
  _number_of_updates          = report._number_of_updates;
  _number_of_divisions        = report._number_of_divisions;
  _birth_time                 = report._birth_time;
  _death_time                 = report._death_time;
  _lifespan                   = report._lifespan;
  _toxicity                   = report._toxicity;
  _inherited_TF_amount        = report._inherited_TF_amount;
  _inherited_E_amount         = report._inherited_E_amount;
  _TF_amount                  = report._TF_amount;
  _E_amount                   = report._E_amount;
  _inherited_metabolic_amount = report._inherited_metabolic_amount;
  _min_metabolic_amount       = report._min_metabolic_amount;
  _metabolic_amount           = report._metabolic_amount;
  _max_metabolic_amount       = report._max_metabolic_amount;
  _metabolic_uptake           = report._metabolic_uptake;
  _metabolic_release          = report._metabolic_release;
  _previous_metabolic_amount  = report._previous_metabolic_amount;
  _min_energy                 = report._min_energy;
  _mean_energy                = report._mean_energy;
  _max_energy                 = report._max_energy;
  _min_score                  = report._min_score;
  _mean_score                 = report._mean_score;
  _max_score                  = report._max_score;
  _metabolic_growth_rate      = report._metabolic_growth_rate;
  _diff_metabolic_growth_rate = report._diff_metabolic_growth_rate;
  _grn_nb_nodes               = report._grn_nb_nodes;
  _grn_nb_edges               = report._grn_nb_edges;
  _metabolic_nb_nodes         = report._metabolic_nb_nodes;
  _metabolic_nb_edges         = report._metabolic_nb_edges;
  _trophic_group              = report._trophic_group;
  _trophic_level              = report._trophic_level;
  
  /*------------------------------------------------------------------ List of mutation events */
  
  _list_of_events.clear();
  for (size_t i = 0; i < report._list_of_events.size(); i++)
  {
    _list_of_events.push_back(new MutationEvent(*report._list_of_events[i]));
  }
  
  /*------------------------------------------------------------------ Point mutations data */
  
  _nb_point_mutations = report._nb_point_mutations;
  
  _nb_NC_point_mutations = report._nb_NC_point_mutations;
  _nb_E_point_mutations  = report._nb_E_point_mutations;
  _nb_TF_point_mutations = report._nb_TF_point_mutations;
  _nb_BS_point_mutations = report._nb_BS_point_mutations;
  _nb_P_point_mutations  = report._nb_P_point_mutations;
  
  _nb_NC_to_E_transitions  = report._nb_NC_to_E_transitions;
  _nb_NC_to_TF_transitions = report._nb_NC_to_TF_transitions;
  _nb_NC_to_BS_transitions = report._nb_NC_to_BS_transitions;
  _nb_NC_to_P_transitions  = report._nb_NC_to_P_transitions;
  
  _nb_E_to_NC_transitions = report._nb_E_to_NC_transitions;
  _nb_E_to_TF_transitions = report._nb_E_to_TF_transitions;
  _nb_E_to_BS_transitions = report._nb_E_to_BS_transitions;
  _nb_E_to_P_transitions  = report._nb_E_to_P_transitions;
  
  _nb_TF_to_NC_transitions = report._nb_TF_to_NC_transitions;
  _nb_TF_to_E_transitions  = report._nb_TF_to_E_transitions;
  _nb_TF_to_BS_transitions = report._nb_TF_to_BS_transitions;
  _nb_TF_to_P_transitions  = report._nb_TF_to_P_transitions;
  
  _nb_BS_to_NC_transitions = report._nb_BS_to_NC_transitions;
  _nb_BS_to_E_transitions  = report._nb_BS_to_E_transitions;
  _nb_BS_to_TF_transitions = report._nb_BS_to_TF_transitions;
  _nb_BS_to_P_transitions  = report._nb_BS_to_P_transitions;
  
  _nb_P_to_NC_transitions = report._nb_P_to_NC_transitions;
  _nb_P_to_E_transitions  = report._nb_P_to_E_transitions;
  _nb_P_to_TF_transitions = report._nb_P_to_TF_transitions;
  _nb_P_to_BS_transitions = report._nb_P_to_BS_transitions;
  
  _mean_s_mutation_size             = report._mean_s_mutation_size;
  _mean_p_mutation_size             = report._mean_p_mutation_size;
  _mean_kcat_mutation_size          = report._mean_kcat_mutation_size;
  _mean_kcat_km_ratio_mutation_size = report._mean_kcat_km_ratio_mutation_size;
  
  _mean_BS_tag_mutation_size  = report._mean_BS_tag_mutation_size;
  _mean_coE_tag_mutation_size = report._mean_coE_tag_mutation_size;
  
  _mean_TF_tag_mutation_size = report._mean_TF_tag_mutation_size;
  
  _mean_basal_expression_level_mutation_size = report._mean_basal_expression_level_mutation_size;
  
  /*------------------------------------------------------------------ HGT data */
  
  _nb_HGT        = report._nb_HGT;
  _mean_HGT_size = report._mean_HGT_size;
  _nb_NC_HGT     = report._nb_NC_HGT;
  _nb_E_HGT      = report._nb_E_HGT;
  _nb_TF_HGT     = report._nb_TF_HGT;
  _nb_BS_HGT     = report._nb_BS_HGT;
  _nb_P_HGT      = report._nb_P_HGT;
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  _nb_rearrangements = report._nb_rearrangements;
  
  _nb_duplicated_NC = report._nb_duplicated_NC;
  _nb_duplicated_E  = report._nb_duplicated_E;
  _nb_duplicated_TF = report._nb_duplicated_TF;
  _nb_duplicated_BS = report._nb_duplicated_BS;
  _nb_duplicated_P  = report._nb_duplicated_P;
  
  _nb_deleted_NC = report._nb_deleted_NC;
  _nb_deleted_E  = report._nb_deleted_E;
  _nb_deleted_TF = report._nb_deleted_TF;
  _nb_deleted_BS = report._nb_deleted_BS;
  _nb_deleted_P  = report._nb_deleted_P;
  
  _nb_duplications   = report._nb_duplications;
  _nb_deletions      = report._nb_deletions;
  _nb_translocations = report._nb_translocations;
  _nb_inversions     = report._nb_inversions;
  
  _mean_rearrangement_size = report._mean_rearrangement_size;
  _mean_duplication_size   = report._mean_duplication_size;
  _mean_deletion_size      = report._mean_deletion_size;
  _mean_translocation_size = report._mean_translocation_size;
  _mean_inversion_size     = report._mean_inversion_size;
  
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
ReplicationReport::~ReplicationReport( void )
{
  for (size_t i = 0; i < _list_of_events.size(); i++)
  {
    delete _list_of_events[i];
    _list_of_events[i] = NULL;
  }
  _list_of_events.clear();
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Add mutation event to the list
 * \details  --
 * \param    MutationEvent* event
 * \return   \e void
 */
void ReplicationReport::add_mutation_event( MutationEvent* event )
{
  /*----------------------------------------*/
  /* 1) save the mutation event in the list */
  /*----------------------------------------*/
  _list_of_events.push_back(event);
  
  /*----------------------------------------*/
  /* 2) update statistical data             */
  /*----------------------------------------*/
  
  /* 2.1) If mutation event is a point mutation */
  
  if (event->get_mutation_type() == POINT_MUTATION)
  {
    genetic_unit* dX = event->get_mutation_vector()->get_dX();
    _nb_point_mutations++;
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* A) point mutations that do not change the type   */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (dX->type == NON_CODING)
    {
      _nb_NC_point_mutations++;
    }
    else if (dX->type == ENZYME)
    {
      _nb_E_point_mutations++;
    }
    else if (dX->type == TRANSCRIPTION_FACTOR)
    {
      _nb_TF_point_mutations++;
    }
    else if (dX->type == BINDING_SITE)
    {
      _nb_BS_point_mutations++;
    }
    else if (dX->type == PROMOTER)
    {
      _nb_P_point_mutations++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* B) point mutations with transitions from NC to X */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (dX->type == NC_TO_E_TRANSITION)
    {
      _nb_NC_to_E_transitions++;
    }
    else if (dX->type == NC_TO_TF_TRANSITION)
    {
      _nb_NC_to_TF_transitions++;
    }
    else if (dX->type == NC_TO_BS_TRANSITION)
    {
      _nb_NC_to_BS_transitions++;
    }
    else if (dX->type == NC_TO_P_TRANSITION)
    {
      _nb_NC_to_P_transitions++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* C) point mutations with transitions from E to X  */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (dX->type == E_TO_NC_TRANSITION)
    {
      _nb_E_to_NC_transitions++;
    }
    else if (dX->type == E_TO_TF_TRANSITION)
    {
      _nb_E_to_TF_transitions++;
    }
    else if (dX->type == E_TO_BS_TRANSITION)
    {
      _nb_E_to_BS_transitions++;
    }
    else if (dX->type == E_TO_P_TRANSITION)
    {
      _nb_E_to_P_transitions++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* D) point mutations with transitions from TF to X */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (dX->type == TF_TO_NC_TRANSITION)
    {
      _nb_TF_to_NC_transitions++;
    }
    else if (dX->type == TF_TO_E_TRANSITION)
    {
      _nb_TF_to_E_transitions++;
    }
    else if (dX->type == TF_TO_BS_TRANSITION)
    {
      _nb_TF_to_BS_transitions++;
    }
    else if (dX->type == TF_TO_P_TRANSITION)
    {
      _nb_TF_to_P_transitions++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* E) point mutations with transitions from BS to X */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (dX->type == BS_TO_NC_TRANSITION)
    {
      _nb_BS_to_NC_transitions++;
    }
    else if (dX->type == BS_TO_E_TRANSITION)
    {
      _nb_BS_to_E_transitions++;
    }
    else if (dX->type == BS_TO_TF_TRANSITION)
    {
      _nb_BS_to_TF_transitions++;
    }
    else if (dX->type == BS_TO_P_TRANSITION)
    {
      _nb_BS_to_P_transitions++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* F) point mutations with transitions from P to X  */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    else if (dX->type == P_TO_NC_TRANSITION)
    {
      _nb_P_to_NC_transitions++;
    }
    else if (dX->type == P_TO_E_TRANSITION)
    {
      _nb_P_to_E_transitions++;
    }
    else if (dX->type == P_TO_TF_TRANSITION)
    {
      _nb_P_to_TF_transitions++;
    }
    else if (dX->type == P_TO_BS_TRANSITION)
    {
      _nb_P_to_BS_transitions++;
    }
    
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    /* G) save only functional attributs mutations      */
    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    if (dX->type == ENZYME)
    {
      _mean_s_mutation_size             += dX->s;
      _mean_p_mutation_size             += dX->p;
      _mean_kcat_mutation_size          += dX->kcat;
      _mean_kcat_km_ratio_mutation_size += dX->kcat_km_ratio;
    }
    else if (dX->type == TRANSCRIPTION_FACTOR)
    {
      _mean_BS_tag_mutation_size  += dX->BS_tag;
      _mean_coE_tag_mutation_size += dX->coE_tag;
    }
    else if (dX->type == BINDING_SITE)
    {
      _mean_TF_tag_mutation_size += dX->TF_tag;
    }
    else if (dX->type == PROMOTER)
    {
      _mean_basal_expression_level_mutation_size += dX->basal_expression_level;
    }
  }
  
  /* 2.2) If mutation event is a HGT */
  
  else if (event->get_mutation_type() == HGT)
  {
    _nb_HGT++;
    _mean_HGT_size += (double)event->get_size();
    _nb_NC_HGT     += event->get_nb_NC();
    _nb_E_HGT      += event->get_nb_E();
    _nb_TF_HGT     += event->get_nb_TF();
    _nb_BS_HGT     += event->get_nb_BS();
    _nb_P_HGT      += event->get_nb_P();
  }
  
  /* 2.3) If mutation event is a duplication */
  
  else if (event->get_mutation_type() == DUPLICATION)
  {
    _nb_rearrangements++;
    _nb_duplications++;
    _nb_duplicated_NC        += event->get_nb_NC();
    _nb_duplicated_E         += event->get_nb_E();
    _nb_duplicated_TF        += event->get_nb_TF();
    _nb_duplicated_BS        += event->get_nb_BS();
    _nb_duplicated_P         += event->get_nb_P();
    _mean_rearrangement_size += event->get_size();
    _mean_duplication_size   += event->get_size();
  }
  
  /* 2.4) If mutation event is a deletion */
  
  else if (event->get_mutation_type() == DELETION)
  {
    _nb_rearrangements++;
    _nb_deletions++;
    _nb_deleted_NC           += event->get_nb_NC();
    _nb_deleted_E            += event->get_nb_E();
    _nb_deleted_TF           += event->get_nb_TF();
    _nb_deleted_BS           += event->get_nb_BS();
    _nb_deleted_P            += event->get_nb_P();
    _mean_rearrangement_size += event->get_size();
    _mean_deletion_size      += event->get_size();
  }
  
  /* 2.5) If mutation event is a translocation */
  
  else if (event->get_mutation_type() == TRANSLOCATION)
  {
    _nb_rearrangements++;
    _nb_translocations++;
    _mean_rearrangement_size += event->get_size();
    _mean_translocation_size += event->get_size();
  }
  
  /* 2.6) If mutation event is an inversion */
  
  else if (event->get_mutation_type() == INVERSION)
  {
    _nb_rearrangements++;
    _nb_inversions++;
    _mean_rearrangement_size += event->get_size();
    _mean_inversion_size     += event->get_size();
  }
}

/**
 * \brief    Compute mean statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void ReplicationReport::compute_mean( void )
{
  if (_nb_E_point_mutations > 0)
  {
    _mean_s_mutation_size             /= _nb_E_point_mutations;
    _mean_p_mutation_size             /= _nb_E_point_mutations;
    _mean_kcat_mutation_size          /= _nb_E_point_mutations;
    _mean_kcat_km_ratio_mutation_size /= _nb_E_point_mutations;
  }
  if (_nb_TF_point_mutations > 0)
  {
    _mean_BS_tag_mutation_size  /= _nb_TF_point_mutations;
    _mean_coE_tag_mutation_size /= _nb_TF_point_mutations;
  }
  if (_nb_BS_point_mutations > 0)
  {
    _mean_TF_tag_mutation_size /= _nb_BS_point_mutations;
  }
  if (_nb_P_point_mutations > 0)
  {
    _mean_basal_expression_level_mutation_size /= _nb_P_point_mutations;
  }
  if (_nb_HGT > 0)
  {
    _mean_HGT_size /= _nb_HGT;
  }
  if (_nb_rearrangements > 0)
  {
    _mean_rearrangement_size /= _nb_rearrangements;
  }
  if (_nb_duplications > 0)
  {
    _mean_duplication_size /= _nb_duplications;
  }
  if (_nb_deletions > 0)
  {
    _mean_deletion_size /= _nb_deletions;
  }
  if (_nb_translocations > 0)
  {
    _mean_translocation_size /= _nb_translocations;
  }
  if (_nb_inversions > 0)
  {
    _mean_inversion_size /= _nb_inversions;
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void ReplicationReport::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ Genome structure */
  
  gzwrite( backup_file, &_old_genome_size,         sizeof(_old_genome_size) );
  gzwrite( backup_file, &_new_genome_size,         sizeof(_new_genome_size) );
  gzwrite( backup_file, &_genome_functional_size,  sizeof(_genome_functional_size) );
  gzwrite( backup_file, &_genome_nb_NC,            sizeof(_genome_nb_NC) );
  gzwrite( backup_file, &_genome_nb_E,             sizeof(_genome_nb_E) );
  gzwrite( backup_file, &_genome_nb_TF,            sizeof(_genome_nb_TF) );
  gzwrite( backup_file, &_genome_nb_BS,            sizeof(_genome_nb_BS) );
  gzwrite( backup_file, &_genome_nb_P,             sizeof(_genome_nb_P) );
  gzwrite( backup_file, &_genome_nb_inner_enzymes, sizeof(_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_genome_nb_inflow_pumps,  sizeof(_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_genome_nb_outflow_pumps, sizeof(_genome_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ GRN data */
  
  gzwrite( backup_file, &_nb_functional_regions,       sizeof(_nb_functional_regions) );
  gzwrite( backup_file, &_nb_enhancers,                sizeof(_nb_enhancers) );
  gzwrite( backup_file, &_nb_operators,                sizeof(_nb_operators) );
  gzwrite( backup_file, &_nb_E_regions,                sizeof(_nb_E_regions) );
  gzwrite( backup_file, &_nb_TF_regions,               sizeof(_nb_TF_regions) );
  gzwrite( backup_file, &_nb_mixed_regions,            sizeof(_nb_mixed_regions) );
  gzwrite( backup_file, &_mean_functional_region_size, sizeof(_mean_functional_region_size) );
  gzwrite( backup_file, &_mean_E_region_size,          sizeof(_mean_E_region_size) );
  gzwrite( backup_file, &_mean_TF_region_size,         sizeof(_mean_TF_region_size) );
  gzwrite( backup_file, &_mean_mixed_region_size,      sizeof(_mean_mixed_region_size) );
  gzwrite( backup_file, &_mean_enhancer_size,          sizeof(_mean_enhancer_size) );
  gzwrite( backup_file, &_mean_operator_size,          sizeof(_mean_operator_size) );
  gzwrite( backup_file, &_mean_operon_size,            sizeof(_mean_operon_size) );
  gzwrite( backup_file, &_mean_E_operon_size,          sizeof(_mean_E_operon_size) );
  gzwrite( backup_file, &_mean_TF_operon_size,         sizeof(_mean_TF_operon_size) );
  gzwrite( backup_file, &_mean_mixed_operon_size,      sizeof(_mean_mixed_operon_size) );
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  gzwrite( backup_file, &_mean_regulation_redundancy, sizeof(_mean_regulation_redundancy) );
  gzwrite( backup_file, &_mean_metabolic_redundancy,  sizeof(_mean_metabolic_redundancy) );
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  gzwrite( backup_file, &_inherited_size,             sizeof(_inherited_size) );
  gzwrite( backup_file, &_inherited_nb_E,             sizeof(_inherited_nb_E) );
  gzwrite( backup_file, &_inherited_nb_TF,            sizeof(_inherited_nb_TF) );
  gzwrite( backup_file, &_inherited_nb_inner_enzymes, sizeof(_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_inherited_nb_inflow_pumps,  sizeof(_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_inherited_nb_outflow_pumps, sizeof(_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ Phenotype */
  
  gzwrite( backup_file, &_id,                         sizeof(_id) );
  gzwrite( backup_file, &_parent_id,                  sizeof(_parent_id) );
  gzwrite( backup_file, &_generation,                 sizeof(_generation) );
  gzwrite( backup_file, &_x,                          sizeof(_x) );
  gzwrite( backup_file, &_y,                          sizeof(_y) );
  gzwrite( backup_file, &_number_of_updates,          sizeof(_number_of_updates) );
  gzwrite( backup_file, &_number_of_divisions,        sizeof(_number_of_divisions) );
  gzwrite( backup_file, &_birth_time,                 sizeof(_birth_time) );
  gzwrite( backup_file, &_death_time,                 sizeof(_death_time) );
  gzwrite( backup_file, &_lifespan,                   sizeof(_lifespan) );
  gzwrite( backup_file, &_toxicity,                   sizeof(_toxicity) );
  gzwrite( backup_file, &_inherited_TF_amount,        sizeof(_inherited_TF_amount) );
  gzwrite( backup_file, &_inherited_E_amount,         sizeof(_inherited_E_amount) );
  gzwrite( backup_file, &_TF_amount,                  sizeof(_TF_amount) );
  gzwrite( backup_file, &_E_amount,                   sizeof(_E_amount) );
  gzwrite( backup_file, &_inherited_metabolic_amount, sizeof(_inherited_metabolic_amount) );
  gzwrite( backup_file, &_min_metabolic_amount,       sizeof(_min_metabolic_amount) );
  gzwrite( backup_file, &_metabolic_amount,           sizeof(_metabolic_amount) );
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
  gzwrite( backup_file, &_trophic_group,              sizeof(_trophic_group) );
  gzwrite( backup_file, &_trophic_level,              sizeof(_trophic_level) );
  
  /*------------------------------------------------------------------ List of mutation events */
  
  size_t n = _list_of_events.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (size_t i = 0; i < n; i++)
  {
    _list_of_events[i]->save( backup_file );
  }
  
  /*------------------------------------------------------------------ Point mutations data */
  
  gzwrite( backup_file, &_nb_point_mutations, sizeof(_nb_point_mutations) );
  
  gzwrite( backup_file, &_nb_NC_point_mutations, sizeof(_nb_NC_point_mutations) );
  gzwrite( backup_file, &_nb_E_point_mutations,  sizeof(_nb_E_point_mutations) );
  gzwrite( backup_file, &_nb_TF_point_mutations, sizeof(_nb_TF_point_mutations) );
  gzwrite( backup_file, &_nb_BS_point_mutations, sizeof(_nb_BS_point_mutations) );
  gzwrite( backup_file, &_nb_P_point_mutations,  sizeof(_nb_P_point_mutations) );
  
  gzwrite( backup_file, &_nb_NC_to_E_transitions,  sizeof(_nb_NC_to_E_transitions) );
  gzwrite( backup_file, &_nb_NC_to_TF_transitions, sizeof(_nb_NC_to_TF_transitions) );
  gzwrite( backup_file, &_nb_NC_to_BS_transitions, sizeof(_nb_NC_to_BS_transitions) );
  gzwrite( backup_file, &_nb_NC_to_P_transitions,  sizeof(_nb_NC_to_P_transitions) );
  
  gzwrite( backup_file, &_nb_E_to_NC_transitions, sizeof(_nb_E_to_NC_transitions) );
  gzwrite( backup_file, &_nb_E_to_TF_transitions, sizeof(_nb_E_to_TF_transitions) );
  gzwrite( backup_file, &_nb_E_to_BS_transitions, sizeof(_nb_E_to_BS_transitions) );
  gzwrite( backup_file, &_nb_E_to_P_transitions,  sizeof(_nb_E_to_P_transitions) );
  
  gzwrite( backup_file, &_nb_TF_to_NC_transitions, sizeof(_nb_TF_to_NC_transitions) );
  gzwrite( backup_file, &_nb_TF_to_E_transitions,  sizeof(_nb_TF_to_E_transitions) );
  gzwrite( backup_file, &_nb_TF_to_BS_transitions, sizeof(_nb_TF_to_BS_transitions) );
  gzwrite( backup_file, &_nb_TF_to_P_transitions,  sizeof(_nb_TF_to_P_transitions) );
  
  gzwrite( backup_file, &_nb_BS_to_NC_transitions, sizeof(_nb_BS_to_NC_transitions) );
  gzwrite( backup_file, &_nb_BS_to_E_transitions,  sizeof(_nb_BS_to_E_transitions) );
  gzwrite( backup_file, &_nb_BS_to_TF_transitions, sizeof(_nb_BS_to_TF_transitions) );
  gzwrite( backup_file, &_nb_BS_to_P_transitions,  sizeof(_nb_BS_to_P_transitions) );
  
  gzwrite( backup_file, &_nb_P_to_NC_transitions, sizeof(_nb_P_to_NC_transitions) );
  gzwrite( backup_file, &_nb_P_to_E_transitions,  sizeof(_nb_P_to_E_transitions) );
  gzwrite( backup_file, &_nb_P_to_TF_transitions, sizeof(_nb_P_to_TF_transitions) );
  gzwrite( backup_file, &_nb_P_to_BS_transitions, sizeof(_nb_P_to_BS_transitions) );
  
  gzwrite( backup_file, &_mean_s_mutation_size,             sizeof(_mean_s_mutation_size) );
  gzwrite( backup_file, &_mean_p_mutation_size,             sizeof(_mean_p_mutation_size) );
  gzwrite( backup_file, &_mean_kcat_mutation_size,          sizeof(_mean_kcat_mutation_size) );
  gzwrite( backup_file, &_mean_kcat_km_ratio_mutation_size, sizeof(_mean_kcat_km_ratio_mutation_size) );
  
  gzwrite( backup_file, &_mean_BS_tag_mutation_size,  sizeof(_mean_BS_tag_mutation_size) );
  gzwrite( backup_file, &_mean_coE_tag_mutation_size, sizeof(_mean_coE_tag_mutation_size) );
  
  gzwrite( backup_file, &_mean_TF_tag_mutation_size, sizeof(_mean_TF_tag_mutation_size) );
  
  gzwrite( backup_file, &_mean_basal_expression_level_mutation_size, sizeof(_mean_basal_expression_level_mutation_size) );
  
  /*------------------------------------------------------------------ HGT data */
  
  gzwrite( backup_file, &_nb_HGT,        sizeof(_nb_HGT) );
  gzwrite( backup_file, &_mean_HGT_size, sizeof(_mean_HGT_size) );
  gzwrite( backup_file, &_nb_NC_HGT,     sizeof(_nb_NC_HGT) );
  gzwrite( backup_file, &_nb_E_HGT,      sizeof(_nb_E_HGT) );
  gzwrite( backup_file, &_nb_TF_HGT,     sizeof(_nb_TF_HGT) );
  gzwrite( backup_file, &_nb_BS_HGT,     sizeof(_nb_BS_HGT) );
  gzwrite( backup_file, &_nb_P_HGT,      sizeof(_nb_P_HGT) );
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  gzwrite( backup_file, &_nb_rearrangements, sizeof(_nb_rearrangements) );
  
  gzwrite( backup_file, &_nb_duplicated_NC, sizeof(_nb_duplicated_NC) );
  gzwrite( backup_file, &_nb_duplicated_E,  sizeof(_nb_duplicated_E) );
  gzwrite( backup_file, &_nb_duplicated_TF, sizeof(_nb_duplicated_TF) );
  gzwrite( backup_file, &_nb_duplicated_BS, sizeof(_nb_duplicated_BS) );
  gzwrite( backup_file, &_nb_duplicated_P,  sizeof(_nb_duplicated_P) );
  
  gzwrite( backup_file, &_nb_deleted_NC, sizeof(_nb_deleted_NC) );
  gzwrite( backup_file, &_nb_deleted_E,  sizeof(_nb_deleted_E) );
  gzwrite( backup_file, &_nb_deleted_TF, sizeof(_nb_deleted_TF) );
  gzwrite( backup_file, &_nb_deleted_BS, sizeof(_nb_deleted_BS) );
  gzwrite( backup_file, &_nb_deleted_P,  sizeof(_nb_deleted_P) );
  
  gzwrite( backup_file, &_nb_duplications,   sizeof(_nb_duplications) );
  gzwrite( backup_file, &_nb_deletions,      sizeof(_nb_deletions) );
  gzwrite( backup_file, &_nb_translocations, sizeof(_nb_translocations) );
  gzwrite( backup_file, &_nb_inversions,     sizeof(_nb_inversions) );
  
  gzwrite( backup_file, &_mean_rearrangement_size, sizeof(_mean_rearrangement_size) );
  gzwrite( backup_file, &_mean_duplication_size,   sizeof(_mean_duplication_size) );
  gzwrite( backup_file, &_mean_deletion_size,      sizeof(_mean_deletion_size) );
  gzwrite( backup_file, &_mean_translocation_size, sizeof(_mean_translocation_size) );
  gzwrite( backup_file, &_mean_inversion_size,     sizeof(_mean_inversion_size) );
}

/**
 * \brief    Clear replication report
 * \details  --
 * \param    void
 * \return   \e void
 */
void ReplicationReport::clear( void )
{
  /*------------------------------------------------------------------ Genome structure */
  
  _old_genome_size         = 0;
  _new_genome_size         = 0;
  _genome_functional_size  = 0;
  _genome_nb_NC            = 0;
  _genome_nb_E             = 0;
  _genome_nb_TF            = 0;
  _genome_nb_BS            = 0;
  _genome_nb_P             = 0;
  _genome_nb_inner_enzymes = 0;
  _genome_nb_inflow_pumps  = 0;
  _genome_nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ GRN data */
  
  _nb_functional_regions       = 0;
  _nb_enhancers                = 0;
  _nb_operators                = 0;
  _nb_E_regions                = 0;
  _nb_TF_regions               = 0;
  _nb_mixed_regions            = 0;
  _mean_functional_region_size = 0.0;
  _mean_E_region_size          = 0.0;
  _mean_TF_region_size         = 0.0;
  _mean_mixed_region_size      = 0.0;
  _mean_enhancer_size          = 0.0;
  _mean_operator_size          = 0.0;
  _mean_operon_size            = 0.0;
  _mean_E_operon_size          = 0.0;
  _mean_TF_operon_size         = 0.0;
  _mean_mixed_operon_size      = 0.0;
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  _mean_regulation_redundancy = 0.0;
  _mean_metabolic_redundancy  = 0.0;
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  _inherited_size             = 0;
  _inherited_nb_E             = 0;
  _inherited_nb_TF            = 0;
  _inherited_nb_inner_enzymes = 0;
  _inherited_nb_inflow_pumps  = 0;
  _inherited_nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ Phenotype */
  
  _id                         = 0;
  _parent_id                  = 0;
  _generation                 = 0;
  _x                          = 0;
  _y                          = 0;
  _number_of_updates          = 0;
  _number_of_divisions        = 0;
  _birth_time                 = 0.0;
  _death_time                 = 0.0;
  _lifespan                   = 0.0;
  _toxicity                   = 0.0;
  _inherited_TF_amount        = 0.0;
  _inherited_E_amount         = 0.0;
  _TF_amount                  = 0.0;
  _E_amount                   = 0.0;
  _inherited_metabolic_amount = 0.0;
  _min_metabolic_amount       = 0.0;
  _metabolic_amount           = 0.0;
  _max_metabolic_amount       = 0.0;
  _metabolic_uptake           = 0.0;
  _metabolic_release          = 0.0;
  _previous_metabolic_amount  = 0.0;
  _min_energy                 = 0.0;
  _mean_energy                = 0.0;
  _max_energy                 = 0.0;
  _min_score                  = 0.0;
  _mean_score                 = 0.0;
  _max_score                  = 0.0;
  _metabolic_growth_rate      = 0.0;
  _diff_metabolic_growth_rate = 0.0;
  _grn_nb_nodes               = 0;
  _grn_nb_edges               = 0;
  _metabolic_nb_nodes         = 0;
  _metabolic_nb_edges         = 0;
  _trophic_group              = 0;
  _trophic_level              = NO_LEVEL;
  
  /*------------------------------------------------------------------ List of mutation events */
  
  for (size_t i = 0; i < _list_of_events.size(); i++)
  {
    delete _list_of_events[i];
    _list_of_events[i] = NULL;
  }
  _list_of_events.clear();
  
  /*------------------------------------------------------------------ Point mutations data */
  
  _nb_point_mutations = 0;
  
  _nb_NC_point_mutations = 0;
  _nb_E_point_mutations  = 0;
  _nb_TF_point_mutations = 0;
  _nb_BS_point_mutations = 0;
  _nb_P_point_mutations  = 0;
  
  _nb_NC_to_E_transitions  = 0;
  _nb_NC_to_TF_transitions = 0;
  _nb_NC_to_BS_transitions = 0;
  _nb_NC_to_P_transitions  = 0;
  
  _nb_E_to_NC_transitions = 0;
  _nb_E_to_TF_transitions = 0;
  _nb_E_to_BS_transitions = 0;
  _nb_E_to_P_transitions  = 0;
  
  _nb_TF_to_NC_transitions = 0;
  _nb_TF_to_E_transitions  = 0;
  _nb_TF_to_BS_transitions = 0;
  _nb_TF_to_P_transitions  = 0;
  
  _nb_BS_to_NC_transitions = 0;
  _nb_BS_to_E_transitions  = 0;
  _nb_BS_to_TF_transitions = 0;
  _nb_BS_to_P_transitions  = 0;
  
  _nb_P_to_NC_transitions = 0;
  _nb_P_to_E_transitions  = 0;
  _nb_P_to_TF_transitions = 0;
  _nb_P_to_BS_transitions = 0;
  
  _mean_s_mutation_size             = 0.0;
  _mean_p_mutation_size             = 0.0;
  _mean_kcat_mutation_size          = 0.0;
  _mean_kcat_km_ratio_mutation_size = 0.0;
  
  _mean_BS_tag_mutation_size  = 0.0;
  _mean_coE_tag_mutation_size = 0.0;
  
  _mean_TF_tag_mutation_size = 0.0;
  
  _mean_basal_expression_level_mutation_size = 0.0;
  
  /*------------------------------------------------------------------ HGT data */
  
  _nb_HGT        = 0;
  _mean_HGT_size = 0.0;
  _nb_NC_HGT     = 0;
  _nb_E_HGT      = 0;
  _nb_TF_HGT     = 0;
  _nb_BS_HGT     = 0;
  _nb_P_HGT      = 0;
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  _nb_rearrangements = 0;
  
  _nb_duplicated_NC = 0;
  _nb_duplicated_E  = 0;
  _nb_duplicated_TF = 0;
  _nb_duplicated_BS = 0;
  _nb_duplicated_P  = 0;
  
  _nb_deleted_NC = 0;
  _nb_deleted_E  = 0;
  _nb_deleted_TF = 0;
  _nb_deleted_BS = 0;
  _nb_deleted_P  = 0;
  
  _nb_duplications   = 0;
  _nb_deletions      = 0;
  _nb_translocations = 0;
  _nb_inversions     = 0;
  
  _mean_rearrangement_size = 0.0;
  _mean_duplication_size   = 0.0;
  _mean_deletion_size      = 0.0;
  _mean_translocation_size = 0.0;
  _mean_inversion_size     = 0.0;
}

/**
 * \brief    Clear event list only
 * \details  --
 * \param    void
 * \return   \e void
 */
void ReplicationReport::clear_event_list( void )
{
  for (size_t i = 0; _list_of_events.size(); i++)
  {
    std::cout << _list_of_events[i]->get_mutation_type() << "\n";
    delete _list_of_events[i];
    _list_of_events[i] = NULL;
  }
  _list_of_events.clear();
}

/**
 * \brief    Write the genome structure header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_genome_structure_header( std::ofstream& filestream )
{
  filestream << "generation" << " " <<
  "old_genome_size" << " " <<
  "genome_size" << " " <<
  "genome_functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "nb_functional_regions" << " " <<
  "nb_enhancers" << " " <<
  "nb_operators" << " " <<
  "nb_E_regions" << " " <<
  "nb_TF_regions" << " " <<
  "nb_mixed_regions" << " " <<
  "mean_functional_region_size" << " " <<
  "mean_E_region_size" << " " <<
  "mean_TF_region_size" << " " <<
  "mean_mixed_region_size" << " " <<
  "mean_enhancer_size" << " " <<
  "mean_operator_size" << " " <<
  "mean_operon_size" << " " <<
  "mean_E_operon_size" << " " <<
  "mean_TF_operon_size" << " " <<
  "mean_mixed_operon_size" << " " <<
  "mean_regulation_redundancy" << " " <<
  "mean_metabolic_redundancy" << "\n";
}

/**
 * \brief    Write the inherited proteins header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_inherited_proteins_header( std::ofstream& filestream )
{
  filestream << "generation" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps" << " " <<
  "inherited_nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write the phenotype header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_phenotype_header( std::ofstream& filestream )
{
  filestream << "generation" << " " <<
  "id" << " " <<
  "parent_id" << " " <<
  "x" << " " <<
  "y" << " " <<
  "number_of_updates" << " " <<
  "number_of_divisions" << " " <<
  "birth_time" << " " <<
  "death_time" << " " <<
  "lifespan" << " " <<
  "toxicity" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "E_amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "min_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "max_metabolic_amount" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "previous_metabolic_amount" << " " <<
  "min_energy" << " " <<
  "mean_energy" << " " <<
  "max_energy" << " " <<
  "min_score" << " " <<
  "mean_score" << " " <<
  "max_score" << " " <<
  "metabolic_growth_rate" << " " <<
  "diff_metabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "trophic_group" << " " <<
  "trophic_level" << "\n";
}

/**
 * \brief    Write the fixed mutations header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_fixed_mutations_header( std::ofstream& filestream )
{
  filestream << "generation" << " " <<
  "old_genome_size" << " " <<
  "nb_point_mutations" << " " <<
  "nb_NC_point_mutations" << " " <<
  "nb_E_point_mutations" << " " <<
  "nb_TF_point_mutations" << " " <<
  "nb_BS_point_mutations" << " " <<
  "nb_P_point_mutations" << " " <<
  "nb_NC_to_E_transitions" << " " <<
  "nb_NC_to_TF_transitions" << " " <<
  "nb_NC_to_BS_transitions" << " " <<
  "nb_NC_to_P_transitions" << " " <<
  "nb_E_to_NC_transitions" << " " <<
  "nb_E_to_TF_transitions" << " " <<
  "nb_E_to_BS_transitions" << " " <<
  "nb_E_to_P_transitions" << " " <<
  "nb_TF_to_NC_transitions" << " " <<
  "nb_TF_to_E_transitions" << " " <<
  "nb_TF_to_BS_transitions" << " " <<
  "nb_TF_to_P_transitions" << " " <<
  "nb_BS_to_NC_transitions" << " " <<
  "nb_BS_to_E_transitions" << " " <<
  "nb_BS_to_TF_transitions" << " " <<
  "nb_BS_to_P_transitions" << " " <<
  "nb_P_to_NC_transitions" << " " <<
  "nb_P_to_E_transitions" << " " <<
  "nb_P_to_TF_transitions" << " " <<
  "nb_P_to_BS_transitions" << " " <<
  "s_mutation_size" << " " <<
  "p_mutation_size" << " " <<
  "kcat_mutation_size" << " " <<
  "kcat_km_ratio_mutation_size" << " " <<
  "BS_tag_mutation_size" << " " <<
  "coE_tag_mutation_size" << " " <<
  "TF_tag_mutation_size" << " " <<
  "beta_mutation_size" << " " <<
  "nb_HGT" << " " <<
  "HGT_size" << " " <<
  "nb_NC_HGT" << " " <<
  "nb_E_HGT" << " " <<
  "nb_TF_HGT" << " " <<
  "nb_BS_HGT" << " " <<
  "nb_P_HGT" << " " <<
  "nb_rearrangements" << " " <<
  "nb_duplicated_NC" << " " <<
  "nb_duplicated_E" << " " <<
  "nb_duplicated_TF" << " " <<
  "nb_duplicated_BS" << " " <<
  "nb_duplicated_P" << " " <<
  "nb_deleted_NC" << " " <<
  "nb_deleted_E" << " " <<
  "nb_deleted_TF" << " " <<
  "nb_deleted_BS" << " " <<
  "nb_deleted_P" << " " <<
  "nb_duplications" << " " <<
  "nb_deletions" << " " <<
  "nb_translocations" << " " <<
  "nb_inversions" << " " <<
  "rearrangement_size" << " " <<
  "duplication_size" << " " <<
  "deletion_size" << " " <<
  "translocation_size" << " " <<
  "inversion_size" << "\n";
}

/**
 * \brief    Write the replication report header
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_replication_report_header( std::ofstream& filestream )
{
  filestream << "generation" << " " <<
  "old_genome_size" << " " <<
  "genome_size" << " " <<
  "genome_functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "nb_functional_regions" << " " <<
  "nb_enhancers" << " " <<
  "nb_operators" << " " <<
  "nb_E_regions" << " " <<
  "nb_TF_regions" << " " <<
  "nb_mixed_regions" << " " <<
  "mean_functional_region_size" << " " <<
  "mean_E_region_size" << " " <<
  "mean_TF_region_size" << " " <<
  "mean_mixed_region_size" << " " <<
  "mean_enhancer_size" << " " <<
  "mean_operator_size" << " " <<
  "mean_operon_size" << " " <<
  "mean_E_operon_size" << " " <<
  "mean_TF_operon_size" << " " <<
  "mean_mixed_operon_size" << " " <<
  "mean_regulation_redundancy" << " " <<
  "mean_metabolic_redundancy" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps" << " " <<
  "inherited_nb_outflow_pumps" << " " <<
  "id" << " " <<
  "parent_id" << " " <<
  "x" << " " <<
  "y" << " " <<
  "number_of_updates" << " " <<
  "number_of_divisions" << " " <<
  "birth_time" << " " <<
  "death_time" << " " <<
  "lifespan" << " " <<
  "toxicity" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "E_amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "min_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "max_metabolic_amount" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "previous_metabolic_amount" << " " <<
  "min_energy" << " " <<
  "mean_energy" << " " <<
  "max_energy" << " " <<
  "min_score" << " " <<
  "mean_score" << " " <<
  "max_score" << " " <<
  "metabolic_growth_rate" << " " <<
  "diff_metabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "trophic_group" << " " <<
  "trophic_level" << " " <<
  "old_genome_size" << " " <<
  "nb_point_mutations" << " " <<
  "nb_NC_point_mutations" << " " <<
  "nb_E_point_mutations" << " " <<
  "nb_TF_point_mutations" << " " <<
  "nb_BS_point_mutations" << " " <<
  "nb_P_point_mutations" << " " <<
  "nb_NC_to_E_transitions" << " " <<
  "nb_NC_to_TF_transitions" << " " <<
  "nb_NC_to_BS_transitions" << " " <<
  "nb_NC_to_P_transitions" << " " <<
  "nb_E_to_NC_transitions" << " " <<
  "nb_E_to_TF_transitions" << " " <<
  "nb_E_to_BS_transitions" << " " <<
  "nb_E_to_P_transitions" << " " <<
  "nb_TF_to_NC_transitions" << " " <<
  "nb_TF_to_E_transitions" << " " <<
  "nb_TF_to_BS_transitions" << " " <<
  "nb_TF_to_P_transitions" << " " <<
  "nb_BS_to_NC_transitions" << " " <<
  "nb_BS_to_E_transitions" << " " <<
  "nb_BS_to_TF_transitions" << " " <<
  "nb_BS_to_P_transitions" << " " <<
  "nb_P_to_NC_transitions" << " " <<
  "nb_P_to_E_transitions" << " " <<
  "nb_P_to_TF_transitions" << " " <<
  "nb_P_to_BS_transitions" << " " <<
  "s_mutation_size" << " " <<
  "p_mutation_size" << " " <<
  "kcat_mutation_size" << " " <<
  "kcat_km_ratio_mutation_size" << " " <<
  "BS_tag_mutation_size" << " " <<
  "coE_tag_mutation_size" << " " <<
  "TF_tag_mutation_size" << " " <<
  "beta_mutation_size" << " " <<
  "nb_HGT" << " " <<
  "HGT_size" << " " <<
  "nb_NC_HGT" << " " <<
  "nb_E_HGT" << " " <<
  "nb_TF_HGT" << " " <<
  "nb_BS_HGT" << " " <<
  "nb_P_HGT" << " " <<
  "nb_rearrangements" << " " <<
  "nb_duplicated_NC" << " " <<
  "nb_duplicated_E" << " " <<
  "nb_duplicated_TF" << " " <<
  "nb_duplicated_BS" << " " <<
  "nb_duplicated_P" << " " <<
  "nb_deleted_NC" << " " <<
  "nb_deleted_E" << " " <<
  "nb_deleted_TF" << " " <<
  "nb_deleted_BS" << " " <<
  "nb_deleted_P" << " " <<
  "nb_duplications" << " " <<
  "nb_deletions" << " " <<
  "nb_translocations" << " " <<
  "nb_inversions" << " " <<
  "rearrangement_size" << " " <<
  "duplication_size" << " " <<
  "deletion_size" << " " <<
  "translocation_size" << " " <<
  "inversion_size" << "\n";
}

/**
 * \brief    Write the genome structure data
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_genome_structure_data( std::ofstream& filestream )
{
  filestream << _generation << " " <<
  _old_genome_size << " " <<
  _new_genome_size << " " <<
  _genome_functional_size << " " <<
  _genome_nb_NC << " " <<
  _genome_nb_E << " " <<
  _genome_nb_TF << " " <<
  _genome_nb_BS << " " <<
  _genome_nb_P << " " <<
  _genome_nb_inner_enzymes << " " <<
  _genome_nb_inflow_pumps << " " <<
  _genome_nb_outflow_pumps << " " <<
  _nb_functional_regions << " " <<
  _nb_enhancers << " " <<
  _nb_operators << " " <<
  _nb_E_regions << " " <<
  _nb_TF_regions << " " <<
  _nb_mixed_regions << " " <<
  _mean_functional_region_size << " " <<
  _mean_E_region_size << " " <<
  _mean_TF_region_size << " " <<
  _mean_mixed_region_size << " " <<
  _mean_enhancer_size << " " <<
  _mean_operator_size << " " <<
  _mean_operon_size << " " <<
  _mean_E_operon_size << " " <<
  _mean_TF_operon_size << " " <<
  _mean_mixed_operon_size << " " <<
  _mean_regulation_redundancy << " " <<
  _mean_metabolic_redundancy << "\n";
}

/**
 * \brief    Write the inherited proteins data
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_inherited_proteins_data( std::ofstream& filestream )
{
  filestream << _generation << " " <<
  _inherited_nb_E << " " <<
  _inherited_nb_TF << " " <<
  _inherited_nb_inner_enzymes << " " <<
  _inherited_nb_inflow_pumps << " " <<
  _inherited_nb_outflow_pumps << "\n";
}

/**
 * \brief    Write the phenotype data
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_phenotype_data( std::ofstream& filestream )
{
  filestream << _generation << " " <<
  _id << " " <<
  _parent_id << " " <<
  _x << " " <<
  _y << " " <<
  _number_of_updates << " " <<
  _number_of_divisions << " " <<
  _birth_time << " " <<
  _death_time << " " <<
  _lifespan << " " <<
  _toxicity << " " <<
  _inherited_TF_amount << " " <<
  _inherited_E_amount << " " <<
  _TF_amount << " " <<
  _E_amount << " " <<
  _inherited_metabolic_amount << " " <<
  _min_metabolic_amount << " " <<
  _metabolic_amount << " " <<
  _max_metabolic_amount << " " <<
  _metabolic_uptake << " " <<
  _metabolic_release << " " <<
  _previous_metabolic_amount << " " <<
  _min_energy << " " <<
  _mean_energy << " " <<
  _max_energy << " " <<
  _min_score << " " <<
  _mean_score << " " <<
  _max_score << " " <<
  _metabolic_growth_rate << " " <<
  _diff_metabolic_growth_rate << " " <<
  _grn_nb_nodes << " " <<
  _grn_nb_edges << " " <<
  _metabolic_nb_nodes << " " <<
  _metabolic_nb_edges << " " <<
  _trophic_group << " " <<
  _trophic_level << "\n";
}

/**
 * \brief    Write the fixed mutations data
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_fixed_mutations_data( std::ofstream& filestream )
{
  filestream << _generation << " " <<
  _old_genome_size << " " <<
  _nb_point_mutations << " " <<
  _nb_NC_point_mutations << " " <<
  _nb_E_point_mutations << " " <<
  _nb_TF_point_mutations << " " <<
  _nb_BS_point_mutations << " " <<
  _nb_P_point_mutations << " " <<
  _nb_NC_to_E_transitions << " " <<
  _nb_NC_to_TF_transitions << " " <<
  _nb_NC_to_BS_transitions << " " <<
  _nb_NC_to_P_transitions << " " <<
  _nb_E_to_NC_transitions << " " <<
  _nb_E_to_TF_transitions << " " <<
  _nb_E_to_BS_transitions << " " <<
  _nb_E_to_P_transitions << " " <<
  _nb_TF_to_NC_transitions << " " <<
  _nb_TF_to_E_transitions << " " <<
  _nb_TF_to_BS_transitions << " " <<
  _nb_TF_to_P_transitions << " " <<
  _nb_BS_to_NC_transitions << " " <<
  _nb_BS_to_E_transitions << " " <<
  _nb_BS_to_TF_transitions << " " <<
  _nb_BS_to_P_transitions << " " <<
  _nb_P_to_NC_transitions << " " <<
  _nb_P_to_E_transitions << " " <<
  _nb_P_to_TF_transitions << " " <<
  _nb_P_to_BS_transitions << " " <<
  _mean_s_mutation_size << " " <<
  _mean_p_mutation_size << " " <<
  _mean_kcat_mutation_size << " " <<
  _mean_kcat_km_ratio_mutation_size << " " <<
  _mean_BS_tag_mutation_size << " " <<
  _mean_coE_tag_mutation_size << " " <<
  _mean_TF_tag_mutation_size << " " <<
  _mean_basal_expression_level_mutation_size << " " <<
  _nb_HGT << " " <<
  _mean_HGT_size << " " <<
  _nb_NC_HGT << " " <<
  _nb_E_HGT << " " <<
  _nb_TF_HGT << " " <<
  _nb_BS_HGT << " " <<
  _nb_P_HGT << " " <<
  _nb_rearrangements << " " <<
  _nb_duplicated_NC << " " <<
  _nb_duplicated_E << " " <<
  _nb_duplicated_TF << " " <<
  _nb_duplicated_BS << " " <<
  _nb_duplicated_P << " " <<
  _nb_deleted_NC << " " <<
  _nb_deleted_E << " " <<
  _nb_deleted_TF << " " <<
  _nb_deleted_BS << " " <<
  _nb_deleted_P << " " <<
  _nb_duplications << " " <<
  _nb_deletions << " " <<
  _nb_translocations << " " <<
  _nb_inversions << " " <<
  _mean_rearrangement_size << " " <<
  _mean_duplication_size << " " <<
  _mean_deletion_size << " " <<
  _mean_translocation_size << " " <<
  _mean_inversion_size << "\n";
}

/**
 * \brief    Write the replication report data
 * \details  --
 * \param    std::ofstream& filestream
 * \return   \e void
 */
void ReplicationReport::write_replication_report_data( std::ofstream& filestream )
{
  filestream << _generation << " " <<
  _old_genome_size << " " <<
  _new_genome_size << " " <<
  _genome_functional_size << " " <<
  _genome_nb_NC << " " <<
  _genome_nb_E << " " <<
  _genome_nb_TF << " " <<
  _genome_nb_BS << " " <<
  _genome_nb_P << " " <<
  _genome_nb_inner_enzymes << " " <<
  _genome_nb_inflow_pumps << " " <<
  _genome_nb_outflow_pumps << " " <<
  _nb_functional_regions << " " <<
  _nb_enhancers << " " <<
  _nb_operators << " " <<
  _nb_E_regions << " " <<
  _nb_TF_regions << " " <<
  _nb_mixed_regions << " " <<
  _mean_functional_region_size << " " <<
  _mean_E_region_size << " " <<
  _mean_TF_region_size << " " <<
  _mean_mixed_region_size << " " <<
  _mean_enhancer_size << " " <<
  _mean_operator_size << " " <<
  _mean_operon_size << " " <<
  _mean_E_operon_size << " " <<
  _mean_TF_operon_size << " " <<
  _mean_mixed_operon_size << " " <<
  _mean_regulation_redundancy << " " <<
  _mean_metabolic_redundancy << " " <<
  _inherited_nb_E << " " <<
  _inherited_nb_TF << " " <<
  _inherited_nb_inner_enzymes << " " <<
  _inherited_nb_inflow_pumps << " " <<
  _inherited_nb_outflow_pumps << " " <<
  _id << " " <<
  _parent_id << " " <<
  _x << " " <<
  _y << " " <<
  _number_of_updates << " " <<
  _number_of_divisions << " " <<
  _birth_time << " " <<
  _death_time << " " <<
  _lifespan << " " <<
  _toxicity << " " <<
  _inherited_TF_amount << " " <<
  _inherited_E_amount << " " <<
  _TF_amount << " " <<
  _E_amount << " " <<
  _inherited_metabolic_amount << " " <<
  _min_metabolic_amount << " " <<
  _metabolic_amount << " " <<
  _max_metabolic_amount << " " <<
  _metabolic_uptake << " " <<
  _metabolic_release << " " <<
  _previous_metabolic_amount << " " <<
  _min_energy << " " <<
  _mean_energy << " " <<
  _max_energy << " " <<
  _min_score << " " <<
  _mean_score << " " <<
  _max_score << " " <<
  _metabolic_growth_rate << " " <<
  _diff_metabolic_growth_rate << " " <<
  _grn_nb_nodes << " " <<
  _grn_nb_edges << " " <<
  _metabolic_nb_nodes << " " <<
  _metabolic_nb_edges << " " <<
  _trophic_group << " " <<
  _trophic_level << " " <<
  _old_genome_size << " " <<
  _nb_point_mutations << " " <<
  _nb_NC_point_mutations << " " <<
  _nb_E_point_mutations << " " <<
  _nb_TF_point_mutations << " " <<
  _nb_BS_point_mutations << " " <<
  _nb_P_point_mutations << " " <<
  _nb_NC_to_E_transitions << " " <<
  _nb_NC_to_TF_transitions << " " <<
  _nb_NC_to_BS_transitions << " " <<
  _nb_NC_to_P_transitions << " " <<
  _nb_E_to_NC_transitions << " " <<
  _nb_E_to_TF_transitions << " " <<
  _nb_E_to_BS_transitions << " " <<
  _nb_E_to_P_transitions << " " <<
  _nb_TF_to_NC_transitions << " " <<
  _nb_TF_to_E_transitions << " " <<
  _nb_TF_to_BS_transitions << " " <<
  _nb_TF_to_P_transitions << " " <<
  _nb_BS_to_NC_transitions << " " <<
  _nb_BS_to_E_transitions << " " <<
  _nb_BS_to_TF_transitions << " " <<
  _nb_BS_to_P_transitions << " " <<
  _nb_P_to_NC_transitions << " " <<
  _nb_P_to_E_transitions << " " <<
  _nb_P_to_TF_transitions << " " <<
  _nb_P_to_BS_transitions << " " <<
  _mean_s_mutation_size << " " <<
  _mean_p_mutation_size << " " <<
  _mean_kcat_mutation_size << " " <<
  _mean_kcat_km_ratio_mutation_size << " " <<
  _mean_BS_tag_mutation_size << " " <<
  _mean_coE_tag_mutation_size << " " <<
  _mean_TF_tag_mutation_size << " " <<
  _mean_basal_expression_level_mutation_size << " " <<
  _nb_HGT << " " <<
  _mean_HGT_size << " " <<
  _nb_NC_HGT << " " <<
  _nb_E_HGT << " " <<
  _nb_TF_HGT << " " <<
  _nb_BS_HGT << " " <<
  _nb_P_HGT << " " <<
  _nb_rearrangements << " " <<
  _nb_duplicated_NC << " " <<
  _nb_duplicated_E << " " <<
  _nb_duplicated_TF << " " <<
  _nb_duplicated_BS << " " <<
  _nb_duplicated_P << " " <<
  _nb_deleted_NC << " " <<
  _nb_deleted_E << " " <<
  _nb_deleted_TF << " " <<
  _nb_deleted_BS << " " <<
  _nb_deleted_P << " " <<
  _nb_duplications << " " <<
  _nb_deletions << " " <<
  _nb_translocations << " " <<
  _nb_inversions << " " <<
  _mean_rearrangement_size << " " <<
  _mean_duplication_size << " " <<
  _mean_deletion_size << " " <<
  _mean_translocation_size << " " <<
  _mean_inversion_size << "\n";
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
