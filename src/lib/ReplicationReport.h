
/**
 * \file      ReplicationReport.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     ReplicationReport class declaration
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

#ifndef __Evo2Sim__ReplicationReport__
#define __Evo2Sim__ReplicationReport__

#include <iostream>
#include <fstream>
#include <vector>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "MutationEvent.h"


class ReplicationReport
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  ReplicationReport( void );
  ReplicationReport( gzFile backup_file );
  ReplicationReport( const ReplicationReport& report );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~ReplicationReport( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ Genome structure */
  
  inline size_t get_old_genome_size( void ) const;
  inline size_t get_new_genome_size( void ) const;
  inline size_t get_genome_functional_size( void ) const;
  inline size_t get_genome_nb_NC( void ) const;
  inline size_t get_genome_nb_E( void ) const;
  inline size_t get_genome_nb_TF( void ) const;
  inline size_t get_genome_nb_BS( void ) const;
  inline size_t get_genome_nb_P( void ) const;
  inline size_t get_genome_nb_inner_enzymes( void ) const;
  inline size_t get_genome_nb_inflow_pumps( void ) const;
  inline size_t get_genome_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ GRN data */
  
  inline size_t get_nb_functional_regions( void ) const;
  inline size_t get_nb_enhancers( void ) const;
  inline size_t get_nb_operators( void ) const;
  inline size_t get_nb_E_regions( void ) const;
  inline size_t get_nb_TF_regions( void ) const;
  inline size_t get_nb_mixed_regions( void ) const;
  inline double get_mean_functional_region_size( void ) const;
  inline double get_mean_E_region_size( void ) const;
  inline double get_mean_TF_region_size( void ) const;
  inline double get_mean_mixed_region_size( void ) const;
  inline double get_mean_enhancer_size( void ) const;
  inline double get_mean_operator_size( void ) const;
  inline double get_mean_operon_size( void ) const;
  inline double get_mean_E_operon_size( void ) const;
  inline double get_mean_TF_operon_size( void ) const;
  inline double get_mean_mixed_operon_size( void ) const;
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  inline double get_mean_regulation_redundancy( void ) const;
  inline double get_mean_metabolic_redundancy( void ) const;
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  inline size_t get_inherited_size( void ) const;
  inline size_t get_inherited_nb_E( void ) const;
  inline size_t get_inherited_nb_TF( void ) const;
  inline size_t get_inherited_nb_inner_enzymes( void ) const;
  inline size_t get_inherited_nb_inflow_pumps( void ) const;
  inline size_t get_inherited_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ Phenotype */
  
  inline unsigned long long int get_id( void ) const;
  inline unsigned long long int get_parent_id( void ) const;
  
  inline size_t get_generation( void ) const;
  inline size_t get_x( void ) const;
  inline size_t get_y( void ) const;
  inline size_t get_number_of_updates( void ) const;
  inline size_t get_number_of_divisions( void ) const;
  inline size_t get_birth_time( void ) const;
  inline size_t get_death_time( void ) const;
  inline size_t get_lifespan( void ) const;
  inline double get_toxicity( void ) const;
  inline double get_inherited_TF_amount( void ) const;
  inline double get_inherited_E_amount( void ) const;
  inline double get_TF_amount( void ) const;
  inline double get_E_amount( void ) const;
  inline double get_inherited_metabolic_amount( void ) const;
  inline double get_min_metabolic_amount( void ) const;
  inline double get_metabolic_amount( void ) const;
  inline double get_max_metabolic_amount( void ) const;
  inline double get_metabolic_uptake( void ) const;
  inline double get_metabolic_release( void ) const;
  inline double get_min_energy( void ) const;
  inline double get_mean_energy( void ) const;
  inline double get_max_energy( void ) const;
  inline double get_min_score( void ) const;
  inline double get_mean_score( void ) const;
  inline double get_max_score( void ) const;
  inline double get_metabolic_growth_rate( void ) const;
  inline double get_Dmetabolic_growth_rate( void ) const;
  inline size_t get_grn_nb_nodes( void ) const;
  inline size_t get_grn_nb_edges( void ) const;
  inline size_t get_metabolic_nb_nodes( void ) const;
  inline size_t get_metabolic_nb_edges( void ) const;
  
  inline unsigned long long int get_trophic_group( void ) const;
  inline trophic_level          get_trophic_level( void ) const;
  
  /*------------------------------------------------------------------ List of mutation events */
  
  inline size_t                       get_number_of_events( void ) const;
  inline std::vector<MutationEvent*>* get_list_of_events( void );
  
  /*------------------------------------------------------------------ Point mutations data */
  
  inline size_t get_nb_point_mutations( void ) const;
  
  inline size_t get_nb_NC_point_mutations( void ) const;
  inline size_t get_nb_E_point_mutations( void ) const;
  inline size_t get_nb_TF_point_mutations( void ) const;
  inline size_t get_nb_BS_point_mutations( void ) const;
  inline size_t get_nb_P_point_mutations( void ) const;
  
  inline size_t get_nb_NC_to_E_transitions( void ) const;
  inline size_t get_nb_NC_to_TF_transitions( void ) const;
  inline size_t get_nb_NC_to_BS_transitions( void ) const;
  inline size_t get_nb_NC_to_P_transitions( void ) const;
  
  inline size_t get_nb_E_to_NC_transitions( void ) const;
  inline size_t get_nb_E_to_TF_transitions( void ) const;
  inline size_t get_nb_E_to_BS_transitions( void ) const;
  inline size_t get_nb_E_to_P_transitions( void ) const;
  
  inline size_t get_nb_TF_to_NC_transitions( void ) const;
  inline size_t get_nb_TF_to_E_transitions( void ) const;
  inline size_t get_nb_TF_to_BS_transitions( void ) const;
  inline size_t get_nb_TF_to_P_transitions( void ) const;
  
  inline size_t get_nb_BS_to_NC_transitions( void ) const;
  inline size_t get_nb_BS_to_E_transitions( void ) const;
  inline size_t get_nb_BS_to_TF_transitions( void ) const;
  inline size_t get_nb_BS_to_P_transitions( void ) const;
  
  inline size_t get_nb_P_to_NC_transitions( void ) const;
  inline size_t get_nb_P_to_E_transitions( void ) const;
  inline size_t get_nb_P_to_TF_transitions( void ) const;
  inline size_t get_nb_P_to_BS_transitions( void ) const;
  
  inline double get_mean_s_mutation_size( void ) const;
  inline double get_mean_p_mutation_size( void ) const;
  inline double get_mean_kcat_mutation_size( void ) const;
  inline double get_mean_kcat_km_ratio_mutation_size( void ) const;
  
  inline double get_mean_BS_tag_mutation_size( void ) const;
  inline double get_mean_coE_tag_mutation_size( void ) const;
  
  inline double get_mean_TF_tag_mutation_size( void ) const;
  
  inline double get_mean_basal_expression_level_mutation_size( void ) const;
  
  /*------------------------------------------------------------------ HGT data */
  
  inline size_t get_nb_HGT( void ) const;
  inline double get_mean_HGT_size( void ) const;
  inline size_t get_nb_NC_HGT( void ) const;
  inline size_t get_nb_E_HGT( void ) const;
  inline size_t get_nb_TF_HGT( void ) const;
  inline size_t get_nb_BS_HGT( void ) const;
  inline size_t get_nb_P_HGT( void ) const;
  
  /*------------------------------------------------------------------ Rearrangements data */
  
  inline size_t get_nb_rearrangements( void ) const;
  
  inline size_t get_nb_duplicated_NC( void ) const;
  inline size_t get_nb_duplicated_E( void ) const;
  inline size_t get_nb_duplicated_TF( void ) const;
  inline size_t get_nb_duplicated_BS( void ) const;
  inline size_t get_nb_duplicated_P( void ) const;

  inline size_t get_nb_deleted_NC( void ) const;
  inline size_t get_nb_deleted_E( void ) const;
  inline size_t get_nb_deleted_TF( void ) const;
  inline size_t get_nb_deleted_BS( void ) const;
  inline size_t get_nb_deleted_P( void ) const;
  
  inline size_t get_nb_duplications( void ) const;
  inline size_t get_nb_deletions( void ) const;
  inline size_t get_nb_translocations( void ) const;
  inline size_t get_nb_inversions( void ) const;
  
  inline double get_mean_rearrangement_size( void ) const;
  inline double get_mean_duplication_size( void ) const;
  inline double get_mean_deletion_size( void ) const;
  inline double get_mean_translocation_size( void ) const;
  inline double get_mean_inversion_size( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ Genome structure */
  
  inline void set_old_genome_size( size_t old_size );
  inline void set_new_genome_size( size_t new_size );
  inline void set_genome_functional_size( size_t functional_size );
  inline void set_genome_nb_NC( size_t nb_NC );
  inline void set_genome_nb_E( size_t nb_E );
  inline void set_genome_nb_TF( size_t nb_TF );
  inline void set_genome_nb_BS( size_t nb_BS );
  inline void set_genome_nb_P( size_t nb_P );
  inline void set_genome_nb_inner_enzymes( size_t nb_inner_enzymes );
  inline void set_genome_nb_inflow_pumps( size_t nb_inflow_pumps );
  inline void set_genome_nb_outflow_pumps( size_t nb_outflow_pumps );
  
  /*------------------------------------------------------------------ GRN data */
  
  inline void set_nb_functional_regions( size_t nb_functional_regions );
  inline void set_nb_enhancers( size_t nb_enhancers );
  inline void set_nb_operators( size_t nb_operators );
  inline void set_nb_E_regions( size_t nb_E_regions );
  inline void set_nb_TF_regions( size_t nb_TF_regions );
  inline void set_nb_mixed_regions( size_t nb_mixed_regions );
  inline void set_mean_functional_region_size( double mean_functional_region_size );
  inline void set_mean_E_region_size( double mean_E_region_size );
  inline void set_mean_TF_region_size( double mean_TF_region_size );
  inline void set_mean_mixed_region_size( double mean_mixed_region_size );
  inline void set_mean_enhancer_size( double mean_enhancer_size );
  inline void set_mean_operator_size( double mean_operator_size );
  inline void set_mean_operon_size( double mean_operon_size );
  inline void set_mean_E_operon_size( double mean_E_operon_size );
  inline void set_mean_TF_operon_size( double mean_TF_operon_size );
  inline void set_mean_mixed_operon_size( double mean_mixed_operon_size );
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  inline void set_mean_regulation_redundancy( double regulation_redundancy );
  inline void set_mean_metabolic_redundancy( double metabolic_redundancy );
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  inline void set_inherited_size( size_t size );
  inline void set_inherited_nb_E( size_t nb_E );
  inline void set_inherited_nb_TF( size_t nb_TF );
  inline void set_inherited_nb_inner_enzymes( size_t nb_inner_enzymes );
  inline void set_inherited_nb_inflow_pumps( size_t nb_inflow_pumps );
  inline void set_inherited_nb_outflow_pumps( size_t nb_outflow_pumps );
  
  /*------------------------------------------------------------------ Phenotype */
  
  inline void set_id( unsigned long long int identifier );
  inline void set_parent_id( unsigned long long int parent_identifier );
  
  inline void set_generation( size_t generation );
  inline void set_x( size_t x );
  inline void set_y( size_t y );
  inline void set_number_of_updates( size_t number_of_updates );
  inline void set_number_of_divisions( size_t number_of_divisions );
  inline void set_birth_time( size_t birth_time );
  inline void set_death_time( size_t death_time );
  inline void set_lifespan( size_t lifespan );
  inline void set_toxicity( double toxicity );
  inline void set_inherited_TF_amount( double inherited_TF_amount );
  inline void set_inherited_E_amount( double inherited_E_amount );
  inline void set_TF_amount( double TF_amount );
  inline void set_E_amount( double E_amount );
  inline void set_inherited_metabolic_amount( double inherited_metabolic_amount );
  inline void set_min_metabolic_amount( double min_metabolic_amount );
  inline void set_metabolic_amount( double metabolic_amount );
  inline void set_max_metabolic_amount( double max_metabolic_amount );
  inline void set_metabolic_uptake( double metabolic_uptake );
  inline void set_metabolic_release( double metabolic_release );
  inline void set_min_energy( double min_energy );
  inline void set_mean_energy( double mean_energy );
  inline void set_max_energy( double max_energy );
  inline void set_min_score( double min_score );
  inline void set_mean_score( double mean_score );
  inline void set_max_score( double max_score );
  inline void set_metabolic_growth_rate( double metabolic_growth_rate );
  inline void set_Dmetabolic_growth_rate( double Dmetabolic_growth_rate );
  inline void set_grn_nb_nodes( size_t grn_nb_nodes );
  inline void set_grn_nb_edges( size_t grn_nb_edges );
  inline void set_metabolic_nb_nodes( size_t metabolic_nb_nodes );
  inline void set_metabolic_nb_edges( size_t metabolic_nb_edges  );
  
  inline void set_trophic_group( unsigned long long int group );
  inline void set_trophic_level( trophic_level level );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void add_mutation_event( MutationEvent* event );
  void compute_mean( void );
  void save( gzFile backup_file );
  void clear( void );
  void clear_event_list( void );
  
  void write_genome_structure_header( std::ofstream& filestream );
  void write_inherited_proteins_header( std::ofstream& filestream );
  void write_phenotype_header( std::ofstream& filestream );
  void write_fixed_mutations_header( std::ofstream& filestream );
  void write_replication_report_header( std::ofstream& filestream );
  
  void write_genome_structure_data( std::ofstream& filestream );
  void write_inherited_proteins_data( std::ofstream& filestream );
  void write_phenotype_data( std::ofstream& filestream );
  void write_fixed_mutations_data( std::ofstream& filestream );
  void write_replication_report_data( std::ofstream& filestream );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ Genome structure (11) */
  
  size_t _old_genome_size;         /*!< Old genome size (before mutations) */
  size_t _new_genome_size;         /*!< New genome size (after mutations)  */
  size_t _genome_functional_size;  /*!< Total size of functional regions   */
  size_t _genome_nb_NC;            /*!< Number of genome NC units          */
  size_t _genome_nb_E;             /*!< Number of genome E units           */
  size_t _genome_nb_TF;            /*!< Number of genome TF units          */
  size_t _genome_nb_BS;            /*!< Number of genome BS units          */
  size_t _genome_nb_P;             /*!< Number of genome NC units          */
  size_t _genome_nb_inner_enzymes; /*!< Number of genome inner enzymes     */
  size_t _genome_nb_inflow_pumps;  /*!< Number of genome inflowing pumps   */
  size_t _genome_nb_outflow_pumps; /*!< Number of genome outflowing pumps  */
  
  /*------------------------------------------------------------------ GRN data (16) */
  
  size_t _nb_functional_regions;       /*!< Number of functional regions                    */
  size_t _nb_enhancers;                /*!< Number of enhancers in functional regions       */
  size_t _nb_operators;                /*!< Number of operators in functional regions       */
  size_t _nb_E_regions;                /*!< Number of functional regions containing only TF */
  size_t _nb_TF_regions;               /*!< Number of functional regions containing only TF */
  size_t _nb_mixed_regions;            /*!< Number of functional regions mixing E and TF    */
  double _mean_functional_region_size; /*!< Mean size of functional regions                 */
  double _mean_E_region_size;          /*!< Mean size of E regions                          */
  double _mean_TF_region_size;         /*!< Mean size of TF regions                         */
  double _mean_mixed_region_size;      /*!< Mean size of mixed regions                      */
  double _mean_enhancer_size;          /*!< Mean size of enhancer sites                     */
  double _mean_operator_size;          /*!< Mean size of operator sites                     */
  double _mean_operon_size;            /*!< Mean size of operons                            */
  double _mean_E_operon_size;          /*!< Mean size of operons containing only E          */
  double _mean_TF_operon_size;         /*!< Mean size of operons containing only TF         */
  double _mean_mixed_operon_size;      /*!< Mean size of operons containing mixing E and TF */
  
  /*------------------------------------------------------------------ Genetic redundancy (2) */
  
  double _mean_regulation_redundancy; /*!< Mean redundancy for TF units */
  double _mean_metabolic_redundancy;  /*!< Mean redundancy for E units  */
  
  /*------------------------------------------------------------------ Inherited proteins structure (6) */
  
  size_t _inherited_size;             /*!< Inherited size                       */
  size_t _inherited_nb_E;             /*!< Number of inherited E units          */
  size_t _inherited_nb_TF;            /*!< Number of inherited TF units         */
  size_t _inherited_nb_inner_enzymes; /*!< Number of inherited inner enzymes    */
  size_t _inherited_nb_inflow_pumps;  /*!< Number of inherited inflowing pumps  */
  size_t _inherited_nb_outflow_pumps; /*!< Number of inherited outflowing pumps */
  
  /*------------------------------------------------------------------ Phenotype (36) */
  
  unsigned long long int _id;            /*!< Identifier                                    */
  unsigned long long int _parent_id;     /*!< Parental identifier                           */
  size_t _generation;                    /*!< Generation of the cell                        */
  size_t _x;                             /*!< Cell's x coordinate                           */
  size_t _y;                             /*!< Cell's y coordinate                           */
  size_t _number_of_updates;             /*!< Number of updates since birth                 */
  size_t _number_of_divisions;           /*!< Number of divisions since birth               */
  size_t _birth_time;                    /*!< Cell's birth time                             */
  size_t _death_time;                    /*!< Cell's death time                             */
  size_t _lifespan;                      /*!< Cell's lifespan                               */
  double _toxicity;                      /*!< Toxicity accumulation                         */
  double _inherited_TF_amount;           /*!< Inherited transcription factors amount        */
  double _inherited_E_amount;            /*!< Inherited enzymes amount                      */
  double _TF_amount;                     /*!< Transcription factors amount                  */
  double _E_amount;                      /*!< Enzymes amount                                */
  double _inherited_metabolic_amount;    /*!< Inherited metabolic amount                    */
  double _min_metabolic_amount;          /*!< Minimum cytoplasmic metabolic amount          */
  double _metabolic_amount;              /*!< Cytoplasmic metabolic amount                  */
  double _max_metabolic_amount;          /*!< Maximum cytoplasmic metabolic amount          */
  double _metabolic_uptake;              /*!< Amount of metabolic uptake                    */
  double _metabolic_release;             /*!< Amount of released metabolites                */
  double _previous_metabolic_amount;     /*!< Previous metabolic amount                     */
  double _min_energy;                    /*!< Minimum energy bilan of the cell through time */
  double _mean_energy;                   /*!< Mean energy bilan of the cell through time    */
  double _max_energy;                    /*!< Maximum energy bilan of the cell through time */
  double _min_score;                     /*!< Minimum score of the cell through time        */
  double _mean_score;                    /*!< Mean score of the cell through time           */
  double _max_score;                     /*!< Maximum score of the cell through time        */
  double _metabolic_growth_rate;         /*!< Cell's metabolic growth rate                  */
  double _diff_metabolic_growth_rate;    /*!< Cell's metabolic growth rate difference       */
  size_t _grn_nb_nodes;                  /*!< Number of nodes in the GRN                    */
  size_t _grn_nb_edges;                  /*!< Number of edges in the GRN                    */
  size_t _metabolic_nb_nodes;            /*!< Number of nodes in the metabolic network      */
  size_t _metabolic_nb_edges;            /*!< Number of edges in the metabolic network      */
  unsigned long long int _trophic_group; /*!< Trophic group index in the trophic network    */
  trophic_level          _trophic_level; /*!< Trophic level                                 */
  
  /*------------------------------------------------------------------ List of mutation events (1) */
  
  std::vector<MutationEvent*> _list_of_events; /*!< List of mutation events */
  
  /*------------------------------------------------------------------ Point mutations data (34) */
  
  size_t _nb_point_mutations; /*!< Number of point mutations */
  
  size_t _nb_NC_point_mutations; /*!< Number of point mutations in NC type */
  size_t _nb_E_point_mutations;  /*!< Number of point mutations in E type  */
  size_t _nb_TF_point_mutations; /*!< Number of point mutations in TF type */
  size_t _nb_BS_point_mutations; /*!< Number of point mutations in BS type */
  size_t _nb_P_point_mutations;  /*!< Number of point mutations in P type  */
  
  size_t _nb_NC_to_E_transitions;  /*!< Number of transitions from NC ro E type  */
  size_t _nb_NC_to_TF_transitions; /*!< Number of transitions from NC ro TF type */
  size_t _nb_NC_to_BS_transitions; /*!< Number of transitions from NC ro BS type */
  size_t _nb_NC_to_P_transitions;  /*!< Number of transitions from NC ro P type  */
  
  size_t _nb_E_to_NC_transitions; /*!< Number of transitions from E ro NC type */
  size_t _nb_E_to_TF_transitions; /*!< Number of transitions from E ro TF type */
  size_t _nb_E_to_BS_transitions; /*!< Number of transitions from E ro BS type */
  size_t _nb_E_to_P_transitions;  /*!< Number of transitions from E ro P type  */
  
  size_t _nb_TF_to_NC_transitions; /*!< Number of transitions from TF ro NC type */
  size_t _nb_TF_to_E_transitions;  /*!< Number of transitions from TF ro E type  */
  size_t _nb_TF_to_BS_transitions; /*!< Number of transitions from TF ro BS type */
  size_t _nb_TF_to_P_transitions;  /*!< Number of transitions from TF ro P type  */
  
  size_t _nb_BS_to_NC_transitions; /*!< Number of transitions from BS ro NC type */
  size_t _nb_BS_to_E_transitions;  /*!< Number of transitions from BS ro E type  */
  size_t _nb_BS_to_TF_transitions; /*!< Number of transitions from BS ro TF type */
  size_t _nb_BS_to_P_transitions;  /*!< Number of transitions from BS ro P type  */
  
  size_t _nb_P_to_NC_transitions; /*!< Number of transitions from P ro NC type */
  size_t _nb_P_to_E_transitions;  /*!< Number of transitions from P ro E type  */
  size_t _nb_P_to_TF_transitions; /*!< Number of transitions from P ro TF type */
  size_t _nb_P_to_BS_transitions; /*!< Number of transitions from P ro BS type */
  
  double _mean_s_mutation_size;             /*!< Mean substrate tag mutation sizes          */
  double _mean_p_mutation_size;             /*!< Mean product tag mutation sizes            */
  double _mean_kcat_mutation_size;          /*!< Mean Kcat constant mutation sizes          */
  double _mean_kcat_km_ratio_mutation_size; /*!< Mean Kcat/Km ratio constant mutation sizes */
  
  double _mean_BS_tag_mutation_size;  /*!< Mean BS tag constant mutation sizes   */
  double _mean_coE_tag_mutation_size; /*!< Mean coE tag constant mutation sizes  */
  
  double _mean_TF_tag_mutation_size; /*!< Mean TF tag constant mutation sizes  */
  
  double _mean_basal_expression_level_mutation_size; /*!< Mean basal expression level mutation size */
  
  /*------------------------------------------------------------------ HGT data (7) */
  
  size_t _nb_HGT;        /*!< Number of HGT         */
  double _mean_HGT_size; /*!< Mean size of HGT      */
  size_t _nb_NC_HGT;     /*!< Number of NC unit HGT */
  size_t _nb_E_HGT;      /*!< Number of E unit HGT  */
  size_t _nb_TF_HGT;     /*!< Number of TF unit HGT */
  size_t _nb_BS_HGT;     /*!< Number of BS unit HGT */
  size_t _nb_P_HGT;      /*!< Number of P unit HGT  */
  
  /*------------------------------------------------------------------ Rearrangements data (20) */
  
  size_t _nb_rearrangements; /*!< Number of rearrangements     */
  
  size_t _nb_duplicated_NC; /*!< Number of duplicated NC type */
  size_t _nb_duplicated_E;  /*!< Number of duplicated E type  */
  size_t _nb_duplicated_TF; /*!< Number of duplicated TF type */
  size_t _nb_duplicated_BS; /*!< Number of duplicated BS type */
  size_t _nb_duplicated_P;  /*!< Number of duplicated P type  */

  size_t _nb_deleted_NC; /*!< Number of deleted NC type */
  size_t _nb_deleted_E;  /*!< Number of deleted E type  */
  size_t _nb_deleted_TF; /*!< Number of deleted TF type */
  size_t _nb_deleted_BS; /*!< Number of deleted BS type */
  size_t _nb_deleted_P;  /*!< Number of deleted P type  */
  
  size_t _nb_duplications;   /*!< Number of duplications   */
  size_t _nb_deletions;      /*!< Number of deletions      */
  size_t _nb_translocations; /*!< Number of translocations */
  size_t _nb_inversions;     /*!< Number of inversions     */
  
  double _mean_rearrangement_size; /*!< Mean size of rearrangements */
  double _mean_duplication_size;   /*!< Mean size of duplications   */
  double _mean_deletion_size;      /*!< Mean size of deletions      */
  double _mean_translocation_size; /*!< Mean size of translocations */
  double _mean_inversion_size;     /*!< Mean size of inversions     */
  
};

/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ Genome structure */

/**
 * \brief    Get the old genome size (before mutations)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_old_genome_size( void ) const
{
  return _old_genome_size;
}

/**
 * \brief    Get the new genome size (after mutations)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_new_genome_size( void ) const
{
  return _new_genome_size;
}

/**
 * \brief    Get the genome functional size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_functional_size( void ) const
{
  return _genome_functional_size;
}

/**
 * \brief    Get the number of NC in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_NC( void ) const
{
  return _genome_nb_NC;
}

/**
 * \brief    Get the number of E in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_E( void ) const
{
  return _genome_nb_E;
}

/**
 * \brief    Get the number of TF in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_TF( void ) const
{
  return _genome_nb_TF;
}

/**
 * \brief    Get the number of BS in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_BS( void ) const
{
  return _genome_nb_BS;
}

/**
 * \brief    Get the number of P in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_P( void ) const
{
  return _genome_nb_P;
}

/**
 * \brief    Get the number of inner enzymes in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_inner_enzymes( void ) const
{
  return _genome_nb_inner_enzymes;
}

/**
 * \brief    Get the number of inflowing pumps in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_inflow_pumps( void ) const
{
  return _genome_nb_inflow_pumps;
}

/**
 * \brief    Get the number of outflowing pumps in the genome
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_genome_nb_outflow_pumps( void ) const
{
  return _genome_nb_outflow_pumps;
}

/*------------------------------------------------------------------ GRN data */

/**
 * \brief    Get the number of functional regions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_functional_regions( void ) const
{
  return _nb_functional_regions;
}

/**
 * \brief    Get the number of enhancers in functional regions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_enhancers( void ) const
{
  return _nb_enhancers;
}

/**
 * \brief    Get the number of operators in functional regions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_operators( void ) const
{
  return _nb_operators;
}

/**
 * \brief    Get the number of functional regions containing only E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_regions( void ) const
{
  return _nb_E_regions;
}

/**
 * \brief    Get the number of functional regions containing only TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_regions( void ) const
{
  return _nb_TF_regions;
}

/**
 * \brief    Get the number of functional regions mixing E and TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_mixed_regions( void ) const
{
  return _nb_mixed_regions;
}

/**
 * \brief    Get the mean size of functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_functional_region_size( void ) const
{
  return _mean_functional_region_size;
}

/**
 * \brief    Get the mean size of functional regions containing only E types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_E_region_size( void ) const
{
  return _mean_E_region_size;
}

/**
 * \brief    Get the mean size of functional regions containing only TF types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_TF_region_size( void ) const
{
  return _mean_TF_region_size;
}

/**
 * \brief    Get the mean size of functional regions mixing E and TF types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_mixed_region_size( void ) const
{
  return _mean_mixed_region_size;
}

/**
 * \brief    Get the mean size of enhancers in functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_enhancer_size( void ) const
{
  return _mean_enhancer_size;
}

/**
 * \brief    Get the mean size of operators in functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_operator_size( void ) const
{
  return _mean_operator_size;
}

/**
 * \brief    Get the mean size of operons
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_operon_size( void ) const
{
  return _mean_operon_size;
}

/**
 * \brief    Get the mean size of operons containing only E types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_E_operon_size( void ) const
{
  return _mean_E_operon_size;
}

/**
 * \brief    Get the mean size of operons containing only TF types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_TF_operon_size( void ) const
{
  return _mean_TF_operon_size;
}

/**
 * \brief    Get the mean size of operons mixing E and TF types
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_mixed_operon_size( void ) const
{
  return _mean_mixed_operon_size;
}

/*------------------------------------------------------------------ Genetic redundancy */

/**
 * \brief    Get mean regulation redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_regulation_redundancy( void ) const
{
  return _mean_regulation_redundancy;
}

/**
 * \brief    Get mean metabolic redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_metabolic_redundancy( void ) const
{
  return _mean_metabolic_redundancy;
}

/*------------------------------------------------------------------ Inherited proteins structure */

/**
 * \brief    Get the inherited proteins size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_size( void ) const
{
  return _inherited_size;
}

/**
 * \brief    Get the number of E in the inherited proteins
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_nb_E( void ) const
{
  return _inherited_nb_E;
}

/**
 * \brief    Get the number of TF in the inherited proteins
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_nb_TF( void ) const
{
  return _inherited_nb_TF;
}

/**
 * \brief    Get the number of inner enzymes in the inherited proteins
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_nb_inner_enzymes( void ) const
{
  return _inherited_nb_inner_enzymes;
}

/**
 * \brief    Get the number of inflowing pumps in the inherited proteins
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_nb_inflow_pumps( void ) const
{
  return _inherited_nb_inflow_pumps;
}

/**
 * \brief    Get the number of outflowing pumps in the inherited proteins
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_inherited_nb_outflow_pumps( void ) const
{
  return _inherited_nb_outflow_pumps;
}

/*------------------------------------------------------------------ Phenotype */

/**
 * \brief    Get cell id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int ReplicationReport::get_id( void ) const
{
  return _id;
}

/*
 * \brief    Get parental cell id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int ReplicationReport::get_parent_id( void ) const
{
  return _parent_id;
}

/**
 * \brief    Get cell generation
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_generation( void ) const
{
  return _generation;
}

/**
 * \brief    Get cell x coordinate
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_x( void ) const
{
  return _x;
}

/**
 * \brief    Get cell y coordinate
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_y( void ) const
{
  return _y;
}

/**
 * \brief    Get the number of updates
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_number_of_updates( void ) const
{
  return _number_of_updates;
}

/**
 * \brief    Get the number of divisions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_number_of_divisions( void ) const
{
  return _number_of_divisions;
}

/**
 * \brief    Get cell birth time
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_birth_time( void ) const
{
  return _birth_time;
}

/**
 * \brief    Get cell death time
 * \details  --
 * \param    size_t
 * \return   \e double
 */
inline size_t ReplicationReport::get_death_time( void ) const
{
  return _death_time;
}

/**
 * \brief    Get cell lifespan
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_lifespan( void ) const
{
  return _lifespan;
}

/**
 * \brief    Get the toxicity accumulation
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_toxicity( void ) const
{
  return _toxicity;
}

/**
 * \brief    Get the inherited transcription factors amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_inherited_TF_amount( void ) const
{
  return _inherited_TF_amount;
}

/**
 * \brief    Get the inherited enzymes amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_inherited_E_amount( void ) const
{
  return _inherited_E_amount;
}

/**
 * \brief    Get the transcription factors amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_TF_amount( void ) const
{
  return _TF_amount;
}

/**
 * \brief    Get the enzymes amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_E_amount( void ) const
{
  return _E_amount;
}

/**
 * \brief    Get the metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_inherited_metabolic_amount( void ) const
{
  return _inherited_metabolic_amount;
}

/**
 * \brief    Get the minimum metabolic amount during cell's life
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_min_metabolic_amount( void ) const
{
  return _min_metabolic_amount;
}

/**
 * \brief    Get the cell's metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_metabolic_amount( void ) const
{
  return _metabolic_amount;
}

/**
 * \brief    Get the maximum metabolic amount during cell's life
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_max_metabolic_amount( void ) const
{
  return _max_metabolic_amount;
}

/**
 * \brief    Get amount of metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_metabolic_uptake( void ) const
{
  return _metabolic_uptake;
}

/**
 * \brief    Get amount of released metabolites
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_metabolic_release( void ) const
{
  return _metabolic_release;
}

/**
 * \brief    Get cell minimum energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_min_energy( void ) const
{
  return _min_energy;
}

/**
 * \brief    Get cell mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_energy( void ) const
{
  return _mean_energy;
}

/**
 * \brief    Get cell maximum energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_max_energy( void ) const
{
  return _max_energy;
}

/**
 * \brief    Get cell minimum score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_min_score( void ) const
{
  return _min_score;
}

/**
 * \brief    Get cell mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_score( void ) const
{
  return _mean_score;
}

/**
 * \brief    Get cell maximum score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_max_score( void ) const
{
  return _max_score;
}

/**
 * \brief    Get the metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_metabolic_growth_rate( void ) const
{
  return _metabolic_growth_rate;
}

/**
 * \brief    Get the metabolic growth rate difference
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_Dmetabolic_growth_rate( void ) const
{
  return _diff_metabolic_growth_rate;
}

/**
 * \brief    Get the number of nodes in the GRN
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_grn_nb_nodes( void ) const
{
  return _grn_nb_nodes;
}

/**
 * \brief    Get the number of edges in the GRN
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_grn_nb_edges( void ) const
{
  return _grn_nb_edges;
}

/**
 * \brief    Get the number of nodes in the metabolic network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_metabolic_nb_nodes( void ) const
{
  return _metabolic_nb_nodes;
}

/**
 * \brief    Get the number of edges in the metabolic network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_metabolic_nb_edges( void ) const
{
  return _metabolic_nb_edges;
}

/**
 * \brief    Get the trophic group index in the trophic network
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int ReplicationReport::get_trophic_group( void ) const
{
  return _trophic_group;
}

/**
 * \brief    Get the trophic level
 * \details  --
 * \param    void
 * \return   \e trophic_level
 */
inline trophic_level ReplicationReport::get_trophic_level( void ) const
{
  return _trophic_level;
}

/*------------------------------------------------------------------ List of mutation events */

/**
 * \brief    Get number of events
 * \details  Get number of events listed in the replication report
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_number_of_events( void ) const
{
  return _list_of_events.size();
}

/**
 * \brief    Get the list of events
 * \details  Get the list of mutational events
 * \param    void
 * \return   \e std::vector<MutationEvent*>*
 */
inline std::vector<MutationEvent*>* ReplicationReport::get_list_of_events( void )
{
  return &_list_of_events;
}

/*------------------------------------------------------------------ Point mutations data */

/**
 * \brief    Get number of point mutations
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_point_mutations( void ) const
{
  return _nb_point_mutations;
}

/**
 * \brief    Get number of point mutations affecting NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_point_mutations( void ) const
{
  return _nb_NC_point_mutations;
}

/**
 * \brief    Get number of point mutations affecting E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_point_mutations( void ) const
{
  return _nb_E_point_mutations;
}

/**
 * \brief    Get number of point mutations affecting TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_point_mutations( void ) const
{
  return _nb_TF_point_mutations;
}

/**
 * \brief    Get number of point mutations affecting BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_point_mutations( void ) const
{
  return _nb_BS_point_mutations;
}

/**
 * \brief    Get number of point mutations affecting P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_point_mutations( void ) const
{
  return _nb_P_point_mutations;
}

/**
 * \brief    Get number of transitions from NC to E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_to_E_transitions( void ) const
{
  return _nb_NC_to_E_transitions;
}

/**
 * \brief    Get number of transitions from NC to TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_to_TF_transitions( void ) const
{
  return _nb_NC_to_TF_transitions;
}

/**
 * \brief    Get number of transitions from NC to BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_to_BS_transitions( void ) const
{
  return _nb_NC_to_BS_transitions;
}

/**
 * \brief    Get number of transitions from NC to P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_to_P_transitions( void ) const
{
  return _nb_NC_to_P_transitions;
}

/**
 * \brief    Get number of transitions from E to NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_to_NC_transitions( void ) const
{
  return _nb_E_to_NC_transitions;
}

/**
 * \brief    Get number of transitions from E to TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_to_TF_transitions( void ) const
{
  return _nb_E_to_TF_transitions;
}

/**
 * \brief    Get number of transitions from E to BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_to_BS_transitions( void ) const
{
  return _nb_E_to_BS_transitions;
}

/**
 * \brief    Get number of transitions from E to P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_to_P_transitions( void ) const
{
  return _nb_E_to_P_transitions;
}

/**
 * \brief    Get number of transitions from TF to NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_to_NC_transitions( void ) const
{
  return _nb_TF_to_NC_transitions;
}

/**
 * \brief    Get number of transitions from TF to E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_to_E_transitions( void ) const
{
  return _nb_TF_to_E_transitions;
}

/**
 * \brief    Get number of transitions from TF to BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_to_BS_transitions( void ) const
{
  return _nb_TF_to_BS_transitions;
}

/**
 * \brief    Get number of transitions from TF to P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_to_P_transitions( void ) const
{
  return _nb_TF_to_P_transitions;
}

/**
 * \brief    Get number of transitions from BS to NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_to_NC_transitions( void ) const
{
  return _nb_BS_to_NC_transitions;
}

/**
 * \brief    Get number of transitions from BS to E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_to_E_transitions( void ) const
{
  return _nb_BS_to_E_transitions;
}

/**
 * \brief    Get number of transitions from BS to TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_to_TF_transitions( void ) const
{
  return _nb_BS_to_TF_transitions;
}

/**
 * \brief    Get number of transitions from BS to P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_to_P_transitions( void ) const
{
  return _nb_BS_to_P_transitions;
}

/**
 * \brief    Get number of transitions from P to NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_to_NC_transitions( void ) const
{
  return _nb_P_to_NC_transitions;
}

/**
 * \brief    Get number of transitions from P to E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_to_E_transitions( void ) const
{
  return _nb_P_to_E_transitions;
}

/**
 * \brief    Get number of transitions from P to TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_to_TF_transitions( void ) const
{
  return _nb_P_to_TF_transitions;
}

/**
 * \brief    Get number of transitions from P to BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_to_BS_transitions( void ) const
{
  return _nb_P_to_BS_transitions;
}

/**
 * \brief    Get the mean size of substrate tag mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_s_mutation_size( void ) const
{
  return _mean_s_mutation_size;
}

/**
 * \brief    Get the mean size of product tag mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_p_mutation_size( void ) const
{
  return _mean_p_mutation_size;
}

/**
 * \brief    Get the mean size of Kcat constant mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_kcat_mutation_size( void ) const
{
  return _mean_kcat_mutation_size;
}

/**
 * \brief    Get the mean size of Kcat/Km ratio constant mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_kcat_km_ratio_mutation_size( void ) const
{
  return _mean_kcat_km_ratio_mutation_size;
}

/**
 * \brief    Get the mean size of BS tag mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_BS_tag_mutation_size( void ) const
{
  return _mean_BS_tag_mutation_size;
}

/**
 * \brief    Get the mean size of coE tag mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_coE_tag_mutation_size( void ) const
{
  return _mean_coE_tag_mutation_size;
}

/**
 * \brief    Get the mean size of TF tag mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_TF_tag_mutation_size( void ) const
{
  return _mean_TF_tag_mutation_size;
}

/**
 * \brief    Get the mean size of basal expression level mutations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_basal_expression_level_mutation_size( void ) const
{
  return _mean_basal_expression_level_mutation_size;
}

/*------------------------------------------------------------------ HGT data */

/**
 * \brief    Get the number of HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_HGT( void ) const
{
  return _nb_HGT;
}

/**
 * \brief    Get the mena size of HGT
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_HGT_size( void ) const
{
  return _mean_HGT_size;
}

/**
 * \brief    Get the number of NC unit HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_NC_HGT( void ) const
{
  return _nb_NC_HGT;
}

/**
 * \brief    Get the number of E unit HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_E_HGT( void ) const
{
  return _nb_E_HGT;
}

/**
 * \brief    Get the number of TF unit HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_TF_HGT( void ) const
{
  return _nb_TF_HGT;
}

/**
 * \brief    Get the number of BS unit HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_BS_HGT( void ) const
{
  return _nb_BS_HGT;
}

/**
 * \brief    Get the number of P unit HGT
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_P_HGT( void ) const
{
  return _nb_P_HGT;
}

/*------------------------------------------------------------------ Rearrangements data */

/**
 * \brief    Get number of rearrangements
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_rearrangements( void ) const
{
  return _nb_rearrangements;
}

/**
 * \brief    Get number of duplicated NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplicated_NC( void ) const
{
  return _nb_duplicated_NC;
}

/**
 * \brief    Get number of duplicated E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplicated_E( void ) const
{
  return _nb_duplicated_E;
}

/**
 * \brief    Get number of duplicated TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplicated_TF( void ) const
{
  return _nb_duplicated_TF;
}

/**
 * \brief    Get number of duplicated BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplicated_BS( void ) const
{
  return _nb_duplicated_BS;
}

/**
 * \brief    Get number of duplicated P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplicated_P( void ) const
{
  return _nb_duplicated_P;
}

/**
 * \brief    Get number of deleted NC types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deleted_NC( void ) const
{
  return _nb_deleted_NC;
}

/**
 * \brief    Get number of deleted E types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deleted_E( void ) const
{
  return _nb_deleted_E;
}

/**
 * \brief    Get number of deleted TF types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deleted_TF( void ) const
{
  return _nb_deleted_TF;
}

/**
 * \brief    Get number of deleted BS types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deleted_BS( void ) const
{
  return _nb_deleted_BS;
}

/**
 * \brief    Get number of deleted P types
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deleted_P( void ) const
{
  return _nb_deleted_P;
}

/**
 * \brief    Get number of duplications
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_duplications( void ) const
{
  return _nb_duplications;
}

/**
 * \brief    Get number of deletions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_deletions( void ) const
{
  return _nb_deletions;
}

/**
 * \brief    Get number of translocations
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_translocations( void ) const
{
  return _nb_translocations;
}

/**
 * \brief    Get number of inversions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t ReplicationReport::get_nb_inversions( void ) const
{
  return _nb_inversions;
}

/**
 * \brief    Get mean size of rearrangements
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_rearrangement_size( void ) const
{
  return _mean_rearrangement_size;
}

/**
 * \brief    Get mean size of duplications
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_duplication_size( void ) const
{
  return _mean_duplication_size;
}

/**
 * \brief    Get mean size of deletions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_deletion_size( void ) const
{
  return _mean_deletion_size;
}

/**
 * \brief    Get mean size of translocations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_translocation_size( void ) const
{
  return _mean_translocation_size;
}

/**
 * \brief    Get mean size of inversions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double ReplicationReport::get_mean_inversion_size( void ) const
{
  return _mean_inversion_size;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*------------------------------------------------------------------ Genome structure */

/**
 * \brief    Set the old genome size (before mutations)
 * \details  --
 * \param    size_t old_genome_size
 * \return   \e void
 */
inline void ReplicationReport::set_old_genome_size( size_t old_genome_size )
{
  _old_genome_size = old_genome_size;
}

/**
 * \brief    Set the new genome size (after mutations)
 * \details  --
 * \param    size_t new_genome_size
 * \return   \e void
 */
inline void ReplicationReport::set_new_genome_size( size_t new_genome_size )
{
  _new_genome_size = new_genome_size;
}

/**
 * \brief    Set the genome functional size
 * \details  --
 * \param    size_t functional_size
 * \return   \e void
 */
inline void ReplicationReport::set_genome_functional_size( size_t functional_size )
{
  _genome_functional_size = functional_size;
}

/**
 * \brief    Set the genome number of NC
 * \details  --
 * \param    size_t nb_NC
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_NC( size_t nb_NC )
{
  _genome_nb_NC = nb_NC;
}

/**
 * \brief    Set the genome number of E
 * \details  --
 * \param    size_t nb_E
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_E( size_t nb_E )
{
  _genome_nb_E = nb_E;
}

/**
 * \brief    Set the genome number of TF
 * \details  --
 * \param    size_t nb_TF
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_TF( size_t nb_TF )
{
  _genome_nb_TF = nb_TF;
}

/**
 * \brief    Set the genome number of BS
 * \details  --
 * \param    size_t nb_BS
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_BS( size_t nb_BS )
{
  _genome_nb_BS = nb_BS;
}

/**
 * \brief    Set the genome number of P
 * \details  --
 * \param    size_t nb_P
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_P( size_t nb_P )
{
  _genome_nb_P = nb_P;
}

/**
 * \brief    Set the genome number of inner enzymes
 * \details  --
 * \param    size_t nb_inner_enzymes
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_inner_enzymes( size_t nb_inner_enzymes )
{
  _genome_nb_inner_enzymes = nb_inner_enzymes;
}

/**
 * \brief    Set the genome number of inflowing pumps
 * \details  --
 * \param    size_t nb_inflow_pumps
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_inflow_pumps( size_t nb_inflow_pumps )
{
  _genome_nb_inflow_pumps = nb_inflow_pumps;
}

/**
 * \brief    Set the genome number of outflowing pumps
 * \details  --
 * \param    size_t nb_outflow_pumps
 * \return   \e void
 */
inline void ReplicationReport::set_genome_nb_outflow_pumps( size_t nb_outflow_pumps )
{
  _genome_nb_outflow_pumps = nb_outflow_pumps;
}

/*------------------------------------------------------------------ GRN data */

/**
 * \brief    Set the number of functional regions
 * \details  --
 * \param    size_t nb_functional_regions
 * \return   \e void
 */
inline void ReplicationReport::set_nb_functional_regions( size_t nb_functional_regions )
{
  _nb_functional_regions = nb_functional_regions;
}

/**
 * \brief    Set the number of enhancers in functional regions
 * \details  --
 * \param    size_t nb_enhancers
 * \return   \e void
 */
inline void ReplicationReport::set_nb_enhancers( size_t nb_enhancers )
{
  _nb_enhancers = nb_enhancers;
}

/**
 * \brief    Set the number of operators in functional regions
 * \details  --
 * \param    size_t nb_operators
 * \return   \e void
 */
inline void ReplicationReport::set_nb_operators( size_t nb_operators )
{
  _nb_operators = nb_operators;
}

/**
 * \brief    Set the number of functional regions containing only E types
 * \details  --
 * \param    size_t nb_E_regions
 * \return   \e void
 */
inline void ReplicationReport::set_nb_E_regions( size_t nb_E_regions )
{
  _nb_E_regions = nb_E_regions;
}

/**
 * \brief    Set the number of functional regions containing only TF types
 * \details  --
 * \param    size_t nb_TF_regions
 * \return   \e void
 */
inline void ReplicationReport::set_nb_TF_regions( size_t nb_TF_regions )
{
  _nb_TF_regions = nb_TF_regions;
}

/**
 * \brief    Set the number of functional regions mixing E and TF types
 * \details  --
 * \param    size_t nb_mixed_regions
 * \return   \e void
 */
inline void ReplicationReport::set_nb_mixed_regions( size_t nb_mixed_regions )
{
  _nb_mixed_regions = nb_mixed_regions;
}

/**
 * \brief    Set the mean size of functional regions
 * \details  --
 * \param    double mean_functional_region_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_functional_region_size( double mean_functional_region_size )
{
  assert(mean_functional_region_size >= 0.0);
  _mean_functional_region_size = mean_functional_region_size;
}

/**
 * \brief    Set the mean size of functional regions containing only E types
 * \details  --
 * \param    double mean_E_region_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_E_region_size( double mean_E_region_size )
{
  assert(mean_E_region_size >= 0.0);
  _mean_E_region_size = mean_E_region_size;
}

/**
 * \brief    Set the mean size of functional regions containing only TF types
 * \details  --
 * \param    double mean_TF_region_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_TF_region_size( double mean_TF_region_size )
{
  assert(mean_TF_region_size >= 0.0);
  _mean_TF_region_size = mean_TF_region_size;
}

/**
 * \brief    Set the mean size of functional regions mixing E and TF types
 * \details  --
 * \param    double mean_mixed_region_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_mixed_region_size( double mean_mixed_region_size )
{
  assert(mean_mixed_region_size >= 0.0);
  _mean_mixed_region_size = mean_mixed_region_size;
}

/**
 * \brief    Set the mean size of enhancers in functional regions
 * \details  --
 * \param    double mean_enhancer_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_enhancer_size( double mean_enhancer_size )
{
  assert(mean_enhancer_size >= 0.0);
  _mean_enhancer_size = mean_enhancer_size;
}

/**
 * \brief    Set the mean size of operators in functional regions
 * \details  --
 * \param    double mean_operator_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_operator_size( double mean_operator_size )
{
  assert(mean_operator_size >= 0.0);
  _mean_operator_size = mean_operator_size;
}

/**
 * \brief    Set the mean size of operons
 * \details  --
 * \param    double mean_operon_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_operon_size( double mean_operon_size )
{
  assert(mean_operon_size >= 0.0);
  _mean_operon_size = mean_operon_size;
}

/**
 * \brief    Set the mean size of operons containing only E types
 * \details  --
 * \param    double mean_E_operon_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_E_operon_size( double mean_E_operon_size )
{
  assert(mean_E_operon_size >= 0.0);
  _mean_E_operon_size = mean_E_operon_size;
}

/**
 * \brief    Set the mean size of operons containing only TF types
 * \details  --
 * \param    double mean_TF_operon_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_TF_operon_size( double mean_TF_operon_size )
{
  assert(mean_TF_operon_size >= 0.0);
  _mean_TF_operon_size = mean_TF_operon_size;
}

/**
 * \brief    Set the mean size of operons mixing E and TF types
 * \details  --
 * \param    double mean_mixed_operon_size
 * \return   \e void
 */
inline void ReplicationReport::set_mean_mixed_operon_size( double mean_mixed_operon_size )
{
  assert(mean_mixed_operon_size >= 0.0);
  _mean_mixed_operon_size = mean_mixed_operon_size;
}

/*------------------------------------------------------------------ Genetic redundancy */

/**
 * \brief    Set the mean regulation redundancy
 * \details  --
 * \param    double regulation_redundancy
 * \return   \e void
 */
inline void ReplicationReport::set_mean_regulation_redundancy( double regulation_redundancy )
{
  assert(regulation_redundancy >= 0.0);
  _mean_regulation_redundancy = regulation_redundancy;
}

/**
 * \brief    Set the mean metabolic redundancy
 * \details  --
 * \param    double metabolic_redundancy
 * \return   \e void
 */
inline void ReplicationReport::set_mean_metabolic_redundancy( double metabolic_redundancy )
{
  assert(metabolic_redundancy >= 0.0);
  _mean_metabolic_redundancy = metabolic_redundancy;
}

/*------------------------------------------------------------------ Inherited proteins structure */

/**
 * \brief    Set the inherited size
 * \details  --
 * \param    size_t size
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_size( size_t size )
{
  _inherited_size = size;
}

/**
 * \brief    Set the inherited number of E
 * \details  --
 * \param    size_t nb_E
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_nb_E( size_t nb_E )
{
  _inherited_nb_E = nb_E;
}

/**
 * \brief    Set the inherited number of TF
 * \details  --
 * \param    size_t nb_TF
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_nb_TF( size_t nb_TF )
{
  _inherited_nb_TF = nb_TF;
}

/**
 * \brief    Set the inherited number of inner enzymes
 * \details  --
 * \param    size_t nb_inner_enzymes
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_nb_inner_enzymes( size_t nb_inner_enzymes )
{
  _inherited_nb_inner_enzymes = nb_inner_enzymes;
}

/**
 * \brief    Set the inherited number of inflowing pumps
 * \details  --
 * \param    size_t nb_inflow_pumps
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_nb_inflow_pumps( size_t nb_inflow_pumps )
{
  _inherited_nb_inflow_pumps = nb_inflow_pumps;
}

/**
 * \brief    Set the inherited number of outflowing pumps
 * \details  --
 * \param    size_t nb_outflow_pumps
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_nb_outflow_pumps( size_t nb_outflow_pumps )
{
  _inherited_nb_outflow_pumps = nb_outflow_pumps;
}

/*------------------------------------------------------------------ Phenotype */

/**
 * \brief    Set cell id
 * \details  --
 * \param    unsigned long long int identifier
 * \return   \e void
 */
inline void ReplicationReport::set_id( unsigned long long int identifier )
{
  _id = identifier;
}

/**
 * \brief    Set parental cell id
 * \details  --
 * \param    unsigned long long int identifier
 * \return   \e void
 */
inline void ReplicationReport::set_parent_id( unsigned long long int identifier )
{
  _parent_id = identifier;
}

/**
 * \brief    Set generation
 * \details  --
 * \param    size_t generation
 * \return   \e void
 */
inline void ReplicationReport::set_generation( size_t generation )
{
  _generation = generation;
}

/**
 * \brief    Set x coordinate
 * \details  --
 * \param    size_t x
 * \return   \e void
 */
inline void ReplicationReport::set_x( size_t x )
{
  _x = x;
}

/**
 * \brief    Set y coordinate
 * \details  --
 * \param    size_t y
 * \return   \e void
 */
inline void ReplicationReport::set_y( size_t y )
{
  _y = y;
}

/**
 * \brief    Set the number of updates
 * \details  --
 * \param    size_t number_of_updates
 * \return   \e void
 */
inline void ReplicationReport::set_number_of_updates( size_t number_of_updates )
{
  _number_of_updates = number_of_updates;
}

/**
 * \brief    Set the number of divisions
 * \details  --
 * \param    size_t number_of_divisions
 * \return   \e void
 */
inline void ReplicationReport::set_number_of_divisions( size_t number_of_divisions )
{
  _number_of_divisions = number_of_divisions;
}

/**
 * \brief    Set cell birth time
 * \details  --
 * \param    size_t birth_time
 * \return   \e void
 */
inline void ReplicationReport::set_birth_time( size_t birth_time )
{
  _birth_time = birth_time;
}

/**
 * \brief    Set cell death time
 * \details  --
 * \param    size_t death_time
 * \return   \e void
 */
inline void ReplicationReport::set_death_time( size_t death_time )
{
  _death_time = death_time;
}

/**
 * \brief    Set lifespan
 * \details  --
 * \param    size_t lifespan
 * \return   \e void
 */
inline void ReplicationReport::set_lifespan( size_t lifespan )
{
  _lifespan = lifespan;
}

/**
 * \brief    Set the toxicity accumulation
 * \details  --
 * \param    double toxicity
 * \return   \e void
 */
inline void ReplicationReport::set_toxicity( double toxicity )
{
  assert(toxicity >= 0.0);
  _toxicity = toxicity;
}

/**
 * \brief    Set inherited TF amount
 * \details  --
 * \param    double inherited_TF_amount
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_TF_amount( double inherited_TF_amount )
{
  assert(inherited_TF_amount >= 0.0);
  _inherited_TF_amount = inherited_TF_amount;
}

/**
 * \brief    Set inherited E amount
 * \details  --
 * \param    double inherited_E_amount
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_E_amount( double inherited_E_amount )
{
  assert(inherited_E_amount >= 0.0);
  _inherited_E_amount = inherited_E_amount;
}

/**
 * \brief    Set TF amount
 * \details  --
 * \param    double TF_amount
 * \return   \e void
 */
inline void ReplicationReport::set_TF_amount( double TF_amount )
{
  assert(TF_amount >= 0.0);
  _TF_amount = TF_amount;
}

/**
 * \brief    Set E amount
 * \details  --
 * \param    double E_amount
 * \return   \e void
 */
inline void ReplicationReport::set_E_amount( double E_amount )
{
  assert(E_amount >= 0.0);
  _E_amount = E_amount;
}

/**
 * \brief    Set inherited metabolic amount
 * \details  --
 * \param    double inherited_metabolic_amount
 * \return   \e void
 */
inline void ReplicationReport::set_inherited_metabolic_amount( double inherited_metabolic_amount )
{
  assert(inherited_metabolic_amount >= 0.0);
  _inherited_metabolic_amount = inherited_metabolic_amount;
}

/**
 * \brief    Set the minimum metabolic amount
 * \details  --
 * \param    double min_metabolic_amount
 * \return   \e void
 */
inline void ReplicationReport::set_min_metabolic_amount( double min_metabolic_amount )
{
  assert(min_metabolic_amount >= 0.0);
  _min_metabolic_amount = min_metabolic_amount;
}

/**
 * \brief    Set the metabolic amount
 * \details  --
 * \param    double metabolic_amount
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_amount( double metabolic_amount )
{
  assert(metabolic_amount >= 0.0);
  _metabolic_amount = metabolic_amount;
}

/**
 * \brief    Set the maximum metabolic amount
 * \details  --
 * \param    double max_metabolic_amount
 * \return   \e void
 */
inline void ReplicationReport::set_max_metabolic_amount( double max_metabolic_amount )
{
  assert(max_metabolic_amount >= 0.0);
  _max_metabolic_amount = max_metabolic_amount;
}

/**
 * \brief    Set the metabolic uptake amount
 * \details  --
 * \param    double metabolic_uptake
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_uptake( double metabolic_uptake )
{
  assert(metabolic_uptake >= 0.0);
  _metabolic_uptake = metabolic_uptake;
}

/**
 * \brief    Set the metabolic release amount
 * \details  --
 * \param    double metabolic_uptake
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_release( double metabolic_release )
{
  assert(metabolic_release >= 0.0);
  _metabolic_release = metabolic_release;
}

/**
 * \brief    Set cell minimum energy
 * \details  --
 * \param    double min_energy
 * \return   \e void
 */
inline void ReplicationReport::set_min_energy( double min_energy )
{
  _min_energy = min_energy;
}

/**
 * \brief    Set cell mean energy
 * \details  --
 * \param    double mean_energy
 * \return   \e void
 */
inline void ReplicationReport::set_mean_energy( double mean_energy )
{
  _mean_energy = mean_energy;
}

/**
 * \brief    Set cell maximum energy
 * \details  --
 * \param    double max_energy
 * \return   \e void
 */
inline void ReplicationReport::set_max_energy( double max_energy )
{
  _max_energy = max_energy;
}

/**
 * \brief    Set cell minimum score
 * \details  --
 * \param    double min_score
 * \return   \e void
 */
inline void ReplicationReport::set_min_score( double min_score )
{
  assert(min_score >= 0.0);
  _min_score = min_score;
}

/**
 * \brief    Set cell mean score
 * \details  --
 * \param    double mean_score
 * \return   \e void
 */
inline void ReplicationReport::set_mean_score( double mean_score )
{
  assert(mean_score >= 0.0);
  _mean_score = mean_score;
}

/**
 * \brief    Set cell maximum score
 * \details  --
 * \param    double max_score
 * \return   \e void
 */
inline void ReplicationReport::set_max_score( double max_score )
{
  assert(max_score >= 0.0);
  _max_score = max_score;
}

/**
 * \brief    Set the metabolic growth rate
 * \details  --
 * \param    double metabolic_growth_rate
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_growth_rate( double metabolic_growth_rate )
{
  _metabolic_growth_rate = metabolic_growth_rate;
}

/**
 * \brief    Set the metabolic growth rate difference
 * \details  --
 * \param    double diff
 * \return   \e void
 */
inline void ReplicationReport::set_Dmetabolic_growth_rate( double diff )
{
  _diff_metabolic_growth_rate = diff;
}

/**
 * \brief    Set the number of nodes in the GRN
 * \details  --
 * \param    size_t nb_nodes
 * \return   \e void
 */
inline void ReplicationReport::set_grn_nb_nodes( size_t nb_nodes )
{
  _grn_nb_nodes = nb_nodes;
}

/**
 * \brief    Set the number of edges in the GRN
 * \details  --
 * \param    size_t nb_edges
 * \return   \e void
 */
inline void ReplicationReport::set_grn_nb_edges( size_t nb_edges )
{
  _grn_nb_edges = nb_edges;
}

/**
 * \brief    Set the number of nodes in the metabolic network
 * \details  --
 * \param    size_t nb_nodes
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_nb_nodes( size_t nb_nodes )
{
  _metabolic_nb_nodes = nb_nodes;
}

/**
 * \brief    Set the number of edges in the metabolic network
 * \details  --
 * \param    size_t nb_edges
 * \return   \e void
 */
inline void ReplicationReport::set_metabolic_nb_edges( size_t nb_edges )
{
  _metabolic_nb_edges = nb_edges;
}

/**
 * \brief    Set the trophic group index in the trophic network
 * \details  --
 * \param    unsigned long long int trophic_group
 * \return   \e void
 */
inline void ReplicationReport::set_trophic_group( unsigned long long int trophic_group )
{
  _trophic_group = trophic_group;
}

/**
 * \brief    Set the trophic level
 * \details  --
 * \param    trophic_level level
 * \return   \e void
 */
inline void ReplicationReport::set_trophic_level( trophic_level level )
{
  _trophic_level = level;
}


#endif /* defined(__Evo2Sim__ReplicationReport__) */
