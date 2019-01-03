
/**
 * \file      TrophicGroup.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      20-10-2015
 * \copyright Copyright (C) 2014-2019 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicGroup class declaration
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

#ifndef __EVOEVO__TrophicGroup__
#define __EVOEVO__TrophicGroup__

#include <iostream>
#include <unordered_map>
#include <vector>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "Cell.h"


class TrophicGroup
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  TrophicGroup( void ) = delete;
  TrophicGroup( unsigned long long int identifier, size_t time, double red, double green, double blue );
  TrophicGroup( gzFile backup_file );
  TrophicGroup( const TrophicGroup& trophic_group );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~TrophicGroup( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ trophic group properties */
  
  inline unsigned long long int               get_identifier( void ) const;
  inline std::vector<unsigned long long int>* get_necrophagy_links( void );
  inline std::vector<unsigned long long int>* get_active_release_links( void );
  inline size_t                               get_appearance_time( void ) const;
  inline size_t                               get_lifespan( void ) const;
  inline size_t                               get_number_of_cells( void ) const;
  inline trophic_level                        get_trophic_level( void ) const;
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  inline double get_mean_generations( void ) const;
  inline double get_mean_inherited_TF_amount( void ) const;
  inline double get_mean_inherited_E_amount( void ) const;
  inline double get_mean_TF_amount( void ) const;
  inline double get_mean_E_amount( void ) const;
  inline double get_mean_inherited_metabolic_amount( void ) const;
  inline double get_mean_metabolic_amount( void ) const;
  inline double get_mean_energy( void ) const;
  inline double get_mean_score( void ) const;
  inline double get_mean_lifespan( void ) const;
  inline double get_mean_number_of_divisions( void ) const;
  inline double get_mean_toxicity( void ) const;
  inline double get_mean_metabolic_uptake( void ) const;
  inline double get_mean_metabolic_release( void ) const;
  inline double get_mean_metabolic_growth_rate( void ) const;
  inline double get_mean_Dmetabolic_growth_rate( void ) const;
  inline double get_mean_grn_nb_nodes( void ) const;
  inline double get_mean_grn_nb_edges( void ) const;
  inline double get_mean_metabolic_nb_nodes( void ) const;
  inline double get_mean_metabolic_nb_edges( void ) const;
  inline double get_mean_regulation_redundancy( void ) const;
  inline double get_mean_metabolic_redundancy( void ) const;
  
  /* GENOME STRUCTURE */
  inline double get_mean_genome_size( void ) const;
  inline double get_mean_functional_size( void ) const;
  inline double get_mean_genome_nb_NC( void ) const;
  inline double get_mean_genome_nb_E( void ) const;
  inline double get_mean_genome_nb_TF( void ) const;
  inline double get_mean_genome_nb_BS( void ) const;
  inline double get_mean_genome_nb_P( void ) const;
  inline double get_mean_genome_nb_inner_enzymes( void ) const;
  inline double get_mean_genome_nb_inflow_pumps( void ) const;
  inline double get_mean_genome_nb_outflow_pumps( void ) const;
  inline double get_mean_genome_nb_functional_regions( void ) const;
  inline double get_mean_genome_nb_enhancers( void ) const;
  inline double get_mean_genome_nb_operators( void ) const;
  inline double get_mean_genome_nb_E_regions( void ) const;
  inline double get_mean_genome_nb_TF_regions( void ) const;
  inline double get_mean_genome_nb_mixed_regions( void ) const;
  inline double get_mean_genome_functional_region_size( void ) const;
  inline double get_mean_genome_E_region_size( void ) const;
  inline double get_mean_genome_TF_region_size( void ) const;
  inline double get_mean_genome_mixed_region_size( void ) const;
  inline double get_mean_genome_enhancer_size( void ) const;
  inline double get_mean_genome_operator_size( void ) const;
  inline double get_mean_genome_operon_size( void ) const;
  inline double get_mean_genome_E_operon_size( void ) const;
  inline double get_mean_genome_TF_operon_size( void ) const;
  inline double get_mean_genome_mixed_operon_size( void ) const;
  
  /* INHERITED STRUCTURE */
  inline double get_mean_inherited_size( void ) const;
  inline double get_mean_inherited_nb_E( void ) const;
  inline double get_mean_inherited_nb_TF( void ) const;
  inline double get_mean_inherited_nb_inner_enzymes( void ) const;
  inline double get_mean_inherited_nb_inflow_pumps( void ) const;
  inline double get_mean_inherited_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ trophic group color */
  
  inline double get_red_color( void ) const;
  inline double get_green_color( void ) const;
  inline double get_blue_color( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ trophic group properties */
  
  inline void set_identifier( unsigned long long int identifier );
  inline void set_necrophagy_links( std::vector<unsigned long long int>* links );
  inline void set_active_release_links( std::vector<unsigned long long int>* links );
  inline void set_appearance_time( size_t appearance_time );
  inline void set_lifespan( size_t lifespan );
  inline void set_number_of_cells( size_t number_of_cells );
  inline void set_trophic_level( trophic_level level );
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  inline void update_lifespan( size_t time );
  inline void add_necrophagy_link( unsigned long long int group_id );
  inline void add_active_release_link( unsigned long long int group_id );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  void add_alive_cell( Parameters* parameters, Cell* cell );
  void compute_mean( void );
  void clear( void );
  void update_profile( size_t N );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ trophic profiles */
  
  std::string _trophic_profile;    /*!< Trophic profile       */
  std::string _production_profile; /*!< Production profile    */
  std::string _uptake_profile;     /*!< Uptake pumps profile  */
  std::string _release_profile;    /*!< Release pumps profile */
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ trophic group properties */
  
  unsigned long long int              _identifier;              /*!< Unique identifier               */
  std::vector<unsigned long long int> _necrophagy_links;        /*!< Necrophagy links                */
  std::vector<unsigned long long int> _active_release_links;    /*!< Active release links            */
  size_t                              _appearance_time;         /*!< Trophic group appearance time   */
  size_t                              _lifespan;                /*!< Trophic group lifespan          */
  size_t                              _number_of_cells;         /*!< Number of cells                 */
  trophic_level                       _trophic_level;           /*!< Trophic level (none, 0, 1 or 2) */
  
  /*------------------------------------------------------------------ trophic group statistics */
  
  /* PHENOTYPE */
  double _mean_generations;                /*!< Number of generations                    */
  double _mean_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _mean_inherited_E_amount;         /*!< Inherited E amount                       */
  double _mean_TF_amount;                  /*!< TF amount                                */
  double _mean_E_amount;                   /*!< E amount                                 */
  double _mean_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _mean_metabolic_amount;           /*!< Metabolic amount                         */
  double _mean_energy;                     /*!< Energy amount                            */
  double _mean_score;                      /*!< Score                                    */
  double _mean_lifespan;                   /*!< Lifespan                                 */
  double _mean_number_of_divisions;        /*!< Number of divisions                      */
  double _mean_toxicity;                   /*!< Toxicity accumulation                    */
  double _mean_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _mean_metabolic_release;          /*!< Amount of metabolic release              */
  double _mean_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _mean_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _mean_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _mean_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _mean_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _mean_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _mean_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _mean_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  
  /* GENOME STRUCTURE */
  double _mean_genome_size;                   /*!< Genome size                                                   */
  double _mean_functional_size;               /*!< Functional size                                               */
  double _mean_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _mean_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _mean_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _mean_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _mean_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _mean_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _mean_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _mean_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _mean_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _mean_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _mean_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _mean_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _mean_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _mean_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _mean_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _mean_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _mean_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _mean_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _mean_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _mean_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _mean_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _mean_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _mean_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _mean_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _mean_inherited_size;             /*!< Number of inherited proteins                                   */
  double _mean_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _mean_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _mean_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _mean_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _mean_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ trophic group color */
  
  double _red_color;   /*!< RGB red   */
  double _green_color; /*!< RGB green */
  double _blue_color;  /*!< RGB blue  */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ trophic group properties */

/**
 * \brief    Get the group identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int TrophicGroup::get_identifier( void ) const
{
  return _identifier;
}

/**
 * \brief    Get necrophagy links
 * \details  --
 * \param    void
 * \return   \e std::vector<unsigned long long int>*
 */
inline std::vector<unsigned long long int>* TrophicGroup::get_necrophagy_links( void )
{
  return &_necrophagy_links;
}

/**
 * \brief    Get active release links
 * \details  --
 * \param    void
 * \return   \e std::vector<unsigned long long int>*
 */
inline std::vector<unsigned long long int>* TrophicGroup::get_active_release_links( void )
{
  return &_active_release_links;
}

/**
 * \brief    Get the appearance time
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicGroup::get_appearance_time( void ) const
{
  return _appearance_time;
}

/**
 * \brief    Get the group lifespan
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicGroup::get_lifespan( void ) const
{
  return _lifespan;
}

/**
 * \brief    Get the number of cells belonging to this group
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicGroup::get_number_of_cells( void ) const
{
  return _number_of_cells;
}

/**
 * \brief    Get the trophic level of this cell
 * \details  --
 * \param    void
 * \return   \e trophic_level
 */
inline trophic_level TrophicGroup::get_trophic_level( void ) const
{
  return _trophic_level;
}

/*------------------------------------------------------------------ trophic group statistics */

/**
 * \brief    Get the mean generation
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_generations( void ) const
{
  return _mean_generations;
}

/**
 * \brief    Get the mean inherited TF amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_TF_amount( void ) const
{
  return _mean_inherited_TF_amount;
}

/**
 * \brief    Get the mean inherited E amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_E_amount( void ) const
{
  return _mean_inherited_E_amount;
}

/**
 * \brief    Get the mean TF amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_TF_amount( void ) const
{
  return _mean_TF_amount;
}

/**
 * \brief    Get the mean E amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_E_amount( void ) const
{
  return _mean_E_amount;
}

/**
 * \brief    Get the mean inherited metabolic amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_metabolic_amount( void ) const
{
  return _mean_inherited_metabolic_amount;
}

/**
 * \brief    Get the mean metabolic amount in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_amount( void ) const
{
  return _mean_metabolic_amount;
}

/**
 * \brief    Get the mean energy in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_energy( void ) const
{
  return _mean_energy;
}

/**
 * \brief    Get the mean score in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_score( void ) const
{
  return _mean_score;
}

/**
 * \brief    Get the mean lifespan in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_lifespan( void ) const
{
  return _mean_lifespan;
}

/**
 * \brief    Get the mean number of division in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_number_of_divisions( void ) const
{
  return _mean_number_of_divisions;
}

/**
 * \brief    Get the mean toxicity in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_toxicity( void ) const
{
  return _mean_toxicity;
}

/**
 * \brief    Get the mean metabolci uptake in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_uptake( void ) const
{
  return _mean_metabolic_uptake;
}

/**
 * \brief    Get the mean metabolic release in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_release( void ) const
{
  return _mean_metabolic_release;
}

/**
 * \brief    Get the mean metabolic growth rate in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_growth_rate( void ) const
{
  return _mean_metabolic_growth_rate;
}

/**
 * \brief    Get the mean metabolic growth rate differential in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_Dmetabolic_growth_rate( void ) const
{
  return _mean_Dmetabolic_growth_rate;
}

/**
 * \brief    Get the mean number of GRN nodes in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_grn_nb_nodes( void ) const
{
  return _mean_grn_nb_nodes;
}

/**
 * \brief    Get the mean number of GRN edges in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_grn_nb_edges( void ) const
{
  return _mean_grn_nb_edges;
}

/**
 * \brief    Get the mean number of metabolic network nodes in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_nb_nodes( void ) const
{
  return _mean_metabolic_nb_nodes;
}

/**
 * \brief    Get the mean number of metabolic network edges in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_nb_edges( void ) const
{
  return _mean_metabolic_nb_edges;
}

/**
 * \brief    Get the mean regulation redundancy in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_regulation_redundancy( void ) const
{
  return _mean_regulation_redundancy;
}

/**
 * \brief    Get the mean metabolic redundancy in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_metabolic_redundancy( void ) const
{
  return _mean_metabolic_redundancy;
}

/**
 * \brief    Get the mean genome size in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_size( void ) const
{
  return _mean_genome_size;
}

/**
 * \brief    Get the mean genome functional size in the group
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_functional_size( void ) const
{
  return _mean_functional_size;
}

/**
 * \brief    Get the mean number of NC units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_NC( void ) const
{
  return _mean_genome_nb_NC;
}

/**
 * \brief    Get the mean number of E units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_E( void ) const
{
  return _mean_genome_nb_E;
}

/**
 * \brief    Get the mean number of TF units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_TF( void ) const
{
  return _mean_genome_nb_TF;
}

/**
 * \brief    Get the mean number of BS units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_BS( void ) const
{
  return _mean_genome_nb_BS;
}

/**
 * \brief    Get the mean number of P units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_P( void ) const
{
  return _mean_genome_nb_P;
}

/**
 * \brief    Get the mean number of inner enzymes
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_inner_enzymes( void ) const
{
  return _mean_genome_nb_inner_enzymes;
}

/**
 * \brief    Get the mean number of inflowing pumps
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_inflow_pumps( void ) const
{
  return _mean_genome_nb_inflow_pumps;
}

/**
 * \brief    Get the mean number of outflowing pumps
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_outflow_pumps( void ) const
{
  return _mean_genome_nb_outflow_pumps;
}

/**
 * \brief    Get the mean number of functional regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_functional_regions( void ) const
{
  return _mean_genome_nb_functional_regions;
}

/**
 * \brief    Get the mean number of enhancers
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_enhancers( void ) const
{
  return _mean_genome_nb_enhancers;
}

/**
 * \brief    Get the mean number of operators
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_operators( void ) const
{
  return _mean_genome_nb_operators;
}

/**
 * \brief    Get the mean number of E regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_E_regions( void ) const
{
  return _mean_genome_nb_E_regions;
}

/**
 * \brief    Get the mean number of TF regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_TF_regions( void ) const
{
  return _mean_genome_nb_TF_regions;
}

/**
 * \brief    Get the mean number of mixed regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_nb_mixed_regions( void ) const
{
  return _mean_genome_nb_mixed_regions;
}

/**
 * \brief    Get the mean size of functional regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_functional_region_size( void ) const
{
  return _mean_genome_functional_region_size;
}

/**
 * \brief    Get the mean size of E regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_E_region_size( void ) const
{
  return _mean_genome_E_region_size;
}

/**
 * \brief    Get the mean size of TF regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_TF_region_size( void ) const
{
  return _mean_genome_TF_region_size;
}

/**
 * \brief    Get the mean size of mixed regions
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_mixed_region_size( void ) const
{
  return _mean_genome_mixed_region_size;
}

/**
 * \brief    Get the mean size of enhancers
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_enhancer_size( void ) const
{
  return _mean_genome_enhancer_size;
}

/**
 * \brief    Get the mean size of operators
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_operator_size( void ) const
{
  return _mean_genome_operator_size;
}

/**
 * \brief    Get the mean size of operons
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_operon_size( void ) const
{
  return _mean_genome_operon_size;
}

/**
 * \brief    Get the mean size of E operons
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_E_operon_size( void ) const
{
  return _mean_genome_E_operon_size;
}

/**
 * \brief    Get the mean size of TF operons
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_TF_operon_size( void ) const
{
  return _mean_genome_TF_operon_size;
}

/**
 * \brief    Get the mean size of mixed operons
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_genome_mixed_operon_size( void ) const
{
  return _mean_genome_mixed_operon_size;
}

/**
 * \brief    Get the mean inherited size
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_size( void ) const
{
  return _mean_inherited_size;
}

/**
 * \brief    Get the mean number of inherited E units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_nb_E( void ) const
{
  return _mean_inherited_nb_E;
}

/**
 * \brief    Get the mean number of inherited TF units
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_nb_TF( void ) const
{
  return _mean_inherited_nb_TF;
}

/**
 * \brief    Get the mean number of inherited inner enzymes
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_nb_inner_enzymes( void ) const
{
  return _mean_inherited_nb_inner_enzymes;
}

/**
 * \brief    Get the mean number of inherited inflowing pumps
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_nb_inflow_pumps( void ) const
{
  return _mean_inherited_nb_inflow_pumps;
}

/**
 * \brief    Get the mean number of inherited outflowing pumps
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_mean_inherited_nb_outflow_pumps( void ) const
{
  return _mean_inherited_nb_outflow_pumps;
}

/*------------------------------------------------------------------ trophic group color */

/**
 * \brief    Get RGB red color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_red_color( void ) const
{
  return _red_color;
}

/**
 * \brief    Get RGB green color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_green_color( void ) const
{
  return _green_color;
}

/**
 * \brief    Get RGB blue color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double TrophicGroup::get_blue_color( void ) const
{
  return _blue_color;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*------------------------------------------------------------------ trophic group properties */

/**
 * \brief    Set the group identifier
 * \details  --
 * \param    unsigned long long int identifier
 * \return   \e void
 */
inline void TrophicGroup::set_identifier( unsigned long long int identifier )
{
  _identifier = identifier;
}

/**
 * \brief    Set necrophagy links
 * \details  --
 * \param    std::vector<unsigned long long int>* links
 * \return   \e void
 */
inline void TrophicGroup::set_necrophagy_links( std::vector<unsigned long long int>* links )
{
  _necrophagy_links.clear();
  for (size_t i = 0; i < links->size(); i++)
  {
    _necrophagy_links.push_back(links->at(i));
  }
}

/**
 * \brief    Set active release links
 * \details  --
 * \param    std::vector<unsigned long long int>* links
 * \return   \e void
 */
inline void TrophicGroup::set_active_release_links( std::vector<unsigned long long int>* links )
{
  _active_release_links.clear();
  for (size_t i = 0; i < links->size(); i++)
  {
    _active_release_links.push_back(links->at(i));
  }
}

/**
 * \brief    Set the appearance time
 * \details  --
 * \param    size_t appearance_time
 * \return   \e void
 */
inline void TrophicGroup::set_appearance_time( size_t appearance_time )
{
  _appearance_time = appearance_time;
}

/**
 * \brief    Set the lifespan
 * \details  --
 * \param    size_t lifespan
 * \return   \e void
 */
inline void TrophicGroup::set_lifespan( size_t lifespan )
{
  _lifespan = lifespan;
}

/**
 * \brief    Set the number of cells belonging to this group
 * \details  --
 * \param    size_t number_of_cells
 * \return   \e void
 */
inline void TrophicGroup::set_number_of_cells( size_t number_of_cells )
{
  _number_of_cells = number_of_cells;
}

/**
 * \brief    Set the trophic level of this cell
 * \details  --
 * \param    trophic_level level
 * \return   \e void
 */
inline void TrophicGroup::set_trophic_level( trophic_level level )
{
  _trophic_level = level;
}

/**
 * \brief    Update the lifespan
 * \details  --
 * \param    size_t time
 * \return   \e void
 */
inline void TrophicGroup::update_lifespan( size_t time )
{
  _lifespan = time-_appearance_time;
}

/**
 * \brief    Add a necrophagy link
 * \details  --
 * \param    unsigned long long int group_id
 * \return   \e void
 */
inline void TrophicGroup::add_necrophagy_link( unsigned long long int group_id )
{
  _necrophagy_links.push_back(group_id);
}

/**
 * \brief    Add an active release link
 * \details  --
 * \param    unsigned long long int group_id
 * \return   \e void
 */
inline void TrophicGroup::add_active_release_link( unsigned long long int group_id )
{
  _active_release_links.push_back(group_id);
}


#endif /* defined(__EVOEVO__TrophicGroup__) */
