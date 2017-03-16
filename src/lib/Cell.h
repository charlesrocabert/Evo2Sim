
/**
 * \file      Cell.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Cell class declaration
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

#ifndef __EVOEVO__Cell__
#define __EVOEVO__Cell__

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Genome.h"
#include "InheritedProteins.h"
#include "ODE.h"
#include "SpeciesList.h"
#include "Environment.h"
#include "Prng.h"


class Cell
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Cell( void ) = delete;
  Cell( Parameters* parameters, Prng* prng );
  Cell( Parameters* parameters, Prng* prng, gzFile backup_file );
  Cell( const Cell& cell ) = delete;
  Cell( Cell& parent, Prng* prng, unsigned long long int new_id, size_t x, size_t y, size_t time );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Cell( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ main cell classes */
  
  inline ReplicationReport* get_replication_report( void );
  inline Genome*            get_genome( void );
  inline InheritedProteins* get_inherited_proteins( void );
  inline SpeciesList*       get_inherited_species_list( void );
  inline SpeciesList*       get_species_list( void );
  
  /*------------------------------------------------------------------ ODE system */
  
  inline ODE* get_ode( void );
  
  /*------------------------------------------------------------------ main cell variables */
  
  inline unsigned long long int get_id( void ) const;
  inline unsigned long long int get_parent_id( void ) const;
  inline size_t                 get_generation( void ) const;
  inline double                 get_amount( void );
  inline double                 get_energy( void ) const;
  inline bool                   get_active( void ) const;
  inline bool                   get_alive( void ) const;
  inline size_t                 get_x( void ) const;
  inline size_t                 get_y( void ) const;
  inline double                 get_score( void ) const;
  inline size_t                 get_number_of_updates( void ) const;
  inline size_t                 get_number_of_divisions( void ) const;
  inline bool                   isActive( void ) const;
  inline bool                   isAlive( void ) const;
  inline bool                   isTagged( void ) const;
  
  /*------------------------------------------------------------------ mutation rates */
  
  inline double get_point_mutation_rate( void ) const;
  inline double get_duplication_rate( void ) const;
  inline double get_deletion_rate( void ) const;
  inline double get_translocation_rate( void ) const;
  inline double get_inversion_rate( void ) const;
  inline double get_transition_rate( void ) const;
  inline double get_breakpoint_rate( void ) const;
  inline double get_substrate_tag_mutation_size( void ) const;
  inline double get_product_tag_mutation_size( void ) const;
  inline double get_kcat_mutation_size( void ) const;
  inline double get_kcat_km_ratio_mutation_size( void ) const;
  inline double get_binding_site_tag_mutation_size( void ) const;
  inline double get_co_enzyme_tag_mutation_size( void ) const;
  inline double get_transcription_factor_tag_mutation_size( void ) const;
  inline double get_basal_expression_level_mutation_size( void ) const;
  
  /*------------------------------------------------------------------ time variables */
  
  inline size_t get_birth_time( void ) const;
  inline size_t get_death_time( void ) const;
  inline size_t get_lifespan( void ) const;
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  inline double            get_toxicity( void ) const;
  inline double            get_inherited_TF_amount( void ) const;
  inline double            get_inherited_E_amount( void ) const;
  inline double            get_TF_amount( void ) const;
  inline double            get_E_amount( void ) const;
  inline double            get_inherited_metabolic_amount( void ) const;
  inline double            get_min_metabolic_amount( void ) const;
  inline double            get_max_metabolic_amount( void ) const;
  inline double            get_metabolic_uptake( void ) const;
  inline double            get_metabolic_release( void ) const;
  inline double            get_min_energy( void ) const;
  inline double            get_mean_energy( void ) const;
  inline double            get_max_energy( void ) const;
  inline double            get_min_score( void ) const;
  inline double            get_mean_score( void ) const;
  inline double            get_max_score( void ) const;
  inline double            get_metabolic_growth_rate( void ) const;
  inline double            get_Dmetabolic_growth_rate( void ) const;
  inline size_t            get_grn_nb_nodes( void ) const;
  inline size_t            get_grn_nb_edges( void ) const;
  inline size_t            get_metabolic_nb_nodes( void ) const;
  inline size_t            get_metabolic_nb_edges( void ) const;
  inline std::vector<int>* get_inflowing_pumps( void );
  inline std::vector<int>* get_outflowing_pumps( void );
  
  inline unsigned long long int get_trophic_group( void ) const;
  inline trophic_level          get_trophic_level( void ) const;
  
  /*------------------------------------------------------------------ cell color */
  
  inline double get_red_color( void ) const;
  inline double get_green_color( void ) const;
  inline double get_blue_color( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  Cell& operator=(const Cell&) = delete;
  
  /*------------------------------------------------------------------ main cell variables */
  
  inline void set_id( unsigned long long int identifier );
  inline void set_parent_id( unsigned long long int identifier );
  inline void set_generation( size_t generation );
  inline void set_energy( double energy );
  inline void set_active( bool active );
  inline void set_alive( bool alive );
  inline void set_x( size_t x );
  inline void set_y( size_t y );
  inline void set_score( double score );
  inline void set_number_of_updates( size_t number_of_updates );
  inline void set_number_of_divisions( size_t number_of_divisions );
  inline void activate( void );
  inline void inactivate( void );
  inline void tag( void );
  inline void untag( void );
  inline void update_number_of_updates( void );
  inline void update_number_of_divisions( void );
  
  /*------------------------------------------------------------------ mutation rates */
  
  inline void set_point_mutation_rate( double point_mutation_rate );
  inline void set_duplication_rate( double duplication_rate );
  inline void set_deletion_rate( double deletion_rate );
  inline void set_translocation_rate( double translocation_rate );
  inline void set_inversion_rate( double inversion_rate );
  inline void set_transition_rate( double transition_rate );
  inline void set_breakpoint_rate( double breakpoint_rate );
  inline void set_substrate_tag_mutation_size( double size );
  inline void set_product_tag_mutation_size( double size );
  inline void set_kcat_mutation_size( double size );
  inline void set_kcat_km_ratio_mutation_size( double size );
  inline void set_binding_site_tag_mutation_size( double size );
  inline void set_co_enzyme_tag_mutation_size( double size );
  inline void set_transcription_factor_tag_mutation_size( double size );
  inline void set_basal_expression_level_mutation_size( double size );
  
  /*------------------------------------------------------------------ time variables */
  
  inline void set_birth_time( size_t birth_time );
  inline void set_death_time( size_t death_time );
  inline void set_lifespan( size_t lifespan );
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  inline void set_toxicity( double toxicity );
  inline void set_min_metabolic_amount( double min_metabolic_amount );
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
  inline void set_Dmetabolic_growth_rate( double diff );
  inline void set_grn_nb_nodes( size_t nb_nodes );
  inline void set_grn_nb_edges( size_t nb_edges );
  inline void set_metabolic_nb_nodes( size_t nb_nodes );
  inline void set_metabolic_nb_edges( size_t nb_edges );
  inline void set_trophic_group( unsigned long long int trophic_group );
  inline void set_trophic_level( trophic_level level );
  
  /*------------------------------------------------------------------ cell color */
  
  inline void set_red_color( double red );
  inline void set_green_color( double green );
  inline void set_blue_color( double blue );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void initialize_inherited_species_list( void );
  void load_genome_in_inherited_proteins( void );
  void load_genome_in_species_lists( void );
  void load_genome_in_ODE_system( Environment* environment, bool from_backup, bool new_individual );
  void synchronize_state_vectors( Environment* environment );
  void update_replication_report_data( void );
  void update( size_t time );
  void compute_score( void );
  void mutate( void );
  void kill( size_t death_time );
  void save( gzFile backup_file );
  void replace_data( Cell* cell );
  
  void write_genome( std::ofstream& filestream );
  void write_inherited_proteins( std::ofstream& filestream );
  void write_genetic_regulation_network( std::ofstream& nodeslist, std::ofstream& edgeslist );
  void write_metabolic_network( std::ofstream& nodeslist, std::ofstream& edgeslist );
  void write_metabolic_amounts( std::ofstream& filestream );
  //void write_current_ODE_system( void );
  void write_state_header( std::ofstream& filestream );
  void write_current_state( std::ofstream& filestream, size_t time, Environment* environment );
  
  void get_pumped_and_produced_metabolites( std::vector<int>& pumped, std::vector<int>& produced, std::vector<int>& pm_pumped, std::vector<int>& pm_produced );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void initialize_mutation_rates( void );
  void save_mutation_rates( gzFile backup_file );
  void load_mutation_rates( gzFile backup_file );
  void copy_mutation_rates( double* mutation_rates );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ simulation parameters */
  
  Parameters* _parameters; /*!< Parameters */
  
  /*------------------------------------------------------------------ prng */
  
  Prng* _prng; /*!< Pseudorandom numbers generator */
  
  /*------------------------------------------------------------------ main cell classes */
  
  ReplicationReport* _replication_report;     /*!< Replication report                  */
  Genome*            _genome;                 /*!< Genome                              */
  InheritedProteins* _inherited_proteins;     /*!< Inherited proteins                  */
  SpeciesList*       _inherited_species_list; /*!< List of inherited metabolic species */
  SpeciesList*       _species_list;           /*!< List of metabolic species           */
  
  /*------------------------------------------------------------------ ODE system */
  
  ODE* _ode; /*!< ODE solver system */
  
  /*------------------------------------------------------------------ main cell variables */
  
  unsigned long long int _id;                  /*!< Identifier                      */
  unsigned long long int _parent_id;           /*!< Parental identifier             */
  size_t                 _generation;          /*!< Generation of the cell          */
  double                 _energy;              /*!< Cell's energy bilan             */
  bool                   _active;              /*!< Defines if the cell is active   */
  bool                   _alive;               /*!< Defines if the cell is alive    */
  size_t                 _x;                   /*!< Cell's x coordinate             */
  size_t                 _y;                   /*!< Cell's y coordinate             */
  double                 _score;               /*!< Cell's score                    */
  size_t                 _number_of_updates;   /*!< Number of updates since birth   */
  size_t                 _number_of_divisions; /*!< Number of divisions since birth */
  bool                   _tagged;              /*!< Cell's tag                      */
  
  /*------------------------------------------------------------------ mutation rates */
  
  double* _mutation_rates; /*!< Vector of mutation rates */
  
  /*------------------------------------------------------------------ time variables */
  
  size_t _birth_time; /*!< Cell's birth time */
  size_t _death_time; /*!< Cell's death time */
  size_t _lifespan;   /*!< Cell's lifespan   */
  
  /*------------------------------------------------------------------ global phenotypic variables */
  
  double _toxicity;                   /*!< Toxicity accumulation                         */
  double _inherited_TF_amount;        /*!< Inherited transcription factors amount        */
  double _inherited_E_amount;         /*!< Inherited enzymes amount                      */
  double _TF_amount;                  /*!< Transcription factors amount                  */
  double _E_amount;                   /*!< Enzymes amount                                */
  double _min_metabolic_amount;       /*!< Minimum cytoplamsic metabolic amount          */
  double _max_metabolic_amount;       /*!< Maximum cytoplamsic metabolic amount          */
  double _metabolic_uptake;           /*!< Amount of metabolic uptake                    */
  double _metabolic_release;          /*!< Amount of released metabolites                */
  double _previous_metabolic_amount;  /*!< Previous metabolic amount                     */
  double _min_energy;                 /*!< Minimum energy bilan of the cell through time */
  double _mean_energy;                /*!< Mean energy bilan of the cell through time    */
  double _max_energy;                 /*!< Maximum energy bilan of the cell through time */
  double _min_score;                  /*!< Minimum score of the cell through time        */
  double _mean_score;                 /*!< Mean score of the cell through time           */
  double _max_score;                  /*!< Maximum score of the cell through time        */
  double _metabolic_growth_rate;      /*!< Cell's metabolic growth rate                  */
  double _diff_metabolic_growth_rate; /*!< Cell's metabolic growth rate difference       */
  size_t _grn_nb_nodes;               /*!< Number of nodes in the GRN                    */
  size_t _grn_nb_edges;               /*!< Number of edges in the GRN                    */
  size_t _metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network      */
  size_t _metabolic_nb_edges;         /*!< Number of edges in the metabolic network      */
  std::vector<int> _inflowing_pumps;  /*!< Types of inflowing pumps                      */
  std::vector<int> _outflowing_pumps; /*!< Types of outflowing pumps                     */
  
  unsigned long long int _trophic_group; /*!< Trophic group index in the trophic network */
  trophic_level          _trophic_level; /*!< Trophic level                              */
  
  /*------------------------------------------------------------------ cell color */
  
  double _red_color;   /*!< RGB red   */
  double _green_color; /*!< RGB green */
  double _blue_color;  /*!< RGB blue  */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ main cell classes */

/**
 * \brief    Get replication report
 * \details  --
 * \param    void
 * \return   \e ReplicationReport*
 */
inline ReplicationReport* Cell::get_replication_report( void )
{
  return _replication_report;
}

/**
 * \brief    Get genome
 * \details  --
 * \param    void
 * \return   \e Genome*
 */
inline Genome* Cell::get_genome( void )
{
  return _genome;
}

/**
 * \brief    Get inherited proteins
 * \details  --
 * \param    void
 * \return   \e InheritedProteins*
 */
inline InheritedProteins* Cell::get_inherited_proteins( void )
{
  return _inherited_proteins;
}

/**
 * \brief    Get inherited species list
 * \details  --
 * \param    void
 * \return   \e SpeciesList*
 */
inline SpeciesList* Cell::get_inherited_species_list( void )
{
  return _species_list;
}

/**
 * \brief    Get species list
 * \details  --
 * \param    void
 * \return   \e SpeciesList*
 */
inline SpeciesList* Cell::get_species_list( void )
{
  return _species_list;
}

/*------------------------------------------------------------------ ODE system */

/**
 * \brief    Get ODE system
 * \details  --
 * \param    void
 * \return   \e ODE*
 */
inline ODE* Cell::get_ode( void )
{
  return _ode;
}

/*------------------------------------------------------------------ main cell variables */

/**
 * \brief    Get cell id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Cell::get_id( void ) const
{
  return _id;
}

/*
 * \brief    Get parental cell id
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Cell::get_parent_id( void ) const
{
  return _parent_id;
}

/**
 * \brief    Get cell generation
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_generation( void ) const
{
  return _generation;
}

/**
 * \brief    Get total amount of metabolites in cell
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_amount( void )
{
  return _species_list->get_amount();
}

/**
 * \brief    Get cell energy bilan
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_energy( void ) const
{
  return _energy;
}

/**
 * \brief    Get cell active boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Cell::get_active( void ) const
{
  return _active;
}

/**
 * \brief    Get cell alive boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Cell::get_alive( void ) const
{
  return _alive;
}

/**
 * \brief    Get cell x coordinate
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_x( void ) const
{
  return _x;
}

/**
 * \brief    Get cell y coordinate
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_y( void ) const
{
  return _y;
}

/**
 * \brief    Get cell score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_score( void ) const
{
  return _score;
}

/**
 * \brief    Get the number of updates
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_number_of_updates( void ) const
{
  return _number_of_updates;
}

/**
 * \brief    Get the number of divisions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_number_of_divisions( void ) const
{
  return _number_of_divisions;
}

/**
 * \brief    Return active state of the cell
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Cell::isActive( void ) const
{
  return _active;
}

/**
 * \brief    Return alive state of the cell
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Cell::isAlive( void ) const
{
  return _alive;
}

/**
 * \brief    Return cell's tag
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Cell::isTagged( void ) const
{
  return _tagged;
}

/*------------------------------------------------------------------ mutation rates */

/**
 * \brief    Get point mutation rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_point_mutation_rate( void ) const
{
  return _mutation_rates[POINT_MUTATION_RATE];
}

/**
 * \brief    Get duplication rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_duplication_rate( void ) const
{
  return _mutation_rates[DUPLICATION_RATE];
}

/**
 * \brief    Get deletion rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_deletion_rate( void ) const
{
  return _mutation_rates[DELETION_RATE];
}

/**
 * \brief    Get translocation rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_translocation_rate( void ) const
{
  return _mutation_rates[TRANSLOCATION_RATE];
}

/**
 * \brief    Get inversion rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_inversion_rate( void ) const
{
  return _mutation_rates[INVERSION_RATE];
}

/**
 * \brief    Get transition rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_transition_rate( void ) const
{
  return _mutation_rates[TRANSITION_RATE];
}

/**
 * \brief    Get breakpoint rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_breakpoint_rate( void ) const
{
  return _mutation_rates[BREAKPOINT_RATE];
}

/**
 * \brief    Get substrate tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_substrate_tag_mutation_size( void ) const
{
  return _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE];
}

/**
 * \brief    Get product tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_product_tag_mutation_size( void ) const
{
  return _mutation_rates[PRODUCT_TAG_MUTATION_SIZE];
}

/**
 * \brief    Get Kcat constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_kcat_mutation_size( void ) const
{
  return _mutation_rates[KCAT_MUTATION_SIZE];
}

/**
 * \brief    Get kcat/Km ratio constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_kcat_km_ratio_mutation_size( void ) const
{
  return _mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE];
}

/**
 * \brief    Get binding site tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_binding_site_tag_mutation_size( void ) const
{
  return _mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE];
}

/**
 * \brief    Get co-enzyme tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_co_enzyme_tag_mutation_size( void ) const
{
  return _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE];
}

/**
 * \brief    Get transcription factor tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_transcription_factor_tag_mutation_size( void ) const
{
  return _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE];
}

/**
 * \brief    Get basal expression level mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_basal_expression_level_mutation_size( void ) const
{
  return _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE];
}

/*------------------------------------------------------------------ time variables */

/**
 * \brief    Get cell birth time
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_birth_time( void ) const
{
  return _birth_time;
}

/**
 * \brief    Get cell death time
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_death_time( void ) const
{
  return _death_time;
}

/**
 * \brief    Get cell lifespan
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_lifespan( void ) const
{
  return _lifespan;
}

/*------------------------------------------------------------------ global phenotypic variables */

/**
 * \brief    Get the toxicity accumulation
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_toxicity( void ) const
{
  return _toxicity;
}

/**
 * \brief    Get the inherited transcription factors amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_inherited_TF_amount( void ) const
{
  return _inherited_TF_amount;
}

/**
 * \brief    Get the inherited enzymes amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_inherited_E_amount( void ) const
{
  return _inherited_E_amount;
}

/**
 * \brief    Get the transcription factors amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_TF_amount( void ) const
{
  return _TF_amount;
}

/**
 * \brief    Get the enzymes amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_E_amount( void ) const
{
  return _E_amount;
}

/**
 * \brief    Get the metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_inherited_metabolic_amount( void ) const
{
  return _inherited_species_list->get_amount();
}

/**
 * \brief    Get the minimum metabolic amount during cell's life
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_min_metabolic_amount( void ) const
{
  return _min_metabolic_amount;
}

/**
 * \brief    Get the maximum metabolic amount during cell's life
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_max_metabolic_amount( void ) const
{
  return _max_metabolic_amount;
}

/**
 * \brief    Get amount of metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_metabolic_uptake( void ) const
{
  return _metabolic_uptake;
}

/**
 * \brief    Get amount of released metabolites
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_metabolic_release( void ) const
{
  return _metabolic_release;
}

/**
 * \brief    Get cell minimum energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_min_energy( void ) const
{
  return _min_energy;
}

/**
 * \brief    Get cell mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_mean_energy( void ) const
{
  return _mean_energy;
}

/**
 * \brief    Get cell maximum energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_max_energy( void ) const
{
  return _max_energy;
}

/**
 * \brief    Get cell minimum score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_min_score( void ) const
{
  return _min_score;
}

/**
 * \brief    Get cell mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_mean_score( void ) const
{
  return _mean_score;
}

/**
 * \brief    Get cell maximum score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_max_score( void ) const
{
  return _max_score;
}

/**
 * \brief    Get the metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_metabolic_growth_rate( void ) const
{
  return _metabolic_growth_rate;
}

/**
 * \brief    Get the metabolic growth rate difference
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Cell::get_Dmetabolic_growth_rate( void ) const
{
  return _diff_metabolic_growth_rate;
}

/**
 * \brief    Get the number of nodes in the GRN
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_grn_nb_nodes( void ) const
{
  return _grn_nb_nodes;
}

/**
 * \brief    Get the number of edges in the GRN
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_grn_nb_edges( void ) const
{
  return _grn_nb_edges;
}

/**
 * \brief    Get the number of nodes in the metabolic network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_metabolic_nb_nodes( void ) const
{
  return _metabolic_nb_nodes;
}

/**
 * \brief    Get the number of edges in the metabolic network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Cell::get_metabolic_nb_edges( void ) const
{
  return _metabolic_nb_edges;
}

/**
 * \brief    Get the different types of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e std::vector<int>*
 */
inline std::vector<int>* Cell::get_inflowing_pumps( void )
{
  return &_inflowing_pumps;
}

/**
 * \brief    Get the different types of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e std::vector<int>*
 */
inline std::vector<int>* Cell::get_outflowing_pumps( void )
{
  return &_outflowing_pumps;
}

/**
 * \brief    Get the trophic group index in the trophic network
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int Cell::get_trophic_group( void ) const
{
  return _trophic_group;
}

/**
 * \brief    Get the trophic level
 * \details  --
 * \param    void
 * \return   \e trophic_level
 */
inline trophic_level Cell::get_trophic_level( void ) const
{
  return _trophic_level;
}

/*------------------------------------------------------------------ cell color */

/**
 * \brief    Get RGB red color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double Cell::get_red_color( void ) const
{
  return _red_color;
}

/**
 * \brief    Get RGB green color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double Cell::get_green_color( void ) const
{
  return _green_color;
}

/**
 * \brief    Get RGB blue color
 * \details  --
 * \param    double
 * \return   \e void
 */
inline double Cell::get_blue_color( void ) const
{
  return _blue_color;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*------------------------------------------------------------------ main cell variables */

/**
 * \brief    Set cell id
 * \details  --
 * \param    unsigned long long int identifier
 * \return   \e void
 */
inline void Cell::set_id( unsigned long long int identifier )
{
  _id = identifier;
}

/**
 * \brief    Set parental cell id
 * \details  --
 * \param    unsigned long long int identifier
 * \return   \e void
 */
inline void Cell::set_parent_id( unsigned long long int identifier )
{
  _parent_id = identifier;
}

/**
 * \brief    Set cell energy bilan
 * \details  --
 * \param    double energy
 * \return   \e void
 */
inline void Cell::set_energy( double energy )
{
  assert(energy >= 0.0);
  _energy = energy;
}

/**
 * \brief    Set generation
 * \details  --
 * \param    size_t generation
 * \return   \e void
 */
inline void Cell::set_generation( size_t generation )
{
  _generation = generation;
}

/**
 * \brief    Set active
 * \details  --
 * \param    bool active
 * \return   \e void
 */
inline void Cell::set_active( bool active )
{
  _active = active;
}

/**
 * \brief    Set alive
 * \details  --
 * \param    bool alive
 * \return   \e void
 */
inline void Cell::set_alive( bool alive )
{
  _alive = alive;
}

/**
 * \brief    Set x coordinate
 * \details  --
 * \param    size_t x
 * \return   \e void
 */
inline void Cell::set_x( size_t x )
{
  assert(x < _parameters->get_width());
  _x = x;
}

/**
 * \brief    Set y coordinate
 * \details  --
 * \param    size_t y
 * \return   \e void
 */
inline void Cell::set_y( size_t y )
{
  assert(y < _parameters->get_height());
  _y = y;
}

/**
 * \brief    Set cell score
 * \details  --
 * \param    double score
 * \return   \e void
 */
inline void Cell::set_score( double score )
{
  assert(score >= 0.0);
  _score = score;
}

/**
 * \brief    Set the number of updates
 * \details  --
 * \param    size_t number_of_updates
 * \return   \e void
 */
inline void Cell::set_number_of_updates( size_t number_of_updates )
{
  _number_of_updates = number_of_updates;
}

/**
 * \brief    Set the number of divisions
 * \details  --
 * \param    size_t number_of_divisions
 * \return   \e void
 */
inline void Cell::set_number_of_divisions( size_t number_of_divisions )
{
  _number_of_divisions = number_of_divisions;
}

/**
 * \brief    Cell is active
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Cell::activate( void )
{
  _active = true;
}

/**
 * \brief    Cell is inactive
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Cell::inactivate( void )
{
  _active = false;
}

/**
 * \brief    Tag the cell
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Cell::tag( void )
{
  _tagged = true;
}

/**
 * \brief    Untag the cell
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Cell::untag( void )
{
  _tagged = false;
}

/**
 * \brief    Update the number of updates
 * \details  Increment the number of updates by one
 * \param    void
 * \return   \e void
 */
inline void Cell::update_number_of_updates( void )
{
  _number_of_updates++;
}

/**
 * \brief    Update the number of divisions
 * \details  Increment the number of divisions by one
 * \param    void
 * \return   \e void
 */
inline void Cell::update_number_of_divisions( void )
{
  _number_of_divisions++;
}

/*------------------------------------------------------------------ mutation rates */

/**
 * \brief    Set point mutation rate
 * \details  --
 * \param    double point_mutation_rate
 * \return   \e void
 */
inline void Cell::set_point_mutation_rate( double point_mutation_rate )
{
  assert(point_mutation_rate >= 0.0);
  assert(point_mutation_rate <= 1.0);
  _mutation_rates[POINT_MUTATION_RATE] = point_mutation_rate;
}

/**
 * \brief    Set duplication rate
 * \details  --
 * \param    double duplication_rate
 * \return   \e void
 */
inline void Cell::set_duplication_rate( double duplication_rate )
{
  assert(duplication_rate >= 0.0);
  assert(duplication_rate <= 1.0);
  _mutation_rates[DUPLICATION_RATE] = duplication_rate;
}

/**
 * \brief    Set deletion rate
 * \details  --
 * \param    double deletion_rate
 * \return   \e void
 */
inline void Cell::set_deletion_rate( double deletion_rate )
{
  assert(deletion_rate >= 0.0);
  assert(deletion_rate <= 1.0);
  _mutation_rates[DELETION_RATE] = deletion_rate;
}

/**
 * \brief    Set translocation rate
 * \details  --
 * \param    double translocation_rate
 * \return   \e void
 */
inline void Cell::set_translocation_rate( double translocation_rate )
{
  assert(translocation_rate >= 0.0);
  assert(translocation_rate <= 1.0);
  _mutation_rates[TRANSLOCATION_RATE] = translocation_rate;
}

/**
 * \brief    Set inversion rate
 * \details  --
 * \param    double inversion_rate
 * \return   \e void
 */
inline void Cell::set_inversion_rate( double inversion_rate )
{
  assert(inversion_rate >= 0.0);
  assert(inversion_rate <= 1.0);
  _mutation_rates[INVERSION_RATE] = inversion_rate;
}

/**
 * \brief    Set transition rate
 * \details  --
 * \param    double transition_rate
 * \return   \e void
 */
inline void Cell::set_transition_rate( double transition_rate )
{
  assert(transition_rate >= 0.0);
  assert(transition_rate <= 1.0);
  _mutation_rates[TRANSITION_RATE] = transition_rate;
}

/**
 * \brief    Set breakpoint rate
 * \details  --
 * \param    double breakpoint_rate
 * \return   \e void
 */
inline void Cell::set_breakpoint_rate( double breakpoint_rate )
{
  assert(breakpoint_rate >= 0.0);
  assert(breakpoint_rate <= 1.0);
  _mutation_rates[BREAKPOINT_RATE] = breakpoint_rate;
}

/**
 * \brief    Set substrate tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_substrate_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE] = size;
}

/**
 * \brief    Set product tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_product_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[PRODUCT_TAG_MUTATION_SIZE] = size;
}

/**
 * \brief    Set Kcat constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_kcat_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[KCAT_MUTATION_SIZE] = size;
}

/**
 * \brief    Set Kcat/Km ratio constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_kcat_km_ratio_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[KCAT_KM_RATIO_MUTATION_SIZE] = size;
}

/**
 * \brief    Set binding site tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_binding_site_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[BINDING_SITE_TAG_MUTATION_SIZE] = size;
}

/**
 * \brief    Set co-enzyme tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_co_enzyme_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE] = size;
}

/**
 * \brief    Set transcription factor tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_transcription_factor_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE] = size;
}

/**
 * \brief    Set basal expression level mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Cell::set_basal_expression_level_mutation_size( double size )
{
  assert(size >= 0.0);
  _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE] = size;
}

/*------------------------------------------------------------------ time variables */

/**
 * \brief    Set cell birth time
 * \details  --
 * \param    size_t birth_time
 * \return   \e void
 */
inline void Cell::set_birth_time( size_t birth_time )
{
  _birth_time = birth_time;
}

/**
 * \brief    Set cell death time
 * \details  --
 * \param    size_t death_time
 * \return   \e void
 */
inline void Cell::set_death_time( size_t death_time )
{
  _death_time = death_time;
}

/**
 * \brief    Set cell lifespan
 * \details  --
 * \param    size_t lifespan
 * \return   \e void
 */
inline void Cell::set_lifespan( size_t lifespan )
{
  _lifespan = lifespan;
}

/*------------------------------------------------------------------ global phenotypic variables */

/**
 * \brief    Set the toxicity accumulation
 * \details  --
 * \param    double toxicity
 * \return   \e void
 */
inline void Cell::set_toxicity( double toxicity )
{
  assert(toxicity >= 0.0);
  _toxicity = toxicity;
}

/**
 * \brief    Set the minimum metabolic amount
 * \details  --
 * \param    double min_metabolic_amount
 * \return   \e void
 */
inline void Cell::set_min_metabolic_amount( double min_metabolic_amount )
{
  assert(min_metabolic_amount >= 0.0);
  _min_metabolic_amount = min_metabolic_amount;
}

/**
 * \brief    Set the maximum metabolic amount
 * \details  --
 * \param    double max_metabolic_amount
 * \return   \e void
 */
inline void Cell::set_max_metabolic_amount( double max_metabolic_amount )
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
inline void Cell::set_metabolic_uptake( double metabolic_uptake )
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
inline void Cell::set_metabolic_release( double metabolic_release )
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
inline void Cell::set_min_energy( double min_energy )
{
  assert(min_energy >= 0.0);
  _min_energy = min_energy;
}

/**
 * \brief    Set cell mean energy
 * \details  --
 * \param    double mean_energy
 * \return   \e void
 */
inline void Cell::set_mean_energy( double mean_energy )
{
  assert(mean_energy >= 0.0);
  _mean_energy = mean_energy;
}

/**
 * \brief    Set cell maximum energy
 * \details  --
 * \param    double max_energy
 * \return   \e void
 */
inline void Cell::set_max_energy( double max_energy )
{
  assert(max_energy >= 0.0);
  _max_energy = max_energy;
}

/**
 * \brief    Set cell minimum score
 * \details  --
 * \param    double min_score
 * \return   \e void
 */
inline void Cell::set_min_score( double min_score )
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
inline void Cell::set_mean_score( double mean_score )
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
inline void Cell::set_max_score( double max_score )
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
inline void Cell::set_metabolic_growth_rate( double metabolic_growth_rate )
{
  assert(metabolic_growth_rate >= 0.0);
  _metabolic_growth_rate = metabolic_growth_rate;
}

/**
 * \brief    Set the metabolic growth rate difference
 * \details  --
 * \param    double diff
 * \return   \e void
 */
inline void Cell::set_Dmetabolic_growth_rate( double diff )
{
  _diff_metabolic_growth_rate = diff;
}

/**
 * \brief    Set the number of nodes in the GRN
 * \details  --
 * \param    size_t nb_nodes
 * \return   \e void
 */
inline void Cell::set_grn_nb_nodes( size_t nb_nodes )
{
  _grn_nb_nodes = nb_nodes;
}

/**
 * \brief    Set the number of edges in the GRN
 * \details  --
 * \param    size_t nb_edges
 * \return   \e void
 */
inline void Cell::set_grn_nb_edges( size_t nb_edges )
{
  _grn_nb_edges = nb_edges;
}

/**
 * \brief    Set the number of nodes in the metabolic network
 * \details  --
 * \param    size_t nb_nodes
 * \return   \e void
 */
inline void Cell::set_metabolic_nb_nodes( size_t nb_nodes )
{
  _metabolic_nb_nodes = nb_nodes;
}

/**
 * \brief    Set the number of edges in the metabolic network
 * \details  --
 * \param    size_t nb_edges
 * \return   \e void
 */
inline void Cell::set_metabolic_nb_edges( size_t nb_edges )
{
  _metabolic_nb_edges = nb_edges;
}

/**
 * \brief    Set the trophic group index in the trophic network
 * \details  --
 * \param    unsigned long long int trophic_group
 * \return   \e void
 */
inline void Cell::set_trophic_group( unsigned long long int trophic_group )
{
  _trophic_group = trophic_group;
}

/**
 * \brief    Set the trophic level
 * \details  --
 * \param    trophic_level level
 * \return   \e void
 */
inline void Cell::set_trophic_level( trophic_level level )
{
  _trophic_level = level;
}

/*------------------------------------------------------------------ cell color */

/**
 * \brief    Set RGB red color
 * \details  --
 * \param    double red
 * \return   \e void
 */
inline void Cell::set_red_color( double red )
{
  assert(red >= 0.0);
  assert(red <= 255.0);
  _red_color = red;
}

/**
 * \brief    Set RGB green color
 * \details  --
 * \param    double green
 * \return   \e void
 */
inline void Cell::set_green_color( double green )
{
  assert(green >= 0.0);
  assert(green <= 255.0);
  _green_color = green;
}

/**
 * \brief    Set RGB blue color
 * \details  --
 * \param    double blue
 * \return   \e void
 */
inline void Cell::set_blue_color( double blue )
{
  assert(blue >= 0.0);
  assert(blue <= 255.0);
  _blue_color = blue;
}

#endif /* defined(__EVOEVO__Cell__) */
